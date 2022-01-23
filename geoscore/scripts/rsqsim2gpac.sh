#!/bin/bash

#ARGS (1) rsqsim fault file (2) rsqsim pore pressure file (3) nstrike (4) ndip (5) rsqsim input file
if [ $# -ne 5 ]; then
   echo "usage: ARGS (1) fault file (2) pore pressure file (3) number of elements in strike direction (4) number of elements in dip direction (5) rsqsim input file"
   exit 1
fi 

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    if [ -d "$GPAC_SCRIPT_PATH" ];then
	echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
    else
	export GPAC_SCRIPT_PATH=/usr/gapps/GEOS/scripts
	if [ -d "$GPAC_SCRIPT_PATH" ];then
	    echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
	else
	    echo "$0 cannot find GPAC_SCRIPT_PATH $GPAC_SCRIPT_PATH"
	    exit 1
	fi
    fi
fi

#####################################################
# EXTRACT FAULT GEOMETRY FROM RSQSIM FILE
#####################################################

#extract the fault (only works for planar, 90-degree dipping faults!)
perl -e 'use strict;
use Math::Trig;
my @x0 = (0,0,0,0,0,0); 
my $ii = 0; 
for(my $i = 0; $i < 3; $i++) 
{ 
  my $dummy = shift; 
  $x0[$ii] = shift; 
  ++$ii; 
  $x0[$ii] = shift; 
  ++$ii;
}

print STDERR "0: $x0[0] 1: $x0[1] 2: $x0[2] 3: $x0[3] 4: $x0[4] 5: $x0[5]\n";
my $nstrike = shift;
my $ndip = shift;

my $xx = $x0[1] - $x0[0]; 
my $yy = $x0[3] - $x0[2]; 
my $rad = atan2($yy, $xx);
my $degree = $rad * 45 / atan(1);
print STDERR "dx: $xx dy: $yy radians: $rad degrees: $degree\n";

$xx *= $xx; 
$yy *= $yy; 
my $length = sqrt($xx + $yy);

my $dr = $length / ($nstrike-1);
my $dx = $dr * 0.5 * cos($rad);
my $dy = $dr * 0.5 * sin($rad);
my $dz = 0.5 * ($x0[4] - $x0[5]) / ($ndip-1);

$x0[0] -= $dx;
$x0[1] += $dx;
$x0[2] -= $dy;
$x0[3] += $dy;
$x0[4] += $dz;
$x0[5] -= $dz;

$dz = ($x0[4] - $x0[5]) / $ndip;
$dr = sqrt(($x0[0] - $x0[1])*($x0[0] - $x0[1]) + ($x0[2] - $x0[3]) * ($x0[2] - $x0[3])) / $nstrike;
print "$dz $dr";

$dz = -0.5 * $ndip + ($x0[4] / $dz);
$dx = $x0[0] / $dr;
$dy = $x0[2] / $dr;
$dr = $nstrike * 0.5;

die "cannot open tmp.py" if(!open(OUT, ">tmp.py"));
print OUT "import cubit
my_cmd = \"create surface rectangle width $nstrike height $ndip yplane\"
cubit.cmd(my_cmd)
my_cmd = \"surface 1 size 1\"
cubit.cmd(my_cmd)
my_cmd = \"surface 1 scheme map\"
cubit.cmd(my_cmd)
my_cmd = \"mesh surface 1\"
cubit.cmd(my_cmd)
my_cmd = \"block 1 volume 1\"
cubit.cmd(my_cmd)
my_cmd = \"move volume 1 x $dr z $dz\"
cubit.cmd(my_cmd)
";
if($degree!=0)
{
  print OUT "my_cmd = \"rotate volume 1 about z angle $degree\"
  cubit.cmd(my_cmd)
";
}
print OUT "my_cmd = \"move volume 1 x $dx y $dy\"
cubit.cmd(my_cmd)
my_cmd = \"export Abaqus \\\"fault_untransformed.inp\\\" Block 1 dimension 3 everything overwrite cubitids\"
cubit.cmd(my_cmd)
my_cmd = \"quit\"
cubit.cmd(my_cmd)\n";
close OUT;' -f `perl $GPAC_SCRIPT_PATH/points2minmax.pl $1 | head -n 1` `perl $GPAC_SCRIPT_PATH/points2minmax.pl $1 | head -n 2 | tail -n 1` `perl $GPAC_SCRIPT_PATH/points2minmax.pl $1 | head -n 3 | tail -n 1` $3 $4 > tmp

#####################################################
# RUN CUBIT
#####################################################

/usr/gapps/cubit/linux64.14.0/cubit -nographics tmp.py

#####################################################
# TRANSFORM FAULT
#####################################################

perl $GPAC_SCRIPT_PATH/abaqus_cylindrical_transform.pl fault_untransformed.inp `cat tmp` > fault.geom
perl $GPAC_SCRIPT_PATH/abaqus2points.pl fault.geom > pts
#exit

#####################################################
# EXTRACT PORE PRESSURE
#####################################################


perl -e 'use strict;
my $filename = shift;
my $tburnin = 0;
die "cannot open $filename\n" if(!open(IN,"<$filename"));
while(<IN>)
{
  if(/tBurnIn/)
  {
    my @arr = split(/\=/);
    $tburnin = $arr[1] * 1.0;
    print STDERR "TBURNIN: $tburnin ($#arr)\n";
  }
}
close IN;
print "$tburnin";
' -f $5 > tburnin

perl -e 'use strict;
my $pp = shift;
my $tburnin = shift;
die "cannot open $pp\n" if(!open(IN,"<$pp"));
my $line = <IN>;
$line =~s/^\s+//;
chomp($line);

my @arr = split(/\s+/,$line);

#21 110 85 30 50 50 80 0 0 -2980
{
  die "cannot open x\n" if(!open(OUT,">x"));
  print STDERR "nx: $arr[1] x0: $arr[7] dx: $arr[4]\n";
  for(my $i = 0; $i < $arr[1]; $i++)
  {  
    my $x = $arr[7] + $i * $arr[4];
    print OUT "$x\n";
  }
  close OUT;
}
{
  die "cannot open y\n" if(!open(OUT,">y"));
  print STDERR "ny: $arr[2] y0: $arr[8] dy: $arr[5]\n";
  for(my $i = 0; $i < $arr[2]; $i++)
  {  
    my $x = $arr[8] + $i * $arr[5];
    print OUT "$x\n";
  }
  close OUT;
}
{
  die "cannot open z\n" if(!open(OUT,">z"));
  print STDERR "nz: $arr[3] z0: $arr[9] dz: $arr[6]\n";
  for(my $i = 0; $i < $arr[3]; $i++)
  {  
    my $x = $arr[9] + $i * $arr[6];
    print OUT "$x\n";
  }
  close OUT;
}

my $nels = $arr[1] * $arr[2] * $arr[3];
print STDERR "nels: $nels\n";
{
  die "cannot open t\n" if(!open(OUT, ">t"));
  die "cannot open voxel\n" if(!open(VOUT, ">voxel"));

  my @matrix = ();
  for(my $ii = 0; $ii < $arr[1]; ++$ii)
  {
    my @yyy = ();
    for(my $jj = 0; $jj < $arr[2]; ++$jj)
    {
      my @zzz = ();
      for(my $kk = 0; $kk < $arr[3]; ++$kk)
      {
        $zzz[$kk] = 0.0;
      }
      push(@yyy,\@zzz);
    }
    push(@matrix,\@yyy);
  } 

  my $repeat = 0; 
  if($tburnin > 0)
  {
    print OUT "0\n";
    $arr[0] += 1;
    print STDERR " --> porepressure: using tburnin = $tburnin\n";
    $repeat = 1;
  }
  print VOUT "$arr[1] $arr[2] $arr[3] $arr[0]\n";
  while(<IN>)
  {
    my $time = $_;
    $time += $tburnin;
    print OUT "$time\n";
    for(my $ii = 0; $ii < $nels; $ii++)
    {
      my $line = <IN>;
      $line=~s/^\s+//;
      chomp($line);
      #die "FIRST LINE: $line\n";
      my @arr2 = split(/\s+/,$line);
      $matrix[$arr2[0]][$arr2[1]][$arr2[2]] = $arr2[3];
    }
    for(my $rr = 0; $rr <= $repeat; $rr++)
    {
      for(my $kk = 0; $kk < $arr[3]; ++$kk)
      {
        for(my $jj = 0; $jj < $arr[2]; ++$jj)
        {
          for(my $ii = 0; $ii < $arr[1]; ++$ii)
          {
            my $ppval = $matrix[$ii][$jj][$kk];
            $ppval *= 1e6;
            print VOUT "$ppval\n";
          }
        }
      }
    }
    $repeat = 0;
  }
  close OUT;
  close VOUT;
}
' -f $2 `cat tburnin`


#####################################################
# EXTRACT INPUT XML SNIPPET FROM RSQSIM INPUT FILE
#####################################################

perl -e 'use strict;
my $filename = shift;
my $nx = shift;
my $ny = shift;
my $nz = shift;
my $tburnin = shift;
my $nel = $nx * $ny * $nz;
my %maxVals = ();
my %minVals = ();
my %otherVals = ();
print STDERR "XML snippet creation is using file $filename: nx=$nx ny=$ny nz=$nz\n";
die "cannot open $filename\n" if(!open(IN,"<$filename"));
while(<IN>)
{
  my $line = $_;
  chomp($line);
  $line=~s/^\s+//;
  my @arr = split(/\=/,$line);
  if($#arr == 1)
  {
    $arr[0]=~s/\s+$//;
    $arr[0]=~s/^\s+//;
    $arr[1]=~s/\s+$//;
    $arr[1]=~s/^\s+//;
    if($arr[0]=~m/Min/)
    {
      $arr[0]=~s/Min//;
      print STDERR "rsqsim: minVals{$arr[0]} = $arr[1]\n";
      $minVals{$arr[0]} = $arr[1];
    }
    elsif($arr[0]=~m/Max/)
    {
      $arr[0]=~s/Max//;
      print STDERR "rsqsim: maxVals{$arr[0]} = $arr[1]\n";
      $maxVals{$arr[0]} = $arr[1];
    }
    else
    {
      print STDERR "rsqsim: otherVals{$arr[0]} = $arr[1]\n";
      $otherVals{$arr[0]} = $arr[1];
    }
  }
}
close IN;

#FIRST, PRINT INPUT FILES
my $key = "";
my @z = ();
die "cannot open z\n" if(!open(IN,"<z"));
while(<IN>)
{
  my $line = $_;
  $line *= 1.0;
  push(@z, $line);
}
close IN;


foreach $key (keys %minVals)
{
  if(exists $maxVals{$key})
  {
    print STDERR "Creating include file: $key\n";
    my $minVal = $minVals{$key};
    my $maxVal = $maxVals{$key};
    my $ddz = 0;
    my $fct = 1.0;
    if($key=~m/sigma/)
    {
      $fct = 1e6;
      if(exists $otherVals{"sigmaZGrad"})
      {
        $ddz = $otherVals{"sigmaZGrad"};
      }
      print STDERR " --> sigma: using dsigmadz = $ddz\n";
    }
    elsif($key=~m/tau/)
    {
      $fct = 1e6;
    }
    die "cannot open $key" if(!open(OUT,">$key"));
    print OUT "$nx $ny $nz\n";
    for(my $k = 0; $k < $nz; $k++)
    {
      my $dval = $ddz * $z[$k];
      for(my $j = 0; $j < $ny; $j++)
      {
        for(my $i = 0; $i < $nx; $i++)
        {
          my $val = $dval + $minVal + ($maxVal - $minVal) * rand();
          $val *= $fct;
          print OUT "$val\n";
        }
      }
    }
    close OUT;
  }
}

#SECOND, PRINT CONSTANT INITIAL VALUES
if(exists $otherVals{"stressOvershootFactor"})
{
  print STDERR "Creating include file: stressOverShootFactor\n";
  my $sof = $otherVals{"stressOvershootFactor"};
  open(OUT,">stressOverShootFactor");
  print OUT "$nx $ny $nz\n";
  for(my $i = 0; $i < $nel; $i++)
  {
    print OUT "$sof\n";
  }
  close OUT;
}
if(exists $otherVals{"fA"})
{
  print STDERR "Creating include file: AReductionFactor\n";
  my $ared = $otherVals{"fA"};
  open(OUT,">AReductionFactor");
  print OUT "$nx $ny $nz\n";
  for(my $i = 0; $i < $nel; $i++)
  {
    print OUT "$ared\n";
  }
  close OUT;
}

#THIRD, PRINT INPUT XML SNIPPET
print STDERR "Creating XML snippet\n";
{
  my $slipRateEQ = 1.0;
  if(exists $otherVals{"ddotEQ1"})
  {
    $slipRateEQ = $otherVals{"ddotEQ1"};
    print STDERR " --> rsqsim: ddotEQ1 = $slipRateEQ\n";
  }
  my $shearmod = 3e10;
  if(exists $otherVals{"lameMu"})
  {
    $shearmod = $otherVals{"lameMu"} * 1e6; 
    print STDERR " --> rsqsim: shearmod = $shearmod\n";
  }
  my $bulkmod = 5e10;
  if(exists $otherVals{"lameLambda"})
  {
    $bulkmod = $otherVals{"lameLambda"} * 1e6; 
    $bulkmod += 2 * $shearmod / 3;
    print STDERR " --> rsqsim: bulkmod = $bulkmod\n";
  }
  my $maxT = 1.5e11 + $tburnin;
  if(exists $otherVals{"maxT"})
  {
    $maxT = $otherVals{"maxT"};
    print STDERR " --> rsqsim: maxT = $maxT\n";
  }
  $maxT /= (365.25 * 24 * 3600);
  open(OUT,">tmp.xml");
  print OUT "  <Solvers>
    <FaultRuptureBEMSolver name=\"solver1\"
            upVector = \"0 0 1\"
            rakeVector = \"0 1 0\"
            maximumHorizontalStressDirection = \"1 0 0\"
            creepFrictionRate = \"0.001\"
            slipRateEQ = \"$slipRateEQ\"
            porePressureTableName = \"pp\">
      <Material ShearModulus = \"$shearmod\" BulkModulus = \"$bulkmod\"/>
    </FaultRuptureBEMSolver>
  </Solvers>

  <SolverApplications>
    <Application name = \"1\" begintime = \"0.0\" endtime = \"$maxT year\">
      <Apply solver = \"solver1\" toregions = \"EB1\"/>
    </Application>
  </SolverApplications>
";
  close OUT;
}
' -f $5 `cat x | wc | awk '{print $1}'` `cat y | wc | awk '{print $1}'` `cat z | wc | awk '{print $1}'` `cat tburnin`
