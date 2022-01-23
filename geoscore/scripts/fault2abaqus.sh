#!/bin/bash

#ARGS (1) strike-degrees (2) dip-degrees (3) strike-length (4) dip-width (5) nstrike (6) ndip (7) x-origin (East) (8) y-origin (North) (9) z-top (Up=+z)
if [ $# -ne 9 ]; then
   echo "usage (output to fault.inp): ARGS (1) strike (2) dip (3) strike-length (4) dip-width (5) nstrike (6) ndip (7) x-origin (8) y-origin (9) z-top : $1"
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
# CREATE CUBIT PYTHON SCRIPT
#####################################################

perl -e 'use strict;
use Math::Trig;

my $strike = shift;
my $dip = shift;

my $nstrike = shift;
my $ndip = shift;

print "import cubit
my_cmd = \"create surface rectangle width $ndip height $nstrike zplane\"
cubit.cmd(my_cmd)
my_cmd = \"surface 1 size 1\"
cubit.cmd(my_cmd)
my_cmd = \"surface 1 scheme map\"
cubit.cmd(my_cmd)
my_cmd = \"mesh surface 1\"
cubit.cmd(my_cmd)
my_cmd = \"block 1 volume 1\"
cubit.cmd(my_cmd)
";
if($dip!=0)
{
    print "my_cmd = \"rotate volume 1 about y angle $dip\"
cubit.cmd(my_cmd)
";
}
if($strike!=0)
{
  print "my_cmd = \"rotate volume 1 about z angle -$strike\"
cubit.cmd(my_cmd)
";
}
print "my_cmd = \"export Abaqus \\\"fault_untransformed.inp\\\" Block 1 dimension 3 everything overwrite cubitids\"
cubit.cmd(my_cmd)
my_cmd = \"quit\"
cubit.cmd(my_cmd)
";' -f  $1 $2 $5 $6 > tmp.py

#####################################################
# RUN CUBIT
#####################################################

/usr/gapps/cubit/linux64.14.0/cubit -nographics tmp.py

#####################################################
# TRANSFORM FAULT
#####################################################

#SCALE--->
perl -e 'use strict;
use Math::Trig;
my $strike = shift;
my $dip = shift;

my $lstrike = shift;
my $ldip = shift;
my $nstrike = shift;
my $ndip = shift;

my @axis = (0,0,1);

die "dip ($dip) must be between 0 and 90\n" if($dip < 0 || $dip > 90);
my $fct = atan(1) / 45.0;
if($dip != 90)
{
  $dip *= $fct;
  $strike *= $fct;

  #my @xstrike = (sin($strike), cos($strike),0);
  $axis[0] = cos($strike)*cos($dip);
  $axis[1] = -sin($strike)*cos($dip);
  $axis[2] = sin($dip);
} else {
  $dip *= $fct;
}

#print "ldip: $ldip\n";
my $ascale = ($ldip / $ndip) * sin($dip);
my $rscale = ($lstrike / $nstrike);

#die "$lstrike $nstrike $ldip $ndip\n";

print "$ascale $rscale $axis[0] $axis[1] $axis[2]";' -f $1 $2 $3 $4 $5 $6 > tmp

echo "TRANSFORMING FAULT USING:"
cat tmp
echo " "

perl $GPAC_SCRIPT_PATH/abaqus_cylindrical_transform.pl fault_untransformed.inp `cat tmp` > tmp.inp

#TRANSLATE--->
perl $GPAC_SCRIPT_PATH/abaqus2points.pl tmp.inp > pts
perl $GPAC_SCRIPT_PATH/points2minmax.pl pts > tmp
perl -e 'use strict;
my $dx = shift;
my $xold = shift;
$dx -= $xold;
my $dy = shift;
my $yold = shift;
$dy -= $yold;
my $dz = shift;
my $zold = shift;
$dz -= $zold;
print "$dx $dy $dz";
' -f $7 `grep ^0 tmp | awk '{print $2}'` $8 `grep ^1 tmp | awk '{print $2}'` $9 `grep ^2 tmp | awk '{print $3}'` > tmp
#note: we just took lowest x, lowest y, and then HIGHEST z

perl $GPAC_SCRIPT_PATH/abaqus_affine_transform.pl tmp.inp `cat tmp` > fault.geom
