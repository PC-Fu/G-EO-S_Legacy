#!/bin/bash

if [ $# -ne 6 ]; then
    echo "usage: $0 <strike maximum stress direction> <strike fault> <dip> <pp at z=0> <pp at z=z0> <z0>"
    exit 1
fi

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpaccore/scripts
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

perl $GPAC_SCRIPT_PATH/stress2voxel.pl sigHhv $1 $2 $3 2>&1 | tee str_out
mv x xs
mv y ys
mv z zs

perl -e 'open(IN,"<str_out");
my $line = "";
while(<IN>)
{
  if(/rakeVector/)
  {
    $line = $_;
    last;
  }
}
close IN;

open(IN,"<run.xml");
while(<IN>)
{
  if(/rakeVector/)
  {
    print $line;
  }
  else
  {
    print $_;
  }
}
close IN;' > tmp
mv tmp run.xml

tail -n `wc zs | awk '{print $1}'` tau > tmp
paste -d ' ' zs tmp > tmp2
tail -n `wc zs | awk '{print $1}'` sigma > tmp
paste -d ' ' tmp2 tmp > z_tau_sig 
rm tmp tmp2

head -n 1 sigma > mu
perl -e 'open(IN,"<z_tau_sig");
my $pp0 = shift;
#4.2e6;
my $pp6000 = shift;
#6.02514e7;
my $z0 = shift;
#-6000
my $dppdz = ($pp6000 - $pp0)/$z0;
while(<IN>)
{
  my $line = $_;
  chomp($line);
  my @arr  = split(/\s+/,$line);
  my $z = $arr[0];
  my $tau = $arr[1];
  my $sigma = $arr[2];
  my $ppcurr = $pp0 + $dppdz * $arr[0];
  my $mu_crit = $tau / ($sigma - $ppcurr);
  my $mu = $mu_crit;
  $mu += 0.1 * rand();
  if($mu < 0.2)
  {
    $mu = 0.2;
  }
  print "$mu\n";
}
' -f $4 $5 $6 >> mu
