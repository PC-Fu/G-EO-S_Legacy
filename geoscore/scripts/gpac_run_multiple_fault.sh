#!/bin/bash

#ARGS (1) Minimum burn-in time in years (2) range of burn-in time perturbation in years (3) post burn-in simulation duration in years (4) number of realizations (5) input file
if [ $# -ne 5 ]; then
   echo "usage $0 : (1) Minimum burn-in time in years (2) range of burn-in time perturbation in years (3) post burn-in simulation duration in years (4) number of realizations (5) input file"
   exit 1
fi 

if [ -e t_yr ]; then
    echo "Using t_yr"
else
    echo "Where is t_yr?"
    exit 1 
fi

#export TBURNMIN=300.0
#export DTBURNRANGE=180.0
#export TSIM=200.0
#export N=10
export TBURNMIN=$1
export DTBURNRANGE=$2
export TSIM=$3
export N=$4

if [ -e $5 ]; then
    echo "Using input file $5"
else
    echo "Cannot find input file: $5"
    exit 1;
fi

if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpaccore/scripts
    if [ -d $GPAC_SCRIPT_PATH ]; then
        echo "$0 setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
    else
	export GPAC_SCRIPT_PATH=~/apps/geos_git/gpaca/scripts
	if [ -d $GPAC_SCRIPT_PATH ]; then
            echo "$0 setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
	else
            echo "$0 EXITING: cannot find GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
            exit 1
	fi
    fi
fi

if [ -z "$GPAC_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_EXEC=~/apps/geos_git/gpaccore/GPAC.opt
    if [ -e $GPAC_EXEC ]; then
        echo "$0 setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
    else
	export GPAC_EXEC=~/apps/geos_git/gpaca/GPAC.opt
	if [ -e $GPAC_EXEC ]; then
            echo "$0 setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
	else
            echo "$0 EXITING: cannot find GPAC_EXEC: $GPAC_EXEC"
            exit 1
	fi
    fi
fi

if [ -e mags ]; then
  rm mags
fi

export III=0
for TCURR in `perl -e 'my $n = shift; my $tburnmin = shift; my $dtburn = shift; for(my $i = 0; $i < $n; $i++) {my $dt = $tburnmin + rand() * $dtburn; print "$dt\n";}' -f $N $TBURNMIN $DTBURNRANGE`; do

    if [ "1" -eq "1" ]; then
    echo "Using tburnin = $TCURR yr"
    perl -e 'use strict;
my $tburnin = shift;
my $dtsim = shift;
my $filename = shift;
my $endtime = $dtsim + $tburnin;
open(IN,"<$filename");
while(<IN>)
{
  if(/begintime/)
  {
    print "    <Application name=\"1\" begintime=\"0.0\" endtime=\"$endtime year\">\n";
  }
  else
  {
    print $_;
  }
}
close IN;
' -f $TCURR $TSIM $5 > runcurr.xml

    perl -e 'use strict;
my $tcurr = shift;
my $fct = 365.25 * 24 * 3600;

$tcurr *= $fct;
print "0\n";
open(IN,"<t_yr");
while(<IN>)
{
  my $t = $_;
  $t *= $fct;
  $t += $tcurr;
  print "$t\n";
}
close IN;' -f $TCURR > t

    rm -f ms_*
    $GPAC_EXEC -i runcurr.xml > cout
    fi
    grep ^eq cout | awk  -v dt=$TCURR '{t = $2 / (365.25 * 24 * 3600); t-=dt; if( $3 > -10) {print t" "$3}}' >> mags
    mv cout cout$III
    tar -czf tmp$III.tar.gz ms_* cout$III
    perl $GPAC_SCRIPT_PATH/events2points.pl cout$III > pts$III
    perl $GPAC_SCRIPT_PATH/points2xmdv.pl pts$III x y z t M E A d > pts$III.okc
    III=$((III + 1))
done
