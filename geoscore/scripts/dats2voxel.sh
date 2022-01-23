#!/bin/bash

#ARGS (1-N) input files
if [ $# -lt 1 ]; then
   echo "usage: ARGS (1-N) input files"
   exit 1
fi 
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpaccore/scripts
    if [ -d $GPAC_SCRIPT_PATH ]; then
	echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
    else
	echo "EXITING: cannot find GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
	exit 1
    fi
fi

if [ -e tmp.pl ]; then
    rm tmp.pl
fi

for cdat in $@
do
    if [ -e tmp.pl ]; then
	tail -n `cat $cdat | wc | awk '{print ($1 - 1)}'` $cdat | awk '{print (6894.7573*$13)}' >> pp
    else
	bash $GPAC_SCRIPT_PATH/dat2voxel.sh $cdat
	perl -e 'my $nt = shift;
my $nx = shift;
my $ny = shift;
my $nz = shift;
print "$nx $ny $nz $nt\n"' -f $# `cat xpp | wc | awk '{print $1}'` `cat ypp | wc | awk '{print $1}'` `cat zpp | wc | awk '{print $1}'` > pp
	cat ppcurr | awk '{print (6894.7573*$13)}' >> pp
    fi
done
