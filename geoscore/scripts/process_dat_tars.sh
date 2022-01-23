#!/bin/bash
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

perl -e 'print "0
0.1
0.2
0.5
1
2
5
10
15
20";' > t_yr

for fcurr in `ls -1 *.tar`
do
    tar -xf $fcurr
    export dcurr=`perl -e 'my $fname=shift; $fname=~s/\.tar//; print "$fname"' -f $fcurr`
    pushd $dcurr
    echo "PROCESSING: $dcurr ..."

    bash $GPAC_SCRIPT_PATH/dats2voxel.sh \
	plot.0yr.dat \
	plot.0.1yr.dat \
	plot.0.2yr.dat \
	plot.0.5yr.dat \
	plot.1yr.dat \
	plot.2yr.dat \
	plot.5yr.dat \
	plot.10yr.dat \
	plot.15yr.dat \
	plot.20yr.dat

    cp ?pp ../
    cp pp `perl -e 'my $cdir = shift; $cdir=~s/\.t13//; print "../pp_$cdir";' -f $dcurr`

    popd
    rm -rf $dcurr
done
