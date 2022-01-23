#!/bin/bash

if [ $# -ne 2 ]; then
   echo "usage: ARGS: (1) BLOCK FILE (2) ABAQUS FILE"
   exit 1
fi 


#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    if [ -d "$GPAC_SCRIPT_PATH" ];then
	echo "spheres2abaqus: Using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
    else
	export GPAC_SCRIPT_PATH=/usr/gapps/GEOS/scripts
	if [ -d "$GPAC_SCRIPT_PATH" ];then
	    echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
	else
	    echo "spheres2abaqus.sh cannot find GPAC_SCRIPT_PATH : $GPAC_SCRIPT_PATH"
	    exit 1
	fi
    fi
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

if [ -z "$QHULL_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export QHULL_EXEC=~/apps/gpac_trunk/scripts/qhull
    echo "setting default for QHULL_EXEC: $QHULL_EXEC"
fi

#REQUIRED SCRIPTS:
#  ldec_geomverts_to_files.pl
#  ldec_geomverts_to_pts.pl
#  points2qhull.pl
#  <exec: qhull>
#  off2abaqus.pl

#split the ldec file into a file per block
perl $GPAC_SCRIPT_PATH/ldec_geomverts_to_files.pl $1 aaaaaa

#for each block convert to OFF via pts then qhull
for f in `ls -1 aaaaaa??????`
do
    perl $GPAC_SCRIPT_PATH/ldec_geomverts_to_pts.pl $f > pts
    perl $GPAC_SCRIPT_PATH/points2qhull.pl pts | $QHULL_EXEC -o > $f.off
done

#combine all OFF files into a single Abaqus file
perl $GPAC_SCRIPT_PATH/off2abaqus.pl aaaaaa??????.off > $2

#clean up
rm pts
