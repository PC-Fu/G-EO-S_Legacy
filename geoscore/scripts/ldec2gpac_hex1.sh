#!/bin/bash

#ARGS: (1) BLOCK FILE (2) ABAQUS FILE (3) ELEMENT SIZE
if [ $# -ne 3 ]; then
   echo "usage: ARGS: (1) BLOCK FILE (2) ABAQUS FILE (3) ELEMENT SIZE"
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

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$QHULL_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export QHULL_EXEC=~/ldecff/trunk/libs/qhull/src/qhull3.1/src/qhull
    echo "setting default for QHULL_EXEC: $QHULL_EXEC"
fi
if [ -z "$CUBIT_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export CUBIT_EXEC=/usr/gapps/cubit/linux64.13.0/cubit
    echo "setting default for CUBIT_EXEC: $CUBIT_EXEC"
fi

#REQUIRED SCRIPTS:
#  ldec_geomverts_to_files.pl
#  ldec_geomverts_to_pts.pl
#  points2qhull.pl
#  <exec: qhull>
#  off2abaqus.pl
#  abaqus2cubit_hex_volume.pl
#  abaqus2tied_node.pl
#  abaqus_combine.pl
# <exec: cubit>

#split the ldec file into a file per block
perl $GPAC_SCRIPT_PATH/ldec_geomverts_to_files.pl $1 aaaaaa

#for each block convert to OFF via pts then qhull
for f in `ls -1 aaaaaa??????`
do
    perl $GPAC_SCRIPT_PATH/ldec_geomverts_to_pts.pl $f > pts
    perl $GPAC_SCRIPT_PATH/points2qhull.pl pts | $QHULL_EXEC -o > $f.off
    perl $GPAC_SCRIPT_PATH/off2abaqus.pl $f.off > $f.inp
    perl $GPAC_SCRIPT_PATH/abaqus2cubit_hex_volume.pl $3 $f.inp tmp.inp > tmp.py
    $CUBIT_EXEC  -nographics tmp.py
    perl $GPAC_SCRIPT_PATH/abaqus2tied_node.pl tmp.inp > $f.inp
done

#combine all OFF files into a single Abaqus file
perl $GPAC_SCRIPT_PATH/abaqus_combine.pl aaaaaa??????.inp > $2

#clean up
rm pts
