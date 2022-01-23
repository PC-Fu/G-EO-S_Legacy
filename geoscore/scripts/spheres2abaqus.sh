#!/bin/bash

#ARGS (1) sphere file (herbold format) (2) unit sphere mesh (3) output ABAQUS file
if [ $# -ne 3 ]; then
   echo "usage: ARGS (1) sphere file (herbold format) (2) unit sphere mesh (3) output ABAQUS file"
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

#REQUIRED SCRIPTS
#(1) abaqus_combine.pl
#(2) abaqus_affine_transform.pl
#(3) affine2randomorientation.pl

rm tmp*.inp
COUNTER=1
while [ $COUNTER -le `head -n 1 $1` ]
do
    echo $COUNTER
    let COUNTER=COUNTER+1
    perl $GPAC_SCRIPT_PATH/abaqus_affine_transform.pl $2 `head -n $COUNTER  $1 | tail -n 1 | awk '{print $3" "$4" "$5" "$2" 0 0 0 "$2" 0 0 0 "$2}' | perl $GPAC_SCRIPT_PATH/affine2randomorientation.pl` > tmp$COUNTER.inp
    #perl $GPAC_SCRIPT_PATH/abaqus_affine_transform.pl $2 `head -n $COUNTER $1 | tail -n 1 | awk '{print $3" "$4" "$5" "$2" 0 0 0 "$2" 0 0 0 "$2}'` > tmp$COUNTER.inp
done
perl $GPAC_SCRIPT_PATH/abaqus_combine.pl tmp*.inp > $3
rm tmp*.inp

exit
