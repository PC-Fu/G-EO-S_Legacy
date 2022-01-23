#!/bin/bash

if [ $# -ne 7 ]; then
    echo "usage: PTS_FILE ZTOP ZBOTTOM NSTRIKE NDIP XORIGIN YORIGIN"
    exit 1
fi

if [ -z $GPAC_SCRIPT_PATH ]; then
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpaca/scripts
    if [ -d $GPAC_SCRIPT_PATH ]; then
	echo "GPAC_SCRIPT_PATH $GPAC_SCRIPT_PATH"
    else
	echo "GPAC_SCRIPT_PATH not found: $GPAC_SCRIPT_PATH"
	exit 1
    fi
fi

export FNAME=$1
export ZTOP=$2
export ZBOTTOM=$3
export NSTRIKE=$4
export NDIP=$5
export X0=$6
export Y0=$7

export WIDTH=$((ZTOP-ZBOTTOM))
echo "Fault width: $WIDTH"
 
awk -v x0=$X0 -v y0=$Y0 '{print ($1 - x0)" "($2 - y0)" "$3}' $FNAME > pts.tmp
#cat pts_greeley.txt | awk '{print ($1 - 260000)" "($2 - 3887000)" "$3}' > fg.txt

export STRIKE=`perl $GPAC_SCRIPT_PATH/points2strike.pl pts.tmp | awk '{print (360+$1)}'`
export ALEN=`perl $GPAC_SCRIPT_PATH/points2strike.pl pts.tmp | awk '{print $2}'`
echo "Strike: $STRIKE  Fault length: $ALEN"

bash $GPAC_SCRIPT_PATH/fault2abaqus.sh $STRIKE 90 $ALEN $WIDTH $NSTRIKE $NDIP 0 0 0
mv fault.geom tmplt.tmp

perl $GPAC_SCRIPT_PATH/abaqus_to_path.pl tmplt.tmp pts.tmp $ZTOP $ZBOTTOM 0 2 > fault.geom
perl $GPAC_SCRIPT_PATH/abaqus2points.pl fault.geom > pts_mesh.tmp
