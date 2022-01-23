#!/usr/bin/bash
if [ $# -ne 3 ]; then
   echo "usage: ARGS (1) file name of abaqus mesh to enclose (2) wall thickness (3) output boundary mesh file name"
   exit 1
fi 

if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi
if [ -z "$CUBIT_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export CUBIT_EXEC=/usr/gapps/cubit/linux64.13.1/cubit
    echo "setting default for CUBIT_EXEC: $CUBIT_EXEC"
fi

#get nodal points of mesh to enclose
perl $GPAC_SCRIPT_PATH/abaqus2points.pl $1 > gpts

#get Cartesian min/max
perl $GPAC_SCRIPT_PATH/points2minmax.pl gpts > gpts_mm

#make cubit input file
#echo "maximum size for $1 is: "
#echo `perl $GPAC_SCRIPT_PATH/abaqus2maximum_face_size.pl $1`
perl $GPAC_SCRIPT_PATH/cubit_wall.pl gpts_mm `perl $GPAC_SCRIPT_PATH/abaqus2maximum_face_size.pl $1` $2 $3 > tmp.py

#run cubit
$CUBIT_EXEC -nographics tmp.py
