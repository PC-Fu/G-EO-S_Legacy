#!/bin/bash

if [ $# -ne 2 ]; then
   echo "usage: ARGS (1) model generation input file (2) sphere output file"
   exit 1
fi 

if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi
if [ -z "$DEM_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export DEM_EXEC=~/apps/dem3d/dem.o3
    echo "setting default for DEM_EXEC: $DEM_EXEC"
fi

#set up CWD
rm -rf sphere_pack
mkdir sphere_pack
cd sphere_pack
cp ../$1 input_spheres

#run packing algorithms
$DEM_EXEC input_spheres > prunout
perl $GPAC_SCRIPT_PATH/dem3d_output2spheres.pl `grep ^write input_spheres | awk '{print $2}'` | awk '{print FNR" "$1" "$2" "$3" "$4}' > tmp
cat tmp | wc -l > $2
cat tmp >> $2
rm tmp
cp $2 ../
cd ../
