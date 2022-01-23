#!/bin/bash

#ARGS (1) input file
if [ $# -ne 1 ]; then
   echo "usage: ARGS (1) input file"
   exit 1
fi 
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

if [ -z "$GPAC_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_EXEC=~/apps/gpac_trunk/src/GPAC.x
    echo "setting default for GPAC_EXEC: $GPAC_EXEC"
fi

#get the mesh
perl $GPAC_SCRIPT_PATH/gpac2fracflowmesh.pl $1 > mesh

#get the boundary conditions
perl $GPAC_SCRIPT_PATH/gpac2fracflowbc.pl $1 > Boundary

#get the prefrac file
export PFPATH=`grep PreFracPath $1 | sed 's/\s*PreFracPath\s*\=\s*//' | sed 's/"//g'`
cp $PFPATH PreFrac

#extract the parameters
perl $GPAC_SCRIPT_PATH/gpac2fracflowparams.pl $1 > Parameters

exit 0
