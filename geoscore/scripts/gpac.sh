#!/bin/bash

#ARGS (1) input file
if [ $# -ne 2 ]; then
   echo "usage: ARGS (1) input file (2) dimensions [2 or 3]"
   exit 1
fi 

if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    #export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    export GPAC_SCRIPT_PATH=~/GPAC/gpac/trunk/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi
if [ -z "$GPAC_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_EXEC=~/apps/gpac_trunk/src/GPAC.x
    echo "setting default for GPAC_EXEC: $GPAC_EXEC"
fi
if [ -z "$FF_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export FF_EXEC=~/apps/fracflow_v1/ff
    echo "setting default for FF_EXEC: $FF_EXEC"
fi

if [ "$2" -eq "3" ]; then
    $GPAC_EXEC -i $1
    exit 0
fi
if [ "$2" -eq "2" ]; then
    bash $GPAC_SCRIPT_PATH/gpac2fracflow.sh $1
    $FF_EXEC
    exit 0
fi

exit 1
