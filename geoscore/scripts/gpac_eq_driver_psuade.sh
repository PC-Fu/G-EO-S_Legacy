#!/bin/bash

#ARGS (1) psuade input file (2) psuade output file (3) index (4) ??? (5) ???
if [ $# -ne 5 ]; then
   echo "usage: ARGS (1) psuade input file (2) psuade output file (3) index (4) ??? (5) ???"
   exit 1
fi
#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_EXEC=~/apps/gpac_trunk/src/GPAC.x
    if [ -a "$GPAC_EXEC" ];then
	echo "$0 using GPAC_EXEC $GPAC_EXEC"
    else
	echo "$0 cannot find GPAC_EXEC $GPAC_EXEC"
	exit 1
    fi
fi
if [ -z "$PSUADE_INPUT_FILE" ]; then
    #NOTE: this only provides a value for the current script
    export PSUADE_INPUT_FILE=$PWD/psuade.in
    if [ -a "$PSUADE_INPUT_FILE" ];then
	echo "$0 using PSUADE_INPUT_FILE $PSUADE_INPUT_FILE"
    else
	echo "$0 cannot find PSUADE_INPUT_FILE $PSUADE_INPUT_FILE"
	exit 1
    fi
fi
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

#########################################################
#
#    EXTRACT VARIABLES FROM PSUADE INPUT
#
#########################################################
perl $GPAC_SCRIPT_PATH/extract_psuade_variables.pl "$PSUADE_INPUT_FILE" $1 > gpac_vars

#########################################################
#
#    GENERATE GPAC INPUT FILES
#
#########################################################
bash $GPAC_SCRIPT_PATH/gpac_eq_driver.sh

#########################################################
#
#    RUN GPAC
#
#########################################################
$GPAC_EXEC -i run.xml | grep ^eq > eqs

#########################################################
#
#    ANALYZE OUTPUT
#
#########################################################
echo "No function yet to process eqs into the expected PSUADE output: $2"

#########################################################
#
#    OUTPUT DUMMY FOR PSUADE
#
#########################################################
echo "-999" > $2
