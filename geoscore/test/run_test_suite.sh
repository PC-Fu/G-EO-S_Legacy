#!/bin/bash

ATS_LEVEL=200
if [ $# -ne 1 ];then
    echo "assuming level 200 and below to test"
else
    ATS_LEVEL=$1
fi

if [ -z "$GPAC_EXEC" ];then
    GPAC_EXEC=$PWD/../GPAC.x
    if [ -e "$GPAC_EXEC" ]; then
	echo "Setting default GPAC path to $GPAC_EXEC"
    else
	echo "Cannot find default GPAC path: $GPAC_EXEC"
	exit 1
    fi
fi

if [ -z "$ATS_EXEC" ];then
    ATS_EXEC=/usr/gapps/ats/chaos_5_x86_64_ib/bin/ats
    if [ -e "$ATS_EXEC" ]; then
	echo "Setting default ATS path to $ATS_EXEC"
    else
	echo "Cannot find default ATS path: $ATS_EXEC"
	exit 1
    fi
#    echo "you must define the path to the ATS executable (under external_libs in local install or /usr/gapps/ats/chaos_5_x86_64_ib/bin/ats on LC)"
#    exit 1
fi

#if [ -z "$GPACRESTARTCOMPARISON_EXEC" ];then
#    echo "you must define the path to the silodiffwrapper.sh script or commensurate comparison executable"
#    exit 1
#fi

$ATS_EXEC --level $ATS_LEVEL testsuite.ats 2>&1 | tee t_out

