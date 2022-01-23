#!/bin/bash

########################################
#ARGS (1) executable path (2) lscratch sub-directory location for working - ignored if LSCRATCH_WD set (3) GPAC trunk path - ignored if GPAC_TRUNK set
########################################
if [ $# -lt 2 ]; then
   echo "usage: ARGS (1) executable path (2) GPAC trunk path - ignored if GPAC_TRUNK set" >&2
   exit 1
fi
if [ -z "$GPAC_TRUNK" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_TRUNK=$2
    if [ -e "$GPAC_TRUNK" ]; then
	echo "setting default for GPAC_TRUNK: $GPAC_TRUNK" >&2
    else
	echo "cannot find default location for GPAC_TRUNK: $GPAC_TRUNK" >&2
	exit 1;
    fi
fi

########################################
#(1) SET NECESSARY ENVIRONMENTAL VARIABLES
########################################

export GPAC_EXEC=$1
if [ -e "$GPAC_EXEC" ]; then
    echo "----> GPAC executable: $GPAC_EXEC" >&2
else
    echo "CANNOT FIND GPAC_EXEC: $GPAC_EXEC" >&2
    exit 1
fi


########################################
#(3) RUN TESTS
########################################
ATS_LEVEL=400
if [ -z "$ATS_EXEC" ];then
    ATS_EXEC=/usr/gapps/ats/chaos_5_x86_64_ib/bin/ats
    if [ -e "$ATS_EXEC" ]; then
	echo "Setting default ATS path to $ATS_EXEC"
    else
	echo "Cannot find default ATS path: $ATS_EXEC"
	exit 1
    fi
fi

echo "#!/bin/csh" > ats.msub
echo "#MSUB -N ats.gpac.job" >> ats.msub
echo "#MSUB -j oe" >> ats.msub
echo "#MSUB -o ats.gpac.job.out" >> ats.msub
echo "#MSUB -q pbatch" >> ats.msub
echo "#MSUB -l nodes=1" >> ats.msub
echo "#MSUB -l partition=cab" >> ats.msub
echo "#MSUB -l walltime=900" >> ats.msub
echo "#MSUB -V" >> ats.msub
echo "# exports all environment var" >> ats.msub
echo "# bank to use" >> ats.msub
echo "#MSUB -A frnet" >> ats.msub
echo "setenv SYS_TYPE chaos_5_x86_64_ib" >> ats.msub
echo "setenv MPI_ARGS \" \"" >> ats.msub
echo "setenv GPAC_TRUNK \"$GPAC_TRUNK\"" >> ats.msub
echo "setenv GPAC_EXEC \"$GPAC_EXEC\"" >> ats.msub
echo "setenv XSD_PATH \"$GPAC_TRUNK/src/schema/gpac.xsd\"" >> ats.msub
echo "date" >> ats.msub
echo "cd $PWD" >> ats.msub
echo "" >> ats.msub
echo "$ATS_EXEC --allInteractive --numNodes=1 --level $ATS_LEVEL testsuite.ats" >> ats.msub
echo "while (1)" >> ats.msub
echo "    sleep 100" >> ats.msub
echo "end" >> ats.msub
echo "date" >> ats.msub
