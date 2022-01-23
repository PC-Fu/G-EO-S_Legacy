#!/bin/bash
RUN_ME=1

########################################
#ARGS (1) executable path (2) lscratch sub-directory location for working - ignored if LSCRATCH_WD set (3) GPAC trunk path - ignored if GPAC_TRUNK set
########################################
if [ $# -lt 3 ]; then
   echo "usage: ARGS (1) executable path (2) lscratch sub-directory location for working - ignored if LSCRATCH_WD set (3) GPAC trunk path - ignored if GPAC_TRUNK set" >&2
   exit 1
fi
if [ -z "$GPAC_TRUNK" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_TRUNK=$3
    if [ -e "$GPAC_TRUNK" ]; then
	echo "setting default for GPAC_TRUNK: $GPAC_TRUNK" >&2
    else
	echo "cannot find default location for GPAC_TRUNK: $GPAC_TRUNK" >&2
	exit 1;
    fi
fi
if [ -z "$SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export SCRIPT_PATH=$GPAC_TRUNK/scripts
    if [ -e "$SCRIPT_PATH" ]; then
	echo "setting default for SCRIPT_PATH: $SCRIPT_PATH" >&2
    else
	export SCRIPT_PATH=/usr/gapps/GEOS/scripts
	if [ -e "$SCRIPT_PATH" ]; then
	    echo "setting default for SCRIPT_PATH: $SCRIPT_PATH" >&2
	else
	    echo "cannot find default location for SCRIPT_PATH: $SCRIPT_PATH" >&2
	    exit 1;
	fi
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

export XMLVALIDATOR_EXEC=$GPAC_TRUNK/external_libs/bin/xmllint
if [ -e "$XMLVALIDATOR_EXEC" ]; then
    echo "----> XML validator executable: $XMLVALIDATOR_EXEC" >&2
else
    export XMLVALIDATOR_EXEC=/usr/gapps/GEOS/external_libs_icc_dbg/bin/xmllint
    if [ -e "$XMLVALIDATOR_EXEC" ]; then
	echo "----> XML validator executable: $XMLVALIDATOR_EXEC" >&2
    else
	echo "CANNOT FIND XMLVALIDATOR_EXEC: $XMLVALIDATOR_EXEC" >&2
	exit 1
    fi
fi

export GPACRESTARTCOMPARISON_EXEC=$GPAC_TRUNK/scripts/silodiffdriver.sh
if [ -e "$GPACRESTARTCOMPARISON_EXEC" ]; then
    echo "----> GPAC restart file comparison utility: $GPACRESTARTCOMPARISON_EXEC"
else
    export GPACRESTARTCOMPARISON_EXEC=/usr/gapps/GEOS/scripts/silodiffdriver.sh
    if [ -e "$GPACRESTARTCOMPARISON_EXEC" ]; then
	echo "----> GPAC restart file comparison utility: $GPACRESTARTCOMPARISON_EXEC"
    else
	echo "CANNOT FIND GPACRESTARTCOMPARISON_EXEC: $GPACRESTARTCOMPARISON_EXEC"
	exit 1
    fi
fi

if [ -e "$XSD_PATH" ]; then
    echo "----> XML schema path: $XSD_PATH"
else
    export XSD_PATH=/usr/gapps/GEOS/schema/gpac.xsd
    if [ -e "$XSD_PATH" ]; then
	echo "----> DEFAULT: XML schema path: $XSD_PATH"
    else
	echo "CANNOT FIND XSD_PATH: $XSD_PATH"
	exit 1
    fi
fi

########################################
#(2) SETUP YOUR WORKING DIRECTORY
########################################
#get the location at which to run (working directory)
#create if necesary
if [ -z "$LSCRATCH_WD" ]; then
    #if the working directory is not explicit, use the best throughput file system on the CZ
    export LSCRATCH_WD=`perl $SCRIPT_PATH/lc_filesystemsummary.pl | tail -n 1 | awk '{print $2}'`"/$USER/$2"
    echo "setting default for LSCRATCH_WD: $LSCRATCH_WD" >&2
fi
echo "----> working directory: $LSCRATCH_WD" >&2

#copy files to the working directory
if [ "$RUN_ME" -eq "1" ]; then
if [ -d $LCRATCH_WD ]; then
    rm -rf $LSCRATCH_WD
fi
mkdir -p $LSCRATCH_WD
cd $LSCRATCH_WD
git clone https://$USER@lc.llnl.gov/stash/scm/gpac/gpaca.git
git clone https://$USER@lc.llnl.gov/stash/scm/gpac/gpacc.git
cd gpaca/test
else
cd $LSCRATCH_WD
cd gpaca/test
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
echo "setenv XSD_PATH \"$XSD_PATH\"" >> ats.msub
echo "setenv XML_SCHEMA_PATH \"$XSD_PATH\"" >> ats.msub
echo "date" >> ats.msub
echo "cd $LSCRATCH_WD/gpaca/test" >> ats.msub
echo "" >> ats.msub
echo "$ATS_EXEC --allInteractive --numNodes=1 --level $ATS_LEVEL testsuite.ats" >> ats.msub
echo "while (1)" >> ats.msub
echo "    sleep 100" >> ats.msub
echo "end" >> ats.msub
echo "date" >> ats.msub

if [ "$RUN_ME" -eq "1" ]; then
msub ats.msub > jnum
fi
