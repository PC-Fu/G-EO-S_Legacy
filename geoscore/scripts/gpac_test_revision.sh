#!/bin/bash
if [ $# -ne 3 ];then
    echo "usage: ARGS (1) revision number (2) input file path (3) geometry file path"
    exit 1
fi

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_EXT_LIB" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_EXT_LIB=~/apps/gpac_trunk/external_libs
    if [ -d "$GPAC_EXT_LIB" ];then
	echo "$0 using GPAC_EXT_LIB $GPAC_EXT_LIB"
    else
	export GPAC_EXT_LIB=/usr/gapps/GEOS/external_libs_icc_dbg
	if [ -d "$GPAC_EXT_LIB" ];then
	    echo "$0 using GPAC_EXT_LIB $GPAC_EXT_LIB"
	else
	    echo "$0 cannot find XMLVALIDATOR_EXEC $XMLVALIDATOR_EXEC"
	    exit 1
	fi
    fi
fi

mkdir rev$1
cd rev$1
svn co -r $1 https://sourceforge.llnl.gov/svn/repos/gpac/trunk/src
cp ../Make.local.defs ./
ln -s $GPAC_EXT_LIB external_libs
export DYLD_LIBRARY_PATH="$GPAC_EXT_LIB/lib/"
. /usr/local/tools/dotkit/init.sh
use ic-12.1.339
cd src
pushd ../../
cp $2 rev$1/src/
cp $3 rev$1/src/
popd
make -j 16
srun -n 1 -ppdebug GPAC.x -i $2 -m $3 > cout
cd ../../

