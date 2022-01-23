#!/bin/bash
export N=16
export MYCOMPILER=icc
if [[ $OSTYPE == darwin* ]]; then
    export MYCOMPILER=mac
    export N=4
elif [[ $OSTYPE == linux* ]]; then
if [ -d /usr/gapps ]; then
. /usr/local/tools/dotkit/init.sh
use ic-13.1.163
use mvapich2-intel-1.7
#  - or - 
use gcc-4.6.1
#mvapich2-gnu-1.7
else
    export MYCOMPILER=ubuntu
    export N=4
fi
fi
echo "COMPILER: $MYCOMPILER"

if [ "1" -eq "1" ]; then
make clean
make -j$N COMPILER=$MYCOMPILER DEBUG=0 $1 2>&1 | tee c_opt
mv GEOS.x GPAC.opt
fi

if [ "1" -eq "1" ]; then
make clean
make -j$N COMPILER=$MYCOMPILER DEBUG=1 $1 2>&1 | tee c_debug
mv GEOS.x GPAC.dbg
fi
