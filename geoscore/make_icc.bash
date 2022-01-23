#!/bin/bash
. /usr/local/tools/dotkit/init.sh
#export PATH=/usr/apps/gnu/4.7.1/bin:$PATH
use gcc-4.6.1
use icc-13.1.163
use mvapich2-intel-1.7

make -j16 COMPILER=icc $1 $2

