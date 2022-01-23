#!/bin/bash
. /usr/local/tools/dotkit/init.sh
use gcc-4.8.2p
use mvapich2-gnu-1.9

make -j2 COMPILER=gcc $1 $2
