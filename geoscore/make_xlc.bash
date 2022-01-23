#!/bin/bash
#export PATH=/usr/apps/gnu/4.7.1/bin:$PATH
make -j12 COMPILER=xlc $1 $2 EXTERNAL_LIBS=/usr/gapps/GEOS/external_libs_VULCAN
