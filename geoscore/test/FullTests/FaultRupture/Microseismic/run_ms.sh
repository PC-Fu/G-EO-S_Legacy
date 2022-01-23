#!/bin/bash
#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GX" ]; then
    export GX=~/apps/geos_git/gpaca/GPAC.x.opt
fi

perl create_toc.pl 0.1 0.5 0.03 0.1
$GX -i ms.xml | tee cout
grep ^Tape cout | awk '{print $2" "$3}' > dfiles/beachpts_betadist_points.dat
grep ^AseismicTape cout | awk '{print $2" "$3}' > dfiles/beachpts_betaaseis_points.dat
perl ~/apps/geos_git/gpaca/scripts/gammadelta2lune.pl betadist betaaseis

