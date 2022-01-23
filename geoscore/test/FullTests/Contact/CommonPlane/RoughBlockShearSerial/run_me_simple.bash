#!/usr/bin/bash
echo "0 0" > statparams
bash gpac_surface.sh statparams 1 1 1.2 5 5 1 1 1
cp mesh.inp mesh5x5.geom
../../../../src/GPAC.x -i mesh.xml -m mesh.geom > prunout
exit
