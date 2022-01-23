#!/usr/bin/bash
echo "0.1 0.1" > statparams
bash gpac_surface.sh statparams 5 5 1.2 5 5 1 0.2 1
cp mesh.inp mesh.geom
../../../../src/GPAC.x -i mesh.xml -m mesh.geom > prunout
exit
