#!/usr/bin/bash

#ARGS (1) STATPARAMFILE (2) n0 (3) n1 (4) hfct (5) x-length (6) y-length (7) z-length (8) block size (9) MATED

#MAKE MESHES MESH

echo "#!/usr/bin/python" > cubit1.py
echo "cubit.cmd('create brick X $5 Y $6 Z $7')" >> cubit1.py
echo "cubit.cmd('volume 1 move 0 0 $7')" >> cubit1.py
echo "cubit.cmd('volume 1 scheme submap')" >> cubit1.py
echo "cubit.cmd('volume 1 size $8')" >> cubit1.py
echo "cubit.cmd('mesh volume 1')" >> cubit1.py
echo "cubit.cmd('Nodeset 1 surface 1')" >> cubit1.py
echo "cubit.cmd('block 1 volume 1')" >> cubit1.py

echo "cubit.cmd('create brick X $5 Y $6 Z $7')" >> cubit1.py
echo "cubit.cmd('volume 2 scheme submap')" >> cubit1.py
echo "cubit.cmd('volume 2 size $8')" >> cubit1.py
echo "cubit.cmd('mesh volume 2')" >> cubit1.py
echo "cubit.cmd('Nodeset 2 surface 8')" >> cubit1.py
echo "cubit.cmd('block 2 volume 2')" >> cubit1.py
echo "cubit.cmd('export Abaqus \"tmp.inp\" block 1 2 overwrite')" >> cubit1.py
echo "cubit.cmd('quit')" >> cubit1.py

/usr/gapps/cubit/linux64.13.0/cubit -nographics cubit1.py

#get points
perl abaqus2points.pl tmp.inp > pts
perl points2minmax.pl pts > mm

#create file1 input for apgen
cat mm | awk '/^0/{print $2}' > file1
cat mm | awk '/^1/{print $2}' >> file1
cat mm | awk '/^0/{print $3}' >> file1
cat mm | awk '/^1/{print $3}' >> file1
echo "$2 $3 $4" >> file1
cat $1 | wc | awk '{print $1}' >> file1
cat $1 | awk '{print $1" "$2}' >> file1

#create file2 input for apgen
cat pts | wc | awk '{print $1}' > file2
cat pts | awk '{print $1" "$2}' >> file2

#make aperture distribution
pushd ../../../../src/SurfaceGeneration/
g++ -O3 -I../ Aperturator.cpp main.cpp -o apgen
popd
../../../../src/SurfaceGeneration/apgen file1 file2 > aps

#alter Abaqus node positions
perl abaqus_scale_nodes_z.pl `grep ^2 mm | awk '{print $2}'` `grep ^2 mm | awk '{print $3}'` aps $9 tmp.inp > mesh.inp
