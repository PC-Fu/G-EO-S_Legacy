#!/bin/bash

#ARGS (1) STATPARAMFILE (2) n0 (3) n1 (4) hfct (5) x-length (6) y-length (7) z-length (8) block size (9) MATED
if [ $# -ne 9 ]; then
   echo "usage: ARGS (1) STATPARAMFILE (2) n0 (3) n1 (4) hfct (5) x-length (6) y-length (7) z-length (8) block size (9) MATED"
   exit 1
fi 

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    if [ -d "$GPAC_SCRIPT_PATH" ];then
	echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
    else
	export GPAC_SCRIPT_PATH=/usr/gapps/GEOS/scripts
	if [ -d "$GPAC_SCRIPT_PATH" ];then
	    echo "$0 using GPAC_SCRIPTH_PATH $GPAC_SCRIPT_PATH"
	else
	    echo "$0 cannot find GPAC_SCRIPT_PATH $GPAC_SCRIPT_PATH"
	    exit 1
	fi
    fi
fi

if [ -z "$CUBIT_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export CUBIT_EXEC=/usr/gapps/cubit/linux64.13.1/cubit
    echo "setting default for CUBIT_EXEC: $CUBIT_EXEC"
fi

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

#flow surface
echo "cubit.cmd('Nodeset 3 surface 7')" >> cubit1.py

#xmin, xmax, ymin, ymax, respectively (on surface 7)
echo "cubit.cmd('Nodeset 4 curve 15')" >> cubit1.py
echo "cubit.cmd('Nodeset 5 curve 13')" >> cubit1.py
echo "cubit.cmd('Nodeset 6 curve 16')" >> cubit1.py
echo "cubit.cmd('Nodeset 7 curve 14')" >> cubit1.py

echo "cubit.cmd('block 2 volume 2')" >> cubit1.py
echo "cubit.cmd('export Abaqus \"tmp.inp\" block 1 2 overwrite')" >> cubit1.py
echo "cubit.cmd('quit')" >> cubit1.py

$CUBIT_EXEC -nographics cubit1.py

#get points
perl $GPAC_SCRIPT_PATH/abaqus2points.pl tmp.inp > pts
perl $GPAC_SCRIPT_PATH/points2minmax.pl pts > mm

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
$GPAC_SCRIPT_PATH/apgen file1 file2 > aps

#alter Abaqus node positions
perl $GPAC_SCRIPT_PATH/abaqus_scale_nodes_z.pl `grep ^2 mm | awk '{print $2}'` `grep ^2 mm | awk '{print $3}'` aps $9 tmp.inp > mesh.inp
