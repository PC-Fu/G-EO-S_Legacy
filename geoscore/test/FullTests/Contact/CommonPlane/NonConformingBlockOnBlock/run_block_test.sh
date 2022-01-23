#!/usr/bin/bash

if [ -e "impact.xml" ];then
    echo "found impact.xml"
else
    echo "CANNOT FIND IMPACT.XML ... EXITING"
    exit 1
fi

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi
if [ -z "$GPAC_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_EXEC=~/apps/gpac_trunk/src/GPAC.x
    echo "setting default for GPAC_EXEC: $GPAC_EXEC"
fi
if [ -z "$VISIT_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export VISIT_EXEC=/Applications/VisIt.app/Contents/MacOS/VisIt
    echo "setting default for VISIT_EXEC: $VISIT_EXEC"
fi

# -v 1.5.5 -nowin -cli -l srun -p pbatch -np 80 -s script.py -verbose -geometry 512x512+0+0 -switch ib -x prism34

#-----------------------
# (1) OBLIQUE IMPACT
rm output_*
cat impact.xml | sed s/AAAAA/\ \ \ \ \ \ \ \ \ \ \ BulkModulus\=\"43.5\"/ > tmp
cat tmp | sed s/BBBBB/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ value\=\"0.0001\ 0.0001\ 0.0\"\>/ > impactoblique.xml
$GPAC_EXEC -i impactoblique.xml -m mesh_2x2_3x3.inp > prunout
grep ^Stress prunout > eobl0
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tobl0
paste -d ' ' tobl0 eobl0 > teobl0

rm output_*
cat impact.xml | sed s/AAAAA/\ \ \ \ \ \ \ \ \ \ \ BulkModulus\=\"43.5\"/ > tmp
cat tmp | sed s/BBBBB/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ value\=\"0.0002\ 0.0002\ 0.0\"\>/ > impactoblique.xml
$GPAC_EXEC -i impactoblique.xml -m mesh_2x2_3x3.inp > prunout
grep ^Stress prunout > eobl1
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tobl1
paste -d ' ' tobl1 eobl1 > teobl1

#add analysis logic here

#-----------------------
# (2) NORMAL IMPACT (RESOLUTION)
cat impact.xml | sed s/AAAAA/\ \ \ \ \ \ \ \ \ \ \ BulkModulus\=\"43.5\"/ > tmp
cat tmp | sed s/BBBBB/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ value\=\"0.0001\ 0.0\ 0.0\"\>/ > impactnormal.xml

rm output_*
$GPAC_EXEC -i impactnormal.xml -m mesh_2x2_3x3.inp > prunout
grep ^Stress prunout > eres0
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tres0
paste -d ' ' tres0 eres0 > teres0
#add analysis logic here


rm output_*
$GPAC_EXEC -i impactnormal.xml -m mesh_4x4_5x5.inp > prunout
grep ^Stress prunout > eres1
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tres1
paste -d ' ' tres1 eres1 > teres1
#add analysis logic here

rm output_*
$GPAC_EXEC -i impactnormal.xml -m mesh_6x6_7x7.inp > prunout
grep ^Stress prunout > eres2
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tres2
paste -d ' ' tres2 eres2 > teres2
#add analysis logic here

rm output_*
$GPAC_EXEC -i impactnormal.xml -m mesh_8x8_10x10.inp > prunout
grep ^Stress prunout > eres3
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tres3
paste -d ' ' tres3 eres3 > teres3
#add analysis logic here

rm output_*
$GPAC_EXEC -i impactnormal.xml -m mesh_18x18_20x20.inp > prunout
grep ^Stress prunout > eres4
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tres4
paste -d ' ' tres4 eres4 > teres4
#add analysis logic here



#-----------------------
# (3) CONTACT FUNCTIONAL
rm tangent_*
$GPAC_EXEC -i movenormalthentangential.xml -m mesh_6x6_7x7.inp > prunout
perl $GPAC_SCRIPT_PATH/gpac_query_general.pl tangential.session "Variable Sum" > tmp.py
$VISIT_EXEC -nowin -cli -s tmp.py > tpath0



#-----------------------
# (4) NORMAL IMPACT (STIFFNESS)
cat impact.xml | sed s/AAAAA/\ \ \ \ \ \ \ \ \ \ \ BulkModulus\=\"43.5e2\"/ > tmp
cat tmp | sed s/BBBBB/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ value\=\"0.01\ 0.0\ 0.0\"\>/ > impactnormal.xml
$GPAC_EXEC -i impactnormal.xml -m mesh_6x6_7x7.inp > prunout
grep ^Stress prunout > estiff0
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tstiff0
paste -d ' ' tstiff0 estiff0 > testiff0
#add analysis logic here

rm output_*
cat impact.xml | sed s/AAAAA/\ \ \ \ \ \ \ \ \ \ \ BulkModulus\=\"43.5e3\"/ > tmp
cat tmp | sed s/BBBBB/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ value\=\"0.01\ 0.0\ 0.0\"\>/ > impactnormal.xml
$GPAC_EXEC -i impactnormal.xml -m mesh_6x6_7x7.inp > prunout
grep ^Stress prunout > estiff1
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tstiff1
paste -d ' ' tstiff1 estiff1 > testiff1
#add analysis logic here

rm output_*
cat impact.xml | sed s/AAAAA/\ \ \ \ \ \ \ \ \ \ \ BulkModulus\=\"43.5e4\"/ > tmp
cat tmp | sed s/BBBBB/\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ value\=\"0.01\ 0.0\ 0.0\"\>/ > impactnormal.xml
$GPAC_EXEC -i impactnormal.xml -m mesh_6x6_7x7.inp > prunout
grep ^Stress prunout > estiff2
cat prunout | awk '/Time=/{split($3,a,"=");split(a[2],b,",");print b[1]}' > tstiff2
paste -d ' ' tstiff2 estiff2 > testiff2
#add analysis logic here

#would be nice to run rigid wall case ... or show it analytically

exit
