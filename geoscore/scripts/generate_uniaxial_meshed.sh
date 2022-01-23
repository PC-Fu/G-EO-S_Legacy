#!/usr/bin/bash

if [ $# -ne 6 ]; then
   echo "usage: ARGS (1) model generation input file (2) unit sphere mesh file (3) wall thickness (4) strain rate (5) combined mesh file name (6) input file name"
   exit 1
fi 

if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi
if [ -z "$DEM_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export DEM_EXEC=~/apps/dem3d/dem.o3
    echo "setting default for DEM_EXEC: $DEM_EXEC"
fi
if [ -z "$CUBIT_EXEC" ]; then
    #NOTE: this only provides a value for the current script
    export CUBIT_EXEC=/usr/gapps/cubit/linux64.13.1/cubit
    echo "setting default for CUBIT_EXEC: $CUBIT_EXEC"
fi

#pack spheres
bash $GPAC_SCRIPT_PATH/generate_uniaxial.sh $1 spheres

#mesh packed spheres
bash $GPAC_SCRIPT_PATH/spheres2abaqus.sh spheres $2 grains.inp

#create boundaries
bash $GPAC_SCRIPT_PATH/generate_uniaxial_bounds.sh grains.inp $3 bounds.inp

#combine meshes
perl $GPAC_SCRIPT_PATH/abaqus_combine.pl grains.inp bounds.inp > $5

#create input file with total strain of 0.3
export TEND=`perl -e 'my $srate = shift; $srate *= 3.0; my $tend = 1.0/$srate; print $tend;' -f $4`
export T0=`perl -e 'my $srate = shift; $srate *= 3.0; my $tend = 1.0/$srate; $tend *= 0.1; print $tend;' -f $4`
export T1=`perl -e 'my $srate = shift; $srate *= 3.0; my $tend = 1.0/$srate; $tend *= 0.9; print $tend;' -f $4`
export DTVIS=`perl -e 'my $dt = $ENV{TEND}; shift; $dt /= 200; print $dt;'`
export ZUPPER=`cat $1 | grep ^upper | awk '{print $4}'`
export ZLOWER=`cat $1 | grep ^lower | awk '{print $4}'`
export VELBOUND=`perl -e 'my $zlower = $ENV{ZLOWER}; my $zupper = $ENV{ZUPPER}; my $srate = shift; my $length = $zupper - $zlower; my $vel = 0.5 * $srate * $length; print $vel;' -f $4`

echo '<?xml version="1.0" ?>' > $6
echo '<Problem>' >> $6
echo '  <Units length="mm" mass="mg" time="us" temperature="K" mole="mol"/>' >> $6
echo '  <Mesh file="combined.geom" units="mm"/>' >> $6
echo '  <Solvers>' >> $6
echo '    <LagrangeExplicitDynamicsSolver name="solver1" courant="0.5"' >> $6
echo '				    tiedNodesFlag="0" tiedNodeNormalRuptureStress="0.1" ' >> $6
echo '				    tiedNodeTolerance="1.0e-3"/>' >> $6
echo '  </Solvers>' >> $6
echo ' ' >> $6
echo '  <SolverApplications>' >> $6
echo "    <Application name=\"1\" begintime=\"0.0\" endtime=\"$TEND\">" >> $6
echo '      <Apply solver="solver1" toregions="PM1 EB5 EB6 EB1 EB2 EB3 EB4"/>' >> $6
echo '    </Application>' >> $6  
echo '  </SolverApplications>' >> $6
echo ' ' >> $6
echo '  <ElementRegions>' >> $6
echo '    <ElementRegion name="PM1" elementtype="poly" BulkModulus="15" ShearModulus="15" Density="2.65" hgStiff="0.01" hgDamp="0.1"/>' >> $6
echo '    <ElementRegion name="EB1" elementtype="uniformstrain" BulkModulus="15" ShearModulus="15" Density="2.65" hgStiff="0.01" hgDamp="0.1"/>' >> $6
echo '    <ElementRegion name="EB2" elementtype="uniformstrain" BulkModulus="15" ShearModulus="15" Density="2.65" hgStiff="0.01" hgDamp="0.1"/>' >> $6
echo '    <ElementRegion name="EB3" elementtype="uniformstrain" BulkModulus="15" ShearModulus="15" Density="2.65" hgStiff="0.01" hgDamp="0.1"/>' >> $6
echo '    <ElementRegion name="EB4" elementtype="uniformstrain" BulkModulus="15" ShearModulus="15" Density="2.65" hgStiff="0.01" hgDamp="0.1"/>' >> $6
echo '    <ElementRegion name="EB5" elementtype="uniformstrain" BulkModulus="15" ShearModulus="15" Density="2.65" hgStiff="0.01" hgDamp="0.1"/>' >> $6
echo '    <ElementRegion name="EB6" elementtype="uniformstrain" BulkModulus="15" ShearModulus="15" Density="2.65" hgStiff="0.01" hgDamp="0.1"/>' >> $6
echo '  </ElementRegions>' >> $6
echo ' ' >> $6
echo ' ' >> $6
echo '  <Contact active="1" ' >> $6
echo '	   aperture="0.005"' >> $6
echo '	   normalApproachYield="0.0049999"' >> $6
echo '	   stressYield="1.5"' >> $6
echo '	   stressSoften="2.0"' >> $6
echo '	   arealStiffnessShear="1.5"' >> $6
echo '	   frictionCoefficient="0.0"' >> $6
echo '	   penetrationTol="0.25"' >> $6
echo '	   cosMinTol="0.5"' >> $6
echo '	   feParentSolnTol="1e-6"' >> $6
echo '	   searchRadiusFactor="0.05"' >> $6
echo '	   searchRadiusVelocityFactor="70.0"' >> $6
echo '	   />' >> $6
echo ' ' >> $6
echo '  <BoundaryConditions>' >> $6
echo ' ' >> $6
for i in `seq 1 6`; do
    echo '    <BoundaryCondition     fieldname="Velocity" ' >> $6
    echo "                           setname=\"NS$i\" " >> $6
    echo '			   component="0"' >> $6
    echo '                           scale="0.0"' >> $6
    echo '			   timetable="ttable2" />' >> $6
    echo '    <BoundaryCondition     fieldname="Velocity" ' >> $6
    echo "                           setname=\"NS$i\" " >> $6
    echo '			   component="1"' >> $6
    echo '                           scale="0.0"' >> $6
    echo '			   timetable="ttable2" />' >> $6
    echo '    <BoundaryCondition     fieldname="Velocity" ' >> $6
    echo "                           setname=\"NS$i\" " >> $6
    echo '			   component="2"' >> $6
    if [ "$i" -gt "4" ]; then
	if [ "$i" -eq "5" ]; then
	    echo "                           scale=\"$VELBOUND\"" >> $6
	else
	    echo "                           scale=\"-$VELBOUND\"" >> $6
	fi
	echo '			   timetable="ttable" />' >> $6    
    else
	echo "                           scale=\"0.0\"" >> $6
	echo '			   timetable="ttable2" />' >> $6    
    fi
done
echo '  </BoundaryConditions>' >> $6
echo ' ' >> $6
echo '  <Tables>' >> $6
echo '    <Table1D name="ttable" ' >> $6
echo "               coord  = \"0.0, $T0, $T1, $TEND\"" >> $6
echo '               value =  "0.0, 1.0,  1.0, 0.0"/>' >> $6
echo '    <Table1D name="ttable2" ' >> $6
echo "               coord  = \"0.0, $TEND\"" >> $6
echo '               value =  "1.0, 1.0"/>' >> $6
echo '  </Tables>' >> $6
echo ' ' >> $6
echo '  <Partition>' >> $6
echo '    <SpatialPartition xpar="1" ypar="1" zpar="1" />' >> $6
echo '  </Partition>' >> $6
echo ' ' >> $6
echo "  <Output  plot_interval=\"$DTVIS\" plotfile_root=\"output\" parallel_silo=\"1\" writeFEMFaces=\"1\" />" >> $6
echo ' ' >> $6
echo '</Problem>' >> $6
