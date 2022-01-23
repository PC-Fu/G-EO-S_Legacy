#!/bin/bash

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

export G=3e10
export K=5e10
export SHEARRATE=3.1709792e-09
export SIGHSTRIKE=135
export TBURN=20

export FAULTSTRIKE=47
export FAULTLENGTH=1690
export FAULTDIP=90
export FAULTHEIGHT=1515
export NDIP=18
export NSTRIKE=20

#####################################################
# CREATE FAULT GEOMETRY
#####################################################

if [ "0" -eq "1" ]; then
bash $GPAC_SCRIPT_PATH/fault2abaqus.sh $FAULTSTRIKE $FAULTDIP $FAULTLENGTH $FAULTHEIGHT $NSTRIKE $NDIP 787 360 -885
echo "-->CREATED ABAQUS FAULT"
fi

#####################################################
# CREATE STRESS FILES
#####################################################

echo "-2500 63088600 38857000 56329107" > sigHhv
echo "-2000 50473600 31067000 45065714" >> sigHhv
echo "-1620 40886200 25146600 36505536" >> sigHhv
echo "-1619 39762640 21338420 36479486" >> sigHhv
echo "-1000 24560000 13180000 22532110" >> sigHhv
echo "    0        0        0        0" >> sigHhv

perl $GPAC_SCRIPT_PATH/stress2voxel.pl sigHhv $SIGHSTRIKE $FAULTSTRIKE $FAULTDIP
mv x xs
mv y ys
mv z zs
echo "-->CREATED STRESS GRADIENT"

export RAKE=`perl -e 'use Math::Trig; my $strike = shift; $strike *= atan(1) / 45.0; my @x = (sin($strike), cos($strike), 0); print "$x[0] $x[1] $x[2]";' -f $FAULTSTRIKE`

#####################################################
# CREATE RANDOM VOXEL FILES
#####################################################

if [ "1" -eq "1" ]; then

echo "A 0.005 0.008" > gpac_vars
echo "B 0.015 0.015" >> gpac_vars
echo "Dc 0.000015 0.000035" >> gpac_vars
echo "frictionCoefficient 0.6 0.9" >> gpac_vars

echo "currentFrictionCoefficient 0.6 0.9" >> gpac_vars
echo "alpha 0.25 0.25" >> gpac_vars
echo "shearRateStar 1e-6 1e-6" >> gpac_vars
echo "shearSlipRateAB 1e-6 1e-6" >> gpac_vars
echo "stateStar 0.1 0.1" >> gpac_vars
echo "state 1e8 1e8" >> gpac_vars
echo "shearRateDrive $SHEARRATE $SHEARRATE" >> gpac_vars
echo "stressRate -1e-18 -1e-18" >> gpac_vars
echo "-->CREATED GPAC_VARS"

fi

seq `cat gpac_vars | wc | awk '{print $1}'` > tmp
if [ -a ic.xml ];then
    rm ic.xml
fi
if [ -a tables.xml ];then
    rm tables.xml
fi
while read line; do
    perl $GPAC_SCRIPT_PATH/random2voxel.pl 4000 4000 -800 0 0 -2490 20 20 20 `echo $line | awk '{print $2" "$3}'`
    export VAR=`echo $line | awk '{print $1}'`
    mv voxel $VAR
    echo "    <InitialConstitutiveValue object=\"Fault_ElementManager\" propertytype=\"$VAR\" tablename=\"$VAR\" toregion=\"EB2\"/>" >> ic.xml
    echo "    <Table3D name=\"$VAR\" x_file=\"x\" y_file=\"y\" z_file=\"z\" voxel_file=\"$VAR\"/>" >> tables.xml
done <gpac_vars

cp frictionCoefficient currentFrictionCoefficient
echo "-->CREATED INITIAL STATE FILES"

#####################################################
# CREATE PORE PRESSURE FILE
#####################################################

#need xpp, ypp, zpp, pp
cp x xpp
cp y ypp
cp z zpp

perl $GPAC_SCRIPT_PATH/injection2voxel.pl xpp ypp zpp $TBURN $FAULTLENGTH `perl -e 'my $l = shift; my $nz = shift; $l /= $nz; print $l;' -f $FAULTHEIGHT $NDIP`
mv voxel pp
export TYEARS=`tail -n 1 t | awk '{x=$1 / (365.25 * 3600 * 24); print x}'`
echo "-->CREATED PORE PRESSURE FILE"

#####################################################
# CREATE INPUT FILE
#####################################################

echo "<?xml version=\"1.0\" ?>" > run.xml
echo "<Problem xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" >> run.xml
echo "	xsi:noNamespaceSchemaLocation=\"gpac.xsd\">" >> run.xml
echo " " >> run.xml
echo "  <Units default_units=\"SI_Units\" />" >> run.xml
echo " " >> run.xml
echo "  <Mesh faultpatch_file=\"fault.geom\" units=\"m\"/>" >> run.xml
echo " " >> run.xml
echo "  <Solvers>" >> run.xml
echo "    <FaultRuptureBEMSolver name=\"solver1\"" >> run.xml
echo "			   upVector=\"0 0 1\"" >> run.xml
echo "			   rakeVector=\"$RAKE\"" >> run.xml
echo "			   maximumHorizontalStressDirection=\"1 0 0\"" >> run.xml
echo "			   dHdu=\"0\" dhdu=\"0\"" >> run.xml
echo "			   slipRateEQ=\"1\" creepFrictionRate=\"0.001\"" >> run.xml
echo "			   porePressureTableName=\"pp\">" >> run.xml
echo "      <Material ShearModulus=\"$G\" BulkModulus=\"$K\"/>" >> run.xml
echo "    </FaultRuptureBEMSolver>" >> run.xml
echo "  </Solvers>" >> run.xml
echo " " >> run.xml
echo "  <SolverApplications>" >> run.xml
echo "    <Application name=\"1\" begintime=\"0.0\" endtime=\"$TYEARS year\">" >> run.xml
echo "      <Apply solver=\"solver1\" toregions=\"EB2\"/>" >> run.xml
echo "    </Application>" >> run.xml
echo "  </SolverApplications>" >> run.xml
echo " " >> run.xml
echo "  <BoundaryConditions>" >> run.xml
echo "    <BoundaryCondition  fieldname=\"shearRateDrive\"" >> run.xml
echo "                        scale=\"$SHEARRATE\"" >> run.xml
echo "                        object=\"FaultElement\"/>" >> run.xml
echo "  </BoundaryConditions>" >> run.xml
echo "  " >> run.xml
echo "  <InitialConstitutiveValues>" >> run.xml
cat ic.xml >> run.xml
echo "    <InitialConstitutiveValue object=\"Fault_ElementManager\" propertytype=\"stressReference\" tablename=\"stressReference\" toregion=\"EB2\"/>" >> run.xml
echo "    <InitialConstitutiveValue object=\"Fault_ElementManager\" propertytype=\"shearStressReference\" tablename=\"shearStressReference\" toregion=\"EB2\"/>" >> run.xml
echo "  </InitialConstitutiveValues>" >> run.xml
echo " " >> run.xml
echo "  <Tables>" >> run.xml
cat tables.xml >> run.xml
echo "    <Table3D name=\"stressReference\" x_file=\"xs\" y_file=\"ys\" z_file=\"zs\" voxel_file=\"sigma\"/>" >> run.xml
echo "    <Table3D name=\"shearStressReference\" x_file=\"xs\" y_file=\"ys\" z_file=\"zs\" voxel_file=\"tau\"/>" >> run.xml
echo "    <Table4D name=\"pp\" t_file=\"t\" x_file=\"xpp\" y_file=\"ypp\" z_file=\"zpp\" time_voxel_file=\"pp\"/>" >> run.xml
echo "  </Tables>" >> run.xml
echo " " >> run.xml
echo "  <Partition>" >> run.xml
echo "    <SpatialPartition xpar=\"1\" ypar=\"1\" zpar=\"1\" />" >> run.xml
echo "  </Partition>" >> run.xml
echo " " >> run.xml
echo "  <Output  writePlot=\"0\" plot_interval=\"1 year\" restart_interval=\"$TYEARS year\" plotfile_root=\"ms\" parallel_silo=\"1\"/>" >> run.xml
echo " " >> run.xml
echo "</Problem>" >> run.xml


echo "-->CREATED INPUT FILE"
echo " "
echo "***** DONE *****"
