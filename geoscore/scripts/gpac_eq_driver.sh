#!/bin/bash

#ASSIGN DEFAULT ENVIRONMENTAL VARIABLES
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpaccore/scripts
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
if [ -z "$PSUADE_INPUT_FILE" ]; then
    #NOTE: this only provides a value for the current script PSUADE file
    export PSUADE_INPUT_FILE=$PWD/rsqsim_fA.in
    #echo "setting default for PSUADE_INPUT_FILE: $PSUADE_INPUT_FILE"
fi


#########################################################################
#########################################################################
# gpac_eq_driver.sh
#       author: Scott Johnson (johnson346@llnl.gov)
#       date:   10/8/2013
#
# usage: 
#       This script produces the input files necessary
#       to run GPAC in the EQ simulation mode for a 
#       single fault case
#
# required:
#       gpac_vars -> file will need to exist in the CWD
#             if it does not, a default will be created:
#             (see default file to see allowable variable names)
#       *     global variable entries are in the form:
#                  "GV <variable name> <variable value>"
#       *     initial condition entries are in the form:
#                  "IC <variable name> <minimum> <maximum>"
#                  (where the variable is sampled in the
#                  uniform distribution with the minimum
#                  and maximum range defined in the entry)
#
#       sigHhv -> file will need to exist in the CWD
#             if it does not, a default will be created:
#       *     entries are in z ascending order with format:
#                  "<z> <max horizontal stress> <min horizontal stress> <vertical stress>"
#       *     maximum horizontal stress in the direction defined by
#             the global variable entry SIGHSTRIKE in "gpac_vars", and
#             the direction is constant with depth
#       *     minimum and maximum horizontal stresses are linearly 
#             interpolated between entries
#
#       fault.geom -> file in ABAQUS ASCII format to describe fault
#             file will need to exist in CWD
#             if it does not, a default will be created using CUBIT:
#       *     entries are in z ascending order with format:
#                  "<z> <max horizontal stress> <min horizontal stress>"
#
#########################################################################
#########################################################################



#####################################################
# CREATE gpac_vars IF IT DOES NOT EXIST
#####################################################

if [ -a gpac_vars ];then

echo "<-- USING PRE-EXISTING gpac_vars"

else

echo "GV G 3e10" > gpac_vars
echo "GV K 5e10" >> gpac_vars
export SHEARRATE=3.1709792e-09
echo "GV SHEARRATE $SHEARRATE" >> gpac_vars

echo "GV SIGHSTRIKE 135" >> gpac_vars
echo "GV TBURN 200" >> gpac_vars

echo "GV FAULTSTRIKE 47" >> gpac_vars
echo "GV FAULTLENGTH 1690" >> gpac_vars
echo "GV FAULTDIP 90" >> gpac_vars
echo "GV FAULTHEIGHT 1515" >> gpac_vars
echo "GV NDIP 70" >> gpac_vars
echo "GV NSTRIKE 80" >> gpac_vars

echo "GV XMIN 0" >> gpac_vars
echo "GV YMIN 0" >> gpac_vars
echo "GV ZMIN -2490" >> gpac_vars

echo "GV XMAX 4000" >> gpac_vars
echo "GV YMAX 4000" >> gpac_vars
echo "GV ZMAX -800" >> gpac_vars

echo "GV XORG 787" >> gpac_vars
echo "GV YORG 360" >> gpac_vars
echo "GV ZTOP -885" >> gpac_vars

echo "GV XINJ 187" >> gpac_vars
echo "GV YINJ 60" >> gpac_vars

echo "GV NVOXX 20" >> gpac_vars
echo "GV NVOXY 20" >> gpac_vars
echo "GV NVOXZ 20" >> gpac_vars

echo "IC Dc 0.000015 0.000035" >> gpac_vars
echo "IC frictionCoefficient 0.6 0.9" >> gpac_vars
echo "IC A 0.005 0.008" >> gpac_vars
echo "IC B 0.015 0.015" >> gpac_vars

echo "IC currentFrictionCoefficient 0.6 0.9" >> gpac_vars
echo "IC alpha 0.25 0.25" >> gpac_vars
echo "IC shearRateStar 1e-6 1e-6" >> gpac_vars
echo "IC shearSlipRateAB 1e-6 1e-6" >> gpac_vars
echo "IC stateStar 0.1 0.1" >> gpac_vars
echo "IC state 1e8 1e8" >> gpac_vars
echo "IC shearRateDrive $SHEARRATE $SHEARRATE" >> gpac_vars
echo "IC stressRate -1e-18 -1e-18" >> gpac_vars

echo "--> CREATED gpac_vars"

fi

grep ^GV gpac_vars | awk '{print $2" "$3}' > gpac_vars_gv

#Set the shear modulus of the reservoir (Pa)
export G=`grep ^G gpac_vars_gv | awk '{print $2}'`
#Set the bulk modulus of the reservoir (Pa)
export K=`grep ^K gpac_vars_gv | awk '{print $2}'`
#Set the shear slip driving rate for the fault (m/s)
export SHEARRATE=`grep ^SHEARRATE gpac_vars_gv | awk '{print $2}'`
#Set the strike of the direction of maximum in situ horizontal stress (degrees)
export SIGHSTRIKE=`grep ^SIGHSTRIKE gpac_vars_gv | awk '{print $2}'`
#Set the time to burn-in (years)
export TBURN=`grep ^TBURN gpac_vars_gv | awk '{print $2}'`

#Set the strike of the direction of the fault (degrees)
export FAULTSTRIKE=`grep ^FAULTSTRIKE gpac_vars_gv | awk '{print $2}'`
#Set the length of the fault along the strike direction (m)
export FAULTLENGTH=`grep ^FAULTLENGTH gpac_vars_gv | awk '{print $2}'`
#Set the dip of the fault (degrees)
export FAULTDIP=`grep ^FAULTDIP gpac_vars_gv | awk '{print $2}'`
#Set the height of the fault along the vertical (z) direction (m)
export FAULTHEIGHT=`grep ^FAULTHEIGHT gpac_vars_gv | awk '{print $2}'`
#Set the number of elements in the dip direction
export NDIP=`grep ^NDIP gpac_vars_gv | awk '{print $2}'`
#Set the number of elements in the strike direction
export NSTRIKE=`grep ^NSTRIKE gpac_vars_gv | awk '{print $2}'`

#Set the minimum extent of the reservoir domain (X)
export XMIN=`grep ^XMIN gpac_vars_gv | awk '{print $2}'`
#Set the minimum extent of the reservoir domain (Y)
export YMIN=`grep ^YMIN gpac_vars_gv | awk '{print $2}'`
#Set the minimum extent of the reservoir domain (Z)
export ZMIN=`grep ^ZMIN gpac_vars_gv | awk '{print $2}'`

#Set the maximum extent of the reservoir domain (X)
export XMAX=`grep ^XMAX gpac_vars_gv | awk '{print $2}'`
#Set the maximum extent of the reservoir domain (Y)
export YMAX=`grep ^YMAX gpac_vars_gv | awk '{print $2}'`
#Set the maximum extent of the reservoir domain (Z)
export ZMAX=`grep ^ZMAX gpac_vars_gv | awk '{print $2}'`

#Set the x-origin (minimum) of the fault (m)
export XORG=`grep ^XORG gpac_vars_gv | awk '{print $2}'`
#Set the y-origin (minimum) of the fault (m)
export YORG=`grep ^YORG gpac_vars_gv | awk '{print $2}'`
#Set the z-top (maximum) of the fault (m)
export ZTOP=`grep ^ZTOP gpac_vars_gv | awk '{print $2}'`

#Set the x-coordinate injector (m)
export XINJ=`grep ^XINJ gpac_vars_gv | awk '{print $2}'`
#Set the y-coordinate injector (m)
export YINJ=`grep ^YINJ gpac_vars_gv | awk '{print $2}'`


#Set number of voxels along each dimension for the randomly generated variables
export NVOXX=`grep ^NVOXX gpac_vars_gv | awk '{print $2}'`
export NVOXY=`grep ^NVOXY gpac_vars_gv | awk '{print $2}'`
export NVOXZ=`grep ^NVOXZ gpac_vars_gv | awk '{print $2}'`

#####################################################
# CREATE FAULT GEOMETRY
#####################################################

if [ -a fault.geom ];then

echo "<-- USING PRE-EXISTING fault.geom"

else

bash $GPAC_SCRIPT_PATH/fault2abaqus.sh $FAULTSTRIKE $FAULTDIP $FAULTLENGTH $FAULTHEIGHT $NSTRIKE $NDIP $XORG $YORG $ZTOP

echo "--> CREATED fault.geom"

fi

#####################################################
# CREATE STRESS FILES
#####################################################

if [ -a sigHhv ];then

echo "<-- USING PRE-EXISTING sigHhv"

else

echo "-2500 63088600 38857000 56329107" > sigHhv
echo "-2000 50473600 31067000 45065714" >> sigHhv
echo "-1620 40886200 25146600 36505536" >> sigHhv
echo "-1619 39762640 21338420 36479486" >> sigHhv
echo "-1000 24560000 13180000 22532110" >> sigHhv
echo "    0        0        0        0" >> sigHhv

echo "--> CREATED sigHhv"

fi

perl $GPAC_SCRIPT_PATH/stress2voxel.pl sigHhv $SIGHSTRIKE $FAULTSTRIKE $FAULTDIP
mv x xs
mv y ys
mv z zs
echo "-->CREATED STRESS GRADIENT"

export RAKE=`perl -e 'use Math::Trig; my $strike = shift; $strike *= atan(1) / 45.0; my @x = (sin($strike), cos($strike), 0); print "$x[0] $x[1] $x[2]";' -f $FAULTSTRIKE`

#####################################################
# CREATE RANDOM VOXEL FILES
#####################################################

if [ -a ic.xml ];then
    rm ic.xml
fi
if [ -a tables.xml ];then
    rm tables.xml
fi
grep ^IC gpac_vars | awk '{print $2" "$3" "$4}' > gpac_vars_ic
while read line; do
    perl $GPAC_SCRIPT_PATH/random2voxel.pl $XMAX $YMAX $ZMAX $XMIN $YMIN $ZMIN $NVOXX $NVOXY $NVOXZ `echo $line | awk '{print $2" "$3}'`
    export VAR=`echo $line | awk '{print $1}'`
    mv voxel $VAR
    echo "    <InitialConstitutiveValue object=\"Fault_ElementManager\" propertytype=\"$VAR\" tablename=\"$VAR\" toregion=\"EB2\"/>" >> ic.xml
    echo "    <Table3D name=\"$VAR\" x_file=\"x\" y_file=\"y\" z_file=\"z\" voxel_file=\"$VAR\"/>" >> tables.xml
done <gpac_vars_ic

cp frictionCoefficient currentFrictionCoefficient
echo "-->CREATED INITIAL STATE FILES"

#####################################################
# CREATE PORE PRESSURE FILE
#####################################################

#need xpp, ypp, zpp, pp
cp x xpp
cp y ypp
cp z zpp

perl $GPAC_SCRIPT_PATH/injection2voxel.pl xpp ypp zpp $TBURN $XINJ $YINJ `perl -e 'my $l = shift; my $nz = shift; $l /= $nz; print $l;' -f $FAULTHEIGHT $NDIP`
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

