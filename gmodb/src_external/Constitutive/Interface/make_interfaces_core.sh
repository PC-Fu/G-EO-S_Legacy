#!/bin/bash
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpacspace/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl InterfaceBase ConstitutiveBase Interface p_r0 ston shearToNormalStiffnessRatio p_r0 mu0 frictionCoefficient p_r0 kappa0 permeabilityInitial p_r0 dkappadsFct permeabilityShearCoefficient p_r0 dkappadsExp permeabilityShearExponent p_r0 dkappadnFct permeabilityNormalCoefficient p_r0 dkappadnExp permeabilityNormalExponent s_r0 ElasticStrainEnergy elasticStrainEnergy s_r0 DissipatedEnergy dissipatedEnergy s_r0 stress stress s_r0 stressShear shearStress s_r0 dxndt normalVelocity s_r0 dxsdt shearVelocity s_r0 xs shear s_r0 normalApproach normalApproach s_r0 dt timestep s_r0 kappa permeability s_r1 stressShearVector shearStressVector b_public ScriptInc/InterfaceBasePublic.cpp x b_protected ScriptInc/InterfaceBaseProtected.cpp x s_public ScriptInc/InterfaceBaseStatePublic.cpp x

#PENALTY DERIVATIVES
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl PenaltyCoulombIntermediate InterfaceBase Interface p_r0 aperture aperture p_r0 normalApproachYield normalApproachYield p_r0 ksoften arealStiffnessSoften p_r0 stressYield stressYield p_r0 stressSoften stressSoften p_r0 kshear arealStiffnessShear p_r0 ktildeAperture ktildeAperture p_r0 kyield arealStiffnessYield p_r0 normalApproachSoften normalApproachSoften b_public ScriptInc/PenaltyCoulombIntermediatePublic.cpp x p_public ScriptInc/PenaltyCoulombIntermediateParameterDataPublic.cpp x p_protected ScriptInc/PenaltyCoulombIntermediateParameterDataProtected.cpp x

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl PenaltyCoulomb PenaltyCoulombIntermediate Interface

#DEM MODELS
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl HertzianIntermediate InterfaceBase Interface s_r0 mu currentFrictionCoefficient s_r0 dissipatedViscousEnergy dissipatedViscousEnergy s_r0 radius radius s_r0 youngs youngsModulus s_r0 mass mass s_r0 hertzCf hertzianCoefficient b_protected ScriptInc/HertzianIntermediateProtected.cpp x s_public ScriptInc/HertzianIntermediateStateDataPublic.cpp x b_public ScriptInc/HertzianIntermediatePublic.cpp x

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl Hertzian HertzianIntermediate Interface

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl Linearized HertzianIntermediate Interface s_r0 a0 a0 s_r0 at at s_r0 ac ac s_int isCohesive cohesionFlag s_r0 dcoh0 cohesionDistance s_r0 dupre dupre s_r0 pullOffForce pullOffForce s_r0 Tk Tk s_int loading loadingFlag s_r0 coefficientOfRestitution2 coefficientOfRestitution2 s_r0 yield yield s_r0 halfRestVel halfRestVel s_r0 ecrit ecrit s_r0 fcrit fcrit s_r0 linearStiffnessElastic linearStiffnessElastic s_r0 vark2 vark2 p_r0 interfaceEnergy interfaceEnergy b_protected ScriptInc/LinearizedProtected.cpp x s_public ScriptInc/LinearizedStateDataPublic.cpp x
