#!/bin/bash
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpacspace/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

#Base class
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl MaterialBase ConstitutiveBase Material p_r0 init_density Density p_r0 E E p_r0 Nu Nu p_r0 init_shearModulus ShearModulus p_r0 Lame Lame s_r0 StressPower StressPower s_r0 DissipatedEnergy DissipatedEnergy s_r0 ElasticStrainEnergy ElasticStrainEnergy s_r0 pressure pressure s_r2s devStress devStress s_r0 BulkModulus BulkModulusCurrent s_r0 ShearModulus ShearModulusCurrent b_public ScriptInc/MaterialBasePublic.cpp x s_public ScriptInc/MaterialBaseStateDataPublic.cpp x

#Materials for FE
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl LinearElasticIntermediate MaterialBase Material p_r0 init_bulkModulus BulkModulus s_r0 density density s_r0 ElasticBulkModulus ElasticBulkModulus s_r0 ElasticShearModulus ElasticShearModulus

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl LinearElastic LinearElasticIntermediate Material p_protected ScriptInc/LinearElasticParameterDataProtected.cpp x b_public ScriptInc/LinearElasticPublic.cpp x b_protected ScriptInc/LinearElasticProtected.cpp x

#Materials for DE
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl LinearElasticDEM MaterialBase Material p_r0 yieldStrength yieldStrength p_r0 cor coefficientOfRestitution p_r0 velHalf velocityOfCoROfHalf p_r0 cement cementBondStrength p_r0 surfaceEnergy surfaceEnergy

#CLEAN UP
echo "You will need to now move declaration of MeanPressureDevStressFromDerived from class file to header before attempting to compile!!"
#head -n 120 MaterialBase.cpp