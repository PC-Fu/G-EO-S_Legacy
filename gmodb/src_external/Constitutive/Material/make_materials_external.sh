#!/bin/bash
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpacspace/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl DamageSphericalPores LinearElasticIntermediate Material p_r0 f1 crackInitiationYield p_r0 f2 crackPropgationYield p_r0 initialVoidRatio initialVoidRatio s_r0 radiusRatio radiusRatio s_r0 fractureNumberDensity fractureNumberDensity s_r0 yAffinity yAffinity s_r0 gAffinity gAffinity s_r0 chordLengthMean chordLengthVariance s_r0 m01 M01 s_r0 m02 M02 s_r0 m11 M11 s_r0 m12 M12 s_r0 m21 M21 s_r0 m22 M22

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl DSIDDamage MaterialBase Material p_r0 alpha alpha p_r0 a1 parameter1 p_r0 a2 parameter2 p_r0 a3 parameter3 p_r0 a4 parameter4 p_r0 c0 initialDamageThreshold p_r0 c1 damageHardeningVariable p_r2s omega0 initialDamage s_r0s fd0 DamageIndicator s_r2s omega currentDamage s_r2s epsilon strain s_r2s epsid IrreversibleStrain s_r2s epsE TotalElasticStrain s_r2s epsel pureElasticStrain b_public ScriptInc/DSIDDamagePublic.cpp x

