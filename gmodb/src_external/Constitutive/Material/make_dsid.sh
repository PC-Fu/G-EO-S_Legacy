#!/bin/bash
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/GPAC/trunk/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl DSIDDamage MaterialBase p_r0 alpha alpha p_r0 a1 parameter1 p_r0 a2 parameter2 p_r0 a3 parameter3 p_r0 a4 parameter4 p_r0 c0 initialDamageThreshold p_r0 c1 damageHardeningVariable p_r2s omega0 initialDamage s_r0s fd0 DamageIndicator s_r2s omega currentDamage s_r2s epsilon strain s_r2s epsid IrreversibleStrain s_r2s epsE TotalElasticStrain s_r2s epsel pureElasticStrain b_public ScriptInc/DSIDDamagePublic.cpp x
