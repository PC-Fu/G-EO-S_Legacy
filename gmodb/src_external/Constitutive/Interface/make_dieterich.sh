#!/bin/bash
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/gpac_trunk/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl Dieterich RateAndStateIntermediate s_int nextState nextState s_int apFail apFail s_int pinned pinnedFlag s_int slowSlip slowSlipFlag s_int neighborInRuptureState neighborInRuptureState s_int currentState currentState s_int Areduced AreducedFlag s_r0 stressReference stressReference s_r0 dstressdtDrive stressRateDrive s_r0 dxsdtDrive shearRateDrive s_r0 stressPin stressPin s_r0 stressShearReference shearStressReference s_r0 dstressShear dstressShear s_r0 dstress dstress s_r0 dstressShearDt shearStressRate s_r0 dstressShearDtDrive shearStressRateDrive s_r0 dxs dshearOverTimestep s_r0 pp porePressure s_r0 dppdt dPorePressureDt s_r0 H H p_r0 tFail tFail p_r0 maxThetaPin maxStatePin s_r0 mu2 frictionCoefficient2 s_r0 mu2High frictionCoefficient2High s_r0 mu2aLow frictionCoefficient2aLow s_r0 mu2Low frictionCoefficient2Low s_r0 mu3High frictionCoefficient3High s_r0 mu3Low frictionCoefficient3Low p_r0 KShearSelf KShearSelf p_r0 KNormalSelf KNormalSelf p_r0 dxsdtAB shearSlipRateAB p_r0 stressShearFail shearStressAtFailure p_r0 AreductionFactor AReductionFactor p_r0 stressOvershootFactor stressOverShootFactor p_r1 rake unitRakeVector b_protected ScriptInc/DieterichProtected.cpp x b_public ScriptInc/DieterichPublic.cpp x
