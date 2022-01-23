#!/bin/bash
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpacspace/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

#ROCK JOINTS
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl JointIntermediate InterfaceBase Interface s_r0 gap gap s_r0 mu currentFrictionCoefficient s_r0 kcurrent currentStiffness s_r0 stressDilation dilationStress s_int ifail failureFlag p_r0 xsResidual shearAtResidual p_r0 muResidual coefficientOfFrictionAtResidual p_r0 cohesion cohesion p_r0 normalApproachYield normalApproachAtYield p_r0 kplasticLoad plasticStiffness p_r0 kplasticUnload unloadStiffness p_r0 kdilation dilationStiffness p_r0 dilationCoefficient0 dilationCoefficient p_r0 xsDilationInit shearAtInitiationOfDilation p_r0 xsDilationLimit shearAtLimitOfDilation p_r0 normalApproachDilationInit normalApproachAtInitiationOfDilation p_r0 stressDilationLimit stressAtLimitOfDilation p_r0 kCoefficientElastic elasticStiffnessCoefficient p_r0 kshear shearStiffness b_public ScriptInc/JointIntermediatePublic.cpp x b_protected ScriptInc/JointIntermediateProtected.cpp x

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl HeuzeJoint JointIntermediate Interface b_protected ScriptInc/HeuzeJointProtected.cpp x b_public ScriptInc/HeuzeJointPublic.cpp

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl ItascaJoint JointIntermediate Interface b_protected ScriptInc/ItascaJointProtected.cpp x


#PENALTY DERIVATIVES
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl RateAndStateIntermediate PenaltyCoulombIntermediate Interface s_r0 theta state s_r0 mu currentFrictionCoefficient s_r0 dstressdt stressRate p_r0 A A p_r0 B B p_r0 vstar shearRateStar p_r0 thetastar stateStar p_r0 Dc Dc p_r0 alpha alpha b_protected ScriptInc/RateAndStateIntermediateProtected.cpp x b_public ScriptInc/RateAndStateIntermediatePublic.cpp x

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl RateAndState RateAndStateIntermediate Interface

perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl Dieterich RateAndStateIntermediate Interface s_int nextState nextState s_int apFail apFail s_int pinned pinnedFlag s_int slowSlip slowSlipFlag s_int neighborInRuptureState neighborInRuptureState s_int currentState currentState s_int Areduced AreducedFlag s_r0 stressReference stressReference s_r0 dstressdtDrive stressRateDrive s_r0 dxsdtDrive shearRateDrive s_r0 stressPin stressPin s_r0 stressShearReference shearStressReference s_r0 dstressShear dstressShear s_r0 dstress dstress s_r0 dstressShearDt shearStressRate s_r0 dstressShearDtDrive shearStressRateDrive s_r0 dxs dshearOverTimestep s_r0 pp porePressure s_r0 dppdt dPorePressureDt s_r0 H H p_r0 tFail tFail p_r0 maxThetaPin maxStatePin s_r0 mu2 frictionCoefficient2 s_r0 mu2High frictionCoefficient2High s_r0 mu2aLow frictionCoefficient2aLow s_r0 mu2Low frictionCoefficient2Low s_r0 mu3High frictionCoefficient3High s_r0 mu3Low frictionCoefficient3Low p_r0 KShearSelf KShearSelf p_r0 KNormalSelf KNormalSelf p_r0 dxsdtAB shearSlipRateAB p_r0 stressShearFail shearStressAtFailure p_r0 AreductionFactor AReductionFactor p_r0 stressOvershootFactor stressOverShootFactor p_r1 rake unitRakeVector b_protected ScriptInc/DieterichProtected.cpp x b_public ScriptInc/DieterichPublic.cpp x
