#!/bin/bash
if [ -z "$GPAC_SCRIPT_PATH" ]; then
    #NOTE: this only provides a value for the current script
    export GPAC_SCRIPT_PATH=~/apps/geos_git/gpacspace/scripts
    echo "setting default for GPAC_SCRIPT_PATH: $GPAC_SCRIPT_PATH"
fi

#Base class
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl CohesiveZoneBase ConstitutiveBase CohesiveZone p_int unloadFlag unloadFlag s_r0 separationCoeff separationCoeff s_r1 traction traction s_r2 stiffness stiffness b_public ScriptInc/CohesiveZoneBasePublic.cpp x

#Materials for DE
perl $GPAC_SCRIPT_PATH/gpac_constitutive_model.pl InitiallyRigidCohesiveZone CohesiveZoneBase CohesiveZone p_r0 failStress failstress p_r0 failGap failgap s_r1 maxGap maxgap s_r0 maxTraction maxtraction b_public ScriptInc/InitiallyRigidCohesiveZonePublic.cpp x

perl -e 'open(IN,"<InitiallyRigidCohesiveZone.cpp");while(<IN>){if(/GeometryUtilities/){print "#include \"Utilities/GeometryUtilities.h\"\n#include \"DataStructures/VectorFields/ElementRegionT.h\"\n#include \"Constitutive/Material/MaterialBase.h\"\n";}else{print $_;}}close IN;' > tmp
mv tmp InitiallyRigidCohesiveZone.cpp
