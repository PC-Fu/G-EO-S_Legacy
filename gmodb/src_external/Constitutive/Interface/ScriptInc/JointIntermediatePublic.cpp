//FUNCTION_BEGIN_PARSE
virtual_void
JointIntermediate::Initialize( const localIndex index,
                               const realT stressNormal,
                               const realT stressShear)
{
  JointIntermediateParameterData& matParams = *this->ParameterData(index);
  JointIntermediateStateData& matState = *this->StateData(index, 0);

  //Normal stress
  matState.stress = stressNormal;

  //Shear stress
  UpdateFriction(matParams, matState);
  matState.stressShear = stressShear;
  const realT strength = ShearStrength( matParams, matState);
  if(strength < fabs(matState.stressShear))
  {
    //At initialization, the system is necessarily at equilibrium, so we
    //need to set the joint at critical stress ... by increasing strength
    matState.mu = matParams.mu0 = fabs(matState.stressShear) / matParams.kshear;
  }
  matState.xs = -matState.stressShear / matParams.kshear;

  //Normal stress
  matState.normalApproach = NormalApproachAtInitialization(stressNormal,
                                                           matParams.kCoefficientElastic,
                                                           matParams.kplasticLoad,
                                                           matParams.kplasticUnload,
                                                           matState.kcurrent,
                                                           matState.gap,
                                                           matParams.normalApproachYield);
}

//FUNCTION_BEGIN_PARSE
virtual_void
JointIntermediate::StrainDrivenUpdate( const localIndex index )
{
  const JointIntermediateParameterData& matParams = *this->ParameterData(index);
  JointIntermediateStateData& matState = *this->StateData(index, 0);

  //(1) evolve the normal stiffness state
  matState.kcurrent = NormalStiffnessRock(matParams.kCoefficientElastic, matParams.kplasticLoad, matParams.kplasticUnload,
            matState.normalApproach, matState.gap, matParams.normalApproachYield);
  //const realT normalApproachNext = matState.normalApproach + matState.dxndt * matState.dt;
  //const realT kNext = normalApproachNext > matState.normalApproach ? NormalStiffnessRock(matParams.kCoefficientElastic, matParams.kplasticLoad, matParams.kplasticUnload,
  //                                                                                   normalApproachNext, matState.gap, matParams.normalApproachYield) : matState.kcurrent;

  //(2) evolve the normal stress
  matState.stress = matState.kcurrent * (matState.normalApproach - matState.gap);

  //(3) evolve the normal stress due to dilation
  const realT dxs = matState.dxsdt * matState.dt;
  matState.stressDilation += DilationalStressIncrement( matParams, matState, dxs);
  matState.stress += matState.stressDilation;

  //(4) evolve the friction
  UpdateFriction(matParams, matState);

  //(5) evolve the shear stress
  matState.stressShear += -matParams.kshear * dxs;
  matState.xs += dxs;

  //(6) determine whether strength exceeded and return to failure surface if necessary
  ThresholdToFailureSurface(matParams, matState, dxs);

  matState.stressShearVector.Normalize();
  matState.stressShearVector *= matState.stressShear;
  return;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
JointIntermediate::StiffnessProjected(const localIndex index)
{
  const JointIntermediateParameterData& matParams = *this->ParameterData(index);
  JointIntermediateStateData& matState = *this->StateData(index, 0);
  const realT normalApproachEstimate = matState.normalApproach +
      (matState.dxndt > 0 ? matState.dxndt * matState.dt : 0.0);
  return NormalStiffnessRock(matParams.kCoefficientElastic, matParams.kplasticLoad, matParams.kplasticUnload,
           normalApproachEstimate, matState.gap, matParams.normalApproachYield);
}
