//FUNCTION_BEGIN_PARSE
virtual_void
JointIntermediate::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                   InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const JointIntermediateParameterData& matParams = static_cast < const JointIntermediateParameterData& > ( matParamsBase );
  JointIntermediateStateData& matState = static_cast < JointIntermediateStateData& > ( matStateBase );

  matState.mu = matParams.mu0;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
JointIntermediate::ShearStrength(const InterfaceBaseParameterData& ,
                                 InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  JointIntermediateStateData& matState = static_cast < JointIntermediateStateData& > ( matStateBase );

  return matState.stress * matState.mu;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
JointIntermediate::DilationalStressIncrement( const InterfaceBaseParameterData&,
                                              InterfaceBaseStateData& ,
                                              const realT ) const
{
  return 0.0;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
JointIntermediate::NormalApproachAtInitialization( const realT stressNormal,
                                                   const realT knCoefficientElastic,
                                                   const realT knPlasticLoading,
                                                   const realT knPlasticUnloading,
                                                   realT& kcurrent,
                                                   realT& normalGap,
                                                   const realT normalApproachNormalYield) const
{
  //Let's make the assumption that the loading monotonically increasing in time before initialization
  //otherwise, must assume that deformation has been elastic in the past ... if not, we're indeterminate
  if(stressNormal <= 0.0)
    return 0.0;

  realT normalApproach = pow(stressNormal / knCoefficientElastic, 2.0/3.0);
  kcurrent = knCoefficientElastic * sqrt(normalApproach);

  if(normalApproach <= normalApproachNormalYield)
    return normalApproach;

  if (knPlasticLoading > knPlasticUnloading)
    throw GPException(
        "Cannot have a kplasticLoad > kplasticUnload without injecting energy!");
  else if (isZero(knPlasticUnloading))
    throw GPException("Cannot have a zero slope plastic unloading curve");

  //get the solution
  if (knPlasticLoading > 0.0)
  {
    normalApproach = stressNormal;
    normalApproach -= knCoefficientElastic * sqrt(normalApproachNormalYield) * normalApproachNormalYield;
    normalApproach /= knPlasticLoading;
    normalApproach += normalApproachNormalYield;
  } //else, normal approach was already calculated correctly!

  const realT fpu = knPlasticUnloading * (normalApproach - normalGap);
  if (fpu > stressNormal)
    normalGap = normalApproach - (stressNormal / knPlasticUnloading);

  kcurrent = knPlasticUnloading;
  return normalApproach;
}

//FUNCTION_BEGIN_PARSE
realT
JointIntermediate::NormalStiffnessRock( const realT knCoefficientElastic,
                                        const realT knPlasticLoading,
                                        const realT knPlasticUnloading,
                                        const realT normalApproach,
                                        realT& normalGap,
                                        const realT normalApproachNormalYield) const
{
  if (normalApproach < 0.0)
    return 0.0;

  const bool noGap = isZero(normalGap);

  //if elastic, give the elastic stiffness, otherwise the plastic unloading stiffness
  if (normalApproach <= normalApproachNormalYield && noGap)
    return knCoefficientElastic * sqrt(normalApproach);

  //adjust the gap appropriately
  if (normalApproach > normalApproachNormalYield)
  {
    //handle case of simple linear contact
    if (isZero(knPlasticLoading - knPlasticUnloading) && noGap)
      return knPlasticLoading;
    //handle case of energy-injecting contact
    else if (knPlasticLoading > knPlasticUnloading)
      throw GPException(
          "Cannot have a kplasticLoad > kplasticUnload without injecting energy!");
    else if (isZero(knPlasticUnloading))
      throw GPException("Cannot have a zero slope plastic unloading curve");

    const realT fpu = knPlasticUnloading * (normalApproach - normalGap);

    realT fpl = 0.0;
    if (knPlasticLoading > 0.0)
    {
      fpl = knPlasticLoading * (normalApproach - normalApproachNormalYield);
      fpl += knCoefficientElastic * sqrt(normalApproachNormalYield) * normalApproachNormalYield;
    }
    else
    {
      fpl = knCoefficientElastic * sqrt(normalApproach);
      if (fpl > knPlasticUnloading)
        fpl = knPlasticUnloading;
      fpl *= normalApproach;
    }

    if (fpu > fpl)
    {
      normalGap = normalApproach - (fpl / knPlasticUnloading);
    }
  }
  return knPlasticUnloading;
}

//FUNCTION_BEGIN_PARSE
realT
JointIntermediate::ShearStrength(const realT stress,
                                 const realT normalStressAtDilationLimit,
                                 const realT tanFrictionCoefficientInitial,
                                 const realT tanFrictionCoefficientResidual,
                                 const realT cohesion,
                                 const int ifail) const
{
  //note: this is the effect of Heuze and Itasca models
  //don't accept if the friction coefficient too close to 90-degrees
  if(tanFrictionCoefficientInitial > 1.0e5)
    return 0.0;

  // NOTE: Otis Walton (Private communication to Joe Morris: 22.FEB.2006) points out that
  // the shear strength should be determined by the SUM of elastic and viscous terms
  // The issue is that at high rate, viscous stresses will be a significant fraction of
  // the total stress and you can no longer argue that viscosity is simply a numerical
  // convenience
  if(isZero(tanFrictionCoefficientResidual - tanFrictionCoefficientInitial))
    return tanFrictionCoefficientInitial * stress + cohesion;

  realT shearStrength = 0.0;
  if (stress > 0.0)
  {
    if(ifail == 1)
    {
      shearStrength += tanFrictionCoefficientResidual * stress;
    }
    else
    {
      // We must figure where on the shear strength
      // envelope we are
      if (stress < normalStressAtDilationLimit)
      {
        // We are below critical normal stress
        // Use initial friction angle
        shearStrength += tanFrictionCoefficientInitial * stress;
      }
      else
      {
        // Initial portion
        shearStrength += tanFrictionCoefficientInitial * normalStressAtDilationLimit;
        // Portion past critical
        shearStrength += tanFrictionCoefficientResidual * (stress - normalStressAtDilationLimit);
      }
    }
  }
  //Add in the contribution due to cohesion
  shearStrength += cohesion;
  return shearStrength;
}
