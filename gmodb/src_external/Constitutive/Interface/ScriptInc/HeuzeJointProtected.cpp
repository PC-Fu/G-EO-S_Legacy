//FUNCTION_BEGIN_PARSE
virtual_void
HeuzeJoint::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                   InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const HeuzeJointParameterData& matParams = static_cast < const HeuzeJointParameterData& > ( matParamsBase );
  HeuzeJointStateData& matState = static_cast < HeuzeJointStateData& > ( matStateBase );

  const realT xsmag = fabs(matState.xs);
  if(xsmag < matParams.xsResidual)
  {
    //not yet failed ...
    const realT slope = xsmag / (matParams.xsResidual-xsmag);
    const realT phi1 = (1.0 - slope) * (matState.mu - matParams.muResidual);
    matState.mu = matParams.muResidual + ( (phi1 < 0.0) ? 0.0 : phi1 );
  }
  else
  {
    //just failed ...
    matState.ifail = 1;
    matState.mu = matParams.muResidual;
  }
}

//FUNCTION_BEGIN_PARSE
virtual_realT
HeuzeJoint::DilationalStressIncrement( const InterfaceBaseParameterData& matParamsBase,
                                       InterfaceBaseStateData& matStateBase,
                                       const realT dxs) const
{
  HeuzeJointStateData& matState = static_cast<HeuzeJointStateData&> (matStateBase);
  const HeuzeJointParameterData& matParams = static_cast<const HeuzeJointParameterData&> (matParamsBase);

  const realT xsmag = fabs(matState.xs);
  const realT dstress = matState.kcurrent * matParams.dilationCoefficient0 * sqrt(xsmag) * dxs;
  const realT stress0 = matState.stress + matState.stressDilation;
  const realT stress1 = stress0 + dstress;

  return stress1 < matParams.stressDilationLimit ? dstress :
      ((matParams.stressDilationLimit - stress0) > 0 ?
          matParams.stressDilationLimit - stress0 :
          0.0);
}

//FUNCTION_BEGIN_PARSE
virtual_realT
HeuzeJoint::NormalApproachAtInitialization( const realT stressNormal,
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
  {
    return normalApproach;
  }

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
static_realT
HeuzeJoint::ZerosFunctionElasticElastic(const realT x, const realT* params )
{
  return params[2] * x * x * x + params[1]  * x + params[0];
}
