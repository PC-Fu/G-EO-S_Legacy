//FUNCTION_BEGIN_PARSE
virtual_void
ItascaJoint::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                             InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const ItascaJointParameterData& matParams = static_cast < const ItascaJointParameterData& > ( matParamsBase );
  ItascaJointStateData& matState = static_cast < ItascaJointStateData& > ( matStateBase );

  if(matState.ifail == 1)
  {
    //once failed, always failed
    matState.mu = matParams.muResidual;
  }
  else
  {
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
}

//FUNCTION_BEGIN_PARSE
virtual_realT
ItascaJoint::DilationalStressIncrement( const InterfaceBaseParameterData& matParamsBase,
                                       InterfaceBaseStateData& matStateBase,
                                       const realT dxs) const
{
  ItascaJointStateData& matState = static_cast<ItascaJointStateData&> (matStateBase);
  const ItascaJointParameterData& matParams = static_cast<const ItascaJointParameterData&> (matParamsBase);

  //check whether dilation is set
  if(matParams.dilationCoefficient0 <= 0)
    return 0.0;

  const realT xsmag = fabs(matState.xs);

  // Check if we have reached relative tangential displacement limit
  if (xsmag < matParams.xsDilationLimit &&
      matState.normalApproach > matParams.normalApproachDilationInit)
  {
    // We have not reached relative tangential displacement limit
    // Dilate using tangentOfDilationAngleAtZeroStress
    // Soften tangentOfDilationAngleAtZeroStress with increasing relative tangential displacement
    realT stressSoftening = matParams.xsDilationLimit > 0.0 ?
        (matParams.xsDilationLimit - xsmag) / matParams.xsDilationLimit :
        1.0;

    // Soften tangentOfDilationAngleAtZeroStress with increasing normal stress
    if (matParams.stressDilationLimit > 0.0)
    {
      stressSoftening *= (1.0 - matState.stress / matParams.stressDilationLimit);
      stressSoftening = stressSoftening < 0 ? 0 : (stressSoftening > 1.0 ? 1.0 : stressSoftening);
    }

    // Replace stiffness with an equilibrium stiffness
    realT kcurr = matParams.kCoefficientElastic * (matParams.kdilation > 0.0 ? matParams.kdilation /
        (matState.kcurrent + matParams.kdilation) : 1.0);
    kcurr *= stressSoftening * matParams.dilationCoefficient0;

    //return change in stress wrt time
    return kcurr * sqrt(xsmag) * dxs;
  }

  // Here we are either in residual or before the onset of dilation
  return 0.0;
}
