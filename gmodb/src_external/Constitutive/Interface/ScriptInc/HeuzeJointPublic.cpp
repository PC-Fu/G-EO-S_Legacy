//FUNCTION_BEGIN_PARSE
virtual_void
HeuzeJoint::Initialize( const localIndex index,
                        const realT stressNormal,
                        const realT stressShear)
{
  //
  //NOTE: BE AWARE OF THE ASSUMPTIONS OF THIS PROCEDURE
  //
  //1) NORMAL DISPLACEMENT OCCURS BEFORE ANY SHEAR DISPLACEMENT
  //2) THE JOINT IS EQUILIBRIUM AT INIT ... THIS IS ENFORCED
  //3) THE JOINT HAS ONLY ELASTIC DILATION ... THIS IS ENFORCED

  HeuzeJointParameterData& matParams = *this->ParameterData(index);
  HeuzeJointStateData& matState = *this->StateData(index, 0);

  matState.stress = stressNormal;
  matState.stressShear = stressShear;

  //set shear strain
  matState.xs = -matState.stressShear / matParams.kshear;

  //set strength parameters, if necessary to keep it in equilibrium at initialization
  const realT xsmag = fabs(matState.xs);
  {
    if(xsmag >= matParams.xsResidual)
      matParams.xsResidual = xsmag * (1.0 + 1e-6);
    matState.mu = matParams.mu0;
    UpdateFriction(matParams, matState);
    const realT strength = ShearStrength( matParams, matState);
    if(strength < fabs(matState.stressShear))
    {
      const realT mu_new = fabs(matState.stressShear) / xsmag;
      matParams.muResidual += mu_new - matState.mu;
      matState.mu = mu_new;
    }
  }

  //now that we know sig_n, sig_s, and x_s, we need x_n
  //for now, let's force elastic dilation ...
  if(stressNormal > matParams.stressDilationLimit)
  {
    matParams.stressDilationLimit = stressNormal;
  }

  //1) Dilation elastic
  {
    const realT f1 = 0.4 * matParams.dilationCoefficient0 * xsmag * xsmag * sqrt(xsmag);
    const realT f3 = 2.0 / 3.0;

    //1a) Is this normal elastic? ... assume normal before
    realT xnelastic = 0;
    {
      //solve the cubic equation for x^1/2 -> elastic
      const realT f0 = -stressNormal / matParams.kCoefficientElastic;
      //solve for x^1/2 -> plastic

      const realT params[] = {f0, f1, f3};
      const realT xnelastic_sqrt = FindRoots::FindRoot(ZerosFunctionElasticElastic,
                                                       params,
                                                       0, 2 * matParams.normalApproachYield,
                                                       1e-6);
      xnelastic = xnelastic_sqrt * xnelastic_sqrt;
    }

    //check ...
    if(xnelastic <= matParams.normalApproachYield)
    {
      //Yes, it is still elastic
      matState.normalApproach = xnelastic;
    }
    else
    {
      //No, it is not elastic, so let's get the plastic part
      matState.normalApproach = matParams.kplasticLoad * f1;
      matState.normalApproach += f3 * matParams.normalApproachYield * sqrt(matParams.normalApproachYield);
      matState.normalApproach /= matParams.kplasticLoad;
      matState.normalApproach += matParams.normalApproachYield;
    }
  }
//  else
//  {
//    //So, now we know we may have a plastic dilation ... more complex
//
//  }
}
