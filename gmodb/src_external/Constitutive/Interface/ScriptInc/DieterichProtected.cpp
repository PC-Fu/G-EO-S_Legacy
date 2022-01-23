//FUNCTION_BEGIN_PARSE
virtual_void
Dieterich::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                           InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const DieterichParameterData& matParams = static_cast < const DieterichParameterData& > ( matParamsBase );
  DieterichStateData& matState = static_cast < DieterichStateData& > ( matStateBase );

  const realT A_ = A(matParams, matState);
  matState.mu = matParams.mu0 + A_ * log( 1.0 + matState.dxsdt / matParams.vstar) +
      matParams.B * log(1.0 + matState.theta / matParams.thetastar);
  const realT dthetadt = 1.0 - matState.theta * matState.dxsdt / matParams.Dc -
      (matParams.alpha * matState.theta * (matState.dstressdt-matState.dppdt)) / (matParams.B * matState.stress);
  matState.theta += dthetadt * matState.dt;
}

//FUNCTION_BEGIN_PARSE
bool
Dieterich::NextTransitionTime( const DieterichParameterData& matParams,
                               DieterichStateData& matState,
                               EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep,
                               const realT dxsdtEQ,
                               const R1Tensor& center,
                               const localIndex i) const
{
  //pinned fault patches are always locked
  if (matState.pinned == 1)
    return false;

  const realT tau = matState.stressShear;

  //set shortcuts
  const realT effectiveStress = matState.stress;
  realT dEffectiveStressDt = (matState.dstressdt-matState.dppdt);
  const EarthquakeSimulation::TransitionState currentStateTrial = EarthquakeSimulation::TransitionStateFromIndex(
      matState.currentState);

  //set failure to not yet
  if (matState.apFail == EarthquakeSimulation::NOW)
    matState.apFail = EarthquakeSimulation::NOTYET;
  const EarthquakeSimulation::FailureState apFailTrial = EarthquakeSimulation::FailureStateFromIndex(
      matState.apFail);

  //LOW_SIGMA
  LowStressCheck(dEffectiveStressDt, effectiveStress, matState.stressPin, matState.dstressShearDt,
                 tau, i, ruptureTimestep);

  //HIGH_THETA
  HighThetaCheck(ruptureTimestep.m_currentTime, matState.theta, matParams.maxThetaPin, i, ruptureTimestep);

  //CHECK OTHER CONDITIONS
  const realT A_ = A(matParams, matState);

  switch (currentStateTrial)
  {
    //THESE ARE HANDLED ABOVE
    case EarthquakeSimulation::POREPRESSURERATECHANGE:
    case EarthquakeSimulation::LOW_SIGMA:
    case EarthquakeSimulation::LOW_TAU:
    case EarthquakeSimulation::HIGH_THETA:
      break;
    case EarthquakeSimulation::LOCK:
      TransitionTimeLocked(dEffectiveStressDt, effectiveStress, matState.dstressShearDt, tau,
                           matParams.vstar, matParams.alpha, A_, matParams.B, matParams.Dc,
                           matState.mu, matState.theta,
                           i, ruptureTimestep);
      TransitionTimeAPriori(ruptureTimestep.m_currentTime, matParams.tFail, apFailTrial, i, ruptureTimestep);
      break;
    case EarthquakeSimulation::NUCLEATE:
      TransitionTimeNucleate(dEffectiveStressDt, effectiveStress, matState.dstressShearDt,
                             tau, matState.dxsdt, matParams.vstar,
                             dxsdtEQ, matParams.alpha, A_, matParams.B, matParams.Dc, matParams.KShearSelf, matState.theta,
                             matState.slowSlip, matState.neighborInRuptureState, matState.H,
                             i, ruptureTimestep);
      TransitionTimeAPriori(ruptureTimestep.m_currentTime, matParams.tFail, apFailTrial,
                            i, ruptureTimestep);
      break;
    case EarthquakeSimulation::RUPTURE:
      TransitionTimeRupture(dEffectiveStressDt, effectiveStress, matState.dstressShearDt,
                            tau, matState.mu2,
                            i, ruptureTimestep);
      break;
    case EarthquakeSimulation::CREEP:
      TransitionTimeCreep(dEffectiveStressDt, effectiveStress, matState.dstressShearDt, tau,
                          matState.mu3Low, matState.mu3High,
                          i, ruptureTimestep);
      break;
    case EarthquakeSimulation::SLOWSLIP_2A:
      TransitionTimeSlowSlipA(dEffectiveStressDt, effectiveStress,
                              matState.dstressShearDt, tau,
                              matState.mu2High, matState.mu2aLow,
                              i, ruptureTimestep);
      break;
    case EarthquakeSimulation::SLOWSLIP_2B:
      TransitionTimeSlowSlipB(dEffectiveStressDt, effectiveStress, matState.dstressShearDt,
                              tau, matState.mu, matState.mu2High,
                              matState.mu2Low, matParams.dxsdtAB, matParams.vstar,
                              dxsdtEQ, A_,
                              i, ruptureTimestep);
      break;
    case EarthquakeSimulation::SLOWSLIP_2C:
      TransitionTimeSlowSlipC(dEffectiveStressDt, effectiveStress,
                              matState.dstressShearDt, tau,
                              matState.mu2High, matState.mu2Low,
                              i, ruptureTimestep);
      break;
    default:
      break;
  }
  return true;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
Dieterich::ShearStrength(const InterfaceBaseParameterData& ,
                         InterfaceBaseStateData& matStateBase) const
{
  DieterichStateData& matState = static_cast<DieterichStateData&> (matStateBase);
  return matState.stress * matState.mu;
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::IncrementalShearDisplacement(const EarthquakeSimulation::TransitionState state,
                                        const realT dt,
                                        const realT dxsdtEQ,
                                        const DieterichParameterData& matParams,
                                        DieterichStateData& matState) const
{
  return IncrementalShearDisplacement(state, dt, matState.stress, (matState.dstressdt-matState.dppdt), matState.stressShear, matState.dstressShearDt,
                                      dxsdtEQ, matParams.Dc, matParams.alpha, matParams.KShearSelf, A(matParams, matState), matParams.B,
                                      matState.dxsdt, matState.H);
}

//FUNCTION_BEGIN_PARSE
static_realT
Dieterich::IncrementalShearDisplacement(const EarthquakeSimulation::TransitionState state,
                                        const realT dt,
                                        const realT sigma0,
                                        const realT sigmaDot0,
                                        const realT tau0,
                                        const realT tauDot0,
                                        const realT slipRateEQ,
                                        const realT Dc,
                                        const realT alpha,
                                        const realT KShearSelf,
                                        const realT A_,
                                        const realT B,
                                        realT& ddot,
                                        realT& H)
{
  if(state != EarthquakeSimulation::NUCLEATE)
    return state == EarthquakeSimulation::LOCK ? 0.0 : ddot * dt;

  const realT dbleps = std::numeric_limits<realT>::min();
  const realT tol=1e-7;

  //set the new ddot
  const realT SDotExt = tauDot0 - ((tau0/sigma0) - alpha)*sigmaDot0;
  {
    realT ddot1r;  /* value to be returned */

    if (fabs((SDotExt/(A_*sigma0))*dt) < tol)
      ddot1r = 1/( (1/ddot)*exp(-(SDotExt/(A_*sigma0))*dt) - (H/A_)*dt);
    else
      ddot1r = 1/( (1/ddot + H*sigma0/SDotExt)*exp(-(SDotExt/(A_*sigma0))*dt) -
                  H*sigma0/SDotExt);

    /* assume that if ddot1r is too big or negative or Inf or NaN that it was supposed
       to go to ddotEQ but missed due to roundoff error.  FIXME */
    if (std::isinf(ddot1r) || ddot1r < 0 || ddot1r > slipRateEQ)
      ddot1r = slipRateEQ;

    ddot = ddot1r;
  }

  const realT f1 = fabs((SDotExt/(A_*sigma0))*dt) < tol ?
      -(ddot*H/A_)*dt :
      (ddot*H*sigma0/SDotExt)*( 1 - exp((SDotExt/(A_*sigma0))*dt) );

  //update H
  H = B/Dc + KShearSelf/sigma0;

  //get the incremental displacement
  const realT dd = fabs(f1) < tol ?
      -(A_/H)*f1 :
      ((f1 <= -1) ? -(A_/H)*log(dbleps) :   /* COMPLETE HACK  FIXME */
      -(A_/H)*log(1 + f1));

  return dd;
}

//FUNCTION_BEGIN_PARSE
static_realT
Dieterich::Theta(const EarthquakeSimulation::TransitionState state,
                 const DieterichParameterData& matParams,
                 DieterichStateData& matState)
{
  return Theta(state, matState.dt, (matState.dstressdt-matState.dppdt), matState.stress, matParams.alpha, matParams.B, matParams.Dc, matState.dxs, matState.theta);
}

//FUNCTION_BEGIN_PARSE
static_realT
Dieterich::Theta(const EarthquakeSimulation::TransitionState state,
                 const realT dt,
                 const realT sigmaDot,
                 const realT sigma0,
                 const realT alpha,
                 const realT B,
                 const realT Dc,
                 const realT dd,
                 const realT theta0)
{
  const realT sigma = sigma0 + dt * sigmaDot;
  const realT theta = (state == EarthquakeSimulation::LOCK) ?
      ThetaLock(dt, sigmaDot, sigma0, alpha, B, theta0) :
      ((state == EarthquakeSimulation::NUCLEATE) ?
          ThetaNucleateAlternate(theta0, dd, Dc, alpha, B, sigma0, sigma) :
          theta0);
  return theta;
}

//FUNCTION_BEGIN_PARSE
static_realT
Dieterich::ThetaLock(const realT dt,
                     const realT deffectiveStressDt,
                     const realT effectiveStress,
                     const realT alpha,
                     const realT B,
                     const realT theta0)
{
  const realT tol = 1e-10;
  const realT sigRatio = deffectiveStressDt*dt/effectiveStress;
  const realT onePlusSigRatio = 1 + sigRatio;
  const realT alphaOverB = alpha/B;

  if(fabs(sigRatio) < tol)
  {
    const realT theta = theta0*(1 - alphaOverB*sigRatio) + dt;
    return theta;
  }
  else if (onePlusSigRatio <= 0)
  {
    return std::numeric_limits<realT>::max();
  }
  else
  {
    const realT onePlusAlphaOverB = 1+alphaOverB;
    const realT p1 = pow(onePlusSigRatio, -alphaOverB);
    const realT p2 = pow(onePlusSigRatio, onePlusAlphaOverB);
    const realT theta = p1 * (theta0 + (p2 - 1) * (effectiveStress/(deffectiveStressDt * onePlusAlphaOverB)));
    return theta;
  }
}

//FUNCTION_BEGIN_PARSE
static_realT
Dieterich::ThetaNucleate(const realT tau,
                         const realT sigma,
                         const realT ddot,
                         const realT mu0,
                         const realT ddotStar,
                         const realT A_,
                         const realT B,
                         const realT Dc)
{
  return((Dc/ddotStar) * exp( (tau - mu0*sigma)/(B*sigma))*pow(ddot/ddotStar, -A_/B));
}

//FUNCTION_BEGIN_PARSE
static_realT
Dieterich::ThetaNucleateAlternate( const realT theta0,
                                   const realT d,
                                   const realT Dc,
                                   const realT alpha,
                                   const realT B,
                                   const realT sigma0,
                                   const realT sigma)
{
  return( theta0*exp(-d/Dc - (alpha/B)*log(sigma/sigma0)) );
}

//FUNCTION_BEGIN_PARSE
static_void
Dieterich::TransitionTimeAPriori(const realT time,
                                 const realT timeFail,
                                 const EarthquakeSimulation::FailureState apFail,
                                 const localIndex local,
                                 EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  const realT dt = timeFail - time;
  if (dt >= 0 && apFail == EarthquakeSimulation::NOTYET)
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::NUCLEATE);
}

///FUNCTION_BEGIN_PARSE
static_void
Dieterich::TransitionTimeLocked(const realT dEffectiveStressDt,
                                const realT effectiveStress,
                                const realT dStressShearDt,
                                const realT stressShear,
                                const realT dShearSlipDtStar,
                                const realT alpha, const realT A_, const realT B,
                                const realT Dc, const realT frictionCoefficient0,
                                const realT theta0, const localIndex local,
                                EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  const realT params[] = {dEffectiveStressDt, effectiveStress,
                         dStressShearDt, stressShear,
                         dShearSlipDtStar, alpha, A_, B, Dc,
                         frictionCoefficient0, theta0};

  /* check sign of f() at dt = 0, if positive, then this patch should
     already have transitioned out of state 0 (this can happen at model
     startup, depending on the initial values of tau and theta)  */
  const realT signf0 = theta0 <= 0 ? 
    ((B-A_) > 0 ? 1 : -1) : 
    (f(0, params) > 0.0 ? 1 : -1);
  if (signf0 > 0)
  {
    ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::NUCLEATE);
    return;
  }

  /* -1 is just a flag so we can tell if it has been reset
                         by any of the if blocks below */
  realT tTrans=-1;

  /* if sigmaDot == 0 and matState.dstressShearDt is negative, then this patch
     will stay in state 0; if solved the equation, might get
     a negative dt, which is invalid because we cannot extrapolate
     tau backwards */
  if (isEqual(dEffectiveStressDt, 0.0) && dStressShearDt <= 0.0)
    tTrans = std::numeric_limits<realT>::max();

  /* if sigmaDot < 0 and tau at -sigma0/sigmaDot is less than zero
     then the patch will stay in state 0 */
  else if (dEffectiveStressDt < 0 && stressShear + dStressShearDt * (-effectiveStress / dEffectiveStressDt) < 0.0)
    tTrans = std::numeric_limits<realT>::max();

  /* next check the sign of f() at the current smallest dt,
     does this still work if sigmaDot != 0 ??? */
  else
  {
    if (ruptureTimestep.m_dt < std::numeric_limits<realT>::max())
    {
      const realT fminDt = f(ruptureTimestep.m_dt, params);
      if (fminDt < 0.0)
      {
        tTrans = std::numeric_limits<realT>::max();
      }
    }
  }

  if (tTrans < 0) /* that is, if we haven't reset it in any of the of blocks above, then root find */
  {
    realT tolTmp = isZero(dEffectiveStressDt) ? std::numeric_limits<realT>::max() : fabs(0.99*effectiveStress/dEffectiveStressDt);
    tolTmp = ruptureTimestep.m_dt < tolTmp ? ruptureTimestep.m_dt : tolTmp;
    tolTmp = 1e30 < tolTmp ? 1e30 : tolTmp;
    tTrans = zbrent(params, 0.0, tolTmp, 1.0e-10);
    //    fprintf(stdout, "S0TransTime: minmin = %f tTrans = %e\n", tolTmp, tTrans);//DEBUG
  }

  //Check to see whether this is controlling
  ruptureTimestep.CheckForTimestepControl(tTrans, local, EarthquakeSimulation::NUCLEATE);
}



///FUNCTION_BEGIN_PARSE
static_void
Dieterich::TransitionTimeNucleate(const realT dEffectiveStressDt,
                                  const realT effectiveStress,
                                  const realT dStressShearDt,
                                  const realT stressShear,
                                  const realT dShearSlipDt,
                                  const realT dShearSlipDtStar,
                                  const realT dShearSlipDtEQ,
                                  const realT alpha,
                                  const realT A_,
                                  const realT B,
                                  const realT Dc,
                                  const realT KShearSelf,
                                  const realT theta0,
                                  const int slowSlip,
                                  const int neighborInRuptureState,
                                  realT& H,
                                  const localIndex local,
                                  EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  const realT tol = 1e-7;
  const realT ddotTrans = (slowSlip == 1 ? dShearSlipDtStar : dShearSlipDtEQ);

  //Check whether we need to immediately lock
  {
    const realT tol_ddot = isZero(theta0) ? std::numeric_limits<realT>::max() : (1-tol) * Dc/theta0;
    if (dShearSlipDt < tol_ddot && dShearSlipDt < ddotTrans && neighborInRuptureState == 0)
    {
      /* then transition immediately back to state 0 */
//      fprintf(stdout, "S1TransTime: dt = 0 nextState = 0\n");//DEBUG
      ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::LOCK);
    }
  }

  //set H
  H = B/Dc + KShearSelf/effectiveStress;

  //set time step from this criterion
  realT tTrans;
  {
    const realT SDotExt = dStressShearDt - ((stressShear/effectiveStress) - alpha) * dEffectiveStressDt;
//    fprintf(stdout, "S1TransTime: SDotExt (%e %e %e %f %e) = %f p->H = %f ", dStressShearDt, stressShear, effectiveStress, alpha, dEffectiveStressDt, SDotExt, H);//DEBUG

    const realT hsigddot = H * effectiveStress * dShearSlipDt;
    if (SDotExt <= -hsigddot)
    {
      tTrans = std::numeric_limits<realT>::max(); /* really should look for transition back to state 0 */
    }
    else if (fabs(SDotExt/ ( H * effectiveStress * ddotTrans) ) < tol)
    {
      realT invddot = 1/dShearSlipDt;
      invddot -= 1/ddotTrans;
      const realT f1 = 1 + (SDotExt/(H * effectiveStress))*invddot;
//      fprintf(stdout, " f1 = %f ddot = %e ddotTrans = %f ", f1, dShearSlipDt, ddotTrans);//DEBUG
      if(fabs(1 - f1) < tol)
      {
        tTrans = (A_ / H) * invddot;
//        fprintf(stdout, " tTrans(a) = %e\n", tTrans);//DEBUG
      }
      else
      {
        tTrans = (A_ * effectiveStress / SDotExt) * log(f1);
//        fprintf(stdout, " tTrans(b) = %e\n", tTrans);//DEBUG
      }
    }
    else
    {
      tTrans = (A_ * effectiveStress / SDotExt) * log( (1 + SDotExt/hsigddot) /
                                        (1 + SDotExt/(H * effectiveStress * ddotTrans) ) );
//      fprintf(stdout, " tTrans(c) = %e\n", tTrans);//DEBUG
    }
  }

  //Check for timestep control
  ruptureTimestep.CheckForTimestepControl(tTrans, local, (slowSlip==1 ? EarthquakeSimulation::SLOWSLIP_2A : EarthquakeSimulation::RUPTURE));
}

///FUNCTION_BEGIN_PARSE
static_void
Dieterich::TransitionTimeRupture(const realT dEffectiveStressDt,
                                 const realT effectiveStress,
                                 const realT dStressShearDt,
                                 const realT stressShear,
                                 const realT frictionCoefficient2,
                                 const localIndex local,
                                 EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  if (dStressShearDt - frictionCoefficient2 * dEffectiveStressDt >= 0.0)
  {
    ruptureTimestep.CheckForTimestepControl(std::numeric_limits<realT>::max(), local, EarthquakeSimulation::RUPTURE);
    return;
  }

  //tau should never be less than muTrans2*sigma, if it is transition out immediately
  if (stressShear < frictionCoefficient2 * effectiveStress)
  {
    ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::LOCK);
  }
  else
  {
    const realT dt = (effectiveStress * frictionCoefficient2 - stressShear)/
        (dStressShearDt - frictionCoefficient2 * dEffectiveStressDt);
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::LOCK);
  }
}

///FUNCTION_BEGIN_PARSE
static_void
Dieterich::TransitionTimeSlowSlipA(const realT dEffectiveStressDt,
                                   const realT effectiveStress,
                                   const realT dStressShearDt,
                                   const realT stressShear,
                                   const realT frictionCoefficient2High,
                                   const realT frictionCoefficient2aLow,
                                   const localIndex local,
                                   EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  //Check 2A against friction coefficients
  {
    const realT mu = stressShear / effectiveStress;
    if (mu > frictionCoefficient2High)
    {
      ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::SLOWSLIP_2B);
      return;
    }
    else if (mu < frictionCoefficient2aLow)
    {
      ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::LOCK);
      return;
    }
  }

  const realT tLow = (frictionCoefficient2aLow * effectiveStress - stressShear)
        / (dStressShearDt - frictionCoefficient2aLow * dEffectiveStressDt);
  const realT tHigh = (frictionCoefficient2High * effectiveStress - stressShear)
      / (dStressShearDt - frictionCoefficient2High * dEffectiveStressDt);
  realT dt = std::numeric_limits<realT>::max();
  if (tLow > 0)
    dt = tLow;
  if (tHigh > 0 && tHigh < dt)
    dt = tHigh;

  if (isEqual(dt,tLow))
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::LOCK);
  else if (isEqual(dt,tHigh))
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::SLOWSLIP_2B);
  else
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::SLOWSLIP_2A);
}

///FUNCTION_BEGIN_PARSE
static_void
Dieterich::TransitionTimeSlowSlipB(const realT dEffectiveStressDt,
                                   const realT effectiveStress,
                                   const realT dStressShearDt,
                                   const realT stressShear,
                                   const realT frictionCoefficient0,
                                   const realT frictionCoefficient2High,
                                   const realT frictionCoefficient2Low,
                                   const realT dShearSlipDtAB,
                                   const realT dShearSlipDtStar,
                                   const realT dShearSlipDtEQ,
                                   const realT A_,
                                   const localIndex local,
                                   EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  //Check friction coefficient limits
  {
    const realT mu = stressShear / effectiveStress;
    if (mu < frictionCoefficient0 + A_ * log(dShearSlipDtAB / dShearSlipDtStar))
    {
      ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::SLOWSLIP_2A);
      return;
    }
    else if (mu > frictionCoefficient0 + A_ * log(dShearSlipDtEQ / dShearSlipDtStar))
    {
      ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::SLOWSLIP_2C);
      return;
    }
  }

  const realT tLow = (frictionCoefficient2Low * effectiveStress - stressShear)
        / (dStressShearDt - frictionCoefficient2Low * dEffectiveStressDt);
  const realT tHigh = (frictionCoefficient2High * effectiveStress - stressShear)
      / (dStressShearDt - frictionCoefficient2High * dEffectiveStressDt);
  realT dt = std::numeric_limits<realT>::max();
  if (tLow > 0)
    dt = tLow;
  if (tHigh > 0 && tHigh < dt)
    dt = tHigh;

  if (isEqual(dt,tLow) && isEqual(frictionCoefficient2Low,frictionCoefficient0 + A_ * log(dShearSlipDtAB / dShearSlipDtStar)))
    ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::SLOWSLIP_2A);
  else if (isEqual(dt, tHigh) && isEqual(frictionCoefficient2High,frictionCoefficient0 + A_ * log(dShearSlipDtEQ / dShearSlipDtStar)))
    ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::SLOWSLIP_2C);
  else
    ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::SLOWSLIP_2B);
}


///FUNCTION_BEGIN_PARSE
static_void
Dieterich::TransitionTimeSlowSlipC(const realT dEffectiveStressDt,
                                   const realT effectiveStress,
                                   const realT dStressShearDt,
                                   const realT stressShear,
                                   const realT frictionCoefficient2High,
                                   const realT frictionCoefficient2Low,
                                   const localIndex local,
                                   EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  //Check friction coefficient limits
  if (stressShear / effectiveStress < frictionCoefficient2Low)
  {
    ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::SLOWSLIP_2B);
    return;
  }

  const realT tLow = (frictionCoefficient2Low * effectiveStress - stressShear)
        / (dStressShearDt - frictionCoefficient2Low * dEffectiveStressDt);
  const realT tHigh = (frictionCoefficient2High * effectiveStress - stressShear)
      / (dStressShearDt - frictionCoefficient2High * dEffectiveStressDt);
  realT dt = std::numeric_limits<realT>::max();
  if (tLow > 0)
    dt = tLow;
  if (tHigh > 0 && tHigh < dt)
    dt = tHigh;

  if (isEqual(dt,tLow))
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::SLOWSLIP_2B);
  else
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::SLOWSLIP_2C);
}

///FUNCTION_BEGIN_PARSE
static_void
Dieterich::TransitionTimeCreep( const realT dEffectiveStressDt,
                                const realT effectiveStress,
                                const realT dStressShearDt,
                                const realT stressShear,
                                const realT frictionCoefficient3Low,
                                const realT frictionCoefficient3High,
                                const localIndex local,
                                EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  const realT tLow = (frictionCoefficient3Low * effectiveStress - stressShear)/
      (dStressShearDt - frictionCoefficient3Low * dEffectiveStressDt);
  const realT tHigh = (frictionCoefficient3High * effectiveStress - stressShear)/
      (dStressShearDt - frictionCoefficient3High * dEffectiveStressDt);
  realT dt = std::numeric_limits<realT>::max();
  if (tLow > 0)
    dt = tLow;
  if (tHigh > 0 && tHigh < dt)
    dt = tHigh;
  ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::CREEP);
}

///FUNCTION_BEGIN_PARSE
static_void
Dieterich::LowStressCheck( const realT dEffectiveStressDt,
                           const realT effectiveStress,
                           const realT stressNormalPin,
                           const realT dStressShearDt,
                           const realT stressShear,
                           const localIndex local,
                           EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  realT dt = std::numeric_limits<realT>::max();
  if (dEffectiveStressDt < 0)
  {
    dt = (stressNormalPin - effectiveStress)/dEffectiveStressDt;
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::LOW_SIGMA);
  }

  if (dStressShearDt < 0 && -stressShear/dStressShearDt < dt)
  {
    dt = -stressShear/dStressShearDt;
    ruptureTimestep.CheckForTimestepControl(dt, local, EarthquakeSimulation::LOW_TAU);
  }
}

///FUNCTION_BEGIN_PARSE
static_void
Dieterich::HighThetaCheck( const realT dt,
                           const realT theta,
                           const realT maxThetaPin,
                           const localIndex local,
                           EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  if (theta > (1e3*dt + maxThetaPin))
    ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::HIGH_THETA);
}

///FUNCTION_BEGIN_PARSE
static_bool
Dieterich::PorePressureChangeCheck( const Table<4, realT>* porePressure,
                                    EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  if(!porePressure)
    return false;

  realT dt = -1.0;
  {
    const std::vector<realT>& times = porePressure->AxisValues(3);
    for(std::vector<realT>::const_iterator it = times.begin(); it != times.end(); ++it)
    {
      if(*it > ruptureTimestep.m_currentTime)
      {
        dt = *it - ruptureTimestep.m_currentTime;
        break;
      }
    }
  }
  if(dt > 0.0)
    return ruptureTimestep.CheckForTimestepControl(
        dt, 0,
        EarthquakeSimulation::POREPRESSURERATECHANGE);
  else
    return false;
}

///FUNCTION_BEGIN_PARSE
static_realT
Dieterich::f(const realT dt,
             const realT* params)
{
  const realT dEffectiveStressDt = params[0];
  const realT effectiveStress = params[1];
  const realT dStressShearDt = params[2];
  const realT stressShear = params[3];
  const realT dShearSlipDtStar = params[4];
  const realT alpha = params[5];
  const realT A_ = params[6];
  const realT B = params[7];
  const realT Dc = params[8];
  const realT frictionCoefficient0 = params[9];
  const realT theta0 = params[10];


  const realT sigma = effectiveStress + dEffectiveStressDt * dt;
  const realT tau = stressShear + dStressShearDt * dt;
  const realT theta = ThetaLock(dt, dEffectiveStressDt, effectiveStress, alpha, B, theta0);
  //const realT thetaDot = 1.0 - (alpha/B)*(dEffectiveStressDt/effectiveStress)*theta;
/*  return( tau -  sigma*(mu0 + (B - A_)*log(theta) + A_*log(fabs(thetaDot))) );
    below seems to work as well on average and simplifies many things */
  return( tau -  sigma*(frictionCoefficient0 + (B - A_) * log(theta * dShearSlipDtStar / Dc)) );
}

//FUNCTION_BEGIN_PARSE
static_realT
Dieterich::zbrent(const realT* params,
                  const realT lower,
                  const realT upper,
                  const realT tol)
{
  return FindRoots::FindRoot(f, params, lower, upper, tol);
}


///FUNCTION_BEGIN_PARSE
static_void
Dieterich::APrioriFail(const realT dxsdtEQ,
                       const DieterichParameterData& matParams,
                       DieterichStateData& matState)
{
  const realT A_ = A(matParams, matState);
  const realT sigma = matState.stress;
  const realT mu0 = matState.mu;
  const realT tauRef = matState.stressShearReference;

  if (matParams.stressShearFail > 0.0) /* set shear stress to tauFail and then set theta for */
  {                   /* immediate failure */
    const realT dtau = matParams.stressShearFail - tauRef;
    const realT tau = dtau + matState.stressShearReference;
    matState.theta = (matParams.Dc/ matParams.vstar) * pow(dxsdtEQ / matParams.vstar, -A_/ matParams.B)*
                      exp((tau/sigma - mu0)/matParams.B);
  }
  else                /* set shear stress to level that should bring patch immediately */
  {                   /* to failure (i.e. a slip speed of ddotEQ) */
    matState.stressShear = sigma*(mu0 + A_ * log(dxsdtEQ / matParams.vstar) +
        matParams.B * log(matState.theta * matParams.vstar / matParams.Dc));
  }

  matState.dxsdt = dxsdtEQ; /* I think this will get overwritten if p is in state 0, but that it is needed if p is in state 1 */
  matState.apFail = EarthquakeSimulation::ALREADY;  /* indicate that this patch has already had its shear
                           stress set at its a priori failure time */
}

///FUNCTION_BEGIN_PARSE
realT
Dieterich::A(const localIndex iRupture) const
{
  const DieterichParameterData& matParams = *ParameterData(iRupture);
  const DieterichStateData& matState = *StateData(iRupture,0);
  return A(matParams, matState);
}

///FUNCTION_BEGIN_PARSE
static_realT
Dieterich::A(const DieterichParameterData& matParams,
             const DieterichStateData& matState)
{
  return matState.Areduced > 0 ?
      matParams.AreductionFactor * matParams.A :
      matParams.A;
}


//FUNCTION_BEGIN_PARSE
void
Dieterich::CheckA(const localIndex i)
{
  DieterichStateData& matState = *StateData(i,0);
  const DieterichParameterData& matParams = *ParameterData(i);
  CheckA(matParams, matState);
}

//FUNCTION_BEGIN_PARSE
void
Dieterich::CheckA(const DieterichParameterData& matParams,
                  DieterichStateData& matState) const
{
  const realT SDotExt = matState.dstressShearDt -
      ((matState.stressShear/matState.stress) -
          matParams.alpha)*(matState.dstressdt-matState.dppdt);

  if (matState.currentState == EarthquakeSimulation::NUCLEATE &&
      SDotExt > 0 &&
      matState.neighborInRuptureState > 0)
  {
    if (matState.Areduced == 0)
      matState.Areduced = 1;
  }
  else if (  matState.Areduced > 0)
  {
    matState.Areduced = 0;
    matState.theta = ThetaNucleate(
        matState.stressShear, matState.stress, matState.dxsdt,
        matParams.mu0,//FIXME: is this right?
        matParams.vstar,
        matParams.A,
        matParams.B,
        matParams.Dc);
  }
}

//FUNCTION_BEGIN_PARSE
bool
Dieterich::UpdateStressDriveRates(const globalIndex giRupture,
                                  const localIndex iCurrent,
                                  const realT ddotDrive,
                                  const EarthquakeSimulation::BoundaryElementDataManagerT& KGD)
{
  const R1Tensor& kgd = KGD(iCurrent, giRupture);
  m_stateData(iCurrent,0).dstressdtDrive += ddotDrive * kgd(0);
  m_stateData(iCurrent,0).dstressShearDtDrive += ddotDrive * kgd(1);
  return true;
}

//FUNCTION_BEGIN_PARSE
bool
Dieterich::UpdateStressRates(const globalIndex giRupture,
                             const localIndex iCurrent,
                             const realT dxsdt,
                             const EarthquakeSimulation::BoundaryElementDataManagerT& KGD)
{
  const R1Tensor& kgd = KGD(iCurrent, giRupture);
  m_stateData(iCurrent,0).dstressdt += dxsdt * kgd(0);
  m_stateData(iCurrent,0).dstressShearDt += dxsdt * kgd(1);
  return true;
}
