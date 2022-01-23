//FUNCTION_BEGIN_PARSE
realT
Dieterich::Advance(EarthquakeSimulation::EarthquakeSimulationTimestep& timestep,
                   const realT dtGlobal,
                   const realT dxsdtEQ,
                   const Array1dT<R1Tensor>& centers,
                   const Table<4, realT>* porePressure)
{
  for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
  {
    DieterichStateData& matState = *StateData(i,0);
    const DieterichParameterData& matParams = *ParameterData(i);

    const EarthquakeSimulation::TransitionState current = EarthquakeSimulation::TransitionStateFromIndex(
        matState.currentState);

    matState.dt = dtGlobal;
    matState.dxs = IncrementalShearDisplacement(current, dtGlobal,
                                                  dxsdtEQ, matParams, matState);
    matState.theta = Theta(current, matParams, matState);

    //advance integration
    matState.xs += matState.dxs;
    matState.dstress += matState.dstressdt * dtGlobal;
    matState.dstressShear += matState.dstressShearDt * dtGlobal + (
        current == EarthquakeSimulation::NUCLEATE ? matParams.KShearSelf * matState.dxs : 0.0);
    SetStress(timestep.m_currentTime, centers[i], porePressure, matState);
  }
  if(timestep.m_local < m_stateData.Dimension(0))
    return StateData(timestep.m_local, 0)->dxsdt;
  else
    return 0.0;
}

//FUNCTION_BEGIN_PARSE
void
Dieterich::SetStressRates(const realT time,
                          const EarthquakeSimulation::BoundaryElementDataManagerT& KGD,
                          const Table<4, realT>* porePressure,
                          const Array1dT<R1Tensor> centers,
                          const bool resetStress)
{
  if(KGD.Dimension(0) != m_stateData.Dimension(0))
    throw GPException("Unsafe condition: KGD length is not the same as the Dieterich state data!");
  for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
  {
    DieterichStateData& matState = *StateData(i,0);
    if(resetStress && matState.dxsdtDrive > 0.0)
    {
      matState.dstressShearDtDrive = 0.0;
      matState.dstressdtDrive = 0.0;
      const R1Tensor* const kgds = KGD[i];
      const Array1dT<R1Tensor>::size_type jMax = KGD.Dimension(1);
      for (Array1dT<R1Tensor>::size_type j = 0; j < jMax; ++j)
      {
        matState.dstressdtDrive -= kgds[j](0);
        matState.dstressShearDtDrive -= kgds[j](1);
      }
      matState.dstressShearDtDrive *= matState.dxsdtDrive;
      matState.dstressdtDrive *= matState.dxsdtDrive;
    }
    realT xyzt[] = {centers[i](0), centers[i](1), centers[i](2), time};
    matState.dstressdt = matState.dstressdtDrive;
    matState.dstressShearDt = matState.dstressShearDtDrive;
    matState.dppdt = porePressure ? porePressure->Gradient(xyzt, 3) : 0;
  }
}

//FUNCTION_BEGIN_PARSE
void
Dieterich::SetStress(const realT time, const R1Tensor& center, const Table<4, realT>* porePressure, DieterichStateData& matState)
{
  const realT xyzt[] = {center(0), center(1), center(2), time};
  matState.pp = porePressure ? porePressure->Lookup(xyzt) : 0;
  matState.dppdt = porePressure ? porePressure->Gradient(xyzt, 3) : 0;
  matState.stress = matState.stressReference + matState.dstress - matState.pp;
  matState.stressShear = matState.stressShearReference + matState.dstressShear;
}

//FUNCTION_BEGIN_PARSE
void
Dieterich::ResetStates(const gArray1d& localToGlobal,
                       const EarthquakeSimulation::BoundaryElementDataManagerT& KGD,
                       const realT time,
                       const realT dxsdtEQ,
                       const globalIndex maxGlobalIndex,
                       const Array1dT<R1Tensor>& centers,
                       const Table<4, realT>* porePressure)
{
  if(KGD.Dimension(0) != m_stateData.Dimension(0))
    throw GPException("Unsafe condition: KGD length is not the same as the Dieterich state data!");

  int creep = 0;

  // reset stressing rates and reset neighbors in rupture state
  for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
  {
    DieterichStateData& matState = *StateData(i,0);
    const DieterichParameterData& matParams = *ParameterData(i);
    SetStress(time, centers[i], porePressure, matState);
    matState.dstressShearDt = matState.dstressShearDtDrive;
    matState.dstressdt = matState.dstressdtDrive;
    matState.neighborInRuptureState = 0;
    if (matState.currentState == EarthquakeSimulation::CREEP)
    {
      matState.stressShear = matState.stressShearReference + matState.dstressShear;
      creep = 1;
      matState.dxsdt = matParams.vstar * exp((matState.stressShear/ matState.stress - matParams.mu0)/(A(i) - matParams.B));
      matState.dxsdt = matState.dxsdt < dxsdtEQ ? matState.dxsdt : dxsdtEQ;
    }
  }

#if GPAC_MPI
  {
    int tmp = creep;
    MPI_Allreduce(&tmp, &creep, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);
  }
#endif

  if (creep==1)
  {
    rArray1d::size_type asize = maxGlobalIndex + 1;
    rArray1d all(asize,0.0);
    iArray1d allStates(asize,0);
#if GPAC_MPI
    int mpisize; MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    if(mpisize > 1)
    {
      rArray1d tmp(asize,0.0);
      iArray1d tmpStates(asize,0.0);
      localIndex i = 0;
      for(gArray1d::const_iterator iter = localToGlobal.begin(); iter != localToGlobal.end(); ++iter, ++i)
      {
        const DieterichStateData& matState = *StateData(i,0);
        tmp[*iter] = matState.dxsdt;
        tmpStates[*iter] = matState.currentState;
      }
      MPI_Allreduce(tmp.data(), all.data(), all.size(),
                 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(tmpStates.data(), allStates.data(), allStates.size(),
                 MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
    }
    else
#endif
    {
      localIndex i = 0;
      for(gArray1d::const_iterator iter = localToGlobal.begin(); iter != localToGlobal.end(); ++iter, ++i)
      {
        const DieterichStateData& matState = *StateData(i,0);
        all[*iter] = matState.dxsdt;
        allStates[*iter] = matState.currentState;
      }
    }

    /* adjust tauDot and sigmaDot for state 3 creeping patches */
    {
      for (localIndex i = 0; i < m_stateData.Dimension(0); i++)
      {
        DieterichStateData& matState = *StateData(i,0);
        const R1Tensor* const kgds = KGD[i];
        const Array1dT<R1Tensor>::size_type jMax = KGD.Dimension(1);
        for (Array1dT<R1Tensor>::size_type j = 0; j < jMax; ++j)
        {
          if (allStates[j] == EarthquakeSimulation::CREEP)
          {
            matState.dstressdt += kgds[j](0) * all[j];
            matState.dstressShearDt += kgds[j](1) * all[j];
          }
        }
      }
    }
  }
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::NextTransitionTime( EarthquakeSimulation::EarthquakeSimulationTimestep& timestep, 
                               const realT dxsdtEQ,
                               const Array1dT<R1Tensor>& centers,
                               const Table<4, realT>* porePressure)
{
  //---------------------------------
  //Pre-Flight
  //TODO: keep this outside of the function
  //rupture.Reset();

  //---------------------------------
  //Flight
  timestep.Reset();
  for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
  {
    DieterichStateData& matState = *StateData(i,0);
    const DieterichParameterData& matParams = *ParameterData(i);
    SetStress(timestep.m_currentTime, centers[i], porePressure, matState);
    NextTransitionTime(matParams, matState, timestep, dxsdtEQ,
                       centers[i], i);
  }
  PorePressureChangeCheck(porePressure, timestep);
  return timestep.m_dt;
}

//FUNCTION_BEGIN_PARSE
bool
Dieterich::PostTimestepSynchronization( EarthquakeSimulation::EarthquakeSimulationTimestep& timestep,
                                        EarthquakeSimulation::TransitionState& current,
                                        realT& dxsdt)
{
  //---------------------------------
  //Post-Flight
  if(timestep.m_local < m_stateData.Dimension(0))
  {
    const localIndex i = timestep.m_local;
    DieterichStateData& matState = *StateData(i,0);
    const DieterichParameterData& matParams = *ParameterData(i);
    if (matParams.tFail >= 0 && timestep.m_currentTime < matParams.tFail)
    {
      if (matState.apFail == EarthquakeSimulation::NOTYET)
        matState.apFail = EarthquakeSimulation::NOW;
      matState.nextState = timestep.m_next;
      timestep.m_dt = 0;
    }
    //this is unconditionally necessary to handle the case of NO transitions
    //in such a case, the lowest global number with dt == std::numeric_limits<realT>::max()
    //will take ownership
    {
      current = EarthquakeSimulation::TransitionStateFromIndex(
          matState.currentState);
      dxsdt = matState.dxsdt;
    }
    return true;
  }

  //---------------------------------
  //Results
  return false;
}

//FUNCTION_BEGIN_PARSE
void
Dieterich::TransitionLockedToNucleate(const localIndex iRupture, const realT dxsdtEQ)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  /*  m->p[patch].matState.dxsdt = (p->Dc/p->theta)*(1 - (p->alpha/p->B)*(p->sigmaDot/p->matState.stress));
   now using the expression below which should give the same answer for a patch which
   actually transitioned 0 -> 1 and should also be right if it immediately transitions
   upon model startup. Note that these expressions for matState.dxsdt have to be consistent with
   what criteria are used for the 0->1 transition time, for example with the currently
   used approximation of matState.stressShear == matState.stress*(matParams.mu0 + (B-A)*ln(theta)), the above should really
   just be (matParams.Dc/theta) */

  realT fpow1 = matState.theta * matParams.vstar / matParams.Dc;
  if(isZero(fpow1))
  {
    matState.dxsdt = 0.0;
  }
  else
  {
    const realT A_ = A(matParams, matState);
    const realT fpow2 = -matParams.B / A_;
    const realT fexp = (matState.stressShear - matParams.mu0 * matState.stress) / (A_ * matState.stress);
    matState.dxsdt = matParams.vstar * exp(fexp);
    if( log(dxsdtEQ / matState.dxsdt) < (fpow2 * log(fpow1)) )
      matState.dxsdt = dxsdtEQ;
    else
      matState.dxsdt *= pow(fpow1, fpow2);
  }

  matState.H = matParams.B / matParams.Dc + matParams.KShearSelf / matState.stress;

  CheckA(matParams, matState); /* will reduce A if appropriate */
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionNucleateToRupture(const localIndex iRupture, const realT dxsdtEQ)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  matState.dstressShearDt += dxsdtEQ * matParams.KShearSelf;
  matState.dstressdt += dxsdtEQ * matParams.KNormalSelf;
  matState.dxsdt = dxsdtEQ;
  matState.theta = matParams.Dc / dxsdtEQ;
  CheckA(matParams, matState);
  const realT A_ = A(matParams, matState);
  matState.mu2 = (matParams.mu0 + (A_ - matParams.B) * log(
      dxsdtEQ / matParams.vstar)) * (1 + matParams.stressOvershootFactor) -
      matParams.stressOvershootFactor * matState.stressShear / matState.stress;
  return matState.dxs;
}

//FUNCTION_BEGIN_PARSE
void
Dieterich::TransitionUpdateContributions(const localIndex iRupture,
                                         const globalIndex giRupture,
                                         const bool isQuiescent,
                                         const EarthquakeSimulation::BoundaryElementDataManagerT& KGD,
                                         const realT dxsdt,
                                         const int checkLevel)
{
  switch (checkLevel)
  {
    case 0:
    {
      for (localIndex i = 0; i < m_stateData.Dimension(0); i++)
      {
        if (i == iRupture)
          continue;
        if(!UpdateStressRates(giRupture, i, dxsdt, KGD))
          continue;
        m_stateData(i,0).neighborInRuptureState = 1;
        CheckA(i);
      }
    }
      break;
    case 1:
    {
      for (localIndex i = 0; i < m_stateData.Dimension(0); i++)
      {
        if (i == iRupture)
          continue;
        if(!UpdateStressRates(giRupture, i, dxsdt, KGD))
          continue;
        if (isQuiescent)
        {
          m_stateData(i,0).neighborInRuptureState = 0;
          CheckA(i);
        }
      }
    }
      break;
    case 2:
    {
      for (localIndex i = 0; i < m_stateData.Dimension(0); i++)
        UpdateStressRates(giRupture, i, dxsdt, KGD);
    }
      break;
    case 3:
    {
      for (localIndex i = 0; i < m_stateData.Dimension(0); i++)
      {
        if(!UpdateStressRates(giRupture, i, dxsdt, KGD))
          continue;
        if (i == iRupture)
          continue;
        m_stateData(i,0).neighborInRuptureState = 1;
        CheckA(i);
      }
    }
      break;
    default:
      throw GPException("Cannot handle state");
  }
}

//FUNCTION_BEGIN_PARSE
void
Dieterich::TransitionNucleateToSlowSlip2A(const localIndex iRupture)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  CheckA(matParams, matState);
  matState.dxsdt = matParams.vstar;
  matState.theta = matParams.Dc / matState.dxsdt;
  matState.mu2aLow = matParams.mu0 *
        ( 1+ matParams.stressOvershootFactor ) -
        matParams.stressOvershootFactor * matState.stressShear / matState.stress;
  matState.mu2 = matParams.mu0 + A(matParams, matState) *
      log(matParams.dxsdtAB / matParams.vstar);
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionRuptureToLocked(const localIndex iRupture, const realT dxsdtEQ)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  matState.dstressShearDt -= dxsdtEQ * matParams.KShearSelf;
  matState.dstressdt -= dxsdtEQ * matParams.KNormalSelf;
  const realT m_ddot = matState.dxsdt;
  matState.dxsdt = 0;
  matState.theta = 1.0 / (dxsdtEQ / matParams.Dc + (matParams.alpha / matParams.B) *
        ((matState.dstressdt - matState.dppdt) / matState.stress));
  return m_ddot;
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionSlowSlip2AToLocked(const localIndex iRupture, const realT dxsdtEQ)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  matState.dstressShearDt -= dxsdtEQ * matParams.KShearSelf;
  matState.dstressdt -= dxsdtEQ * matParams.KNormalSelf;
  const realT m_ddot = matState.dxsdt;
  matState.dxsdt = 0;
  matState.theta = matParams.Dc / matParams.vstar;
  return m_ddot;
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionSlowSlip2AToSlowSlip2B(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  const realT A_ = A(matParams, matState);
  realT m_ddot = matState.dxsdt;
  matState.dxsdt = matParams.vstar * exp((matState.stressShear /
      matState.stress - matParams.mu0) / A_);
  matState.dxsdt = dxsdtEQ < matState.dxsdt ? dxsdtEQ : matState.dxsdt;

  matState.theta = matParams.Dc / matState.dxsdt;
  m_ddot = matState.dxsdt - m_ddot;
  matState.mu2Low = std::max(matState.stressShear / matState.stress -
                             dMuCreep, matParams.mu0 + A_ * log(matParams.dxsdtAB / matParams.vstar));
  matState.mu2High = std::min(matState.stressShear / matState.stress + dMuCreep,
                          matParams.mu0 + A_ * log(dxsdtEQ / matParams.vstar));
  matState.mu2aLow = std::min(matState.mu2aLow, matParams.mu0 -
                          matParams.stressOvershootFactor * (matState.stressShear /
                              matState.stress - matParams.mu0));
  return m_ddot;
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionSlowSlip2BToSlowSlip2B(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  const realT A_ = A(matParams, matState);
  realT m_ddot = matState.dxsdt;
  matState.dxsdt = matParams.vstar * exp((matState.stressShear / matState.stress - matParams.mu0) /
                                         A_);
  matState.dxsdt = dxsdtEQ < matState.dxsdt ? dxsdtEQ : matState.dxsdt;

  matState.theta = matParams.Dc / matState.dxsdt;
  m_ddot = matState.dxsdt - m_ddot;
  matState.mu2Low = std::max(matState.stressShear / matState.stress - dMuCreep, matParams.mu0 +
                             A_ * log((matParams.B - A_) / A_));
  matState.mu2High = std::min(matState.stressShear / matState.stress + dMuCreep,
                            matParams.mu0 + A_ * log(dxsdtEQ / matParams.vstar));
  matState.mu2aLow = std::min(matState.mu2aLow, matParams.mu0 - matParams.stressOvershootFactor * (matState.stressShear / matState.stress - matParams.mu0));
  return m_ddot;
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionSlowSlip2BToSlowSlip2A(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  const realT A_ = A(matParams, matState);
  realT m_ddot = matState.dxsdt;
  matState.dxsdt = matParams.vstar * exp( (matState.stressShear/matState.stress -
      matParams.mu0)/A_ );
  matState.theta = matParams.Dc / matState.dxsdt;

  m_ddot = matState.dxsdt - m_ddot;

  matState.mu2aLow = std::min(matState.mu2aLow, matParams.mu0 - matParams.stressOvershootFactor *
                              ( matState.stressShear / matState.stress - matParams.mu0));
  matState.mu2High = matParams.mu0 + A_ * log(matParams.dxsdtAB / matParams.vstar);
  return m_ddot;
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionSlowSlip2BToSlowSlip2C(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  realT m_ddot = matState.dxsdt;
  const realT A_ = A(matParams, matState);

  matState.dxsdt = dxsdtEQ;

  matState.theta = matParams.Dc / matState.dxsdt;
  m_ddot = matState.dxsdt - m_ddot;

  matState.mu2Low = matParams.mu0 + A_ * log(dxsdtEQ / matParams.vstar);
  matState.mu2High = std::numeric_limits<realT>::max();
  matState.mu2aLow = std::min(matState.mu2aLow, matParams.mu0 -
                                matParams.stressOvershootFactor * ( matState.stressShear /
                                    matState.stress - matParams.mu0));

  return m_ddot;
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionSlowSlip2CToSlowSlip2B(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  realT m_ddot = matState.dxsdt;
  const realT A_ = A(matParams, matState);
  matState.dxsdt = matParams.vstar * exp( (matState.stressShear/matState.stress - matParams.mu0)/A_ );
  matState.dxsdt = dxsdtEQ < matState.dxsdt ? dxsdtEQ : matState.dxsdt;

  matState.theta = matParams.Dc / matState.dxsdt;
  m_ddot = matState.dxsdt - m_ddot;

  matState.mu2Low = std::max(matState.stressShear / matState.stress - dMuCreep,
                         matParams.mu0 + A_ * log(matParams.dxsdtAB /
                                                                        matParams.vstar));
  matState.mu2High = std::min(matState.stressShear / matState.stress + dMuCreep,
                          matParams.mu0 + A_ * log(dxsdtEQ / matParams.vstar));
  matState.mu2aLow = std::min(matState.mu2aLow,
                          matParams.mu0 - matParams.stressOvershootFactor * ( matState.stressShear / matState.stress - matParams.mu0));
  return m_ddot;
}

//FUNCTION_BEGIN_PARSE
realT
Dieterich::TransitionCreepToCreep(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep)
{
  DieterichStateData& matState = *StateData(iRupture,0);
  const DieterichParameterData& matParams = *ParameterData(iRupture);

  realT m_ddot = matState.dxsdt;
  const realT A_ = A(matParams, matState);
  matState.dxsdt = matParams.vstar * exp((matState.stressShear / matState.stress - matParams.mu0) /
                                           (A_ - matParams.B));
  matState.dxsdt = dxsdtEQ < matState.dxsdt ? dxsdtEQ : matState.dxsdt;

  matState.theta = matParams.Dc / matState.dxsdt;
  m_ddot = matState.dxsdt - m_ddot;

  matState.mu3Low = matState.stressShear / matState.stress - dMuCreep;
  matState.mu3High = matState.stressShear / matState.stress + dMuCreep;

  return m_ddot;
}

//FUNCTION_BEGIN_PARSE
void
Dieterich::TransitionLowSigma(const localIndex iRupture,
                              const globalIndex giRupture,
                              const EarthquakeSimulation::BoundaryElementDataManagerT& KGD)
{
  if(StateData(iRupture,0)->dxsdtDrive > 0)
  {
    for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
    {
      DieterichStateData& matState = *StateData(i,0);
      matState.pinned = 1;
      //adjust backslip stress rates for the fact that this patch is now locked
      UpdateStressDriveRates(giRupture, i, StateData(iRupture,0)->dxsdtDrive, KGD);
    }
  }
  else
  {
    for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
    {
      StateData(i,0)->pinned = 1;
    }
  }
}


