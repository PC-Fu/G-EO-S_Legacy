//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2014, Lawrence Livermore National Security, LLC.
//
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//
//  Randolph Settgast		Stuart Walsh
//  Scott Johnson		Pengcheng Fu
//  Joshua White
//
//  LLNL-CODE-656690
//  GMOD-B, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GMOD-B. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
//  Please also read "Additional BSD Notice" below.
//
//  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the 
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Additional BSD Notice
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 * Dieterich.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"

#include "Dieterich.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "Constitutive/Interface/InterfaceFactory.h"
#include <typeinfo>
#include <assert.h>

Dieterich::Dieterich( ):
RateAndStateIntermediate( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

Dieterich::~Dieterich()
{

}

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

void
Dieterich::SetStress(const realT time, const R1Tensor& center, const Table<4, realT>* porePressure, DieterichStateData& matState)
{
  const realT xyzt[] = {center(0), center(1), center(2), time};
  matState.pp = porePressure ? porePressure->Lookup(xyzt) : 0;
  matState.dppdt = porePressure ? porePressure->Gradient(xyzt, 3) : 0;
  matState.stress = matState.stressReference + matState.dstress - matState.pp;
  matState.stressShear = matState.stressShearReference + matState.dstressShear;
}

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


void
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

realT
Dieterich::ShearStrength(const InterfaceBaseParameterData& ,
                         InterfaceBaseStateData& matStateBase) const
{
  DieterichStateData& matState = static_cast<DieterichStateData&> (matStateBase);
  return matState.stress * matState.mu;
}

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

realT
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

realT
Dieterich::Theta(const EarthquakeSimulation::TransitionState state,
                 const DieterichParameterData& matParams,
                 DieterichStateData& matState)
{
  return Theta(state, matState.dt, (matState.dstressdt-matState.dppdt), matState.stress, matParams.alpha, matParams.B, matParams.Dc, matState.dxs, matState.theta);
}

realT
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

realT
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

realT
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

realT
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

void
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

void
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



void
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

void
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

void
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

void
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


void
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

void
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

void
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

void
Dieterich::HighThetaCheck( const realT dt,
                           const realT theta,
                           const realT maxThetaPin,
                           const localIndex local,
                           EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep)
{
  if (theta > (1e3*dt + maxThetaPin))
    ruptureTimestep.CheckForTimestepControl(0.0, local, EarthquakeSimulation::HIGH_THETA);
}

bool
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

realT
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

realT
Dieterich::zbrent(const realT* params,
                  const realT lower,
                  const realT upper,
                  const realT tol)
{
  return FindRoots::FindRoot(f, params, lower, upper, tol);
}


void
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

realT
Dieterich::A(const localIndex iRupture) const
{
  const DieterichParameterData& matParams = *ParameterData(iRupture);
  const DieterichStateData& matState = *StateData(iRupture,0);
  return A(matParams, matState);
}

realT
Dieterich::A(const DieterichParameterData& matParams,
             const DieterichStateData& matState)
{
  return matState.Areduced > 0 ?
      matParams.AreductionFactor * matParams.A :
      matParams.A;
}


void
Dieterich::CheckA(const localIndex i)
{
  DieterichStateData& matState = *StateData(i,0);
  const DieterichParameterData& matParams = *ParameterData(i);
  CheckA(matParams, matState);
}

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

/// Register class in the class factory
REGISTER_INTERFACE( Dieterich )
