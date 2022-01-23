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

#include "FaultRuptureDieterich.h"
#include "FaultRupture.h"
#include "Constitutive/Interface/Dieterich.h"

namespace EarthquakeSimulation
{
  FaultRuptureDieterich::FaultRuptureDieterich():
    m_global(std::numeric_limits<globalIndex>::max()),
    m_ddot(0.0),
    m_current(LOCK), m_onNeighbor(false),
    m_isOn(false), m_isSerial(true),
    m_nRuptures(0),
    dxsdtEQ(1.0), dMuCreep(0.0),
    m_timestep(), m_rupture()
  {
    m_timestep.Initialize();
  }

  FaultRuptureDieterich::~FaultRuptureDieterich()
  {
  }

  void FaultRuptureDieterich::ReadXML( TICPP::HierarchicalDataNode*  hdn)
  {
    dxsdtEQ = hdn->GetAttributeOrDefault<realT>("slipRateEQ",1.0);
    dMuCreep = hdn->GetAttributeOrDefault<realT>("creepFrictionRate",0.0);
  }

  void
  FaultRuptureDieterich::WriteSilo( SiloFile& siloFile )
  {
    std::string subDirectory =   "FaultRuptureDInfo";
    if ( DBMkDir(siloFile.m_dbFilePtr, subDirectory.c_str()) != 0)
      throw GPException("Cannot make directory " + subDirectory);
    DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());

    //data members
    siloFile.DBWriteWrapper("m_global", m_global );
    siloFile.DBWriteWrapper("m_ddot" , m_ddot );

    int tmp = m_current; siloFile.DBWriteWrapper("m_current", tmp);
    tmp = m_onNeighbor ? 1 : 0; siloFile.DBWriteWrapper("m_onNeighbor", tmp);
    tmp = m_isOn ? 1 : 0; siloFile.DBWriteWrapper("m_isOn", tmp);
    tmp = m_isSerial ? 1 : 0; siloFile.DBWriteWrapper("m_isSerial", tmp );

    siloFile.DBWriteWrapper("m_nRuptures", m_nRuptures);
    siloFile.DBWriteWrapper("dxsdtEQ", dxsdtEQ);
    siloFile.DBWriteWrapper("dMuCreep", dMuCreep);

    //m_timestep;
    siloFile.DBWriteWrapper("m_local", m_timestep.m_local);
    siloFile.DBWriteWrapper("m_dt", m_timestep.m_dt);
    siloFile.DBWriteWrapper("m_currentTime" , m_timestep.m_currentTime);
    tmp = m_timestep.m_next; siloFile.DBWriteWrapper("m_next", tmp );
    tmp = m_timestep.m_isPorePressureRateChange; siloFile.DBWriteWrapper("m_isPorePressureRateChange", tmp );

    //m_rupture
    m_rupture.WriteSilo(siloFile);

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }

  void
  FaultRuptureDieterich::ReadSilo( const SiloFile& siloFile )
  {
    std::string subDirectory =   "FaultRuptureDInfo";
    if ( DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str()) == 0)
    {
      //data members
      siloFile.DBReadWrapper("m_global", m_global );
      siloFile.DBReadWrapper("m_ddot" , m_ddot );

      int tmp; siloFile.DBReadWrapper("m_current", tmp);
      m_current = EarthquakeSimulation::TransitionStateFromIndex(tmp);
      siloFile.DBReadWrapper("m_onNeighbor", tmp); m_onNeighbor = tmp == 1;
      siloFile.DBReadWrapper("m_isOn", tmp); m_isOn = tmp == 1;
      siloFile.DBReadWrapper("m_isSerial", tmp ); m_isSerial = tmp == 1;

      siloFile.DBReadWrapper("m_nRuptures", m_nRuptures );
      siloFile.DBReadWrapper("dxsdtEQ", dxsdtEQ);
      siloFile.DBReadWrapper("dMuCreep", dMuCreep);

      //m_timestep;
      siloFile.DBReadWrapper("m_local", m_timestep.m_local);
      siloFile.DBReadWrapper("m_dt", m_timestep.m_dt);
      siloFile.DBReadWrapper("m_currentTime" , m_timestep.m_currentTime);
      siloFile.DBReadWrapper("m_next", tmp );
      m_timestep.m_next = EarthquakeSimulation::TransitionStateFromIndex(tmp);
      siloFile.DBReadWrapper("m_isPorePressureRateChange", tmp );
      m_timestep.m_isPorePressureRateChange = tmp == 1;

      //m_rupture
      m_rupture.ReadSilo(siloFile);

      DBSetDir(siloFile.m_dbFilePtr, "..");
    }
  }


  bool
  FaultRuptureDieterich::Transition(Dieterich& m_manager,
                                    const BoundaryElementDataManagerT& KG)
  {
    switch(Transition())
    {
      case LOCK_TO_NUCLEATE:
        return TransitionLockedToNucleate(m_manager);
      case NUCLEATE_TO_LOCK:
        return TransitionNucleateToLocked(m_manager);
      case NUCLEATE_TO_RUPTURE:
        return TransitionNucleateToRupture(m_manager, KG);
      case NUCLEATE_TO_SLOWSLIP_2A:
        return TransitionNucleateToSlowSlip2A(m_manager, KG);
      case RUPTURE_TO_LOCK:
        return TransitionRuptureToLocked(m_manager, KG);
      case CREEP_TO_CREEP:
        return TransitionCreepToCreep(m_manager, KG);
      case SLOWSLIP_2A_TO_LOCK:
        return TransitionSlowSlip2AToLocked(m_manager, KG);
      case SLOWSLIP_2A_TO_SLOWSLIP_2B:
        return TransitionSlowSlip2AToSlowSlip2B(m_manager, KG);
      case SLOWSLIP_2B_TO_SLOWSLIP_2A:
        return TransitionSlowSlip2BToSlowSlip2A(m_manager, KG);
      case SLOWSLIP_2B_TO_SLOWSLIP_2B:
        return TransitionSlowSlip2BToSlowSlip2B(m_manager, KG);
      case SLOWSLIP_2B_TO_SLOWSLIP_2C:
        return TransitionSlowSlip2BToSlowSlip2C(m_manager, KG);
      case SLOWSLIP_2C_TO_SLOWSLIP_2B:
        return TransitionSlowSlip2CToSlowSlip2B(m_manager, KG);
      default:
        return false;
    }
  }

  bool
  FaultRuptureDieterich::TransitionLockedToNucleate(Dieterich& m_manager)
  {
    if (!OnNeighbor())
    {
      SetCurrentState(m_manager, NUCLEATE);
      m_manager.TransitionLockedToNucleate(m_timestep.m_local, dxsdtEQ);
    }
    return false;
  }

  bool
  FaultRuptureDieterich::TransitionNucleateToLocked(Dieterich& m_manager)
  {
    if(!OnNeighbor())
      SetCurrentState(m_manager, LOCK);
    return false;
  }

  void
  FaultRuptureDieterich::TransitionUpdateContributions(Dieterich& m_manager,
                                                       const BoundaryElementDataManagerT& KG,
                                                       const realT dxsdt,
                                                       const int checkLevel)
  {
    m_manager.TransitionUpdateContributions(m_timestep.m_local, m_global, this->IsQuiescent(),
                                            KG, dxsdt, checkLevel);
  }

  bool
  FaultRuptureDieterich::TransitionNucleateToRupture(Dieterich& m_manager,
                                                     const BoundaryElementDataManagerT& KG)
  {
    /* first make changes to transitioning element */
    if(!OnNeighbor())
    {
      SetCurrentState(m_manager, RUPTURE);
      m_manager.TransitionNucleateToRupture(m_timestep.m_local, dxsdtEQ);

      const realT shearSlip0 = m_manager.StateData(m_timestep.m_local,0)->xs;
      AddRupture(m_timestep.m_local, m_timestep.m_currentTime, shearSlip0);
      m_ddot = m_manager.StateData(m_timestep.m_local,0)->dxsdt;
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, dxsdtEQ, 0);
    return false;
  }

  bool
  FaultRuptureDieterich::TransitionNucleateToSlowSlip2A(Dieterich& m_manager,
                                                        const BoundaryElementDataManagerT& KG)
  {
    //If this is mine, update
    if(!OnNeighbor())
    {
      SetCurrentState(m_manager, SLOWSLIP_2A);
      m_manager.TransitionNucleateToSlowSlip2A(m_timestep.m_local);

      const realT shearSlip0 = m_manager.StateData(m_timestep.m_local,0)->xs;
      AddRupture(m_timestep.m_local, m_timestep.m_currentTime, shearSlip0);
      m_ddot = m_manager.StateData(m_timestep.m_local,0)->dxsdt;
    }
  #if GPAC_MPI
    //Make sure all processes have the new matState.dxsdt
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, m_ddot, 3);
    return true;
  }

  bool
  FaultRuptureDieterich::TransitionRuptureToLocked(Dieterich& m_manager,
                                                   const BoundaryElementDataManagerT& KG)
  {
    if(!OnNeighbor())
    {
      SetCurrentState(m_manager, LOCK);
      EndRupture();
      m_ddot = m_manager.TransitionRuptureToLocked(m_timestep.m_local, dxsdtEQ);
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions( m_manager, KG, -dxsdtEQ, 1);
    return true;
  }

  bool
  FaultRuptureDieterich::TransitionSlowSlip2AToLocked(Dieterich& m_manager,
                                                      const BoundaryElementDataManagerT& KG)
  {
    if(!OnNeighbor())
    {
      SetCurrentState(m_manager, LOCK);
      EndRupture();
      m_ddot = m_manager.TransitionSlowSlip2AToLocked(m_timestep.m_local, dxsdtEQ);
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, m_ddot, 1);
    return true;
  }

  bool
  FaultRuptureDieterich::TransitionSlowSlip2AToSlowSlip2B(Dieterich& m_manager,
                                                          const BoundaryElementDataManagerT& KG)
  {
    if (!OnNeighbor())
    {
      SetCurrentState(m_manager, SLOWSLIP_2B);
      m_ddot = m_manager.TransitionSlowSlip2AToSlowSlip2B(m_timestep.m_local, dxsdtEQ, dMuCreep);
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, m_ddot, 2);
    return true;
  }

  bool
  FaultRuptureDieterich::TransitionSlowSlip2BToSlowSlip2B(Dieterich& m_manager,
                                                          const BoundaryElementDataManagerT& KG)
  {
    if (!OnNeighbor())
    {
      SetCurrentState(m_manager, SLOWSLIP_2B);
      m_ddot = m_manager.TransitionSlowSlip2BToSlowSlip2B(m_timestep.m_local, dxsdtEQ, dMuCreep);
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, m_ddot, 2);
    return true;
  }

  bool
  FaultRuptureDieterich::TransitionSlowSlip2BToSlowSlip2A(Dieterich& m_manager,
                                                          const BoundaryElementDataManagerT& KG)
  {
    if (!OnNeighbor())
    {
      SetCurrentState(m_manager, SLOWSLIP_2A);
      m_ddot = m_manager.TransitionSlowSlip2BToSlowSlip2A(m_timestep.m_local, dxsdtEQ, dMuCreep);
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, m_ddot, 2);
    return true;
  }

  bool
  FaultRuptureDieterich::TransitionSlowSlip2BToSlowSlip2C(Dieterich& m_manager,
                                                          const BoundaryElementDataManagerT& KG)
  {
    if (!OnNeighbor())
    {
      SetCurrentState(m_manager, SLOWSLIP_2C);
      m_ddot = m_manager.TransitionSlowSlip2BToSlowSlip2C(m_timestep.m_local, dxsdtEQ, dMuCreep);
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, m_ddot, 2);
    return true;
  }

  bool
  FaultRuptureDieterich::TransitionSlowSlip2CToSlowSlip2B(Dieterich& m_manager,
                                                          const BoundaryElementDataManagerT& KG)
  {
    if (!OnNeighbor())
    {
      SetCurrentState(m_manager, SLOWSLIP_2B);
      m_ddot = m_manager.TransitionSlowSlip2CToSlowSlip2B(m_timestep.m_local, dxsdtEQ, dMuCreep);
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, m_ddot, 2);
    return true;
  }

  bool
  FaultRuptureDieterich::TransitionCreepToCreep(Dieterich& m_manager,
                                                const BoundaryElementDataManagerT& KG)
  {
    if (!OnNeighbor())
    {
      SetCurrentState(m_manager, CREEP);
      m_ddot = m_manager.TransitionCreepToCreep(m_timestep.m_local, dxsdtEQ, dMuCreep);
    }
  #if GPAC_MPI
    SynchronizeState();
  #endif
    TransitionUpdateContributions(m_manager, KG, m_ddot, 2);
    return true;
  }

  /* for now, we are just locking such patches permanently */
  bool
  FaultRuptureDieterich::TransitionLowSigma(Dieterich& m_manager,
                                            const BoundaryElementDataManagerT& KG,
                                            const realT dxsdt)
  {
    if(OnNeighbor())
      return false;
    if (m_timestep.m_next == RUPTURE)
      m_manager.TransitionRuptureToLocked(m_timestep.m_local, dxsdt);
    else
      SetCurrentState(m_manager, LOCK);

    m_manager.TransitionLowSigma(m_timestep.m_local, m_global, KG);
    return true;
  }

  void
  FaultRuptureDieterich::SetCurrentState(Dieterich& m_manager,
                                         const TransitionState state)
  {
    m_manager.StateData(m_timestep.m_local, 0)->currentState = state;
    m_current = state;
  }

  TransitionState
  FaultRuptureDieterich::CurrentState(const Dieterich& m_manager) const
  {
    return TransitionStateFromIndex(
        m_manager.StateData(m_timestep.m_local, 0)->currentState);
  }

  void
  FaultRuptureDieterich::Reset()
  {
    m_global = std::numeric_limits<globalIndex>::max();
    m_ddot = 0.0;
    m_current = LOCK;
    m_onNeighbor = false;
    m_timestep.Reset();
  }

  void
  FaultRuptureDieterich::ClearState()
  {
    Reset();

    ////THESE ARE SET AT PROGRAM START ... DO NOT CLEAR
    //m_isSerial;
    //dxsdtEQ, dMuCreep;

    m_timestep.Initialize();

    m_isOn = false;
    m_nRuptures = 0;
    m_rupture.Reset();
  }

  #if GPAC_MPI
  void FaultRuptureDieterich::SynchronizeState()
  {
    if(m_isSerial)
      return;

    //TODO: make sure that SOMEONE takes ownership
    int mpisize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    int mpirank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    if (mpisize < 2)
    {
      return;
    }
    else
    {
      const size_t sz = sizeof(FaultRuptureDieterich);
      const int sendsz = mpisize - 1;
      if (!m_onNeighbor)
      {
        //std::cout << mpirank << ": sending to neighbors -- " << sendsz << std::endl;
        std::vector<MPI_Request> requests(sendsz);
        std::vector<MPI_Status> statuses(sendsz);
        Array1dT<Array1dT<char> > ddots(sendsz);
        localIndex ii = 0;
        for (int i = 0; i < mpisize; i++)
        {
          if(i==mpirank)
            continue;
          ddots[ii].resize(sz);
          memcpy(ddots[ii].data(), this, sz);
          //std::cout << mpirank << ": sending to ... " << i << std::endl;
          MPI_Isend(ddots[ii].data(), sz, MPI_CHAR, i, i, MPI_COMM_WORLD, &requests[ii]);
          ii++;
        }
        MPI_Waitall(sendsz, requests.data(), statuses.data());
      }
      else
      {
        //std::cout << mpirank << ": receiving from neighbors" << std::endl;

        Array1dT<char> cArray(sz);
        MPI_Status status;
        MPI_Recv(cArray.data(), sz, MPI_CHAR, MPI_ANY_SOURCE, mpirank, MPI_COMM_WORLD, &status);

        localIndex tmp = this->m_timestep.m_local;
        memcpy(this, cArray.data(), sz);

        //don't overwrite local index or onNeighbor
        this->m_timestep.m_local = tmp;
        this->m_onNeighbor = true;
      }
    }
    //std::cout << mpirank << ": done with SYNC" << std::endl;
  }
  #endif

  bool
  FaultRuptureDieterich::Synchronize(const realT globalDt,
                                     const gArray1d& localToGlobalMap,
                                     const std::map<globalIndex, localIndex>& globalToLocalMap)
  {
  #if GPAC_MPI
    if(!m_isSerial)
    {
      //NOTE: THIS ASSUMES THAT MIN(GLOBAL NUMBERS) >= 1
      this->m_global = isEqual(this->m_timestep.m_dt, globalDt) ? localToGlobalMap[this->m_timestep.m_local] : 0;
      {
        globalIndex localGlobal = this->m_global;
        MPI_Allreduce(&localGlobal, &this->m_global, 1,
                        MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
      }

      const localIndex* const a = stlMapLookupPointer(globalToLocalMap,
                                                      this->m_global);
      int mpisize;
      MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
      this->m_onNeighbor = !(mpisize < 2 || a);

      //synchronize the state of the ruptured element across processes
      this->SynchronizeState();
      return this->m_onNeighbor;
    }
  #endif
    // (1-2) Update the ruptured element
    this->m_onNeighbor = false;
    return this->m_onNeighbor;
  }

  bool FaultRuptureDieterich::SkipTransition() const
  {
    return this->m_timestep.m_isPorePressureRateChange;
  }

  bool FaultRuptureDieterich::TransitionLowSigma() const
  {
    return this->m_timestep.m_next == LOW_SIGMA || this->m_timestep.m_next == LOW_TAU || this->m_timestep.m_next == HIGH_THETA;
  }

  int FaultRuptureDieterich::Transition() const
  {
    if(SkipTransition())
      return ANY_TO_POREPRESSURERATECHANGE;
    switch(this->m_current)
    {
      case LOCK:
        if (this->m_timestep.m_next == NUCLEATE)
          return LOCK_TO_NUCLEATE;//0to1
        else
          throw GPException("MakeTransition(): illegal State0 transition: LOCK->" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        break;
      case NUCLEATE:
        switch(this->m_timestep.m_next)
        {
          case LOCK:
            return NUCLEATE_TO_LOCK; //1to0
          case RUPTURE:
            return NUCLEATE_TO_RUPTURE;//1to2
          case SLOWSLIP_2A:
            return NUCLEATE_TO_SLOWSLIP_2A;//1to2a
          case SLOWSLIP_2B:
          case SLOWSLIP_2C:
          case NUCLEATE:
          case POREPRESSURERATECHANGE:
          case LOW_TAU:
          case LOW_SIGMA:
          case HIGH_THETA:
          case CREEP:
          default:
            throw GPException("MakeTransition(): illegal State1 transition: NUCLEATE->" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        }
        break;
      case RUPTURE:
        if (this->m_timestep.m_next == LOCK)
          return RUPTURE_TO_LOCK;//2to0
        else
          throw GPException("MakeTransition(): illegal State2 transition: RUPTURE->" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        break;
      case CREEP:
        if (this->m_timestep.m_next == CREEP)
          return CREEP_TO_CREEP;//3to3
        else
          throw GPException("MakeTransition(): illegal State3 transition: CREEP->" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        break;
      case SLOWSLIP_2A:
        if (this->m_timestep.m_next == LOCK)
          return SLOWSLIP_2A_TO_LOCK;//2ato0
        else if (this->m_timestep.m_next == SLOWSLIP_2B)
          return SLOWSLIP_2A_TO_SLOWSLIP_2B;//2ato2b
        else
          throw GPException("MakeTransition(): illegal State2A transition: SLOWSLIP_2A->" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        break;
      case SLOWSLIP_2B:
        switch(this->m_timestep.m_next)
        {
          case SLOWSLIP_2A:
            return SLOWSLIP_2B_TO_SLOWSLIP_2A;//2bto2a
          case SLOWSLIP_2B:
            return SLOWSLIP_2B_TO_SLOWSLIP_2B;//2bto2b
          case SLOWSLIP_2C:
            return SLOWSLIP_2B_TO_SLOWSLIP_2C;//2bto2c
          case RUPTURE:
          case NUCLEATE:
          case LOCK:
          case POREPRESSURERATECHANGE:
          case LOW_TAU:
          case LOW_SIGMA:
          case HIGH_THETA:
          case CREEP:
          default:
            throw GPException("MakeTransition(): illegal State2B transition: SLOWSLIP_2B->" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        }
        break;
      case SLOWSLIP_2C:
        if (this->m_timestep.m_next == SLOWSLIP_2B)
          return SLOWSLIP_2C_TO_SLOWSLIP_2B;
        else
          throw GPException("MakeTransition(): illegal State2C transition: SLOWSLIP_2C->" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        break;
      case LOW_SIGMA:
      case LOW_TAU:
      case HIGH_THETA:
      case POREPRESSURERATECHANGE:
        throw GPException("MakeTransition(): illegal State transition - state not a transitional state: next=" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        break;
      default:
        throw GPException("MakeTransition(): illegal State transition - state not recognized: next=" + EarthquakeSimulation::TransitionStateToString(m_timestep.m_next));
        break;
    }
  }

  //std::ostream& operator<<(std::ostream& out, const FaultRuptureDieterich& rupture)
  //{
  //  out << " local/global: " << m_timestep.m_local << " / " << m_global << "\n dt = " << m_timestep.m_dt << "\n matState.dxsdt = " << m_ddot << "\n";
  //  //out << " current state: " << FaultRuptureDieterich::TransitionStateToString(current) << "\n";
  //  out << " next state: " << FaultRuptureDieterich::TransitionStateToString(m_timestep.m_next) << "\n";
  //  return out;
  //}

  void FaultRuptureDieterich::AddRupture(const localIndex a,
                                         const realT time,
                                         const realT shear0)
  {
    if(!m_isOn)
    {
      m_rupture.Initialize(a, time, shear0);
      m_isOn = true;
      m_nRuptures = 0;
    }
    m_rupture.AddRupture(a, time, shear0);
    ++m_nRuptures;
  }

  realT
  FaultRuptureDieterich::FinalizeJustFinishedEarthquake(const realT G,
                                                        const rArray1d& currentShear,
                                                        const Array1dT<R1Tensor>& centers,
                                                        const rArray1d& areas)
  {
    //TODO: calculate rise time and magnitude for each sub-event: use m_timestep.m_currentTime for rise
    const realT magnitude = m_rupture.UpdateTotalMagnitude(G, currentShear, centers, areas, m_isSerial);
    this->m_isOn = false;
    return magnitude;
  }

  void
  FaultRuptureDieterich::FinalizeJustFinishedEarthquake2()
  {
    m_rupture.UpdateSubEvents(m_timestep.m_currentTime);
  }
}
