/*
 * FaultRuptureDieterich.h
 *
 *  Created on: Jun 6, 2012
 *      Author: johnson346
 */

#ifndef FAULTRUPTUREDIETERICH_H_
#define FAULTRUPTUREDIETERICH_H_

#include "Common/typedefs.h"
#include "Utilities/Utilities.h"
#include "Common/GPException.h"
#include "EarthquakeSimulationCommon.h"
#include "FaultRupture.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "IO/silo/SiloFile.h"
#include "BoundaryElementDataManagerT.h"
#include <string>
#include <iostream>

class Dieterich;

#if GPAC_MPI
#include "mpi.h"
#endif

namespace EarthquakeSimulation
{

  class FaultRuptureDieterich
  {
  public:

    FaultRuptureDieterich();
    virtual ~FaultRuptureDieterich();

    void ReadXML(TICPP::HierarchicalDataNode* hdn);

    void SetCurrentState(const TransitionState state);

//    inline void SetTransitionState(const int i, const realT ddot)
//    {
//      i_transition = i;
//      ddot_transition = ddot;
//    };

    bool Transition(Dieterich& m_manager,
                    const BoundaryElementDataManagerT& KG);

    void WriteSilo( SiloFile& siloFile );

    void ReadSilo( const SiloFile& siloFile );

  protected:
    bool TransitionLockedToNucleate(Dieterich& m_manager);
    bool TransitionNucleateToLocked(Dieterich& m_manager);
    void TransitionUpdateContributions(Dieterich& m_manager,
                                       const BoundaryElementDataManagerT& KG,
                                       const realT dxsdt,
                                       const int checkLevel);
    bool TransitionNucleateToRupture(Dieterich& m_manager,
                                     const BoundaryElementDataManagerT& KG);
    bool TransitionNucleateToSlowSlip2A(Dieterich& m_manager,
                                        const BoundaryElementDataManagerT& KG);
    bool TransitionRuptureToLocked(Dieterich& m_manager,
                                   const BoundaryElementDataManagerT& KG);
    bool TransitionSlowSlip2AToLocked(Dieterich& m_manager,
                                      const BoundaryElementDataManagerT& KG);
    bool TransitionSlowSlip2AToSlowSlip2B(Dieterich& m_manager,
                                          const BoundaryElementDataManagerT& KG);
    bool TransitionSlowSlip2BToSlowSlip2A(Dieterich& m_manager,
                                          const BoundaryElementDataManagerT& KG);
    bool TransitionSlowSlip2BToSlowSlip2B(Dieterich& m_manager,
                                          const BoundaryElementDataManagerT& KG);
    bool TransitionSlowSlip2BToSlowSlip2C(Dieterich& m_manager,
                                          const BoundaryElementDataManagerT& KG);
    bool TransitionSlowSlip2CToSlowSlip2B(Dieterich& m_manager,
                                          const BoundaryElementDataManagerT& KG);
    bool TransitionCreepToCreep(Dieterich& m_manager,
                                const BoundaryElementDataManagerT& KG);
    bool TransitionLowSigma(Dieterich& m_manager,
                            const BoundaryElementDataManagerT& KG,
                            const realT dxsdt);

    void SetCurrentState(Dieterich& m_manager, const TransitionState state);
    TransitionState CurrentState(const Dieterich& m_manager) const;

  public:
    void Reset();

    void ClearState();

  #if GPAC_MPI
    void SynchronizeState();
  #endif
    bool Synchronize(const realT globalDt,
                     const gArray1d& localToGlobalMap,
                     const std::map<globalIndex, localIndex>& globalToLocalMap);

    bool SkipTransition() const;
    bool TransitionLowSigma() const;
    int Transition() const;

    void AddRupture(const localIndex a,
                    const realT time,
                    const realT shear0);

    void CurrentRuptureFaultElementIndex(lArray1d& index) const { m_rupture.CurrentRuptureFaultElementIndex(index); }

  public:
    inline bool EarthquakeJustFinished() const { return m_isOn && m_nRuptures == 0; }

    realT FinalizeJustFinishedEarthquake(const realT G,
                                         const rArray1d& currentShear,
                                         const Array1dT<R1Tensor>& centers,
                                         const rArray1d& areas);

    void FinalizeJustFinishedEarthquake2();

    inline void EndRupture()
    {
      if(this->m_nRuptures > 0)
        --this->m_nRuptures;
    };
    inline bool IsQuiescent() const
    {
      return this->m_nRuptures == 0;
    };
    inline bool OnNeighbor() const
    {
      return this->m_onNeighbor;
    };

    //track the current rupture
    globalIndex m_global;
    realT m_ddot;
    TransitionState m_current;
    bool m_onNeighbor;
    bool m_isOn;
    bool m_isSerial;
    globalIndex m_nRuptures;

    realT dxsdtEQ, dMuCreep;

    EarthquakeSimulationTimestep m_timestep;
    FaultRupture m_rupture;
  };
}
//std::ostream& operator<<(std::ostream& out, const FaultRuptureDieterich& rupture);

#endif /* FAULTRUPTUREDIETERICH_H_ */
