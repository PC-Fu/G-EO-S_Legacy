/*
 * FaultRupture.h
 *
 *  Created on: Mar 9, 2013
 *      Author: johnson346
 */

#ifndef FAULTRUPTURE_H_
#define FAULTRUPTURE_H_

#include "Common/typedefs.h"
#include "Utilities/Utilities.h"
#include "Common/GPException.h"
#include "FrequencyMagnitude.h"
#include "IO/silo/SiloFile.h"
#include <string>
#include <iostream>

#ifdef GPAC_MPI
#include <mpi.h>
#endif

namespace EarthquakeSimulation
{
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // G-R Data
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class GutenbergRichter : public FrequencyMagnitude
  {
  public:
    GutenbergRichter();
    virtual ~GutenbergRichter();

    void AddEvent(const realT magnitude);

    realT BValue(const int imin, const int imax, const realT duration) const;

    realT AValue(const int imin, const int imax, const realT duration, realT& b_value) const;

    static bool Compare(const GutenbergRichter& gr0, const int imin0, const int imax0, const realT duration0,
                        const GutenbergRichter& gr1, const int imin1, const int imax1, const realT duration1,
                        const realT tol = 1e-3, const bool avalue = false, const bool bvalue = true);
  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Sub Data
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class FaultRuptureData
  {
  public:
    FaultRuptureData();
    virtual ~FaultRuptureData();

    void Reset();

    void Set(const localIndex elementIndex, const realT time, const realT riseTime,
             const realT slip, const realT area, const realT magnitude,
             const realT moment, const R1Tensor& hypocenter);

    void WriteSilo( SiloFile& siloFile );

    void ReadSilo( const SiloFile& siloFile );

    inline FaultRuptureData&
    operator=(const FaultRuptureData& rupture)
    {
      m_elementIndex = rupture.m_elementIndex;
      m_time = rupture.m_time;
      m_riseTime = rupture.m_riseTime;
      m_slip = rupture.m_slip;
      m_area = rupture.m_area;
      m_magnitude = rupture.m_magnitude;
      m_moment = rupture.m_moment;
      m_hypocenter = rupture.m_hypocenter;
      return *this;
    }

    localIndex m_elementIndex;
    realT m_time;
    realT m_riseTime;
    realT m_slip;
    realT m_area;
    realT m_magnitude;
    realT m_moment;
    R1Tensor m_hypocenter;
  };

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Full Data
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  class FaultRupture
  {
  public:
    FaultRupture();

    FaultRupture( const FaultRupture& fr );

    virtual ~FaultRupture();

    void Reset();

    void AddRupture(const localIndex index, const realT time, const realT slip);

    void Initialize(const localIndex index, const realT time, const realT slip);

    inline void CurrentRuptureFaultElementIndex(lArray1d& index) const
    {
      index.resize(m_events.size());
      localIndex i = 0;
      for(Array1dT<FaultRuptureData>::const_iterator it = m_events.begin(); it != m_events.end(); ++it, ++i)
        index[i] = it->m_elementIndex;
    }

    realT UpdateTotalMagnitude(const realT G,
                               const rArray1d& currentShear,
                               const Array1dT<R1Tensor>& centers,
                               const rArray1d& areas,
                               const bool isSerial = true);

    void UpdateSubEvents(const realT time);

    //in dyne*cm the eqn is the familiar (2/3)*log10(M0) - 10.73 -> converting to SI yields: (threshold to 1e-20 -> Mw=-20)
    static realT Magnitude(const realT moment) { return moment > 1e-20 ? (log10(moment) - 9.1)/1.5 : -20.0; }

    void WriteSilo( SiloFile& siloFile );

    void ReadSilo( const SiloFile& siloFile );

    inline FaultRupture&
    operator=(const FaultRupture& rupture)
    {
      m_main = rupture.m_main;
      m_events.clear();
      m_events = rupture.m_events;
      return *this;
    }

    FaultRuptureData m_main;
    Array1dT<FaultRuptureData> m_events;
  };
}

#endif /* FAULTRUPTURE_H_ */
