/*
 * FaultRupture.cpp
 *
 *  Created on: Mar 9, 2013
 *      Author: johnson346
 */

#include "FaultRupture.h"

namespace EarthquakeSimulation
{
  GutenbergRichter::GutenbergRichter() : FrequencyMagnitude() {}

  GutenbergRichter::~GutenbergRichter(){}

  void
  GutenbergRichter::AddEvent(const realT magnitude)
  {
    realT current = m_minMagnitude;
    for(iArray1d::iterator it = m_frequencies.begin(); it != m_frequencies.end(); ++it, current+=m_binMagnitude)
    {
      if(magnitude >= current)
        ++(*it);
    }
    ++m_total;
  }

  realT
  GutenbergRichter::BValue(const int imin, const int imax, const realT duration) const
  {
    realT b_value;
    AValue(imin, imax, duration, b_value);
    return b_value;
  }

  realT
  GutenbergRichter::AValue(const int imin, const int imax, const realT duration, realT& b_value) const
  {
    //single linear regression line
    realT xy_sum = 0;
    realT x_sum  = 0;
    realT x2_sum = 0;
    realT y_sum  = 0;

    realT current = m_minMagnitude + imin * m_binMagnitude;
    int ii = 0;
    for (int i = imin; i <= imax; current += m_binMagnitude, ++i)
    {
      //calculate N per annum
      if(m_frequencies[i] > 0)
      {
        const realT nnorm = m_frequencies[i] / duration;
        const realT val = log10(nnorm);
        xy_sum += val * current;
        x_sum += current;
        y_sum += val;
        x2_sum += current * current;
        ++ii;
      }
    }
    if(ii == 0)
      throw GPException("GutenbergRichter::AValue : no frequencies to compute A-Value!");

    const realT inv_n = 1.0 / ii; //m_frequencies.size();
    b_value = -(xy_sum - (inv_n * x_sum * y_sum)) / (x2_sum - (inv_n * x_sum * x_sum));
    const realT a_value = inv_n * (y_sum + b_value * x_sum);

    return a_value;
  }

  bool
  GutenbergRichter::Compare(const GutenbergRichter& gr0, const int imin0, const int imax0, const realT duration0,
                            const GutenbergRichter& gr1, const int imin1, const int imax1, const realT duration1,
                            const realT tol, const bool avalue, const bool bvalue)
  {
    realT b0;
    const realT a0 = gr0.AValue(imin0, imax0, duration0, b0);

    realT b1;
    const realT a1 = gr1.AValue(imin1, imax1, duration1, b1);

    bool ret = true;
    if(avalue)
    {
      const realT diff = 2 * fabs(a0 - a1) / (a0 + a1);
      ret = ret && diff < tol;
      std::cout << "-->G-R comparison: t0=" << duration0 << " t1=" << duration1 << " a0=" << a0 << " a1=" << a1 << " diff=" << diff << " tol=" << tol << std::endl;
    }
    if(bvalue)
    {
      const realT diff = 2 * fabs(b0 - b1) / (b0 + b1);
      ret = ret && diff < tol;
      std::cout << "-->G-R comparison: t0=" << duration0 << " t1=" << duration1 << " b0=" << b0 << " b1=" << b1 << " diff=" << diff << " tol=" << tol << std::endl;
    }
    return ret;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Sub Data
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  FaultRuptureData::FaultRuptureData() : m_elementIndex(0), m_time(0), m_riseTime(0),
      m_slip(0), m_area(0), m_magnitude(0), m_moment(0), m_hypocenter(0)
  {
  }

  FaultRuptureData::~FaultRuptureData(){}

  void
  FaultRuptureData::Reset()
  {
    m_elementIndex = 0;
    m_time = 0;
    m_riseTime = 0;
    m_slip = 0;
    m_area = 0;
    m_magnitude = 0;
    m_moment = 0;
    m_hypocenter = 0;
  }

  void
  FaultRuptureData::Set(const localIndex elementIndex, const realT time, const realT riseTime,
                        const realT slip, const realT area, const realT magnitude,
                        const realT moment, const R1Tensor& hypocenter)
  {
    m_elementIndex = elementIndex;
    m_time = time;
    m_riseTime = riseTime;
    m_slip = slip;
    m_area = area;
    m_magnitude = magnitude;
    m_moment = moment;
    m_hypocenter = hypocenter;
  }

  void
  FaultRuptureData::WriteSilo( SiloFile& siloFile )
  {
    siloFile.DBWriteWrapper("m_elementIndex", m_elementIndex);
    siloFile.DBWriteWrapper("m_time", m_time);
    siloFile.DBWriteWrapper("m_riseTime", m_riseTime);
    siloFile.DBWriteWrapper("m_slip", m_slip);
    siloFile.DBWriteWrapper("m_area", m_area);
    siloFile.DBWriteWrapper("m_magnitude", m_magnitude);
    siloFile.DBWriteWrapper("m_moment", m_moment);
    siloFile.DBWriteWrapper("m_hypocenter", m_hypocenter);
  }

  void
  FaultRuptureData::ReadSilo( const SiloFile& siloFile )
  {
    siloFile.DBReadWrapper("m_elementIndex", m_elementIndex);
    siloFile.DBReadWrapper("m_time", m_time);
    siloFile.DBReadWrapper("m_riseTime", m_riseTime);
    siloFile.DBReadWrapper("m_slip", m_slip);
    siloFile.DBReadWrapper("m_area", m_area);
    siloFile.DBReadWrapper("m_magnitude", m_magnitude);
    siloFile.DBReadWrapper("m_moment", m_moment);
    siloFile.DBReadWrapper("m_hypocenter", m_hypocenter);
  }



  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Full Data
  /////////////////////////////////////////////////////////////////////////////////////////////////////////

  FaultRupture::FaultRupture() : m_main(), m_events() {}

  FaultRupture::FaultRupture( const FaultRupture& fr )
  {
    m_main = fr.m_main;
    m_events = fr.m_events;
  }

  FaultRupture::~FaultRupture(){}

  void
  FaultRupture::Reset()
  {
    m_main.Reset();
    m_events.clear();
  }

  void
  FaultRupture::AddRupture(const localIndex index, const realT time, const realT slip)
  {
    FaultRuptureData data;
    data.m_elementIndex = index;
    data.m_time = time;
    data.m_slip = slip;

    //Following are calculated during UpdateTotalMagnitude ...
    //  realT m_area;
    //  realT m_magnitude;
    //  realT m_moment;
    //  R1Tensor m_hypocenter;

    //TODO: what do we do about rise time?
    //realT m_riseTime;

    m_events.push_back(data);
  }

  void
  FaultRupture::Initialize(const localIndex index, const realT time, const realT slip)
  {
    Reset();
    m_main.m_elementIndex = index;
    m_main.m_time = time;
    m_main.m_slip = slip;
  }

  void
  FaultRupture::UpdateSubEvents(const realT time)
  {
    for(Array1dT<FaultRuptureData>::iterator it = m_events.begin(); it != m_events.end(); ++it)
    {
      FaultRuptureData& d = *it;
      d.m_magnitude = Magnitude(d.m_moment);
      d.m_riseTime = time - d.m_time;
    }
  }

  realT
  FaultRupture::UpdateTotalMagnitude(const realT G,
                                     const rArray1d& currentShear,
                                     const Array1dT<R1Tensor>& centers,
                                     const rArray1d& areas,
                                     const bool isSerial)
  {
    localIndex i = 0;
    m_main.m_moment = 0.0;
    m_main.m_area = 0.0;
    m_main.m_hypocenter = 0.0;
    R1Tensor tmp;
    for(Array1dT<FaultRuptureData>::iterator it = m_events.begin(); it != m_events.end(); ++it, ++i)
    {
      FaultRuptureData& d = *it;
      d.m_slip = currentShear[i] - d.m_slip;
      d.m_area = areas[d.m_elementIndex];
      d.m_hypocenter = centers[d.m_elementIndex];
      d.m_moment =  G * d.m_area * d.m_slip;
      m_main.m_moment += d.m_moment;
      m_main.m_area += d.m_area;
      tmp = d.m_hypocenter;
      tmp *= d.m_area;
      m_main.m_hypocenter += tmp;
    }
#if GPAC_MPI
    if(!isSerial)
    {
      realT Mlocal[] = {m_main.m_moment, m_main.m_area,
                        m_main.m_hypocenter(0),
                        m_main.m_hypocenter(1),
                        m_main.m_hypocenter(2)};
      realT Mglobal[] = {m_main.m_moment, m_main.m_area,
                         m_main.m_hypocenter(0),
                         m_main.m_hypocenter(1),
                         m_main.m_hypocenter(2)};
      MPI_Allreduce(&Mlocal, &Mglobal, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      m_main.m_moment = Mglobal[0];
      m_main.m_area = Mglobal[1];
      m_main.m_hypocenter(0) = Mglobal[2];
      m_main.m_hypocenter(1) = Mglobal[3];
      m_main.m_hypocenter(2) = Mglobal[4];
    }
#endif
    //NOTE: ASSUMES QUANTITIES ARE SI!!!
    m_main.m_hypocenter *= 1.0/m_main.m_area;
    m_main.m_magnitude = Magnitude(m_main.m_moment);
    m_main.m_slip = m_main.m_moment / (G * m_main.m_area);
    return m_main.m_magnitude;
  }

  void
  FaultRupture::WriteSilo(SiloFile& siloFile)
  {
    std::string s;
    std::string subDirectory =   "FaultRuptureInfo";
    if ( DBMkDir(siloFile.m_dbFilePtr, subDirectory.c_str()) != 0)
      throw GPException("Cannot make directory " + subDirectory);
    DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());

    //write main
    m_main.WriteSilo(siloFile);

    //write sub-events
    int nentries = m_events.size();
    int i = 0;
    siloFile.DBWriteWrapper("nentries", nentries);
    for(Array1dT<FaultRuptureData>::iterator it = m_events.begin(); it != m_events.end(); ++it, ++i)
    {
      std::stringstream ss;
      ss << i;
      s = ss.str();
      if ( DBMkDir(siloFile.m_dbFilePtr, s.c_str()) != 0)
        throw GPException("Cannot make directory for sub-events under FaultRuptureInfo");

      DBSetDir(siloFile.m_dbFilePtr, s.c_str());
      it->WriteSilo(siloFile);
      DBSetDir(siloFile.m_dbFilePtr, "..");
    }

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }

  void
  FaultRupture::ReadSilo(const SiloFile& siloFile)
  {
    std::string s;
    std::string subDirectory =   "FaultRuptureInfo";
    if ( DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str()) == 0)
    {
      //read main
      m_main.ReadSilo(siloFile);

      int nentries;
      siloFile.DBReadWrapper("nentries", nentries);
      m_events.clear();
      m_events.reserve(nentries);

      //read sub-events
      for(int i = 0; i < nentries; ++i)
      {
        std::stringstream ss;
        ss << i;
        s = ss.str();
        if ( DBSetDir(siloFile.m_dbFilePtr, s.c_str()) != 0)
          throw GPException("Cannot open expected sub-directory under FaultRuptureInfo: " + s);
        FaultRuptureData d;
        d.ReadSilo(siloFile);
        m_events.push_back(d);
        DBSetDir(siloFile.m_dbFilePtr, "..");
      }

      DBSetDir(siloFile.m_dbFilePtr, "..");
    }
  }
}
