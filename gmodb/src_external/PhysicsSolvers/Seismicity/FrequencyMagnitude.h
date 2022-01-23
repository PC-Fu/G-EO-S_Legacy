/*
 * FrequencyMagnitude.h
 *
 *  Created on: Mar 9, 2013
 *      Author: johnson346
 */

#ifndef FREQUENCY_MAGNITUDE_H_
#define FREQUENCY_MAGNITUDE_H_

#include "Common/typedefs.h"
#include "Utilities/Utilities.h"
#include "Common/GPException.h"

namespace EarthquakeSimulation
{
  class FrequencyMagnitude
  {
  public:
    FrequencyMagnitude() : m_range(3.0), m_minMagnitude(1.0), m_binMagnitude(0.3), m_frequencies(), m_total(0) {}

    ~FrequencyMagnitude(){}

    size_t SetRange(const realT min, const realT max, const realT dbin)
    {
      m_minMagnitude = min;
      m_range = max - min;
      m_binMagnitude = dbin;

      int num = m_range / m_binMagnitude;
      if((m_range - (num * m_binMagnitude)) > (0.9 * m_binMagnitude))
        ++num;
      m_binMagnitude = m_range / num;
      m_total = 0;
      m_frequencies.clear();
      m_frequencies.resize(num+1, 0);
      return m_frequencies.size();
    }

    inline localIndex BinIndex(const realT magnitude) const
    {
      const realT current = (magnitude - m_minMagnitude) / m_frequencies.size();
      localIndex ret = (localIndex) current;
      ++ret;
      if(ret >= m_frequencies.size())
        --ret;
      return ret;
    }

    void ResetFrequencies()
    {
      m_total = 0;
      for(iArray1d::iterator it = m_frequencies.begin(); it != m_frequencies.end(); ++it)
        *it = 0;
    }

    inline localIndex TotalEvents() const { return m_total; }

    inline FrequencyMagnitude&
    operator=(const FrequencyMagnitude& rupture)
    {
      m_range = rupture.m_range;
      m_minMagnitude = rupture.m_minMagnitude;
      m_binMagnitude = rupture.m_binMagnitude;

      m_frequencies.clear();
      m_frequencies = rupture.m_frequencies;
      return *this;
    }

    realT m_range;
    realT m_minMagnitude;
    realT m_binMagnitude;
    iArray1d m_frequencies;
    localIndex m_total;
  };
}

#endif /* FREQUENCY_MAGNITUDE_H_ */
