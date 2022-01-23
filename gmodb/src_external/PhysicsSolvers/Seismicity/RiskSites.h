/*
 * RiskSites.h
 *
 *  Created on: Apr 8, 2013
 *      Author: johnson346
 */

#ifndef RISK_SITES_H_
#define RISK_SITES_H_

#include "Common/typedefs.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "SurfaceGeneration/StatisticalDistributionBaseT.h"

#include "IO/ticpp/HierarchicalDataNode.h"

#include "GroundMotion1D.h"

namespace EarthquakeSimulation
{
  class RiskSites
  {

  public:
    RiskSites();
    virtual ~RiskSites();

    void ReadXML(TICPP::HierarchicalDataNode* hdn, const realT timeAleatoric,
                 const realT timeLongTermSeismicity);

    void
    IncrementHazard(const realT time,
                    const R1Tensor& location,
                    const realT magnitude);

//    void
//    UpdateGroundMotionHazard(const realT time,
//                             const R1Tensor& location,
//                             const realT magnitude,
//                             const int nbins = 10,
//                             const int nsamples = 5);

    void FinalizeAleatoricGroundMotionHazard();
    void PrintGroundMotionHazardAndRisk(const realT weight) const;

//    bool IsGroundMotionHazardFilled(const int nmin);

    void RiskBins(rArray1d& riskBins) const;

    inline realT LongTermSeismicityTime() const
    {
      if(m_times.size() < 1)
        throw GPException("Have not defined LTS time");
      return m_times[0];
    }

  private:

    inline localIndex TimeWindowIndex(const realT time)
    {
      localIndex i = 0;
      for(i = 0; i < m_times.size(); i++)
        if(m_times[i] >= time)
          break;
      return i;
    }

    void AddSite(const std::string& name, const std::string& gmName, const R1Tensor& location, const realT weight);
    void SetHazardBins(const realT minLogPeakAcceleration,
                       const realT maxLogPeakAcceleration,
                       const realT dLogPeakAcceleration);
    void Initialize();

    std::map<std::string, GroundMotion1D> m_gmMap;
    Array1dT<GroundMotion1D*> m_gms;

    rArray1d m_hazardBins;
    rArray1d m_times;

    Array2dT<rArray1d> m_siteHazardsFinal;//dim of sites and times, respectively (inner array for each hazard bin)

    //site-wise aleatoric variables
    Array2dT<rArray1d> m_siteHazards;//dim of sites and times, respectively (inner array for each hazard bin)
    //Array2dT<iArray1d> m_siteHazardCounts;//dim of sites and times, respectively (inner array for each hazard bin)

    //site-wise parameters
    sArray1d m_siteNames, m_siteGroundMotionModelNames;
    Array1dT<R1Tensor> m_siteLocations;
    rArray1d m_siteWeights;

    //nuisance
    realT m_logAccelerationAnchor;
    realT m_nuisanceBeta;
  };
}

#endif /* RISK_SITES_H_ */
