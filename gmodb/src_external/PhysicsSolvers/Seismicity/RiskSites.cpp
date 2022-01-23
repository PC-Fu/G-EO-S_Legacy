/*
 * RiskSites.cpp
 *
 *  Created on: May 13, 2013
 *      Author: johnson346
 */

#include "RiskSites.h"
#include "Common/GPException.h"

namespace EarthquakeSimulation
{
  RiskSites::RiskSites() :
      m_gmMap(),
      m_gms(),
      m_hazardBins(),
      m_times(),
      m_siteHazardsFinal(),
      m_siteHazards(),
      //m_siteHazardCounts(),
      m_siteNames(),
      m_siteGroundMotionModelNames(),
      m_siteLocations(),
      m_siteWeights(),
      m_logAccelerationAnchor(0.0),
      m_nuisanceBeta(0.8)
  {

  }

  RiskSites::~RiskSites()
  {
  }

  void RiskSites::ReadXML(TICPP::HierarchicalDataNode* hdn, const realT timeAleatoric,
                          const realT timeLongTermSeismicity)
  {
    //set hazard bins
    const realT logmin = hdn->GetAttributeOrDefault<realT>("logmin", 0.0);
    const realT logmax = hdn->GetAttributeOrDefault<realT>("logmax", 0.0);
    const realT logsize = hdn->GetAttributeOrDefault<realT>("logsize", 0.0);

    //get time windows
    m_times.clear();
    {
      m_times.push_back(timeLongTermSeismicity);
      TICPP::HierarchicalDataNode* timeNode = hdn->GetChild("TimeWindows");
      if (timeNode)
      {
        for (TICPP::HierarchicalDataNode* child = timeNode->Next(true); child;
            child = timeNode->Next())
        {
          const realT time = timeLongTermSeismicity + child->GetAttributeValue<realT>("endTime");
          m_times.push_back(time);
        }
      }
      if (m_times.back() < (timeLongTermSeismicity + timeAleatoric))
        m_times.push_back(timeLongTermSeismicity + timeAleatoric);
    }

    //set the logarithmically ticked bins
    SetHazardBins(logmin, logmax, logsize);

    //set the nuisance variables
    //NOTE: must be put into units of log(cm/s/s)
    m_logAccelerationAnchor = log(
        100 * hdn->GetAttributeOrDefault<realT>("nuisanceAnchorAcceleration", 0.1));
    m_nuisanceBeta = hdn->GetAttributeOrDefault<realT>("nuisanceStandardDeviation", 0.8);

    //set sites
    for (TICPP::HierarchicalDataNode* child = hdn->Next(true); child; child = hdn->Next())
    {
      if (streq(child->Heading(), "Site"))
      {
        const std::string name = child->GetAttributeString("name");
        const std::string gmName = child->GetAttributeString("groundMotionModelName");
        const R1Tensor location = child->GetAttributeTensor("location");
        const realT weight = child->GetAttributeOrDefault<realT>("weight", 0.0);
        AddSite(name, gmName, location, weight);
      }
    }

    //set ground motion models
    TICPP::HierarchicalDataNode* GMNode = hdn->GetChild("GroundMotionModels");
    if (GMNode)
    {
      for (TICPP::HierarchicalDataNode* gmNode = GMNode->Next(true); gmNode;
          gmNode = GMNode->Next())
      {
        const std::string str = gmNode->GetAttributeString("name");
        if (str.length() == 0)
          throw GPException(
              "RiskSites::ReadXML : cannot define a ground motion model that has no name !");
        m_gmMap[gmNode->GetAttributeString("name")].ReadXML(gmNode);
      }
    }

    //validate that the sites have referenced valid ground motion models
    for (sArray1d::const_iterator it = m_siteGroundMotionModelNames.begin();
        it != m_siteGroundMotionModelNames.end(); ++it)
    {
      std::map<std::string, GroundMotion1D>::iterator itm = m_gmMap.find(*it);
      if (itm == m_gmMap.end())
        throw GPException(
            "RiskSites::ReadXML : could not find referenced ground motion model " + (*it));
      else
        m_gms.push_back(&itm->second);
    }

    //finish up
    Initialize();
  }

  void RiskSites::AddSite(const std::string& name, const std::string& gmName,
                          const R1Tensor& location, const realT weight)
  {
    if (m_hazardBins.size() == 0)
      throw GPException(
          "RiskSites: You must define the hazard bins before attempting to add a site");

    m_siteNames.push_back(name);
    m_siteGroundMotionModelNames.push_back(gmName);
    m_siteWeights.push_back(weight);
    m_siteLocations.push_back(location);

    m_siteHazards.resize2(m_siteNames.size(), m_times.size());
    m_siteHazardsFinal.resize2(m_siteNames.size(), m_times.size()); //don't need LTS
    //m_siteHazardCounts.resize2(m_siteNames.size(), m_times.size());

    const localIndex i = m_siteHazards.Dimension(0) - 1;
    for (localIndex j = 0; j < m_times.size(); j++)
    {
      m_siteHazards(i, j).resize(m_hazardBins.size(), 0.0);
      m_siteHazardsFinal(i, j).resize(m_hazardBins.size(), 0.0);
    }
  }

  void RiskSites::SetHazardBins(const realT minLogPeakAcceleration,
                                const realT maxLogPeakAcceleration,
                                const realT dLogPeakAcceleration)
  {
    const realT range = maxLogPeakAcceleration - minLogPeakAcceleration;
    if (range <= 0)
      throw GPException("RiskSites: Must have a positive range of log peak acceleration");
    if (dLogPeakAcceleration <= 0)
      throw GPException("RiskSites: Must have a positive bin size of log peak acceleration");
    int nbins = (int) (range / dLogPeakAcceleration);
    if (nbins == 0)
      ++nbins;
    const realT dbin = range / nbins;
    m_hazardBins.resize(nbins, 0.0);
    for (int i = 0; i < nbins; i++)
      m_hazardBins[i] = minLogPeakAcceleration + i * dbin;
  }

  void RiskSites::Initialize()
  {
    //normalize weights
    {
      realT isum = 0.0;
      for (rArray1d::const_iterator it = m_siteWeights.begin(); it != m_siteWeights.end(); ++it)
        isum += *it;
      if (isZero(isum))
        return;
      isum = 1.0 / isum;
      for (rArray1d::iterator it = m_siteWeights.begin(); it != m_siteWeights.end(); ++it)
        *it *= isum;
    }
  }

  void RiskSites::IncrementHazard(const realT time,
                                  const R1Tensor& location,
                                  const realT magnitude)
  {
    const localIndex j = TimeWindowIndex(time);
    if(j >= m_times.size())
      return; //time exceeds all windows

    R1Tensor dx;
    for (localIndex i = 0; i < m_siteNames.size(); ++i)
    {
      dx = location;
      dx -= m_siteLocations[i];
      const realT distance = dx.L2_Norm();
      m_gms[i]->IncrementHazard(distance, magnitude, m_hazardBins, m_siteHazards(i, j));
    }
  }

//  void RiskSites::UpdateGroundMotionHazard(const realT time, const R1Tensor& location,
//                                           const realT magnitude, const int nbins,
//                                           const int nsamples)
//  {
//    const localIndex j = TimeWindowIndex(time);
//
//    R1Tensor dx;
//    for (localIndex i = 0; i < m_siteNames.size(); ++i)
//    {
//      dx = location;
//      dx -= m_siteLocations[i];
//      const realT distance = dx.L2_Norm();
//      m_gms[i]->GroundMotionHazard(distance, magnitude, m_hazardBins, m_siteHazards(i, j),
//                                   m_siteHazardCounts(i, j), nbins, nsamples);
//    }
//  }

  void RiskSites::FinalizeAleatoricGroundMotionHazard()
  {
    const realT yrtos = 365.25 * 24 * 3600;

    for (localIndex i = 0; i < m_siteNames.size(); ++i)
    {
      //calculate the total hazard (time index 0 = LTS)
      realT tLast = 0;//m_times[0];
      for (localIndex j = 0; j < m_times.size(); j++)
      {
        const realT idt = yrtos / (m_times[j] - tLast);
        tLast = m_times[j];

        rArray1d& hbins = m_siteHazards(i, j);
        rArray1d& hbinsFinal = m_siteHazardsFinal(i, j);
        rArray1d::iterator ihf = hbinsFinal.begin();

        for (rArray1d::iterator ih = hbins.begin(); ih != hbins.end();
             ++ih, ++ihf)//++ihc0, ++ih0, ++ii)
        {
          //calculate the total per annum hazard for the site for the hazard bin
          realT& hazard = *ih;
          hazard *= idt;

          //increment the final hazard by the relative per annum hazard
          realT& hazardFinal = *ihf;
          hazardFinal += hazard;
        } //for each hazard bin
      } //for each time after LTS
    } //for each site

      //reset the values of the aleatoric hazards and counts
    for (localIndex i = 0; i < m_siteHazards.Dimension(0); ++i)
    {
      for (localIndex j = 0; j < m_siteHazards.Dimension(1); ++j)
      {
        m_siteHazards(i, j) = 0;
      }
    }
  }

  void RiskSites::RiskBins(rArray1d& riskBins) const
  {
    //get the evaluation of the Gaussian CDF for each value in the hazard bin
    riskBins.clear();
    riskBins.resize(m_hazardBins.size(), 0.0);
    std::copy(m_hazardBins.begin(), m_hazardBins.end(), riskBins.begin());
    for (rArray1d::iterator it = riskBins.begin(); it != riskBins.end(); ++it)
      *it = StatisticalDistributionBaseT::PhiNormal(
          (*it - m_logAccelerationAnchor) / m_nuisanceBeta);
  }

  void RiskSites::PrintGroundMotionHazardAndRisk(const realT weight) const
  {
    rArray2d hazardFinal;
    hazardFinal.resize2(m_siteHazardsFinal.Dimension(1), m_hazardBins.size());

    for (localIndex i = 0; i < m_siteHazardsFinal.Dimension(0); ++i)
    {
      const realT currWeight = weight * m_siteWeights[i];
      for (localIndex j = 0; j < m_siteHazardsFinal.Dimension(1); j++)
      {
        const rArray1d& haz = m_siteHazardsFinal(i, j);
        localIndex ii = 0;
        for (rArray1d::const_iterator it = haz.begin(); it != haz.end(); ++it, ++ii)
        {
          hazardFinal(j, ii) += currWeight * (*it);
        } //for each hazard bin
      } //for each time window
    } //for each site

    //PRINT THE HAZARD AND RISK
    rArray1d riskBins;
    RiskBins(riskBins);
    const realT stoyr = 1.0 / (365.25 * 24 * 3600);
//    realT tbegin = 0;
    for (localIndex i = 0; i < hazardFinal.Dimension(0); ++i)
    {
      const realT t0 = i > 0 ?
          (stoyr * (m_times[i-1] - m_times[0])) :
          (-stoyr * m_times[i]);
      const realT t1 = stoyr * (m_times[i] - m_times[0]);
      for (localIndex j = 0; j < hazardFinal.Dimension(1); j++)
      {
        std::cout << "HAZARD(year_start,year_end,log(hazard_x),log(risk_x),risk):"
            << " " << t0
            << " " << t1
            << " " << m_hazardBins[j]
            << " " << riskBins[j]
            << " " << hazardFinal(i, j)
            << std::endl;
      }
    }
  }

//  bool RiskSites::IsGroundMotionHazardFilled(const int nmin)
//  {
//    for (localIndex i = 0; i < m_siteHazardCounts.Dimension(0); i++)
//    {
//      for (localIndex j = 0; j < m_siteHazardCounts.Dimension(1); j++)
//      {
//        const iArray1d& itmp = m_siteHazardCounts(i, j);
//        for (iArray1d::const_iterator it = itmp.begin(); it != itmp.end(); ++it)
//          if (*it < nmin)
//            return false;
//      }
//    }
//    return true;
//  }
}
