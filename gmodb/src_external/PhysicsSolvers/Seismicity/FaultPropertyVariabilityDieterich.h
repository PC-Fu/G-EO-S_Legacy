/*
 * FaultPropertyVariabilityDieterich.h
 *
 *  Created on: May 9, 2013
 *      Author: johnson346
 */

#ifndef FAULTPROPERTYVARIABILITYDIETERICH_H_
#define FAULTPROPERTYVARIABILITYDIETERICH_H_

#include "SurfaceGeneration/StatisticalDistributionBaseT.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include <map>

class Dieterich;

namespace EarthquakeSimulation
{
  class FaultPropertyVariabilityDieterichData
  {
  public:
    std::string m_name;
    StatisticalDistributionBaseT m_epistemic, m_aleatoric;
    bool m_isParameter;

    void ReadXML(TICPP::HierarchicalDataNode* hdn)
    {
      //set hazard bins
      m_name = hdn->GetAttributeString("name");
      if(m_name.length() == 0)
        throw GPException("Must have a valid Dieterich parameter or state name");

      const realT emin = hdn->GetAttributeOrDefault<realT>("epistemicMin", 0.0);
      const realT emax = hdn->GetAttributeOrDefault<realT>("epistemicMax", 0.0);

      if((emax - emin) <= 0)
        throw GPException("FaultPropertyVariabilityDieterichData::ReadXML : Must have a positive definite range for the epistemic value of " + m_name);

      m_epistemic.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, emin);
      m_epistemic.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, emax);
    }
  };



  class FaultPropertyVariabilityDieterich
  {
  public:
    FaultPropertyVariabilityDieterich();
    virtual ~FaultPropertyVariabilityDieterich();

    void ReadXML(TICPP::HierarchicalDataNode* );

    void NextEpistemic();
    void NextAleatoric(Dieterich& contact);

    inline size_t NumberOfAleatoric() const { return m_numberOfAleatoric;}
    inline size_t NumberOfEpistemic() const { return m_numberOfEpistemic;}

  private:
    bool m_filled;
    size_t m_numberOfAleatoric, m_numberOfEpistemic;
    Array1dT<FaultPropertyVariabilityDieterichData> m_entries;
    std::map<std::string, size_t> m_parameterOffsets;
    std::map<std::string, size_t> m_stateOffsets;

    void FillOffsets(Dieterich& contact, const bool forceRefill = false);
  };

} /* namespace EarthquakeSimulation */
#endif /* FAULTPROPERTYVARIABILITYDIETERICH_H_ */
