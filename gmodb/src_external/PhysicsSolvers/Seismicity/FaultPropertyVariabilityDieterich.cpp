/*
 * FaultPropertyVariabilityDieterich.cpp
 *
 *  Created on: May 9, 2013
 *      Author: johnson346
 */

#include "FaultPropertyVariabilityDieterich.h"
#include "Constitutive/Interface/Dieterich.h"

namespace EarthquakeSimulation
{

  FaultPropertyVariabilityDieterich::FaultPropertyVariabilityDieterich() : m_filled(false),
      m_numberOfAleatoric(10), m_numberOfEpistemic(5), m_entries(), m_parameterOffsets(), m_stateOffsets()
  {
    sArray1d intNames, paramNames, stateNames, R1TensorNames, R2TensorNames, R2SymTensorNames;
    DieterichParameterData::GetVariableNames(intNames, paramNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    DieterichStateData::GetVariableNames(intNames, stateNames, R1TensorNames, R2TensorNames, R2SymTensorNames);

    for(sArray1d::const_iterator it = paramNames.begin(); it != paramNames.end(); ++it)
    {
      m_parameterOffsets[*it] = std::numeric_limits<size_t>::max();
    }
    for(sArray1d::const_iterator it = stateNames.begin(); it != stateNames.end(); ++it)
    {
      m_stateOffsets[*it] = std::numeric_limits<size_t>::max();
    }
  }

  FaultPropertyVariabilityDieterich::~FaultPropertyVariabilityDieterich()
  {
    // TODO Auto-generated destructor stub
  }

  void
  FaultPropertyVariabilityDieterich::ReadXML(TICPP::HierarchicalDataNode* hdn)
  {
    m_numberOfAleatoric = hdn->GetAttributeOrDefault<size_t>("numberOfAleatoricIterations", 10);
    m_numberOfEpistemic = hdn->GetAttributeOrDefault<size_t>("numberOfEpistemicIterations", 10);

    //Read variable properties
    for(TICPP::HierarchicalDataNode* varNode = hdn->Next(true); varNode; varNode = hdn->Next())
    {
      m_entries.push_back(FaultPropertyVariabilityDieterichData());
      FaultPropertyVariabilityDieterichData& datum = m_entries.back();
      datum.ReadXML(varNode);
      if(m_stateOffsets.find(datum.m_name) == m_stateOffsets.end())
      {
        datum.m_isParameter = m_parameterOffsets.find(datum.m_name) != m_parameterOffsets.end();
        if(!datum.m_isParameter)
          throw GPException("FaultPropertyVariabilityDieterich::ReadXML : Cannot add FaultPropertyVariabilityDieterichData with a name that is not in the parameters or states : " + datum.m_name);
      }
    }
  }

  void
  FaultPropertyVariabilityDieterich::NextEpistemic()
  {
    for(Array1dT<FaultPropertyVariabilityDieterichData>::iterator it = m_entries.begin(); it != m_entries.end(); ++it)
    {
      const realT* epistemicMin = it->m_epistemic.GetParameter(StatisticalDistributionBaseT::MINIMUM_VALUE);
      if(!epistemicMin)
        throw GPException("FaultPropertyVariabilityDieterich::NextEpistemic : Must define a minimum value for the range of values for " + it->m_name);

      const realT* epistemicMax = it->m_epistemic.GetParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE);
      if(!epistemicMax)
        throw GPException("FaultPropertyVariabilityDieterich::NextEpistemic : Must define a maximum value for the range of values for " + it->m_name);

      const realT aleatoricMean = StatisticalDistributionBaseT::UniformSample(*epistemicMin, *epistemicMax);
      const realT epistemicMean = (*epistemicMin) + 0.5 * ((*epistemicMax) - (*epistemicMin));

      realT halfRange = (aleatoricMean > epistemicMean) ? ((*epistemicMax) - aleatoricMean) : (aleatoricMean - (*epistemicMin));
      halfRange *= StatisticalDistributionBaseT::UniformSample(0.0, 1.0);

      const realT aleatoricMin = aleatoricMean - halfRange;
      const realT aleatoricMax = aleatoricMean + halfRange;

      it->m_aleatoric.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, aleatoricMin);
      it->m_aleatoric.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, aleatoricMax);
    }
  }

  void
  FaultPropertyVariabilityDieterich::NextAleatoric(Dieterich& contact)
  {
    //if offsets are not yet filled, then fill them
    FillOffsets(contact, true);

    for(Array1dT<FaultPropertyVariabilityDieterichData>::const_iterator it = m_entries.begin(); it != m_entries.end(); ++it)
    {
      if(!it->m_isParameter)
      {
        //STATES
        const std::map<std::string, size_t>::const_iterator istate = m_stateOffsets.find(it->m_name);
        if(istate == m_stateOffsets.end())
          throw GPException("Requested non-existent state!");
        const size_t sz = istate->second;
        for(localIndex i = 0; i < contact.NumStateIndex0(); i++)
        {
          const realT value = it->m_aleatoric.UniformSample();
          for(localIndex j = 0; j < contact.NumStateIndex1(); j++)
            memcpy(((char*)contact.StateData(i, j) + sz), &value, sizeof(realT));
        }
      }
      else
      {
        //PARAMETERS
        const std::map<std::string, size_t>::const_iterator iparam = m_parameterOffsets.find(it->m_name);
        if(iparam == m_parameterOffsets.end())
          throw GPException("Requested non-existent state!");
        const size_t sz = iparam->second;
        for(localIndex i = 0; i < contact.NumStateIndex0(); i++)
        {
          const realT value = it->m_aleatoric.UniformSample();
          memcpy(((char*)contact.ParameterData(i) + sz), &value, sizeof(realT));
        }
      }
    }
  }

  void FaultPropertyVariabilityDieterich::FillOffsets(Dieterich& contact, const bool forceRefill)
  {
    std::map<std::string, size_t> intOffsets;
    std::map<std::string, size_t> realOffsets;
    std::map<std::string, size_t> R1TensorOffsets;
    std::map<std::string, size_t> R2TensorOffsets;
    std::map<std::string, size_t> R2SymTensorOffsets;

    if(forceRefill)
      m_filled = false;

    if(!m_filled && contact.NumStateIndex0() > 0)
    {
      //STATES
      {
        DieterichStateData& state = *contact.StateData(0,0);
        state.GetVariableOffsets(intOffsets,
                                 realOffsets,
                                 R1TensorOffsets,
                                 R2TensorOffsets,
                                 R2SymTensorOffsets);
        for(std::map<std::string, size_t>::iterator it = m_stateOffsets.begin(); it != m_stateOffsets.end(); ++it)
        {
          std::map<std::string, size_t>::const_iterator ir0 = realOffsets.find(it->first);
          if(ir0 != realOffsets.end())
            it->second = ir0->second;
          else
            throw GPException("FaultPropertyVariabilityDieterich::FillOffsets : Cannot find state field " + it->first);
        }
      }

      intOffsets.clear();
      realOffsets.clear();
      R1TensorOffsets.clear();
      R2TensorOffsets.clear();
      R2SymTensorOffsets.clear();

      //PARAMETERS
      {
        DieterichParameterData& param = *contact.ParameterData(0);
        param.GetVariableOffsets(intOffsets,
                                 realOffsets,
                                 R1TensorOffsets,
                                 R2TensorOffsets,
                                 R2SymTensorOffsets);
        for(std::map<std::string, size_t>::iterator it = m_parameterOffsets.begin(); it != m_parameterOffsets.end(); ++it)
        {
          std::map<std::string, size_t>::const_iterator ir0 = realOffsets.find(it->first);
          if(ir0 != realOffsets.end())
            it->second = ir0->second;
          else
            throw GPException("FaultPropertyVariabilityDieterich::FillOffsets : Cannot find parameter field " + it->first);
        }
      }
      m_filled = true;
    }
  }
} /* namespace EarthquakeSimulation */
