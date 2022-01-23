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
/**
 * @file SeismicRiskSolver.cpp
 * @author Scott Johnson
 * @date created on March 9, 2013
 */

#include "SeismicRiskSolver.h"
#include "FaultRupture.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Utilities/Utilities.h"
#include "Utilities/Kinematics.h"
#include "Utilities/StringUtilities.h"

#include "Common/Common.h"

// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

#include "ObjectManagers/TableManager.h"
#include <algorithm>

#include "RiskSites.h"

#if GPAC_MPI
#include <mpi.h>
#endif

using namespace BoundaryConditionFunctions;

SeismicRiskSolver::SeismicRiskSolver(  const std::string& name,
                                       ProblemManagerT* const pm ):
FaultRuptureBEMSolver(name,pm), m_sites(), m_faultVariability(),
m_numberEventsWindowBurnInDetermination(1000),
m_minM(1.0), m_maxM(4.0), m_dM(0.3), m_grTol(1.0e-2), m_timeBurnIn(-1.0),
m_catalogFileName("")
{
}

void
SeismicRiskSolver::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  //domain.m_faultElementManager.ReadXML(hdn); //do this outside of this function: ProblemManagerT
  m_catalogFileName = hdn->GetAttributeString("catalogFileName");
  const bool bypassPhysicalSimulation = !m_catalogFileName.empty();
  if(!bypassPhysicalSimulation)
    FaultRuptureBEMSolver::ReadXML(hdn);

  m_dM = hdn->GetAttributeOrDefault<realT>("magnitudeBinSize", 0.3);
  m_minM = hdn->GetAttributeOrDefault<realT>("minimumMagnitude", 1.0);
  m_maxM = hdn->GetAttributeOrDefault<realT>("maximumMagnitude", 4.0);
  m_timeBurnIn = hdn->GetAttributeOrDefault<realT>("burnInTime", -1.0);

  m_timeAleatoric = hdn->GetAttributeOrDefault<realT>("aleatoricRunTime", -1.0);
  if(!bypassPhysicalSimulation && m_timeAleatoric <= 0)
  {
    throw GPException("SeismicRiskSolver::ReadXML : Must have a simulation time defined for the aleatoric loops");
  }
  m_grTol = hdn->GetAttributeOrDefault<realT>("gutenbergRichterComparisonTolerance", 1.0e-2);

  //sites
  {
    TICPP::HierarchicalDataNode* sitesNode = hdn->GetChild("RiskSites");
    if(sitesNode)
    {
      //get the long-term seismicity window, then offset
      const realT tltsDefault = 1000 * 365.25 * 24 * 3600;
      const realT dtlts = hdn->GetAttributeOrDefault<realT>("longTermSeismicityTimePeriod",
                                                            tltsDefault);
      m_sites.ReadXML(sitesNode, m_timeAleatoric, dtlts);
    }
  }

  //variability
  if(!bypassPhysicalSimulation)
  {
    TICPP::HierarchicalDataNode* varsNode = hdn->GetChild("FaultVariables");
    if(varsNode)
      m_faultVariability.ReadXML(varsNode);
  }

  //material
  {
    //LOOK IN FAULT ELEMENT MANAGER INITIALIZATION!!
  }
}

double SeismicRiskSolver::TimeStep(const realT& time,
                                 const realT& dt,
                                 const int cycleNumber,
                                 PhysicalDomainT& domain,
                                 const sArray1d& namesOfSolverRegions ,
                                 SpatialPartition& partition ,
                                 FractunatorBase* const fractunator)
{
  //todo: add parallelization

  ////////////////////////////////////////////////////////////////
  //     LOOP THROUGH ALL SIMULATIONS
  ////////////////////////////////////////////////////////////////
  const bool bypass = !m_catalogFileName.empty();
  if(bypass)
  {
    //you have an empirical catalog, so skip the sampling
    BypassPhysicalSimulationAndReadCatalog();
  }
  else
  {
    for(localIndex i = 0; i < m_faultVariability.NumberOfEpistemic(); i++)
    {
      std::cout << "EPISTEMIC LOOP: " << (i+1) << " of " << m_faultVariability.NumberOfEpistemic() << std::endl;
      EpistemicLoop(domain);
    }
  }

  ////////////////////////////////////////////////////////////////
  //     CLEAN UP
  ////////////////////////////////////////////////////////////////
  {
    const realT weight = bypass ? 1.0 : 1.0 / (m_faultVariability.NumberOfEpistemic() * m_faultVariability.NumberOfAleatoric());
    m_sites.PrintGroundMotionHazardAndRisk(weight);
  }
  return dt;
}

void SeismicRiskSolver::SetMaxStableTimeStep( PhysicalDomainT& domain,
                                              const sArray1d& namesOfSolverRegions ,
                                              SpatialPartition& partition )
{
  m_stabledt.m_maxdt = 0.9 * std::numeric_limits<realT>::max();
}

void SeismicRiskSolver::EpistemicLoop(PhysicalDomainT& domain)
{
  //todo: add changes to fault geometry
  m_faultVariability.NextEpistemic();
  for(localIndex i = 0; i < m_faultVariability.NumberOfAleatoric(); i++)
  {
    std::cout << "ALEATORIC LOOP: " << (i+1) << " of " << m_faultVariability.NumberOfAleatoric() << std::endl;
    AleatoricLoop(domain, i==0);
  }
}

bool SeismicRiskSolver::TimeStepFaultRuptureBEM(realT& ct, PhysicalDomainT& domain, const bool steadyState)
{
  m_stabledt.m_maxdt = domain.m_faultElementManager.CalculateTimestep(steadyState);
  const bool endOfEQ = domain.m_faultElementManager.TimeStep(m_stabledt.m_maxdt, steadyState);
  ct += m_stabledt.m_maxdt;
  return endOfEQ;
}

void SeismicRiskSolver::BypassPhysicalSimulationAndReadCatalog()
{
  R1Tensor x;
  realT m, t;
  std::ifstream input;
  input.open(m_catalogFileName.c_str());
  const bool exist = input.is_open();
  if(!exist)
    throw GPException("Could not open seismic catalog file!");
  while(!input.eof())
  {
    input >> t >> m >> x(0) >> x(1) >> x(2);
    m_sites.IncrementHazard(t, x, m);
  }
  m_sites.FinalizeAleatoricGroundMotionHazard();
  input.close();
}

void SeismicRiskSolver::AleatoricLoop(PhysicalDomainT& domain, const bool resetStress)
{
  ////////////////////////////////////////////////////////////////
  //     INITIALIZE
  ////////////////////////////////////////////////////////////////
  domain.m_faultElementManager.SetPorePressureSteadyState();

  domain.m_faultElementManager.ResetStatesAndParameters(resetStress);

  //get the next aleatoric realization of fault parameters ... override any
  //initial or boundary conditions (which were set in ResetStatesAndParameters)
  //NOTE: the aleatoric values WILL be overridden by any _boundary_ conditions
  //during successive timesteps.
  m_faultVariability.NextAleatoric(domain.m_faultElementManager.DieterichObject());

  //run until burn-in
  RunUntilBurnIn(domain);

  //run post-burn-in; only update hazard with events after burn-in
  realT ct = 0.0;
  const realT ctEnd = m_timeAleatoric + m_sites.LongTermSeismicityTime();
  bool keepGoing = true;
  while(keepGoing)
  {
    //determine whether to reset the time stamp in the BEM solver
    //based on whether you have reached the LTS limit
    //i.e., the pore pressure evaluation should see the reset time
    const bool steadyState = ct < m_sites.LongTermSeismicityTime();

    ////////////////////////////////////////////////////////////////
    //     TIMESTEP
    ////////////////////////////////////////////////////////////////
    const bool endOfEQ = TimeStepFaultRuptureBEM(ct, domain, steadyState);
    if(endOfEQ)
    {
      //todo: determine whether to do the ground motion calculation based on whether
      //you've reached your magnitude quota for this time window
      m_sites.IncrementHazard(ct,
                              domain.m_faultElementManager.LastEvent().m_hypocenter,
                              domain.m_faultElementManager.LastEvent().m_magnitude);
    }
    if(ct > ctEnd)
      break;
  }

  ////////////////////////////////////////////////////////////////
  //     CLEAN UP
  ////////////////////////////////////////////////////////////////
  m_sites.FinalizeAleatoricGroundMotionHazard();
}

realT SeismicRiskSolver::RunUntilBurnIn(PhysicalDomainT& domain)
{
  if(m_timeBurnIn < 0.0 && m_numberEventsWindowBurnInDetermination <= 20)
  {
    throw GPException(" SeismicRiskSolver::RunUntilBurnIn : burn-in time not specified and number of events to use in burn-in auto-determination are too few!");
  }

  //todo: set G-R parameters from UI
  EarthquakeSimulation::GutenbergRichter gr0, gr1;
  size_t n0 = gr0.SetRange(2.0, 5.0, 0.25);
  gr1.SetRange(2.0, 5.0, 0.25);
  int imin0 = 0, imax0 = n0 - 1;
  int imin1 = imin0, imax1 = imax0;

  //define states
  bool filled0 = false, filled1 = false, set = false;
  realT ct = 0.0;
  realT begin0 = 0.0, end0 = 0.0, begin1 = 0.0, end1 = 0.0;
  realT m_min = std::numeric_limits<realT>::max();
  realT m_max = -m_min;

  //loop until burn-in is assessed true
  while(true)
  {
    /////////////////////////
    //timestep
    /////////////////////////

    const bool endOfEQ = TimeStepFaultRuptureBEM(ct, domain, true);

    /////////////////////////
    //check burn-in condition
    /////////////////////////

    //If we've set a hard burn-in time limit, use that ...
    if(m_timeBurnIn > 0.0)
    {
      if(ct > m_timeBurnIn)
        return ct;
    }
    //otherwise, check whether we've reached burn-in by comparing consecutive windows
    else if(endOfEQ)
    {
      const realT m = domain.m_faultElementManager.LastMagnitude();
      if(!filled0)
      {
        gr0.AddEvent(m);
        filled0 = gr0.TotalEvents() >= m_numberEventsWindowBurnInDetermination;
        end0 = ct;

        //during the first pass, gather the range of the magnitudes
        //use these to refine the G-R range to compare
        if(!set)
        {
          if(filled0)
          {
            m_min *= m_min < 0 ? 1.1 : 0.9;
            m_max *= m_max < 0 ? 0.9 : 1.1;
            const realT dm = 0.05 * (m_max - m_min);
            n0 = gr0.SetRange(m_min, m_max, dm);
            gr1.SetRange(m_min, m_max, dm);
            imax0 = n0 - 1;
            imax1 = imax0;
            set = true;
            filled0 = false;
            std::cout << "-->burn-in: m_min=" << m_min << " m_max=" << m_max << " dm=" << dm << " n=" << n0 << std::endl;
          }
          else
          {
            m_min = m < m_min ? m : m_min;
            m_max = m > m_max ? m : m_max;
          }
        }


      }
      else
      {
        if(filled0 && filled1)
        {
          const realT duration0 = end0 - begin0;
          const realT duration1 = end1 - begin1;
          if(EarthquakeSimulation::GutenbergRichter::Compare(gr0, imin0, imax0, duration0,
                                                             gr1, imin1, imax1, duration1,
                                                             m_grTol))
          {
            const realT stoyr = 1.0 / (365.25 * 24 * 3600);
            std::cout << "-->burn-in converged at: " << (stoyr * ct) << " years (" << ct << " seconds)" << std::endl;
            return ct;
          }
          gr0 = gr1;
          begin0 = begin1;
          end0 = end1;
          gr1.ResetFrequencies();
          begin1 = ct;
        }
        gr1.AddEvent(m);
        filled1 = gr1.TotalEvents() >= m_numberEventsWindowBurnInDetermination;
        end1 = ct;
      }
    }
  }
  return -1.0;
}

/// Register solver in the solver factory
REGISTER_SOLVER( SeismicRiskSolver )
