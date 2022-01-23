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
 * @file SeismicRiskSolver.h
 * @author Scott Johnson
 * @date created on March 9, 2013
 */

#ifndef SEISMICRISKSOLVER_H_
#define SEISMICRISKSOLVER_H_

#include "FaultRuptureBEMSolver.h"
#include "RiskSites.h"
#include "FaultPropertyVariabilityDieterich.h"

class SeismicRiskSolver: public FaultRuptureBEMSolver
{
public:
  SeismicRiskSolver( const std::string& name,
                     ProblemManagerT* const pm );

  ~SeismicRiskSolver(){};

  double
  TimeStep(const realT& time, const realT& dt,
           const int cycleNumber,
           PhysicalDomainT& domain,
           const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
           FractunatorBase* const fractunator);

  virtual void SetMaxStableTimeStep( PhysicalDomainT& domain,
                                     const sArray1d& namesOfSolverRegions ,
                                     SpatialPartition& partition );

  virtual void
  ReadXML( TICPP::HierarchicalDataNode* const hdn );

  static const char*
  SolverName() {return "SeismicRiskSolver";};

protected:
  void EpistemicLoop(PhysicalDomainT& domain);
  void AleatoricLoop(PhysicalDomainT& domain, const bool resetStress = false);
  realT RunUntilBurnIn(PhysicalDomainT& domain);
  bool TimeStepFaultRuptureBEM(realT& ct, PhysicalDomainT& domain, const bool steadyState = false);
  void BypassPhysicalSimulationAndReadCatalog();

  EarthquakeSimulation::RiskSites m_sites;
  EarthquakeSimulation::FaultPropertyVariabilityDieterich m_faultVariability;

  localIndex m_numberEventsWindowBurnInDetermination;
  realT m_minM, m_maxM, m_dM, m_grTol, m_timeBurnIn, m_timeAleatoric;
  std::string m_catalogFileName;
};

#endif /* SEISMICRISKSOLVER_H_ */

//1) Epistemic:
//The user defines the number of epistemic loops. There is currently no adaptivity.
//Also, each epistemic iteration is weighted equally in the finally hazard calculation.
//The epistemic loop functions to define several variable ranges to be used in RSQSim
//(e.g., A-value ranges, B-value ranges, etc.) as well as the fault trace to be interpreted
//by RSQSim to construct a fault mesh. All ranges are determined using a uniform random distribution.
//The epistemic loop also creates a per annum long-term seismicity recurrence relationship to use in
//weighting the results of the simulations; this is currently formed using a simple truncated exponential
//function parameterized by alpha, beta, and lambda (provided by the user).
//
//{
//A) Catalog Generation:
//Catalogs are generated for a user-defined number (kCalib) of magnitude bins; for each magnitude bin a
//number of (aleatorically varying according to a uniform random distribution) catalogs are generated.
//The user defines the number of aleatoric samples for all magnitude bins. The catalogs are generated over
//a user-specified length of time with a user-defined "burn-in" period required for RSQSim to reach stationary
//statistics. All simulations are performed before any processing for weighting. Therefore, the whole catalog is
//saved (indeed all of the element history is saved as well), though, later only the calibration bin in which you
//are interested is extracted from this dataset. This approach is quite unnecessarily memory intensive, since an
//array of dimension: (number_of_samples) x (number_of_calibration_bins) x (number_of_magnitude_bins) x (number_of_states)
//is saved for processing later.
//
//B) Weighting:
//To determine the catalog weights, the simple case of a single calibration bin would yield a weight for the epistemic iteration of:
//per annum long-term seismicity / average per annum recurrence for the given calibration as a composite of all aleatoric simulations
//
//For multiple calibration bins, a linear system of equations is solved to get the average weighting. Details are in WeightCat.f90.
//For the point source solution as well as the EMPSYN ground motion calculation in SIMRISK, the weights are all reset to 0.1, so be careful.
//
//C) Ground motion
//Ground motion is calculated for each event in each calibration bin's ensemble of cached event details (including element-wise history).
//}
//
//After all epistemic loops are finished, the equally-weighted results are then combined into a composite set of hazard curves
//(by percentile probability of exceedence).
