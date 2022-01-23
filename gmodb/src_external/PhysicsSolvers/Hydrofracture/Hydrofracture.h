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
 * @file LagrangeDynamicsParallelPlateFlowExplicit.h
 * @author settgast1
 * @date Feb 28, 2011
 */

#ifndef HYDROFRACTURE_H_
#define HYDROFRACTURE_H_

#include "PhysicsSolvers/SolverBase.h"
#include "PhysicsSolvers/Lagrange/LagrangeSmallStrainLinearElastic.h"
#include "PhysicsSolvers/Flow/ParallelPlateFlowSolverFV.h"

// Explicit
#include "PhysicsSolvers/Lagrange/LagrangeLargeStrain.h"
#include "PhysicsSolvers/ParallelPlateFlowSolver.h"
#include "PhysicsSolvers/HydroStaticParallelPlateFlowSolver.h"
#include "PhysicsSolvers/FractureFlow/MatrixFlowSolver.h"

class ParallelPlateProppantSolver;

class Hydrofracture: public SolverBase
{
public:

  enum TimeIntegrationEnum
  {
    FirstTimeIntegrationEnum = 0,
    Implicit = 0,
    Explicit = 1,
    ImplicitHydroExplicitMechanics = 2,
    numTimeIntegrationEnums = 3
  };

  TimeIntegrationEnum IntToTimeIntegration( const int ival )
  {
    switch( ival )
    {
      case 0:
        return Implicit;
      case 1:
        return Explicit;
    }
    throw GPException("Invalid value input into Hydrofracture::IntToTimeIntegration()");
    return numTimeIntegrationEnums;
  }


  //int m_enableTimeIntegrationOption[numTimeIntegrationEnums];
  TimeIntegrationEnum m_timeIntegrationOption;


  Hydrofracture(  const std::string& name,
                  ProblemManagerT* const pm );

  virtual ~Hydrofracture();

  double TimeStep( const realT& time,
                 const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions,
                 SpatialPartition& partition,
                 FractunatorBase* const fractunator );

  void TimeStepImplicit( const realT& time,
                         const realT& dt,
                         const int cycleNumber,
                         PhysicalDomainT& domain,
                         const sArray1d& namesOfSolverRegions,
                         SpatialPartition& partition,
                         FractunatorBase* const fractunator );

  void TimeStepExplicit( const realT& time,
                         const realT& dt,
                         const int cycleNumber,
                         PhysicalDomainT& domain,
                         const sArray1d& namesOfSolverRegions,
                         SpatialPartition& partition,
                         FractunatorBase* const fractunator );


  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const sArray1d& namesOfSolverRegions);

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;

  virtual void RegisterFields( PhysicalDomainT& domain );

  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition );

  void InitializeCommunications( PartitionBase& partition );

  void SetupSystem(PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);
  realT Assemble(PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt);



  /// name of solver class
  /// NB: Name() function of parent class returns name of specific solver objects
  static const char* SolverName(){return "Hydrofracture";};

  virtual realT CheckSolutionBlock( const PhysicalDomainT& domain);


  virtual void PropagateSolutionBlock( const realT scalingFactor,
                                       PhysicalDomainT& domain);


  virtual void PostSyncConsistency( PhysicalDomainT& domain, SpatialPartition& partition );

  realT AssembleCouplingTerms( const PhysicalDomainT& domain,
                               const Array1dT< rArray1d >& dwdu,
                               const rArray1d& dwdw );

  realT SetGlobalMaxTimeStep(realT local_dt);


private:

  void CalculateContactStress(PhysicalDomainT& domain,
                              const realT& time,
                              const realT& dt,
                              localIndex& kf,
                              const localIndex faceIndex[],
                              realT& stressPen,
                              R1Tensor& stressShear);


  LagrangeSolverBase* m_ldSolve;
  ParallelPlateFlowSolverBase* m_ppSolve;
  MatrixFlowSolver* m_mfSolve;
  ParallelPlateProppantSolver* m_propSolve;

  std::string m_lgSolverName;
  std::string m_ppSolverName;
  std::string m_mfSolverName;
  std::string m_propSolverName;

  int m_deformationVariablesLocalRows;
  int m_deformationVariablesGlobalRows;

  int m_flowVariablesLocalRows;
  int m_flowVariablesGlobalRows;

  Epetra_System  m_epetraBlockSystem;

  const Epetra_MpiComm & epetra_comm;

  realT m_kJn, m_kJs;
  realT m_COFJ;
  realT m_jointCohesion;
  realT m_fLockedInSIF;
  realT m_faceStrengthRandomFactor;
//  bool m_COFJSetByUser;
//  bool m_jointCohesionSetByUser;
  bool m_COFJSetByInitialCondition;
  bool m_jointCohesionSetByInitialCondition;

  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> m_syncedFields;

  std::map<std::string,SolverBase*>& m_solvers;

  realT m_dtMatrixFlow;
  realT m_tNextMatrixFlowUpdate;

  // proppant
  realT m_oldPPStepTime;
  realT m_newPPStepTime;
  realT m_pp_dt;
  int m_nSubsteps;

};




#endif /* HYDROFRACTURE_H_ */
