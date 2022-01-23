//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2013, Lawrence Livermore National Security, LLC.
//
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)     Stuart Walsh(walsh24@llnl.gov)
//  Scott Johnson (johnson346@llnl.gov)        Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)           
//
//  All rights reserved.
//
//  This file is part of GPAC.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file SteadyStateParallelPlateFlowSolver_TwoD.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef SSPARALLELPLATEFLOWSOLVER_MPI_TWOD_H_
#define SSPARALLELPLATEFLOWSOLVER_MPI_TWOD_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"
#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include <set>

#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"


class SteadyStateParallelPlateFlowSolver_TwoD : public SolverBase
{
public:
  SteadyStateParallelPlateFlowSolver_TwoD(const std::string& name,
                                          ProblemManagerT* const pm );

  virtual ~SteadyStateParallelPlateFlowSolver_TwoD();
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );


  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 
  void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain );
  void CalculateSteadyStatePressureDistribution(PhysicalDomainT& domain, SpatialPartition& partition,realT time );
  void UpdateEOS( const realT dt, PhysicalDomainT& domain );
  void UpdateFlux( const realT time, const realT dt, PhysicalDomainT& domain);
  
  // Pressure controlled boundary condition
  void PressureBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
         BoundaryConditionBase* bc, const lSet& set, realT time);
                                  
  void PressureBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
         BoundaryConditionBase* bc, const lSet& set, realT time);
                                                 
  // Net Flux boundary condition
  void FluxControlledBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
    BoundaryConditionBase* bc, const lSet& set,realT time);
                                     
  void FluxControlledBoundaryCondition_Setup(PhysicalDomainT& domain, int& localCount, const realT& time);       
  void FluxControlledBoundaryCondition_Count(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bc, const lSet& set,realT time);
                                                                          
  void FluxControlledBoundaryCondition_Sparsity(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
    BoundaryConditionBase* bc, const lSet& set,realT time);                            
  void FluxControlledBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
    BoundaryConditionBase* bc, const lSet& set,realT time); 
  

  void SetConcentrationField(std::string& concentrationFieldName,std::string reactionRateFieldName,std::string rrDerivFieldName);

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "SteadyStateParallelPlateFlowSolver_TwoD";};

private:


  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition);

  //void GetStiffnessMatrix(PhysicalDomainT& domain);
  void UpdateAperture(PhysicalDomainT&  domain);
  
  const lSet* m_faceSet;
  std::string m_faceSetName;
  localIndex m_numFaces;
  std::map<localIndex,localIndex> m_faceDofMap; 
  
  std::map<localIndex,lArray1d> m_edgesToFaces;
  
  realT m_tol; // Solver convergence criteria
  int m_maxSolverIters; // Maximum number of solver iterations
  
  // Flags
  bool doDataWrite;
  bool verboseFlag;
  bool updateApertureFlag;
  
  realT m_SHP_FCT;
  realT m_mu;
  realT m_bulk_modulus;
  realT m_min_aperture;
  realT m_max_aperture;
  
  // Boundary conditions
  unsigned int m_flux_bc_counter;
  iArray1d m_flux_bc_trilinos_host;
  lArray1d m_flux_bc_trilinos_index;
  rArray1d m_flux_bc_pressure;
  
  // MPI
  const int this_mpi_process;
  const int n_mpi_processes;
  
  #if GPAC_MPI
    const Epetra_MpiComm & m_epetra_comm;
  #else
    const Epetra_SerialComm & m_epetra_comm;
  #endif
  
  Teuchos::RCP<Epetra_Map>         row_map;
  Teuchos::RCP<Epetra_FECrsGraph>  sparsity;
  Teuchos::RCP<Epetra_FECrsMatrix> matrix;
  Teuchos::RCP<Epetra_FEVector>    solution;
  Teuchos::RCP<Epetra_FEVector>    rhs;
  
  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;
  
  std::string m_TrilinosIndexStr;
  static int m_instances;
  
  struct
  {
    double krylov_tol;
  }
  numerics;
  
};

int SteadyStateParallelPlateFlowSolver_TwoD::m_instances = 0;

#endif /* SSPARALLELPLATEFLOWSOLVER_H_ */
