//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2013, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)     Stuart Walsh(walsh24@llnl.gov)
//  Scott Johnson (johnson346@llnl.gov)        Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//
//  LLNL-CODE-6182322
//  GPAC, Version 2.0
//
//  All rights reserved.
//  This file is part of GPAC.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file ParallelPlateFlowSolverSteadyState.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef SSPARALLELPLATEFLOWSOLVER_H_
#define SSPARALLELPLATEFLOWSOLVER_H_

#include "PhysicsSolvers/ParallelPlateFlowSolverBase.h"
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


class ParallelPlateFlowSolverSteadyState : public ParallelPlateFlowSolverBase
{
public:
  ParallelPlateFlowSolverSteadyState( const std::string& name , ProblemManagerT* const problemManager);
  virtual ~ParallelPlateFlowSolverSteadyState();
  
  void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 
  void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain, realT time, realT dt);
  void CalculateSteadyStatePressureDistribution(PhysicalDomainT& domain, SpatialPartition& partition,realT time );
  using ParallelPlateFlowSolverBase::UpdateEOS;
  void UpdateEOS( const realT dt, PhysicalDomainT& domain );
  using ParallelPlateFlowSolverBase::UpdateFlux;
  void UpdateFlux( const realT time, const realT dt, PhysicalDomainT& domain, SpatialPartition& partition);
  
  // Pressure controlled boundary condition
  void PressureBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
         BoundaryConditionBase* bc, const lSet& set, realT time);
                                  
  void PressureBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
         BoundaryConditionBase* bc, const lSet& set, realT time);


  // Single partition periodic boundary condition
  void SinglePartitionBC_NeighborUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                        BoundaryConditionBase* bc, realT time );

  void SinglePartitionBC_Sparsity(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                  BoundaryConditionBase* bc, realT time );

  void SinglePartitionBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                         BoundaryConditionBase* bc, realT time );

  void SinglePartitionBC_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                        BoundaryConditionBase* bc, realT time );
                                                 
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
  static const char* SolverName(){return "ParallelPlateFlowSolverSteadyState";};

private:

  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);

  using ParallelPlateFlowSolverBase::Assemble;
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition);
  
  const lSet* m_faceSet;

  localIndex m_numFaces;
  std::map<localIndex,localIndex> m_faceDofMap; 
  
  std::map<localIndex,lArray1d> m_edgesToFaces;
  
  // Flags
  bool m_doApertureUpdate;
  bool m_useMLPrecond;
  
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
  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFieldsB;
  
  std::string m_TrilinosIndexStr;
  static int m_instances;
  
  struct
  {
    realT m_tol; // Solver convergence criteria
    int m_maxIters; // Maximum number of solver iterations
    bool m_verboseFlag;
  }
  m_numerics;
  
};

int ParallelPlateFlowSolverSteadyState::m_instances = 0;

#endif /* SSPARALLELPLATEFLOWSOLVER_H_ */
