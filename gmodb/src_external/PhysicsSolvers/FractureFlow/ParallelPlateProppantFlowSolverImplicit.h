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
 * @file ParallelPlateProppantFlowSolverImplicit.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef PARALLELPLATEPROPPANTFLOWSOLVER_IM_H_
#define PARALLELPLATEPROPPANTFLOWSOLVER_IM_H_

#include "PhysicsSolvers/ParallelPlateFlowSolverBase.h"

#include "ProppantModels.h"
#include "LeakoffModels.h"

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

// forward declarations
class ParallelPlateProppantSolver;


class ParallelPlateProppantFlowSolverImplicit : public ParallelPlateFlowSolverBase
{
public:
	ParallelPlateProppantFlowSolverImplicit( const std::string& name, ProblemManagerT* const problemManager);
  virtual ~ParallelPlateProppantFlowSolverImplicit();
  
  void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 
  void DefineFlowSets( PhysicalDomainT& domain );

  void CalculateMassRate(PhysicalDomainT& domain, SpatialPartition& partition,realT time,realT dt );
  using ParallelPlateFlowSolverBase::UpdateEOS;
  void UpdateEOS( PhysicalDomainT& domain,const realT time,const realT dt ); // const bool updateMass
  using ParallelPlateFlowSolverBase::UpdateFlux;
  void UpdateFlux( const realT time,
                   const realT dt, PhysicalDomainT& domain, SpatialPartition& partition);
  
  // Pressure controlled boundary condition
  void PressureBoundaryCondition( PhysicalDomainT& domain,
                                  ObjectDataStructureBaseT& object,
                                  BoundaryConditionBase* const bc,
                                  const lSet& set,
                                  const realT time,
                                  const realT dt );
                                  
  void PressureBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
         BoundaryConditionBase* bc, const lSet& set, realT time);

  // Face Pressure BC
  void FixedFacePressureBoundaryCondition( PhysicalDomainT& domain,ObjectDataStructureBaseT& object,
          BoundaryConditionBase* const bc, const lSet& set, const realT time, const realT dt );

  // Face Aperture BC
  void FixedApertureBoundaryCondition( PhysicalDomainT& domain,ObjectDataStructureBaseT& object,
          BoundaryConditionBase* const bc, const lSet& set, const realT time, const realT dt );

  // Velocity boundary condition
  void FixedFluxBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
    BoundaryConditionBase* bc, const lSet& set,const realT time, const realT dt);

  void FixedFluxBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
         BoundaryConditionBase* bc, const lSet& set, realT time);

  // Net flux boundary condition
  void FixedNetFluxBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
    BoundaryConditionBase* bc, const lSet& set,const realT time, const realT dt);

  void FixedNetFluxBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
         BoundaryConditionBase* bc, const lSet& set, realT time);

  void FixedNetFluxBoundaryConditionB(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
    BoundaryConditionBase* bc, const lSet& set,const realT time, const realT dt);

  void FixedNetFluxBoundaryConditionB_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
         BoundaryConditionBase* bc, const lSet& set, realT time);
                                                 
  void SetConcentrationField(std::string& concentrationFieldName,std::string reactionRateFieldName,std::string rrDerivFieldName);

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "ParallelPlateProppantFlowSolverImplicit";};

  void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,
                                                 realT time,
                                                 realT dt);

  void GenerateSlurryParameters( PhysicalDomainT& domain );

  void OverwriteOldGeometricQuantities( PhysicalDomainT& domain);

  realT ProppantPackPermeability(realT h, realT w,realT la,const realT& max_proppant_vf, const realT& fluid_mu);
  realT PartiallyFilledEdgePermeability(realT packVf, realT kappaOpen, realT kappaPack);

  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const sArray1d& namesOfSolverRegions);

  // Flags
  bool m_doApertureUpdate;

protected:

  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);
  using ParallelPlateFlowSolverBase::Assemble;
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition);


  void InitializeDensity( PhysicalDomainT& domain);
  void UpdateAperture(PhysicalDomainT&  domain);

  realT TwoFacePermeability(const Array1dT<R1Tensor>& edgeCenters,
                            const rArray1d& edgeLengths,
                            const Array1dT<R1Tensor>& faceCenters,
                            const rArray1d& apertures,
                            const rArray1d& packVfs,
                            localIndex eg,localIndex kf, localIndex kfb);


  realT TwoFacePermeability_PowerLaw(const Array1dT<R1Tensor>& edgeCenters,
                                     const rArray1d& edgeLengths,
                                     const Array1dT<R1Tensor>& faceCenters,
                                     const rArray1d& apertures,
                                     const Array1dT<R1Tensor>& fluidFlow,
                                     const rArray1d& packVfs,
                                     localIndex eg,localIndex kf, localIndex kfb);


  realT TwoFacePermeability_HerschelBulkley(const Array1dT<R1Tensor>& edgeCenters,
                                     const rArray1d& edgeLengths,
                                     const Array1dT<R1Tensor>& faceCenters,
                                     const rArray1d& apertures,
                                     const Array1dT<R1Tensor>& fluidFlow,
                                     const rArray1d& packVfs,
                                     localIndex eg,localIndex kf, localIndex kfb);



  realT OneFacePermeability(const Array1dT<R1Tensor>& edgeCenters,
                            const rArray1d& edgeLengths,
                            const Array1dT<R1Tensor>& faceCenters,
                            const rArray1d& apertures,
                            const rArray1d& packVfs,
                            localIndex eg,localIndex kf);


  realT OneFacePermeability_PowerLaw(const Array1dT<R1Tensor>& edgeCenters,
                                     const rArray1d& edgeLengths,
                                     const Array1dT<R1Tensor>& faceCenters,
                                     const rArray1d& apertures,
                                     const Array1dT<R1Tensor>& fluidFlow,
                                     const rArray1d& packVfs,
                                     localIndex eg,localIndex kf);



  realT OneFacePermeability_HerschelBulkley(const Array1dT<R1Tensor>& edgeCenters,
                                     const rArray1d& edgeLengths,
                                     const Array1dT<R1Tensor>& faceCenters,
                                     const rArray1d& apertures,
                                     const Array1dT<R1Tensor>& fluidFlow,
                                     const rArray1d& packVfs,
                                     localIndex eg,localIndex kf);

  ProppantBase m_proppantData;
  LeakoffBase m_leakoffModel;
  
  realT CalculateFacePermeability(realT app,realT qMag,realT concentration, realT packVf);

private:

  lSet m_faceSet;
  std::string m_faceSetName;
  localIndex m_numFaces;
  std::map<localIndex,localIndex> m_faceDofMap; 
  
  std::map<localIndex,lArray1d> m_edgesToFaces;
  
  
  realT m_phi; // Mixed Euler parameter (0 = forward difference, 0.5 = central difference, 1.0 = backward difference)
  
  // Barton joint model
  realT m_bBarton;
  realT m_aBarton;
  realT m_wZeroStress;

  
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
  }
  m_numerics;

  int m_cycleNumber;

  // flags
  bool m_verboseFlag;
  bool m_useHBPowerlawFluid;
  bool m_doErosion;


  realT m_dt; // for boundary conditions
  realT m_maxdt;
  realT m_apertureRelaxationTime;

  // boundary condition data
  realT m_fixedFlux_pInj;

  friend class ParallelPlateProppantSolver; // only so proppant data can be shared - need to fix
};


#endif /* PARALLELPLATEFLOWSOLVER_IM_H_ */
