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
//  LLNL-CODE-656616
//  GEOS-CORE, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GEOS-CORE. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
//
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
 * @file ImplicitMechanicsSolver.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef LAGRANGESOLVERBASE_H_
#define LAGRANGESOLVERBASE_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"
#include "SurfaceGeneration/FractunatorBase.h"

#ifdef SRC_EXTERNAL
#include "Contact/CommonPlaneContact.h"
#endif
/*
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
*/

class MaterialBaseParameterData;


class LagrangeSolverBase: public SolverBase
{
public:
  LagrangeSolverBase(  const std::string& name,
                       ProblemManagerT* const pm );
  virtual ~LagrangeSolverBase();

  double TimeStep(const realT& time,
                const realT& dt,
                const int cycleNumber,
                PhysicalDomainT& domain,
                const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                FractunatorBase* const fractunator);


  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition );

  void InitializeCommunications( PartitionBase& partition );

  virtual void RegisterFields(PhysicalDomainT& domain);

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn );

  void RegisterTemporaryFields( PhysicalDomainT& domain );
  void DeregisterTemporaryFields( PhysicalDomainT& domain );
  void FillTemporaryFields( PhysicalDomainT& domain );
  void OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain );

  void PostProcess (PhysicalDomainT& domain,
                    SpatialPartition& partition,
                    const sArray1d& namesOfSolverRegions);

  enum TimeIntegrationEnum
  {
    FirstTimeIntegrationEnum = 0,
    QuasiStatic = 0,
    ImplicitDynamic = 1,
    ExplicitDynamic = 2,
    numTimeIntegrationEnums = 3
  };

  enum TwoDOptions
  {
    PlaneStrain = 0,
    PlaneStress = 1,
    Axisym = 2,
    num2dOptions = 3
  };

  TimeIntegrationEnum IntToTimeIntegration( const int ival )
  {
    switch( ival )
    {
      case 0:
        return QuasiStatic;
      case 1:
        return ImplicitDynamic;
      case 2:
        return ExplicitDynamic;
    }
    throw GPException("Invalid value input into LagrangeSolverBase::IntToTimeIntegration()");
    return numTimeIntegrationEnums;
  }

  TwoDOptions IntToTwoDOptions( const int ival )
  {
    switch( ival )
    {
      case 0:
        return PlaneStrain;
      case 1:
        return PlaneStress;
      case 2:
        return Axisym;
    }
    throw GPException("Invalid value input into LagrangeSolverBase::IntToTwoDOptions()");
    return num2dOptions;
  }

  void SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                     SpatialPartition& partition,
                                     int& numLocalRows,
                                     int& numGlobalRows );

  typedef std::map<ElementManagerT::RegKeyType, ElementRegionT > RegionMap;



  virtual realT Assemble    ( PhysicalDomainT& domain,
                             const SpatialPartition& partition,
                             const realT time );

  /*
  realT AssembleFluidPressureContributions( PhysicalDomainT& domain,
                                           const iArray1d& deformationTrilinosIndex,
                                           const iArray1d& flowTrilinosIndex,
                                           const Array1dT< rArray1d >& dwdu,
                                           const rArray1d& dwdw,
                                           const int flowDofOffset );
*/


  std::string m_trilinosIndexStr;
  static int m_instances;

  realT DampingM() const { return m_dampingM; }



  int m_enableTimeIntegrationOption[numTimeIntegrationEnums];
  TimeIntegrationEnum m_timeIntegrationOption;

protected:

  virtual void InsertGlobalIndices( PhysicalDomainT& domain){}

  virtual void UpdateContactDataStructures( PhysicalDomainT& domain, const realT time){}

  virtual void GetContactStiffnessContribution( PhysicalDomainT& domain){}

  virtual void UpdatePlasticStrainsForConvergedIteration( PhysicalDomainT& domain){}

  virtual void StoreHistoryVariablesForCurrentLoadStepAndResetTheField( PhysicalDomainT& domain){}

  int dim;
  const unsigned this_mpi_process;
  const unsigned n_mpi_processes;
  const bool     verbose;

  bool m_writeNodalStress;


public:
  realT m_dampingM;
  realT m_dampingK;

  realT m_timeToSnapshotDisp;
  // Fu: If we specify an initial stress field, the stress field itself might not be balanced.
  // We need to store the displacement after the model has settled.  We use the difference between the total disp and this ref disp for the Displace transformation in VisIt.

  realT m_refTemperature;
  int m_useNodalTemperature;
  sArray1d m_thermalRegionNames;
  R1Tensor m_gravityVector;
protected:
  int m_tiedNodesFlag;
  Array1dT<lSet> m_KinematicConstraintNodes;
  realT m_tiedNodeNormalRuptureStress;
  realT m_tiedNodeShearRuptureStress;
  realT m_tiedNodeTolerance;

  iArray1d dummyDof;
  

  realT m_cfl;


  TwoDOptions m_2dOption;

  realT m_bulkQLinear;
  realT m_bulkQQuadratic;
  

private:
  virtual void SetupSystem ( PhysicalDomainT& domain, SpatialPartition& partition, const realT time, const realT& dt);

  virtual void TimeStepExplicitDynamic( const realT& time,
                                        const realT& dt ,
                                        PhysicalDomainT& domain,
                                        const sArray1d& namesOfSolverRegions,
                                        SpatialPartition& partition );

public:  // Fu note: We need to access this from the new hydrofrac solver, which is not a derived class.
  void ApplyForcesFromContact(PhysicalDomainT& domain,
                              StableTimeStep& timeStep,
                              const realT dt );


  virtual void ApplyThermalStress( ElementRegionT& elemRegion,
                                   const NodeManagerT& nodeManager,
                                   const localIndex& elementID,
                                   Epetra_SerialDenseVector& rhs,
                                   const bool updateTemperatureOnly){}

  virtual void ProcessElementRegions( NodeManagerT& nodeManager,
                                      ElementManagerT& elemManager,
                                      const sArray1d& namesOfSolverRegions,
                                      const realT dt  ) ;

  void SnapshotNodalDisplacement( NodeManagerT& nodeManager) ;
private:
  virtual void ProcessElementRegion( NodeManagerT& nodeManager,
                                     ElementRegionT& elemRegion,
                                     const realT dt ) = 0;

  virtual void ProcessCohesiveZones( NodeManagerT& nodeManager,
                                     FaceManagerT& faceManager,
                                     const realT dt );

  //virtual void RegisterFields_Derived(PhysicalDomainT& domain);
  //virtual void ReadXML_Derived(TICPP::HierarchicalDataNode* hdn);


  virtual realT CalculateElementResidualAndDerivative( const MaterialBaseParameterData& matParams,
                                                      const FiniteElementBase& fe,
                                                      const Array2dT<R1Tensor>& dNdX,
                                                      const realT* const detJ,
                                                      const Epetra_SerialDenseVector& dof_np1,
                                                      Epetra_SerialDenseMatrix& dRdU,
                                                      Epetra_SerialDenseVector& R ) = 0;

  virtual void IntegrateCohesiveZoneContributions( const int numNodesInFace,
                                                   const realT area,
                                                   const R1Tensor& traction,
                                                   const R2Tensor& stiffness,
                                                   Epetra_SerialDenseMatrix& face_matrix,
                                                   Epetra_SerialDenseVector& face_rhs ) ;
  virtual void ApplyGravity( NodeManagerT& nodeManager,
                                       ElementRegionT& elemRegion,
                                       const realT dt ) = 0;



private:
  virtual void SetupMLPreconditioner( const PhysicalDomainT& domain, ML_Epetra::MultiLevelPreconditioner* MLPrec );


  // boundary conditions
  virtual void TractionBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time){}

  virtual void PressureBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                          BoundaryConditionBase* bc, const lSet& set, realT time){}

  virtual void DisplacementBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                              BoundaryConditionBase* bc, const lSet& set, realT time){}

};



#endif /* ImplicitMechanicsSolver_H_ */
