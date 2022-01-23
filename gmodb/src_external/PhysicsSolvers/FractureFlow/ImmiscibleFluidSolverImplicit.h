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
 * @file ImmiscibleFluidSolverImplicit.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef PARALLELPLATEFLOWSOLVER_IM_H_
#define PARALLELPLATEFLOWSOLVER_IM_H_

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


class ImmiscibleFluidSolverImplicit : public ParallelPlateFlowSolverBase
{
public:
  ImmiscibleFluidSolverImplicit( const std::string& name,
                                 ProblemManagerT* const pm );
  virtual ~ImmiscibleFluidSolverImplicit();
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );

  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition );

  void InitializeCommunications( PartitionBase& partition );

  virtual void SetMaxStableTimeStep( PhysicalDomainT& domain,
                                     const sArray1d& namesOfSolverRegions,
                                     SpatialPartition& partition);

  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 

  void CalculateMassRate(PhysicalDomainT& domain, SpatialPartition& partition,realT time,realT dt );

  using ParallelPlateFlowSolverBase::UpdateEOS;
  void UpdateEOS( PhysicalDomainT& domain,const realT dt);
  void EquilibratePressure(realT& rho_red, realT& rho_blue, realT& phi, realT& Pnew);

  using ParallelPlateFlowSolverBase::UpdateFlux;
  void UpdateFlux( const realT time, SpatialPartition& partition,PhysicalDomainT& domain);

  void UpdateColorFunction( PhysicalDomainT& domain,SpatialPartition& partition,const realT time, const realT dt  );
  
  // Pressure controlled boundary condition
  void PressureBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, 
         BoundaryConditionBase* bc, const lSet& set, realT time);
                                  
  void PressureBoundaryCondition_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
         BoundaryConditionBase* bc, const lSet& set, realT time);

  void PressureBoundaryCondition_ColorFunctionUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object,
                                                BoundaryConditionBase* bc, const lSet& set, realT time);

  void OutflowPressureBoundaryCondition_ColorFunctionUpdate(PhysicalDomainT& domain,
                                                            ObjectDataStructureBaseT& object,
                                                            BoundaryConditionBase* bc, const lSet& set, realT time);
                                                 
  void SetConcentrationField(std::string& concentrationFieldName,std::string reactionRateFieldName,std::string rrDerivFieldName);

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "ImmiscibleFluidSolverImplicit";};

private:


  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);

  using ParallelPlateFlowSolverBase::Assemble;
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition);

  using ParallelPlateFlowSolverBase::GenerateParallelPlateGeometricQuantities;
  void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain, SpatialPartition* partition, realT dt );
  void InitializeDensity( PhysicalDomainT& domain);
  void UpdateAperture(PhysicalDomainT&  domain);
  void OverwriteOldGeometricQuantities( PhysicalDomainT& domain);

  realT TwoFacePermeability(const Array1dT<R1Tensor>& edgeCenters,
                            const rArray1d& edgeLengths,
                            const Array1dT<R1Tensor>& faceCenters,
                            const rArray1d& apertures,
                            localIndex eg,localIndex kf, localIndex kfb);

  std::string m_flowFaceSetName;
  lSet* m_nodeSet;
  const lSet* m_faceSet;
  
  localIndex m_numFaces;
  std::map<localIndex,localIndex> m_faceDofMap; 
  
  std::map<localIndex,lArray1d> m_edgesToFaces;
  
  // Flags
  bool m_doApertureUpdate;
  bool m_verboseFlag;
  bool m_useMLPrecond; // fixme - ultimately want to replace this with generic preconditioner object
  
  realT m_phi; // Mixed Euler parameter (0 = forward difference, 0.5 = central difference, 1.0 = backward difference)
  
  // Boundary conditions
  realT m_dt; // current dt

  // Timestep range
  realT m_min_dt; // initial dt
  realT m_max_dt;
  
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
  
  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;  // synced fields for the flow model
  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedImmiscibleFields; // synced fields for the immiscible model
  
  std::string m_TrilinosIndexStr;
  static int m_instances;
  
  struct NumericsStruct
  {
    NumericsStruct():
      m_tol(1e-10),
      m_maxIters(1000){};

    realT m_tol; // Solver convergence criteria
    int m_maxIters; // Maximum number of solver iterations
  }
  m_numerics;

  //////////////////////


  // fluid properties
  struct FluidData{
    realT m_mu;
    realT m_bulk_modulus;
    realT m_rho_o;
    R1Tensor m_bodyForce;
    void ReadXML( TICPP::HierarchicalDataNode* hdn ){
      m_mu = hdn->GetAttributeOrDefault("mu","1.0e-3 N.s/m^2");
      m_bulk_modulus = hdn->GetAttributeOrDefault(PS_STR::BulkModulusStr,"2.0e9 Pa");
      m_rho_o = hdn->GetAttributeOrDefault("rho_o","1 kg/L");
      m_bodyForce =  hdn->GetAttributeOrDefault<R1Tensor>("body_force",R1Tensor());
     /* std::cout << "        FluidData: mu " << m_mu << std::endl;
      std::cout << "        FluidData: K_o " << m_bulk_modulus << std::endl;
      std::cout << "        FluidData: rho_o " << m_rho_o << std::endl;
      std::cout << "        FluidData: body_force " << m_bodyForce << std::endl;*/
    }
  };
  FluidData m_redFluidData;
  FluidData m_blueFluidData;


  // Permeability
  struct {
    realT m_SHP_FCT;
    realT m_min_aperture;
    realT m_max_aperture;
    void ReadXML( TICPP::HierarchicalDataNode* hdn ){
      m_SHP_FCT = hdn->GetAttributeOrDefault<realT>("shapefactor",1.0);
      m_min_aperture = hdn->GetAttributeOrDefault(PS_STR::MinimumApertureStr,"1 um");
      m_max_aperture = hdn->GetAttributeOrDefault(PS_STR::MaximumApertureStr,"4 mm");
    }
  } m_permeabilityData;

  // surface tension
  struct {
    realT m_sigma;   //
    bool m_calculateInPlaneCurvature;
    bool m_useConvergenceAngle;
    realT m_alpha;   // wetting angle
    void ReadXML( TICPP::HierarchicalDataNode* hdn ){
      m_sigma = hdn->GetAttributeOrDefault("surface_tension",0.1);
      m_calculateInPlaneCurvature = hdn->GetAttributeOrDefault("calculate_in_plane_curvature",false);
      m_useConvergenceAngle = false; // FIXME - not implemented
      realT w = hdn->GetAttributeOrDefault("wetting_angle",30.0); // wetting angle in degrees
      m_alpha = w*( 3.141592653589793 / 180.0);  // wetting angle in radians
     /* std::cout << "        SurfaceTension: sigma " << m_sigma << std::endl;
      std::cout << "        SurfaceTension: alpha " << m_alpha << std::endl;*/
    };
    realT CalculateNormalRadiusOfCurvature(realT h, realT beta){
      return 0.5*h/cos(m_alpha + beta);
    };

  } m_surfaceTensionData;

  // color function
  struct {
    realT m_D;
    realT m_beta;
    void ReadXML( TICPP::HierarchicalDataNode* hdn ){
      m_D = hdn->GetAttributeOrDefault("diffusivity",0.1);
      m_beta = hdn->GetAttributeOrDefault("beta",0.9);
     /* std::cout << "        ColorFunction: diffusivity " << m_D << std::endl;
      std::cout << "        ColorFunction: beta " << m_beta << std::endl;*/
    };
  } m_colorFunctionData;

  
};

int ImmiscibleFluidSolverImplicit::m_instances = 0;

#endif /* PARALLELPLATEFLOWSOLVER_IM_H_ */
