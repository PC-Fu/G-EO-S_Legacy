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
#ifndef TWODADRSOLVER_PROPPANT_MPI_H_
#define TWODADRSOLVER_PROPPANT_MPI_H_

/*
 * ParallelPlateProppantSolver.h
 *
 *  Created on: Jun 7, 2011
 *      Author: walsh24
 */

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/FunctionManager.h"
#include "ObjectManagers/EdgeManagerT.h"
#include "Common/Common.h"

#include "ElementLibrary/FiniteElement.h"

#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include "ProppantModels.h"

#include <map>
#include <string>
#include <vector>

namespace TDSSADR_STR{
  const std::string ReactionRateDerivStr = "ReactionRateDerivative";
  const std::string PreviousConcentrationStr = "PreviousConcentration";
}

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

class ParallelPlateProppantFlowSolverImplicit;

class ParallelPlateProppantSolver : public SolverBase
{
public:
  ParallelPlateProppantSolver(const std::string& name, ProblemManagerT* pm);
  ~ParallelPlateProppantSolver();
  
  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn);
  
  virtual void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );


  virtual double TimeStep( const realT& time,
                         const realT& dt,
                         const int cycleNumber,
                         PhysicalDomainT& domain,
                         const sArray1d& namesOfSolverRegions,
                         SpatialPartition& partition,
                         FractunatorBase* const fractunator );
                         
  virtual void RegisterFields( PhysicalDomainT& domain );
  
  virtual void DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bc, const lSet& set,realT time);
  virtual void OutflowBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bc, const lSet& set,realT time);

  void SetConcentrationField(std::string& concentrationFieldName,std::string reactionRateFieldName ="",std::string RRDerivFieldName ="");

  static const char* SolverName(){return "ParallelPlateProppantSolver";};
  
  realT GetDiffusivity(void) const{return m_D.Trace()/3.0; };
  
protected:

  void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition);
  void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time,const realT& dt);
  void Solve       (PhysicalDomainT& domain, SpatialPartition& partition);
  
  void DefineFlowSets( PhysicalDomainT& domain );

  virtual void GenerateFaceCenters(NodeManagerT& nodeManager,FaceManagerT& faceManager);
  
  inline void CalculateSlipVelocity(const realT volumeFraction, const R1Tensor& fractureNormal, const R1Tensor& mixtureVelocity, R1Tensor& slipVelocity);
  
  void CalculateScreenout(PhysicalDomainT& domain, SpatialPartition& partition, const realT& time,const realT& dt);
  void CalculateErosionRate(PhysicalDomainT& domain, SpatialPartition& partition, const realT& time,const realT& dt);

  //
  localIndex m_numNodes;
  localIndex m_numFaces;
  localIndex m_numLocalFaces;
  //localIndex m_cumNodeCount;
  
  std::map<localIndex,localIndex> m_faceDofMap; 
  std::map<localIndex,localIndex> m_nodeDofMap; // maps dof for each node to node number. 
  
  std::map<localIndex,lArray1d> m_edgesToFaces;
  
  lSet* m_faceSet;
  lSet* m_nodeSet;
  std::string m_faceSetName;
  
  // Flags
  //bool doDataWrite;
  bool m_verboseFlag;
  bool m_doSteadyState;

  bool m_doScreenOut;
  bool m_doErosion;  // calculate erosion rate and update erosion aperture
    
  rArray1d m_weakBCdiag;
  
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
    bool useMLPrecond;
  }
  numerics;
  
  ParallelPlateProppantFlowSolverImplicit* m_flowSolverPtr;
  ProppantBase*  m_proppantDataPtr;
  SlurryBase* m_slurryModelPtr;

  ScourErosionBase m_ErosionModel; // fixme wrong place for this
  MatrixFluxErosionBase m_MatrixErosionModel; // fixme wrong place for this

  realT m_maxFluxCoefficient;

  /*
  struct ProppantData
  {
    realT m_density;
    realT m_fluidDensity;
    realT m_fluidViscosity;
    realT m_diameter;
    realT m_gravity;
    realT m_maxVf;  // maximum volume fraction
  } m_proppantData;

  struct SlurryModel {
    realT m_singleParticleSettlingSpeed;
    realT m_hinderedSettlingCoefficient;
    realT m_maxVf;

    void Initialize(ProppantData proppantData,std::string singleParticleModelType){
    	realT g = proppantData.m_gravity;
    	realT mu =  proppantData.m_fluidViscosity;
        realT d = proppantData.m_diameter;
        realT rhos = proppantData.m_density;
        realT rhol = proppantData.m_fluidDensity;
    	m_hinderedSettlingCoefficient = -5.9;
        m_maxVf = proppantData.m_maxVf;
    	if( streq(singleParticleModelType,"stokes") ){
    	  m_singleParticleSettlingSpeed = g*(rhos-rhol)*d*d/18.0*mu;
    	} else if( streq(singleParticleModelType,"newton") ) {

      	  m_singleParticleSettlingSpeed = 1.74*sqrt(d)*sqrt(g*(rhos-rhol)/rhol);
    	} else if( streq(singleParticleModelType,"allen") ) {
      	  m_singleParticleSettlingSpeed = 0.2*pow(d,1.18)*pow( (g*(rhos-rhol)/rhol),0.72)/pow(mu/rhol,0.45);

    	} else {
          throw GPException("Error Slurry Model: Unsupported model type:" +singleParticleModelType +".");
    	}
    };

    realT HinderedSettlingSpeed(realT phi){
    	   return m_singleParticleSettlingSpeed*exp(m_hinderedSettlingCoefficient*phi);
    };
  } m_slurryModel;
  */

private: 

  // Problem specific flags
  bool hasReaction;
  
  // Element data
  R2Tensor m_D;   // diffusivity
  realT m_alpha_L; // longitudinal dispersivity
  realT m_alpha_T; // transverse dispersivity
    
  // Field names
  std::string concentrationFieldName_;
  std::string reactionRateFieldName_;
  std::string rrDerivFieldName_; 

  sArray1d m_concentrationFieldNames;
  sArray1d m_reactionRateFieldNames;
  
  bool m_useCollisionalSlipVelocity;
  realT SuperadvectiveTransportFactor(const realT vf);

  realT m_dt; // delta t - needed for boundary conditions
  realT m_init_dt;

  // solver map - provide link to flow solver
  std::map<std::string,SolverBase*>* m_solverMapPtr;
  std::string m_FlowSolverName;

  R1Tensor m_upVector;


};


#endif /*TWODADRSOLVER_PROPPANT_H_*/
