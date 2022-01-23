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
 * @file FractureFlowSolver.h
 * @author hao1
 * @date Oct. 21, 2013
 */

#ifndef FRACTUREFLOWSOLVER_H_
#define FRACTUREFLOWSOLVER_H_

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

using namespace std;

class FractureFlowSolver : public SolverBase
{

public:

  typedef Array1dT<rArray1d> rSArray2d; 
  typedef Array1dT<iArray1d*> iArrayPtrs;
  typedef Array1dT<rArray1d*> rArrayPtrs;

  FractureFlowSolver(const std::string &name, ProblemManagerT* const problemManager);
  virtual ~FractureFlowSolver();
  
  void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, 
                 SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 
  static const char* SolverName(){return "FractureFlowSolver";};

private:

  void Assemble(PhysicalDomainT& domain, SpatialPartition& partition, const realT& time);
  realT Solve(PhysicalDomainT& domain, SpatialPartition& partition);
  
  void MakeOldAcc(PhysicalDomainT& domain);

  void SetDirichletNodeBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time);

  void DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set, realT time);

  void SetDirichletBCs(PhysicalDomainT& domain);

  void UpdateDirichletBCs(PhysicalDomainT& domain, const realT &time);

  void SetSrcFluxBCs(PhysicalDomainT& domain);

  void SetSrcFluxBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time);

  void UpdateSrcFluxBCs(PhysicalDomainT& domain, const realT &time);

  void SrcFluxBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time);

  void UpdateFluidProps(PhysicalDomainT& domain);

  void UpdateRockProps(PhysicalDomainT& domain);

  void MakeFractureParameters(PhysicalDomainT& domain);

  void MakeMatrixParameters(PhysicalDomainT& domain);

  void CalculateInterpCoefs(PhysicalDomainT& domain);

  // Flags
  bool m_doApertureUpdate;
  bool m_doRockPropUpdate;
  bool m_useMLPrecond;
  bool m_doMFCoupling;
  
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
    realT m_tol; // Solver convergence criteria
    int m_maxIters; // Maximum number of solver iterations

  } m_numerics;

  realT m_dt;


  // NR iteration convergence criteria

  struct 
  {
    realT m_tol;
    realT m_maxIters;

  } m_NRNumerics;


  struct 
  {
    //liquid phase properties

    realT m_mu; 
    realT m_rho_o; 
    realT m_press_o;
    realT m_compress;

    realT m_fluidCond;
    realT m_fluidHeatCap;

    //rock properties

    rArray1d m_rockCompress; //for multiple matrix element regions
    rArray1d m_rockPermCoef; //for multiple matrix element regions

    realT  m_maxPorosity;
    realT  m_minPorosity;

    realT m_maxAperture;
    realT m_minAperture;

    /*
        The following parameters are modified to
       account for multiple matrix element regions.
       last modified by Pratanu Roy (roy24) Nov 2014
    */

    //realT m_rockHeatCap;
    //realT m_rockCond;
    //realT m_rockDensity;
    rArray1d m_rockHeatCap;
    rArray1d m_rockCond;
    rArray1d m_rockDensity;


  } m_fluidRockProperty;

  struct {

    int ncomp;
    int nvar;
    bool thermal;

  } m_eqt;
  
  realT m_gravFactor;
  R1Tensor m_downVector;


  // parameters used to define matrix/fracture element connections

  /* 
     In order to save memory cost some of the following parameters 
     may be defined or calculated during time stepping instead of 
     being stored permanently. 

     We will test fracture-flow solver efficiency later, in particular 
     re-examining/organizing its data structure for better 
     performance. 

  */

  pArray1d m_elemIndex; 
  iArray1d m_accumIndex;

  rSArray2d m_volume; 
  rSArray2d m_rockHeatCap;
 
  Array1dT<rSArray2d> m_femCoefs;

  rArray1d m_density;
  rArray1d m_elev;

  rArray1d m_oldDensity;

  //trilinos index

  iArrayPtrs m_elem_is_ghost;

  //global variables

  sArray1d m_fieldName;

  iArray1d  m_elemRegionIndex;
  Array1dT<ElementRegionT *>  m_elemRegion;

  sArray1d m_fracture_regionName;
  sArray1d m_matrix_regionName;

  localIndex m_fracture_regionNumber;
  localIndex m_matrix_regionNumber;

  Array1dT<rSArray2d> m_femCoefs2;
  rArray1d m_oldTemperature;

  sArray1d m_thermalProductionBCSetName;

  realT m_dtMin;
  realT m_dtMax;

  Array1dT<rArray2d> m_interpCoefs;

};

int FractureFlowSolver::m_instances = 0;

#endif /* FRACTUREFLOWSOLVER_H_ */
