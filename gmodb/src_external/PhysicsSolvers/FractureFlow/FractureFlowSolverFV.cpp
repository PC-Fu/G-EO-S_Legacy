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
 * @file FractureFlowSolver.cpp
 * @author hao1
 * @date Oct. 21, 2013
 */

#include "FractureFlowSolver.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"


// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

using namespace BoundaryConditionFunctions;
using namespace PS_STR;

FractureFlowSolver::FractureFlowSolver(  const std::string& name,
                                         ProblemManagerT* const pm ):
          SolverBase(name,pm),
          m_useMLPrecond(false),
          this_mpi_process(pm->m_epetraComm.MyPID()),
          n_mpi_processes(pm->m_epetraComm.NumProc()),
          m_epetra_comm(pm->m_epetraComm),
          row_map(),
          sparsity(),
          matrix(),
          solution(),
          rhs(),
          syncedFields(),
          m_TrilinosIndexStr(),
          m_dt(0.0)
{
  ++m_instances; 
  m_TrilinosIndexStr = "FFS_SS_" +  toString<int>(m_instances) + "_GlobalDof";
}

FractureFlowSolver::~FractureFlowSolver() {}

/*

      <FractureFlowSolver name="fracMatrixFlowSolver"
                          coupleMatrixFractureFlow="1"
                          fractureRegionName="EB200" 
                          updateAperture="0" 
                          fieldName="Pressure"
                          rockCompress="3.5e-10 1/Pa"
                          downVector="0.0 0.0 0.0"
      />


*/

void FractureFlowSolver::ReadXML( TICPP::HierarchicalDataNode* hdn )
{ 

  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("updateAperture",false);
  m_doMFCoupling = hdn->GetAttributeOrDefault<bool>("coupleMatrixFractureFlow",false);

  // Linear solver parameters
  m_numerics.m_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-5);
  m_numerics.m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);
  
  m_useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);

  // Nonlinear solver parameters
  m_NRNumerics.m_tol = hdn->GetAttributeOrDefault<realT>("nrTol", 1e-5);
  m_NRNumerics.m_maxIters = hdn->GetAttributeOrDefault<realT>("nrMaxIterations", 10);

  // equations
  m_eqt.ncomp = 1;
  m_eqt.thermal = hdn->GetAttributeOrDefault<bool>("thermal",false);
  m_eqt.nvar = m_eqt.thermal ? m_eqt.ncomp + 1 : m_eqt.ncomp;


  // fluid and rock properties

  m_fluidRockProperty.m_compress = hdn->GetAttributeOrDefault("fluidCompress", "5.0e-10 1/Pa");

  m_fluidRockProperty.m_mu = hdn->GetAttributeOrDefault("fluidViscosity", "0.001 Pa.s");

  m_fluidRockProperty.m_rho_o = hdn->GetAttributeOrDefault("rhoRef","1000.0 kg/m^3");

  m_fluidRockProperty.m_press_o = hdn->GetAttributeOrDefault("pressRef","1.0e5 Pa");

  m_fluidRockProperty.m_rockCompress = hdn->GetAttributeOrDefault("rockCompress", "3.5e-10 1/Pa");

  m_fluidRockProperty.m_minAperture = hdn->GetAttributeOrDefault(MinimumApertureStr,"1 um");

  m_fluidRockProperty.m_maxAperture = hdn->GetAttributeOrDefault(MaximumApertureStr,"4 mm");

  m_fluidRockProperty.m_maxPorosity = hdn->GetAttributeOrDefault("maxrockPorosity","0.999");

  m_fluidRockProperty.m_minPorosity = hdn->GetAttributeOrDefault("minRockPorosity","0.001");

  m_gravFactor = hdn->GetAttributeOrDefault("gravFactor", "9.81 m^2/s");

  m_downVector = hdn->GetAttributeTensorOrDefault("downVector", R1Tensor(0.0, 0.0, -1.0));

  m_fracture_regionName = hdn->GetStringVector("fractureRegionName");

  m_fracture_regionNumber = m_fracture_regionName.size();

  if(m_fracture_regionName.empty())
      throw GPException("ERROR: 'fractureRegionName' is not set");

  m_fieldName = hdn->GetStringVector("fieldName");

  if(m_fieldName.empty())
      throw GPException("ERROR: fieldName is not set!");

}

void FractureFlowSolver::RegisterFields( PhysicalDomainT& domain )
{

  ElementManagerT& elementManager = domain.m_feElementManager;

  int num;

  num = (elementManager.m_ElementRegions).size();

  m_elemRegionIndex.resize(num);
  m_elemRegionIndex = -1;

  for(int i =0; i < m_fracture_regionNumber; ++i){

    std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.find(m_fracture_regionName[i]);

    if(it == elementManager.m_ElementRegions.end())
      throw GPException("Fracture region " + m_fracture_regionName[i] + " is not found!");

    ElementRegionT& elemRegion = it->second;
    m_elemRegionIndex[elemRegion.m_regionNumber] = i;

    m_elemRegion.push_back(&elemRegion);

    elemRegion.AddKeylessDataField<realT>(ApertureStr,true,true);

    elemRegion.AddKeyedDataField<FieldInfo::pressure>();    

    elemRegion.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);

  }

  if(m_doMFCoupling) {

    int index;

    m_matrix_regionNumber = 0;

    std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.begin();

    for(localIndex i =0; i < elementManager.m_ElementRegions.size(); i++, it++){
      ElementRegionT& elemRegion = it->second;

      if(m_elemRegionIndex[elemRegion.m_regionNumber] == -1) {

	index = m_matrix_regionNumber + m_fracture_regionNumber;
	m_elemRegionIndex[elemRegion.m_regionNumber] = index;

	m_elemRegion.push_back(&elemRegion);

	elemRegion.AddKeylessDataField<realT>(PermeabilityStr,true,true);
	elemRegion.AddKeylessDataField<realT>(PorosityStr,true,true);

	elemRegion.AddKeyedDataField<FieldInfo::pressure>();    

	elemRegion.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);

	m_matrix_regionNumber++;

      }
    }
  }

}

void FractureFlowSolver::InitializeCommunications( PartitionBase& partition )
{

    syncedFields.clear();

    syncedFields[PhysicalDomainT::FiniteElementElementManager].push_back(Field<FieldInfo::pressure>::Name());

    if(m_eqt.thermal)
      syncedFields[PhysicalDomainT::FiniteElementElementManager].push_back("Temperature");

}


void FractureFlowSolver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{

  FaceManagerT& faceManager = domain.m_feFaceManager;
  EdgeManagerT& edgeManager = domain.m_feEdgeManager;

  int num, index;
  int m, n, p, q, num1, num2;

  int nvar = m_eqt.nvar;
			    
  //define m_trilinos_index  

  // allocate memory for parameters used to defined element connection
  // again some of them can be calculated during time stepping 
  // for memory saving

  num = m_elemRegion.size();

  m_elem_is_ghost.resize(num);
  m_trilinos_index.resize(num);
  m_pressure.resize(num);

  m_accumIndex.resize(num);
  m_density.resize(num);
  m_volume.resize(num);

  m_oldAcc.resize(m_eqt.nvar);

  for(int i = 0; i < m_eqt.nvar; i++)
    m_oldAcc[i].resize(num);

  index = 0;

  for(int i = 0; i < num; i++) {

    m_elem_is_ghost[i] = &(m_elemRegion[i]->GetFieldData<FieldInfo::ghostRank>());
    m_trilinos_index[i] = &(m_elemRegion[i]->GetFieldData<int>(m_TrilinosIndexStr));
    m_pressure[i] = &(m_elemRegion[i]->GetFieldData<FieldInfo::pressure>());

    m_density[i].resize(m_elemRegion[i]->m_numElems);

    for(int j = 0; j < m_eqt.nvar; j++)
      m_oldAcc[j][i].resize(m_elemRegion[i]->m_numElems);

    m_volume[i].resize(m_elemRegion[i]->m_numElems);

    m_accumIndex[i] = index;
    index += m_elemRegion[i]->m_numElems;

    for(localIndex j = 0; j < m_elemRegion[i]->m_numElems; j++) {
      m_elemIndex.push_back(std::make_pair(i, j));
    }
  }

  for(localIndex i =0; i < m_trilinos_index.size(); i++) 
    (*m_trilinos_index[i]) = 0;

  m_fractureVolume.resize(m_fracture_regionNumber);

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    m_fractureVolume[i].resize(m_elemRegion[i]->m_numElems);

  }

  if(m_doMFCoupling) {

    m_porosity.resize(m_matrix_regionNumber);

    for(int i = 0; i < m_matrix_regionNumber; i++) {

      index = m_fracture_regionNumber + i;
      ElementRegionT* elemRegion = m_elemRegion[index];
      m_porosity[i] = &(elemRegion->GetFieldData<realT>(PorosityStr));

    }

  }

  //define connenctions

  iArray1d& edge_is_ghost       = edgeManager.GetFieldData<FieldInfo::ghostRank>();

  iArray1d& face_is_ghost       = faceManager.GetFieldData<FieldInfo::ghostRank>();

  //fracture elements 

  localIndex kf, numEdges, eg;
  std::map<localIndex, localIndex>  fracture_faceMap; 

  index = 0;

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    const ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {

      kf = elemRegion->m_toFacesRelation[j][0];
      fracture_faceMap[kf] = i * m_accumIndex[i] + j;

      numEdges = faceManager.m_toEdgesRelation[kf].size();

      for(localIndex a =0; a < numEdges; ++a){

	eg = faceManager.m_toEdgesRelation[kf][a];

	if(edge_is_ghost[eg] < 0) {

	  std::map<localIndex,localIndex>::iterator itr = m_edgeIndex.find(eg);

	  if(itr == m_edgeIndex.end()) {

	    m_edgeIndex[eg] = index++;

	  }
	}
      }
    }
  }


  m_edgeToFace.resize(index); 

  std::map<localIndex,localIndex>::const_iterator edge = m_edgeIndex.begin();

  for(int i = 0; i < index; ++i, ++edge) {

    eg = edge->first;

    lSet& edgeFaces = edgeManager.m_toFacesRelation[eg];
      
    for( lSet::iterator edgeFace=edgeFaces.begin() ; edgeFace!=edgeFaces.end() ; ++edgeFace ){

      if(isMember(*edgeFace, fracture_faceMap)){

	m_edgeToFace[edge->second].push_back(fracture_faceMap[*edgeFace]);
      
      }
    }
  }

  // initialize Dirichlet boundary conditions

  SetDirichletBCs(domain);

  if(m_doMFCoupling) {

    m_faceToElem.resize(faceManager.m_numFaces);

    for(localIndex i = 0; i < faceManager.m_numFaces; i++) {

      if(face_is_ghost[i] < 0) {

	// a face has either one or two elements 
	num = faceManager.m_toElementsRelation[i].size();

	for(int j=0; j < num; j++) {

	  const std::pair< ElementRegionT*, localIndex >& elem = faceManager.m_toElementsRelation[i][j];

	  m = m_elemRegionIndex[(elem.first)->m_regionNumber];

	  //skip the fracure element that is on an external face
	  if(m < m_fracture_regionNumber) 
	    continue;

	  n = elem.second;

	  index = m * m_accumIndex[m] + n;
    
	  m_faceToElem[i].push_back(index);

	}

	if(isMember(i,fracture_faceMap)) {

	  m_faceToElem[i].push_back(fracture_faceMap[i]);

	}

      }

    }

  }

  MakeConnections(domain); // only called once


  num = m_edgeToFace.size();

  m_edgeToFaceAL.resize(num);
  m_edge_trans.resize(num);
  m_edge_dz.resize(num);
  m_edge_wt1.resize(num);
  m_edge_wt2.resize(num);

  for(localIndex i = 0; i < m_edgeToFace.size(); i++) {

    num = m_edgeToFace[i].size(); // an edge could be joined by  more than two faces 
    m_edgeToFaceAL[i].resize(num);

    num = m_edge_connection[i].size();
    m_edge_trans[i].resize(num);
    m_edge_dz[i].resize(num);
    m_edge_wt1[i].resize(num);
    m_edge_wt2[i].resize(num);

  }

  if(m_doMFCoupling) {

    num = m_faceToElem.size();
    m_faceToElemAL.resize(num);
    m_face_trans.resize(num);
    m_face_wt.resize(num);
    m_face_dz.resize(num);

    for(localIndex i = 0; i < m_faceToElem.size(); i++) {

      num = m_faceToElem[i].size();
      m_faceToElemAL[i].resize(num);

      num = m_face_connection[i].size();
      m_face_trans[i].resize(num);
      m_face_wt[i].resize(num);
      m_face_dz[i].resize(num);

    }

  }

  // calculate some geometric parameters 
  // (e.g. element volume, connection area)

  MakeGeometryParameters(domain);

  MakeTransParameters(domain);


  // define trilinos index

  // we need to examine the effect of element order on trilinos linear solver performance later !!

  // local rows

  int n_local_rows = 0;

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    const ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex  j = 0; j < elemRegion->m_numElems; j++) {

      if((*m_elem_is_ghost[i])[j] < 0 && (*m_trilinos_index[i])[j] != -1) {

	n_local_rows++;

      }
    }

  }

  if(m_doMFCoupling) {

    for(int i = 0; i < m_matrix_regionNumber; i++) {

      index = i + m_fracture_regionNumber;

      const ElementRegionT* elemRegion = m_elemRegion[index];

      for(localIndex  j = 0; j < elemRegion->m_numElems; j++) {

	if((*m_elem_is_ghost[index])[j] < 0 && (*m_trilinos_index[index])[j] != -1) {

	  n_local_rows++;

	}
      }

    }
  }
			    
  std::vector<int> gather(n_mpi_processes);

  m_epetra_comm.GatherAll(&n_local_rows,
                        &gather.front(),
                        1);

  int first_local_row = 0;
  int n_global_rows = 0;

  for(int proc=0; proc<n_mpi_processes; ++proc)
  {
    n_global_rows += gather[proc];
    if(proc<this_mpi_process)
      first_local_row += gather[proc];
  }
  
  unsigned local_count = 0;

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    const ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j = 0; j < elemRegion->m_numElems; j++) {
    
      if((*m_elem_is_ghost[i])[j] < 0) {

	if((*m_trilinos_index[i])[j] != -1) {
	  (*m_trilinos_index[i])[j] = first_local_row + local_count;
	  local_count++;
	}

      }
      else {
	
	(*m_trilinos_index[i])[j] = -INT_MAX;	
	
      }
    }
  }
			    
  if(m_doMFCoupling) {

    for(int i = 0; i < m_matrix_regionNumber; i++) {

      index = i + m_fracture_regionNumber;

      const ElementRegionT* elemRegion = m_elemRegion[index];

      for(localIndex  j = 0; j < elemRegion->m_numElems; j++) {
    
	if((*m_elem_is_ghost[index])[j] < 0) {

	  if((*m_trilinos_index[index])[j] != -1) {

	    (*m_trilinos_index[index])[j] = first_local_row + local_count;
	    local_count++;

	  }
	}
	else {
	
	  (*m_trilinos_index[index])[j] = -INT_MAX;	
	
	}
      }
    }
  }

  assert(static_cast<int>(local_count) == n_local_rows);

  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedIndexFields;
  syncedIndexFields[PhysicalDomainT::FiniteElementElementManager].push_back(m_TrilinosIndexStr);

  partition.SetBufferSizes(syncedIndexFields, CommRegistry::fractureFlowSolver);
  partition.SynchronizeFields(syncedIndexFields, CommRegistry::fractureFlowSolver);


  SetSrcFluxBCs(domain);


  // create epetra map

  row_map = Teuchos::rcp(new Epetra_Map(n_global_rows * nvar, n_local_rows * nvar,0,m_epetra_comm));

  // set up sparsity graph

  sparsity = Teuchos::rcp(new Epetra_FECrsGraph(Copy,*row_map,0));
	
  // loop over fracture-fracture connections

  std::vector<int> dofIndex(2*nvar);

  for(localIndex i = 0; i < m_edge_connection.size(); i++) {

    if(m_edge_is_bc[i] < 0) {

      for(localIndex j = 0; j < m_edge_connection[i].size(); ++j) {

	num1 = m_edge_connection[i][j].first;
	num2 = m_edge_connection[i][j].second;

	if(num1 != num2) {

	  m = m_elemIndex[num1].first;
	  n = m_elemIndex[num1].second;

	  p = m_elemIndex[num2].first;
	  q = m_elemIndex[num2].second;

	  if((*m_trilinos_index[m])[n] != -1 && (*m_trilinos_index[p])[q] != -1) {
	    for(int k = 0; k < m_eqt.nvar; k++) {

	      dofIndex[k] = (*m_trilinos_index[m])[n] * m_eqt.nvar + k;
	      dofIndex[k+1] = (*m_trilinos_index[p])[q] * m_eqt.nvar + k;

	      sparsity->InsertGlobalIndices(dofIndex.size(),
					    &dofIndex.front(),
					    dofIndex.size(),
					    &dofIndex.front());

	    }

	  }

	}

      }

    }

  }

  if(m_doMFCoupling) {

    // loop over matrix-fracture-matrix connections
    for(localIndex i = 0; i < faceManager.m_numFaces; i++) {

      if(face_is_ghost[i] < 0) {

	num = m_faceToElem[i].size();

	if(num > 1) {

	  for(localIndex j = 0; j < m_face_connection[i].size(); j++) {

	    num1 = m_face_connection[i][j].first;
	    num2 = m_face_connection[i][j].second;

	    m = m_elemIndex[num1].first;
	    n = m_elemIndex[num1].second;

	    p = m_elemIndex[num2].first;
	    q = m_elemIndex[num2].second;

	    if((*m_trilinos_index[m])[n] != -1 && (*m_trilinos_index[p])[q] != -1) {

	      for(int k = 0; k < m_eqt.nvar; k++) {

		dofIndex[k] = (*m_trilinos_index[m])[n] * m_eqt.nvar + k;
		dofIndex[k+1] = (*m_trilinos_index[p])[q] * m_eqt.nvar + k;

		sparsity->InsertGlobalIndices(dofIndex.size(),
					      &dofIndex.front(),
					      dofIndex.size(),
					      &dofIndex.front());

	      }

	    }
	  }

	}

      }

    }

  }

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  matrix   = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*sparsity));
  solution = Teuchos::rcp(new Epetra_FEVector(*row_map));
  rhs      = Teuchos::rcp(new Epetra_FEVector(*row_map));

  partition.SetBufferSizes(syncedFields, CommRegistry::fractureFlowSolver);

}

void FractureFlowSolver::SetDirichletElementBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time) {

  if(object.GetObjectType() == ObjectDataStructureBaseT::ElementRegion && bcBase->m_objectKey == PhysicalDomainT::FiniteElementElementRegion) {

    MultiVarDirichletBoundaryCondition* bc = dynamic_cast<MultiVarDirichletBoundaryCondition*> (bcBase);

    if(bc) {

      bc->checkVars(m_fieldName);

      const std::map<std::string, ElementRegionT>::const_iterator it = domain.m_feElementManager.m_ElementRegions.find(bc->m_regionName);

      if(it == domain.m_feElementManager.m_ElementRegions.end())
	throw GPException("BC region " + bc->m_regionName + " is not found!");

      localIndex m =m_elemRegionIndex[(it->second).m_regionNumber];

      lSet::const_iterator kf=set.begin();
    
      for(localIndex i =0; i < set.size() ; ++i, ++kf) {

	(*m_trilinos_index[m])[*kf] = -1;

      }

    }

  }

}

void FractureFlowSolver::SetDirichletEdgeBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time) {

  if(object.GetObjectType() == ObjectDataStructureBaseT::EdgeManager && bcBase->m_objectKey == PhysicalDomainT::FiniteElementEdgeManager) {

    MultiVarDirichletBoundaryCondition* bc = dynamic_cast<MultiVarDirichletBoundaryCondition*> (bcBase);

    if(bc) {
      bc->checkVars(m_fieldName);

      iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

      lSet::const_iterator kf=set.begin();
    
      for(localIndex i =0; i < set.size() ; ++i, ++kf) {

	if(edge_is_ghost[*kf] < 0) {

	  if(isMember(*kf, m_edgeIndex)) {

	    m_edge_is_bc[m_edgeIndex[*kf] ] = m_bcEdgeCount++;

	  }
	  else  
	    throw GPException("Boundary edge is not joined by a fracture face");
	}

      }

    }
  }

}


void FractureFlowSolver::SetSrcFluxBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time)
{

  realT value = 0.0;
  realT sum, w;
  int num, index;
  localIndex m, n, k;

  MultiVarSrcFluxBoundaryCondition* bc = dynamic_cast<MultiVarSrcFluxBoundaryCondition*> (bcBase);

  if(bc) {

      bc->checkVars(m_fieldName);

      if(bc->isAllocedByWeight()) {

	if(bc->m_objectKey == PhysicalDomainT::FiniteElementFaceManager && m_doMFCoupling) {

	  const iArray1d& face_is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

	  lSet::const_iterator kf=set.begin();

	  for( localIndex i =0; i < set.size() ; ++i, ++kf ){

	    if(face_is_ghost[*kf] < 0) {

	      num = m_face_connection[*kf].size();

	      if(num == 1) {

		k = m_face_connection[*kf][0].first;

		m = m_elemIndex[k].first;
		n = m_elemIndex[k].second;

		if((*m_trilinos_index[m])[n] >= 0) {

		  w = domain.m_feFaceManager.SurfaceArea(domain.m_feNodeManager, *kf);
		  value += w;

		  bc->allocFac().push_back(w);

		}

	      }
	      else  
		throw GPException("Error: SrcFlux BC face is not on external boundary");

	    }

	  }

	}
	else if(bc->m_objectKey == PhysicalDomainT::FiniteElementEdgeManager) {

	  const iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

	  lSet::const_iterator eg=set.begin();

	  for( localIndex i =0; i < set.size() ; ++i, ++eg ){

	    if(edge_is_ghost[*eg] < 0) {

	      if(isMember(*eg, m_edgeIndex)) {

		k = m_edgeIndex[*eg];
		num = m_edge_connection[k].size();

		if(num == 1) {

		  index = m_edge_connection[k][0].first;
		  m = m_elemIndex[index].first;
		  n = m_elemIndex[index].second;

		  if((*m_trilinos_index[m])[n] >= 0) {

		    w = domain.m_feEdgeManager.EdgeLength(domain.m_feNodeManager, *eg);
		    value += w;

		    bc->allocFac().push_back(w);

		  }

		}
		else 
		  throw GPException("Error: SrcFlux BC fracture edge is not on external boundary");

	      }
	      else 
		throw GPException("Error: SrcFlux BC edge is not joined by a fracture face");

	    }

	  }

	}
	else if(bc->m_objectKey == PhysicalDomainT::FiniteElementElementRegion) {

	  const std::map<std::string, ElementRegionT>::const_iterator it = domain.m_feElementManager.m_ElementRegions.find(bc->m_regionName);

	  if(it == domain.m_feElementManager.m_ElementRegions.end())
	    throw GPException("BC region " + bc->m_regionName + " is not found!");

	  m =m_elemRegionIndex[(it->second).m_regionNumber];

	  lSet::const_iterator elem=set.begin();

	  for( localIndex i =0; i < set.size() ; ++i, ++elem ){

	    if((*m_elem_is_ghost[m])[*elem] < 0) {

	      if((*m_trilinos_index[m])[*elem] >= 0) {

		value += m_volume[m][*elem];
		bc->allocFac().push_back(m_volume[m][*elem]);

	      }

	    }

	  }

	}

	MPI_Allreduce(&value, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	for(localIndex i = 0; i < bc->allocFac().size(); ++i) {

	  bc->allocFac()[i] /= sum;

	}

      }

  }

}

void FractureFlowSolver::SrcFluxBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time)
{

  realT q;
 
  localIndex m, n, k, p, index;

  int nvar = m_eqt.nvar;
  int ncomp = m_eqt.ncomp;

  Epetra_IntSerialDenseVector  elem_index(nvar);
  Epetra_SerialDenseVector     elem_rhs(nvar);


  MultiVarSrcFluxBoundaryCondition* bc = dynamic_cast<MultiVarSrcFluxBoundaryCondition*> (bcBase);

  if(bc) {

    const rArray1d &fluxValue = bc->GetValues(time);

    if(bc->m_objectKey == PhysicalDomainT::FiniteElementFaceManager && m_doMFCoupling) {

      const iArray1d& face_is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

      lSet::const_iterator kf=set.begin();

      index = 0;

      for( localIndex i =0; i < set.size() ; ++i, ++kf ){

	if(face_is_ghost[*kf] < 0) {

	  k = m_face_connection[*kf][0].first;

	  m = m_elemIndex[k].first;
	  n = m_elemIndex[k].second;

	  if((*m_trilinos_index[m])[n] >= 0) {

	    elem_index[0] = (*m_trilinos_index[m])[n] * nvar;
	    q = fluxValue[0];

	    if(bc->isAllocedByWeight())
	      q = fluxValue[0] * bc->allocFac()[index++];

	    elem_rhs[0] = -m_dt * q;

	    if(m_eqt.thermal) {

	      elem_index[ncomp] = (*m_trilinos_index[m])[n] * nvar + ncomp;
	      elem_rhs[ncomp] = -m_dt * q * fluxValue[ncomp];

	    }
	  
	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	  }

	}

      }

    }
    else if(bc->m_objectKey == PhysicalDomainT::FiniteElementEdgeManager) {

      const iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

      lSet::const_iterator eg=set.begin();

      index = 0;

      for( localIndex i =0; i < set.size() ; ++i, ++eg ){

	if(edge_is_ghost[*eg] < 0) {

	  k = m_edgeIndex[*eg];

	  p = m_edge_connection[k][0].first;
	  m = m_elemIndex[p].first;
	  n = m_elemIndex[p].second;

	  if((*m_trilinos_index[m])[n] >= 0) {

	    elem_index[0] = (*m_trilinos_index[m])[n] * nvar;
	    q = fluxValue[0];

	    if(bc->isAllocedByWeight())
	      q = fluxValue[0] * bc->allocFac()[index++];

	    elem_rhs[0] = -m_dt * q;

	    if(m_eqt.thermal) {

	      elem_index[ncomp] = (*m_trilinos_index[m])[n] * nvar + ncomp;
	      elem_rhs[ncomp] = -m_dt * q * fluxValue[ncomp];

	    }

	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	  }
	}

      }

    }
    else if(bc->m_objectKey == PhysicalDomainT::FiniteElementElementRegion) {

      const std::map<std::string, ElementRegionT>::const_iterator it = domain.m_feElementManager.m_ElementRegions.find(bc->m_regionName);

      m =m_elemRegionIndex[(it->second).m_regionNumber];

      lSet::const_iterator elem=set.begin();

      index = 0;

      for( localIndex i =0; i < set.size() ; ++i, ++elem ){

	if((*m_elem_is_ghost[m])[*elem] < 0) {

	  if((*m_trilinos_index[m])[*elem] >= 0) {

	    elem_index[0] = (*m_trilinos_index[m])[*elem] * nvar;
	    q = fluxValue[0];

	    if(bc->isAllocedByWeight())
	      q = fluxValue[0] * bc->allocFac()[index++];

	    elem_rhs[0] = -m_dt * q;

	    if(m_eqt.thermal) {

	      elem_index[ncomp] = (*m_trilinos_index[m])[*elem] * nvar + ncomp;
	      elem_rhs[ncomp] = -m_dt * q * fluxValue[ncomp];

	    }

	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	  }
	}

      }

    }

  }

}

void FractureFlowSolver::DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set, realT time){
 
  int nvar = m_eqt.nvar;
  int ncomp = m_eqt.ncomp;

  MultiVarDirichletBoundaryCondition* bc = dynamic_cast<MultiVarDirichletBoundaryCondition*> (bcBase);

  if(bc) {

    if(!bc->isClamped()) {

      const rArray1d &value = bc->GetValues(time);

      if(object.GetObjectType() == ObjectDataStructureBaseT::EdgeManager && bc->m_objectKey == PhysicalDomainT::FiniteElementEdgeManager) {

	const iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

	lSet::const_iterator eg=set.begin();
    
	for(localIndex i =0; i < set.size() ; ++i, ++eg) {

	  if(edge_is_ghost[*eg] < 0) {

	    for(int j = 0; j < nvar; j++) {
	      m_bcEdgeValue(j, m_bcEdgeCount) = value[j];
	    }

	    m_bcEdgeDensityValue[m_bcEdgeCount] = m_fluidRockProperty.m_rho_o * exp(m_fluidRockProperty.m_compress * (m_bcEdgeValue(0, m_bcEdgeCount) - m_fluidRockProperty.m_press_o));

	    m_bcEdgeCount++;

	  }
	}
      }
      else if(object.GetObjectType() == ObjectDataStructureBaseT::ElementRegion && bc->m_objectKey == PhysicalDomainT::FiniteElementElementRegion) {

	const std::map<std::string, ElementRegionT>::const_iterator it = domain.m_feElementManager.m_ElementRegions.find(bc->m_regionName);

	localIndex m =m_elemRegionIndex[(it->second).m_regionNumber];

	lSet::const_iterator kf=set.begin();
    
	for(localIndex i =0; i < set.size() ; ++i, ++kf) {

	  (*m_pressure[m])[*kf] = value[0];

	  if(m_eqt.thermal)   
	    (*m_temperature[m])[*kf] = value[ncomp];

	}

      }

    }

  }

}


void FractureFlowSolver::Assemble(PhysicalDomainT&  domain,
				  SpatialPartition& partition,
				  const realT& time)
{
  
  matrix->PutScalar(0.0);
  rhs->PutScalar(0.0);
  solution->PutScalar(0.0);

  int nvar = m_eqt.nvar;

  Epetra_IntSerialDenseVector  elem_index (nvar);
  Epetra_SerialDenseVector     elem_rhs     (nvar);
  Epetra_SerialDenseMatrix     elem_matrix  (nvar,nvar);

  int num = 2*nvar;

  Epetra_IntSerialDenseVector  face_index(num);
  Epetra_SerialDenseVector     face_rhs(num);
  Epetra_SerialDenseMatrix     face_matrix(num, num);

  int m, n, p, q, n1, n2;
  int index;
  realT pv, rho, rhoav;
  realT term0, term1, term2, term3, term4, term5;

  realT m_compress = m_fluidRockProperty.m_compress;
  realT m_mu = m_fluidRockProperty.m_mu;

  // loop over elements (accumulation term)

  // fracture elements

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    for(localIndex j=0; j< m_elemRegion[i]->m_numElems; j++) {

      if((*m_elem_is_ghost[i])[j] < 0 && (*m_trilinos_index[i])[j] >= 0) {

	elem_matrix(0,0) = m_compress * m_density[i][j] * m_fractureVolume[i][j];
	elem_index(0) = (*m_trilinos_index[i])[j] * nvar;

	elem_rhs(0) =  (m_density[i][j] - (m_oldAcc[0])[i][j]) * m_fractureVolume[i][j]; 

	matrix->SumIntoGlobalValues(elem_index, elem_matrix);
	rhs->SumIntoGlobalValues(elem_index, elem_rhs);

      }

    }

  }

  // loop over fracture connections 

  for(localIndex i = 0; i < m_edge_connection.size(); ++i) {

    index = m_edge_is_bc[i];

    if(!(m_edge_connection[i][0].first == m_edge_connection[i][0].second && index < 0)) {

      for(localIndex j = 0; j < m_edge_connection[i].size(); ++j) {

	n1 = m_edge_connection[i][j].first;

	m = m_elemIndex[n1].first;
	n = m_elemIndex[n1].second;

	if(index >= 0) {

	  pv = m_bcEdgeValue[0][index];
	  rho = m_bcEdgeDensityValue[index];

	}
	else {

	  n2 = m_edge_connection[i][j].second;
	  p = m_elemIndex[n2].first;
	  q = m_elemIndex[n2].second;

	  pv = (*m_pressure[p])[q];
	  rho = m_density[p][q];

	}

	rhoav = m_edge_wt1[i][j] * m_density[m][n] + m_edge_wt2[i][j] * rho;


	term0 = (*m_pressure[m])[n] - pv - rhoav * m_edge_dz[i][j];

	term1 = rhoav * m_edge_trans[i][j] / m_mu * m_dt;
 
	term2 = term0 * m_edge_wt1[i][j] * m_density[m][n] * m_compress * m_edge_trans[i][j] / m_mu * m_dt;
	term4 = m_edge_dz[i][j] * m_edge_wt1[i][j] * m_density[m][n] * m_compress;

	if(index < 0) {

	  term3 = term0 * m_edge_wt2[i][j] * m_density[p][q] * m_compress * m_edge_trans[i][j] / m_mu * m_dt;

	  term5 = m_edge_dz[i][j] * m_edge_wt1[i][j] * m_density[p][q] * m_compress;

	  if((*m_trilinos_index[m])[n] >= 0 && (*m_trilinos_index[p])[q] >= 0) {

	    face_matrix(0,0) = (1.0 - term4) * term1 + term2;
	    face_matrix(1,1) = (1.0 + term5) * term1 - term3;
  	    	 
	    face_matrix(0,1) = (-1.0 - term5) * term1 + term3;
	    face_matrix(1,0) = (-1.0 + term4) * term1 - term2;
        
	    face_index[0] = (*m_trilinos_index[m])[n] * nvar;
	    face_index[1] = (*m_trilinos_index[p])[q] * nvar;
  	 
	    matrix->SumIntoGlobalValues(face_index, face_matrix);

	    elem_index(0) = (*m_trilinos_index[m])[n];
	    elem_rhs(0) =  term1 * term0;
	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	    elem_index(0) = (*m_trilinos_index[p])[q];
	    elem_rhs(0) = -term1 * term0;
	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	  } 
	  else if((*m_trilinos_index[m])[n] >= 0) {

	    elem_matrix(0,0) = (1.0 - term4) * term1 + term2;
	    elem_index[0] = (*m_trilinos_index[m])[n];

	    matrix->SumIntoGlobalValues(elem_index, elem_matrix);

	    elem_rhs(0) =  term1 * term0;
	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);


	  } 
	  else if((*m_trilinos_index[p])[q] >= 0) {

	    elem_matrix(0,0) = (1.0 + term5) * term1 - term3;
	    elem_index[0] = (*m_trilinos_index[p])[q];

	    matrix->SumIntoGlobalValues(elem_index, elem_matrix);

	    elem_rhs(0) =  -term1 * term0;
	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	  }
	}
	else if((*m_trilinos_index[m])[n] >= 0) {

	    elem_matrix(0,0) = (1.0 - term4) * term1 + term2;
	    elem_index[0] = (*m_trilinos_index[m])[n];

	    matrix->SumIntoGlobalValues(elem_index, elem_matrix);

	    elem_rhs(0) =  term1 * term0;
	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	}
      }

    }

  }

  if(m_doMFCoupling) {

    // matrix elements

    for(int i = 0; i < m_matrix_regionNumber; i++) {

      index = i + m_fracture_regionNumber;

      for(localIndex j=0; j< m_elemRegion[index]->m_numElems; j++) {

	if((*m_elem_is_ghost[index])[j] < 0 && (*m_trilinos_index[index])[j] >= 0) {

	  elem_matrix(0,0) = (m_fluidRockProperty.m_rockCompress + m_compress )* m_density[index][j] * m_volume[index][j] * (*m_porosity[i])[j];

	  elem_index(0) = (*m_trilinos_index[index])[j];

	  matrix->SumIntoGlobalValues(elem_index, elem_matrix);

	  elem_rhs(0) =  (m_density[index][j] * (*m_porosity[i])[j] - (m_oldAcc[0])[index][j]) * m_volume[index][j]; 

	  rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	}

      }
    }

    const iArray1d& face_is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

    for(localIndex i = 0; i < domain.m_feFaceManager.m_numFaces; i++) {

      if(face_is_ghost[i] < 0) {

	if(m_face_connection[i][0].first != m_face_connection[i][0].second) {

	  for(localIndex j = 0; j < m_face_connection[i].size(); j++) {

	    n1 = m_face_connection[i][j].first;
	    n2 = m_face_connection[i][j].second;

	    m = m_elemIndex[n1].first;
	    n = m_elemIndex[n1].second;

	    p = m_elemIndex[n2].first;
	    q = m_elemIndex[n2].second;

	    rhoav = m_face_wt[i][j] * m_density[m][n] + (1.0 - m_face_wt[i][j]) * m_density[p][q];

	    term0 = (*m_pressure[m])[n] - (*m_pressure[p])[q] - rhoav * m_face_dz[i][j];

	    term1 = rhoav * m_face_trans[i][j] / m_mu * m_dt;
 
	    term2 = term0 * m_face_wt[i][j] * m_density[m][n]  * m_compress * m_face_trans[i][j] / m_mu * m_dt;
	    term3 = term0 * (1.0 - m_face_wt[i][j]) * m_density[p][q]  * m_compress * m_face_trans[i][j] / m_mu * m_dt;

	    term4 = m_face_dz[i][j] * m_face_wt[i][j] * m_density[m][n] * m_compress;
	    term5 = m_face_dz[i][j] * (1.0 - m_face_wt[i][j]) * m_density[p][q] * m_compress;

	    if((*m_trilinos_index[m])[n] >= 0 && (*m_trilinos_index[p])[q] >= 0) {

	      face_matrix(0,0) = (1.0 - term4) * term1 + term2;
	      face_matrix(1,1) = (1.0 + term5) * term1 - term3;
  	    	 
	      face_matrix(0,1) = (-1.0 - term5) * term1 + term3;
	      face_matrix(1,0) = (-1.0 + term4) * term1 - term2;
  	  
	      face_index[0] = (*m_trilinos_index[m])[n];
	      face_index[1] = (*m_trilinos_index[p])[q];
 
	      matrix->SumIntoGlobalValues(face_index, face_matrix);

	      elem_index(0) = (*m_trilinos_index[m])[n];
	      elem_rhs(0) =  term1 * term0; 
	      rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	      elem_index(0) = (*m_trilinos_index[p])[q];
	      elem_rhs(0) = -term1 * term0;
	      rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	    }
	    else if((*m_trilinos_index[m])[n] >= 0) {

	      elem_index(0) = (*m_trilinos_index[m])[n];
	      elem_matrix(0,0) = (1.0 - term4) * term1 + term2;

	      matrix->SumIntoGlobalValues(elem_index, elem_matrix);

	      elem_rhs(0) =  term1 * term0; 
	      rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	    }
	    else if((*m_trilinos_index[p])[q] >= 0) {

	      elem_index(0) = (*m_trilinos_index[p])[q];
	      elem_matrix(0,0) = (1.0 + term5) * term1 - term3;

	      matrix->SumIntoGlobalValues(elem_index, elem_matrix);

	      elem_rhs(0) =  -term1 * term0; 
	      rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	    }

	  }

	}

      }

    }

  }

  //update srcflux boundary conditions
  UpdateSrcFluxBCs(domain, time);

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

}

realT FractureFlowSolver::Solve(PhysicalDomainT&  domain,
				SpatialPartition& partition)
{

  realT out, tmp;

  int index;
  double* local_solution = NULL;

  solution->ExtractView(&local_solution,&index);


  Epetra_LinearProblem problem(&(*matrix),
                               &(*solution),
                               &(*rhs));

  Teuchos::ParameterList MLList;

  ML_Epetra::SetDefaults("SA",MLList);

  ML_Epetra::MultiLevelPreconditioner* MLPrec = NULL;
  if(m_useMLPrecond){
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*matrix, MLList);
  }


  // Later we need to test performance of different solvers/preconditioners 
  // for fracture flow problems 

  AztecOO solver(problem);

  if(m_useMLPrecond) solver.SetPrecOperator(MLPrec);

  solver.SetAztecOption(AZ_solver,AZ_bicgstab);
  solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
  solver.SetAztecOption(AZ_conv,AZ_rhs);

  solver.Iterate(m_numerics.m_maxIters,m_numerics.m_tol);

  if(m_useMLPrecond){
     delete MLPrec;
  }
  
  solution->MinValue(&out);
  if(out < 0.0) out = -out;

  solution->MaxValue(&tmp);
  if(tmp < 0.0) tmp = -tmp;

  if(out < tmp)
    out = tmp;

  long long int rowTmp;
  int lid;

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    for(localIndex j = 0; j < m_elemRegion[i]->m_numElems; j++) {

      if((*m_elem_is_ghost[i])[j] < 0 && (*m_trilinos_index[i])[j] >= 0) {

	rowTmp = static_cast<long long int>((*m_trilinos_index[i])[j]);
	lid = row_map->LID(rowTmp);
	(*m_pressure[i])[j] -= local_solution[lid];

      }

    }
  }

  if(m_doMFCoupling) {

    //update matrix porosity
    for(int i = 0; i < m_matrix_regionNumber; i++) {

      index = i + m_fracture_regionNumber;

      for(localIndex j = 0; j < m_elemRegion[index]->m_numElems; j++) {

	if((*m_elem_is_ghost[index])[j] < 0 && (*m_trilinos_index[index])[j] >= 0) {

	  rowTmp = static_cast<long long int>((*m_trilinos_index[index])[j]);
	  lid = row_map->LID(rowTmp);
	  (*m_pressure[index])[j] -= local_solution[lid];

	  tmp = (*m_porosity[i])[j] * (1.0 - m_fluidRockProperty.m_rockCompress * local_solution[lid]);

	  if(tmp > m_fluidRockProperty.m_maxPorosity)
	    tmp = m_fluidRockProperty.m_maxPorosity;

	  if(tmp < m_fluidRockProperty.m_minPorosity)
	    tmp = m_fluidRockProperty.m_minPorosity;

	  (*m_porosity[i])[j] = tmp;

	}

      }
    }

  }

  // re-sync ghost nodes
  partition.SynchronizeFields(syncedFields, CommRegistry::fractureFlowSolver);

  return out;

}

void FractureFlowSolver::TimeStep(const realT& time,
				  const realT& dt,
				  const int cycleNumber,
				  PhysicalDomainT& domain,
				  const sArray1d& namesOfSolverRegions,
				  SpatialPartition& partition,
				  FractunatorBase* const fractunator)
{

  m_dt = dt;

  if(m_doApertureUpdate)
    MakeTransParameters(domain);

  UpdateDirichletBCs(domain, time);

  UpdateFluidRockProps(domain);

  MakeOldAcc();

  bool converge = 0;
  realT out;

  for(int i = 0; i < m_NRNumerics.m_maxIters; i++) {

    Assemble(domain,partition, time);
    out = Solve(domain,partition);
    if(out < m_NRNumerics.m_tol) {
      converge = 1;
      break;
    }

    UpdateFluidRockProps(domain);

  }

  if(!converge) {

    printf("Max. Residual = %e\n", out);
    throw GPException("Error: hit max. NR iterations -- need timestep cutback (not implemented: ");
  }

}

void FractureFlowSolver::MakeConnections(PhysicalDomainT& domain)
{
  int num;
  iArray1d a;

  m_edge_connection.resize(m_edgeToFace.size());

  for(localIndex i = 0; i < m_edgeToFace.size(); ++i) {

    num = m_edgeToFace[i].size();
    a.resize(num);

    for(int j = 0; j < num; ++j) {

      a[j] = m_edgeToFace[i][j];

    }

    if(m_edge_is_bc[i] >= 0 || num == 1) {

      for(int j = 0; j < num; ++j) {

	m_edge_connection[i].push_back(std::make_pair(a[j], a[j]));

      }

    }
    else if(num > 1) {

      for(int j = 0; j < num -1; ++j)
	for(int k = j + 1; k < num; ++k) {

	  m_edge_connection[i].push_back(std::make_pair(a[j], a[k]));

	}
    }

  }

  if(m_doMFCoupling) {

    // matrix elements
    FaceManagerT& faceManager = domain.m_feFaceManager;

    iArray1d& face_is_ghost       = faceManager.GetFieldData<FieldInfo::ghostRank>();

    m_face_connection.resize(faceManager.m_numFaces);

    for(localIndex i = 0; i < faceManager.m_numFaces; i++) {

      if(face_is_ghost[i] < 0) {

	num = m_faceToElem[i].size();
	a.resize(num);

	for(int j = 0; j < num; ++j) {

	  a[j] = m_faceToElem[i][j];
   
	}

	if(num == 1) {

	  m_face_connection[i].push_back(std::make_pair(a[0], a[0]));

	}
	else {

	  for(int j = 0; j < num - 1; j++) {

	    m_face_connection[i].push_back(std::make_pair(a[j], a[num-1]));

	  }

	}

      }

    }

  }

}

void FractureFlowSolver::MakeTransParameters(PhysicalDomainT& domain)
{

  int m, n, index;
  R1Tensor la, lb;
  realT app;

  Array1dT<rArray1d *> aperture(m_fracture_regionNumber);

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];
    aperture[i] = &(elemRegion->GetFieldData<realT>(ApertureStr));

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {

      if((*aperture[i])[j] > m_fluidRockProperty.m_maxAperture)  
	(*aperture[i])[j] = m_fluidRockProperty.m_maxAperture;

      if((*aperture[i])[j] < m_fluidRockProperty.m_minAperture)  
	(*aperture[i])[j] = m_fluidRockProperty.m_minAperture;

    }

  }

  // fracture elements

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {

      if((*m_elem_is_ghost[i])[j] < 0) {

	m_fractureVolume[i][j] = m_volume[i][j] * (*aperture[i])[j];
       
      }
    }

  }

  rArray1d a, b;
  realT tmp, sum;
  int num;

  for(localIndex i = 0; i < m_edge_connection.size(); ++i) {

    num = m_edgeToFace[i].size();
    a.resize(num);
    b.resize(num);

    for(int j = 0; j < num; ++j) {

      index = m_edgeToFace[i][j];
      m = m_elemIndex[index].first;
      n = m_elemIndex[index].second;

      app = (*aperture[m])[n];

      a[j] =  m_edgeToFaceAL[i][j] * app * app * app;
      b[j] =  m_edgeToFaceAL[i][j];

    }

    if(m_edge_is_bc[i] >= 0 || num == 1) {

      for(int j = 0; j < num; ++j) {

	m_edge_trans[i][j] = a[j];
	m_edge_wt1[i][j] = 0.5;
	m_edge_wt2[i][j] = 0.5;

      }

    }
    else if(num > 1) {

      tmp = 0.0;
      sum = 0.0;

      for(int j = 0; j < num; ++j) {
	tmp += a[j];
	sum += b[j];
      }

      index = 0;
      for(int j = 0; j < num - 1; ++j)
	for(int k = j + 1; k < num; ++k) {

	  m_edge_trans[i][index] = a[j] * a[k] / tmp;
	  m_edge_wt1[i][index] = b[j] / sum;
	  m_edge_wt2[i][index] = b[k] / sum;
	  index++;

	}
    }

  }

  if(m_doMFCoupling) {

    FaceManagerT& faceManager = domain.m_feFaceManager;

    // matrix elements

    Array1dT<rArray1d *> permeability(m_matrix_regionNumber);

    for(int i = 0; i < m_matrix_regionNumber; i++) {

      index = i + m_fracture_regionNumber;
      ElementRegionT* elemRegion = m_elemRegion[index];
      permeability[i] = &(elemRegion->GetFieldData<realT>(PermeabilityStr));

    }

    const iArray1d& face_is_ghost       = faceManager.GetFieldData<FieldInfo::ghostRank>();

    for(localIndex i = 0; i < faceManager.m_numFaces; i++) {

      if(face_is_ghost[i] < 0) {

	num = m_faceToElem[i].size();
	a.resize(num);
	b.resize(num);

	for(int j = 0; j < num; ++j) {

	  index = m_faceToElem[i][j];
   
	  m = m_elemIndex[index].first;
	  n = m_elemIndex[index].second;

	  if(m < m_fracture_regionNumber) {

	    app = (*aperture[m])[n];
	    a[j] =  m_faceToElemAL[i][j] * app;
            b[j] =  m_faceToElemAL[i][j] * 12.0 / app;
	  }
	  else {

	    a[j] =  m_faceToElemAL[i][j] * (*permeability[m - m_fracture_regionNumber])[j];
	    b[j] =  m_faceToElemAL[i][j];
	  }

	}

	if(num == 1) {

	  m_face_trans[i][0] = a[0];
	  m_face_wt[i][0] = 1.0;

	}
	else {

	  for(int j = 0; j < num - 1; j++) {

	    m_face_trans[i][j] = a[j] * a[num-1] / (a[j] + a[num-1]);
	    m_face_wt[i][j] = b[j] / (b[j] + b[num-1]);

	  }

	}

      }

    }

  }

}

void FractureFlowSolver::MakeGeometryParameters(PhysicalDomainT& domain) {

  R1Tensor la, lb;
  localIndex index;
  realT w;
  int m, n, kf, eg;
  localIndex num;

  FaceManagerT& faceManager = domain.m_feFaceManager;
  NodeManagerT& nodeManager = domain.m_feNodeManager;
  EdgeManagerT& edgeManager = domain.m_feEdgeManager;

  Array1dT<R1Tensor> faceCenter(faceManager.m_numFaces);

  // fracture elements
  for(int i = 0; i < m_fracture_regionNumber; i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j=0; j < elemRegion->m_numElems; j++) {

      kf = elemRegion->m_toFacesRelation[j][0];
      faceManager.FaceCenter(nodeManager, kf, faceCenter[kf]);

      if((*m_elem_is_ghost[i])[j] < 0) {

	m_volume[i][j] = faceManager.SurfaceArea(nodeManager, kf);
       
      }
    }
  }

  rArray1d a;

  for(std::map<localIndex,localIndex>::const_iterator edge = m_edgeIndex.begin(); edge != m_edgeIndex.end(); ++edge) {

    eg = edge->first;

    index = edge->second;

    edgeManager.EdgeCenter(nodeManager, eg, la);
    w = domain.m_feEdgeManager.EdgeLength(nodeManager, eg);

    num = m_edgeToFace[index].size();
    a.resize(num);

    for(localIndex i = 0; i < num; ++i) {

      m = m_elemIndex[m_edgeToFace[index][i] ].first;
      n = m_elemIndex[m_edgeToFace[index][i] ].second;

      kf = m_elemRegion[m]->m_toFacesRelation[n][0];

      lb = faceCenter[kf] - la;

      m_edgeToFaceAL[index][i] =  w / 12.0 / lb.L2_Norm(); // ?? * cos(theta)

      a[i] =  Dot(lb, m_downVector) * m_gravFactor;


    }

    if(m_edge_is_bc[index] >= 0 || num == 1) {

      for(localIndex j = 0; j < num; ++j) {

	m_edge_dz[index][j] = a[j];

      }

    }
    else if(num > 1) {

      m = 0;      

      for(localIndex j = 0; j < num -1; ++j)
	for(localIndex k = j + 1; k < num; ++k) {

	  m_edge_dz[index][m] = a[j] - a[k];
	  m++;

	}
    }

  }

  if(m_doMFCoupling) {

    // matrix elements
    const iArray1d& face_is_ghost       = faceManager.GetFieldData<FieldInfo::ghostRank>();

    const Array1dT<R1Tensor>& nodePos = nodeManager.GetFieldData<FieldInfo::referencePosition>();

    Array1dT<Array1dT<R1Tensor> > elemCenter(m_matrix_regionNumber);
    localIndex faceIndex, nnodes;
    R1Tensor dummy;

    for(int i = 0; i < m_matrix_regionNumber; i++) {

      index = i + m_fracture_regionNumber;

      ElementRegionT* elemRegion = m_elemRegion[index];

      for(localIndex j=0; j < elemRegion->m_numElems; j++) {

        elemCenter[i].push_back(elemRegion->GetElementCenter(j, nodeManager));

      }

      for(localIndex j=0; j< elemRegion->m_numElems; j++) {

	if((*m_elem_is_ghost[index])[j] < 0) {

	  m_volume[index][j] = 0.0;

	  for(localIndex k = 0; k < elemRegion->m_toFacesRelation.Dimension(1); k++) {
	    faceIndex = elemRegion->m_toFacesRelation(j, k);
	    nnodes = faceManager.m_toNodesRelation[faceIndex].size();

	    R1Tensor x0(nodePos[faceManager.m_toNodesRelation[faceIndex][0]]);

	    for (localIndex l = 2; l < nnodes; l++) {
	      R1Tensor x1(nodePos[faceManager.m_toNodesRelation[faceIndex][l - 1]]);
	      R1Tensor x2(nodePos[faceManager.m_toNodesRelation[faceIndex][l]]);

	      m_volume[index][j] += GeometryUtilities::CentroidAndVolume_3DTetrahedron(elemCenter[i][j], x0, x1, x2, dummy);

	    }

	  }

	}

      }

    }


    // loop over faces

    realT faceArea;
    R1Tensor fCenter, fNormal;

    for(localIndex i = 0; i < faceManager.m_numFaces; i++) {

      if(face_is_ghost[i] < 0) {

	num = faceManager.m_toElementsRelation[i].size();
	faceArea = faceManager.FaceCenterAndNormal(nodeManager, i, fCenter, fNormal);

	a.resize(m_faceToElem[i].size());

	index = 0;

	for(localIndex j=0; j < num; j++) {



	  m = m_elemRegionIndex[(elem.first)->m_regionNumber] - m_fracture_regionNumber;

	  if(m < 0)
	    continue;

	  n = elem.second;

	  la = elemCenter[m][n] - fCenter;

	  m_faceToElemAL[i][index] =  fabs(Dot(la.UnitVector(), fNormal)) / la.L2_Norm();
	  a[index] = Dot(la, m_downVector) * m_gravFactor;

	  index++;
	}

	if(m_faceToElem[i].size() > index) {

	  //the face is also a fracture

	  m_faceToElemAL[i][index] =  faceArea / 6.0; //A/12.0/(d / 2.0)

	  a[index] = 0.0;
	}

	num = a.size();

	if(num == 1) {

	  m_face_dz[i][0] = a[0];

	}
	else {

	  for(localIndex j = 0; j < num - 1; j++) {

	    m_face_dz[i][j] = a[j] - a[num-1];

	  }

	}
      }

    }

  }
}


void FractureFlowSolver::UpdateFluidRockProps(PhysicalDomainT& domain) {

  for(localIndex i = 0; i < m_elemRegion.size(); ++i)

      for(localIndex j = 0; j < m_elemRegion[i]->m_numElems; ++j) {

	m_density[i][j] = m_fluidRockProperty.m_rho_o * exp(m_fluidRockProperty.m_compress * ((*m_pressure[i])[j] - m_fluidRockProperty.m_press_o));

      }

}

void FractureFlowSolver::SetDirichletBCs(PhysicalDomainT& domain)
{

  for(localIndex i =0; i < m_elemRegion.size(); i++) {

    ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SetDirichletElementBoundaryCondition, domain, (*m_elemRegion[i]), "MultiVars", 0.0);

  }

  // The following approach is not a good way to handle edge boundary 
  // conditions, which needs to be replaced by a new object.

  m_edge_is_bc.resize(m_edgeToFace.size());
  m_edge_is_bc = -1;

  m_bcEdgeCount = 0;

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SetDirichletEdgeBoundaryCondition, domain, domain.m_feEdgeManager, "MultiVars", 0.0);

  m_bcEdgeValue.resize2(m_eqt.nvar, m_bcEdgeCount);
  m_bcEdgeDensityValue.resize(m_bcEdgeCount);

}


void FractureFlowSolver::SetSrcFluxBCs(PhysicalDomainT& domain)
{

  // flux only into fracture elements (edge-length-based)

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SetSrcFluxBoundaryCondition, domain, domain.m_feEdgeManager, "MultiVarSrcFlux", 0.0);                                       
  // flux only into boundary matrix elements (face-area-based)

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SetSrcFluxBoundaryCondition, domain, domain.m_feFaceManager, "MultiVarSrcFlux", 0.0);                                        
  for(localIndex i =0; i < m_elemRegion.size(); i++) {

    // flux into elements (element-volume-based) including fracture elements
 
    ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SetSrcFluxBoundaryCondition, domain, (*m_elemRegion[i]), "MultiVarSrcFlux", 0.0);                                        
  }

}


void FractureFlowSolver::UpdateSrcFluxBCs(PhysicalDomainT& domain, const realT &time)
{

  // loop over all possible srcflux boundary conditions

  // flux only into fracture elements (edge-length-based)

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SrcFluxBoundaryCondition, domain, domain.m_feEdgeManager, "MultiVarSrcFlux", 0.0);                                       
  // flux only into boundary matrix elements (face-area-based)

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SrcFluxBoundaryCondition, domain, domain.m_feFaceManager, "MultiVarSrcFlux", 0.0);                                        
  for(localIndex i =0; i < m_elemRegion.size(); i++) {

    // flux into elements (element-volume-based) including fracture elements
 
    ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SrcFluxBoundaryCondition, domain, (*m_elemRegion[i]), "MultiVarSrcFlux", 0.0);                                        

  }

}

void FractureFlowSolver::UpdateDirichletBCs(PhysicalDomainT& domain, const realT &time)
{

  for(localIndex i =0; i < m_elemRegion.size(); i++) {

    ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::DirichletBoundaryCondition, domain, (*m_elemRegion[i]), "MultiVars", 0.0);

  }

  // The following approach is not a good way to handle edge boundary 
  // conditions, which needs to be replaced by a new object.

  m_bcEdgeCount = 0;

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::DirichletBoundaryCondition, domain, domain.m_feEdgeManager, "MultiVars", 0.0);

}


void FractureFlowSolver::MakeOldAcc() {

  for(int i =0; i < m_fracture_regionNumber; ++i){

    for(localIndex j =0; j < m_elemRegion[i]->m_numElems; ++j){

      (m_oldAcc[0])[i][j] = m_density[i][j];

    }

  }

  if(m_doMFCoupling) {

    int index;

    for(int i =0; i < m_matrix_regionNumber; ++i){

      index = m_fracture_regionNumber + i;

      for(localIndex j =0; j < m_elemRegion[index]->m_numElems; ++j) {

	(m_oldAcc[0])[index][j] = m_density[index][j] * (*m_porosity[i])[j];

      }

    }

  }

}

/// Register solver in the solver factory
REGISTER_SOLVER( FractureFlowSolver )
