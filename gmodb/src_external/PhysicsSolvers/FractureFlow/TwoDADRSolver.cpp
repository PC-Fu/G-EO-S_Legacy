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
 * @file TwoDADRSolver.cpp
 * @author walsh24
 * @date June 1, 2011 
 */

#include "TwoDADRSolver.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"

#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/FunctionManager.h"
#include "ObjectManagers/FaceManagerT.h"
#include "ObjectManagers/EdgeManagerT.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"


// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

// Flow Solver
#include "PhysicsSolvers/ParallelPlateFlowSolverBase.h"

using namespace BoundaryConditionFunctions;
using namespace PS_STR;
using namespace TDSSADR_STR;
namespace{
  std::string NodeWeightsStr =  "NodeWeights";
  //std::string FaceAreaStr = "FaceAreas";
  //std::string FaceNormalStr = "FaceNormals";
}



namespace{

 R2Tensor Eye3d(1.0, 0.0,  0.0,
                0.0, 1.0,  0.0,
                0.0, 0.0,  1.0);
 
  const realT TINY = 1e-64;
  const realT SMALL = 1e-8;
}


/// Least-Squares Finite-element advection diffusion reaction solver
ParallelPlateADRSolver::ParallelPlateADRSolver(const std::string& name,
                                               ProblemManagerT* const pm ):
SolverBase(name,pm),
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
numerics(),
hasReaction(),
m_D(),
concentrationFieldName_(),
reactionRateFieldName_(),
rrDerivFieldName_()
{
  m_TrilinosIndexStr = "TwoDPPADRS_" +  toString<int>(m_instances) + "_GlobalDof";  
  ++m_instances;
}

ParallelPlateADRSolver::~ParallelPlateADRSolver()
{
  // TODO Auto-generated destructor stub
}

void ParallelPlateADRSolver::ReadXML(TICPP::HierarchicalDataNode* const hdn)
{
  SolverBase::ReadXML( hdn );

  realT D = hdn->GetAttributeOrDefault(DiffusivityStr,"1e-6m^2/s");
  m_D = Eye3d; m_D *= D;
  
  //m_tol = hdn->GetAttributeOrDefault("tol",1e-10);
 // m_maxSolverIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);
  //doDataWrite = hdn->GetAttributeOrDefault<bool>("doDataWrite",false);
  m_verboseFlag = hdn->GetAttributeOrDefault<bool>(VerboseStr,true);

  m_doSteadyState = hdn->GetAttributeOrDefault<bool>("steadyState",false);
  
  // face set
  m_faceSetName = hdn->GetAttributeString("faceset");
  
  // Fieldnames
  //concentrationFieldName_ = hdn->GetAttributeStringOrDefault("species",ConcentrationStr);
  //reactionRateFieldName_ = hdn->GetAttributeStringOrDefault("reactionRateField","");

  concentrationFieldName_ = ConcentrationStr;
  reactionRateFieldName_ = "";

  m_concentrationFieldNames = hdn->GetStringVector("species"); 
  m_reactionRateFieldNames = hdn->GetStringVector("reactionRateField");
  Trim(m_concentrationFieldNames," \t\n\r");
  Trim(m_reactionRateFieldNames," \t\n\r");

  if(m_concentrationFieldNames.size() > 0){
    concentrationFieldName_ = m_concentrationFieldNames[0];
    
    if(m_reactionRateFieldNames.size() == 0) 
        m_reactionRateFieldNames.resize(m_concentrationFieldNames.size() );
     

    if(m_reactionRateFieldNames.size() != m_concentrationFieldNames.size()) 
        throw GPException("Error ParallelPlateADRSolver: Number of reaction rate fields does not match number of species fields.");
    reactionRateFieldName_ = m_reactionRateFieldNames[0];
  } else {
    m_concentrationFieldNames.push_back(concentrationFieldName_);
    m_reactionRateFieldNames.push_back(reactionRateFieldName_);
  }
  
  if(reactionRateFieldName_!=""){
  	hasReaction = true;
    rrDerivFieldName_ = hdn->GetAttributeStringOrDefault("reactionRateDerivField",ReactionRateDerivStr);
  } else {
  	hasReaction = false;
  }
  
  numerics.krylov_tol     = hdn->GetAttributeOrDefault("tol",1.0e-10);

  numerics.useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);
  
  this->m_courant  = hdn->GetAttributeOrDefault("courant",1.0);

  // Fluid data - needed for Triple point advection
  // should really get these from flow solver
  m_mu = hdn->GetAttributeOrDefault("mu","1.0e-3 N.s/m^2");
  m_SHP_FCT = hdn->GetAttributeOrDefault<realT>("shapefactor",1.0);


}


void ParallelPlateADRSolver::RegisterFields( PhysicalDomainT& domain )
{
  
  for(unsigned i =0; i < m_concentrationFieldNames.size();++i){
    domain.m_feNodeManager.AddKeylessDataField<realT>( m_concentrationFieldNames[i], true, true );
    domain.m_feFaceManager.AddKeylessDataField<realT>( m_concentrationFieldNames[i], true, true );

    if(!(m_reactionRateFieldNames[i].empty())){
      domain.m_feFaceManager.AddKeylessDataField<realT>( m_reactionRateFieldNames[i], true, true );
      domain.m_feFaceManager.AddKeylessDataField<realT>( rrDerivFieldName_, true, true );
    }
  }
  
  domain.m_feNodeManager.AddKeylessDataField<realT>( concentrationFieldName_, true, true );
  domain.m_feFaceManager.AddKeylessDataField<realT>( concentrationFieldName_, true, true );
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>( FaceCenterStr, true, true );
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>( FaceNormalStr, true, true );
  domain.m_feFaceManager.AddKeylessDataField<realT>( ApertureStr, true, true );
  domain.m_feFaceManager.AddKeylessDataField<realT>( FaceAreaStr, true, true );
  
  domain.m_feEdgeManager.AddKeylessDataField<realT>( VolumetricFluxStr, true, true );
  
  
  if(hasReaction){
    domain.m_feFaceManager.AddKeylessDataField<realT>( reactionRateFieldName_, true, true );
    domain.m_feFaceManager.AddKeylessDataField<realT>( rrDerivFieldName_, true, true );
  }
  
  domain.m_feFaceManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);
  domain.m_feNodeManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);
     
}

void ParallelPlateADRSolver::SetConcentrationField(std::string& concentrationFieldName,std::string reactionRateFieldName,std::string rrDerivFieldName){
  concentrationFieldName_ = concentrationFieldName;
  if(reactionRateFieldName != ""){
  	hasReaction = true;
  	reactionRateFieldName_ = reactionRateFieldName;
  	rrDerivFieldName_ = rrDerivFieldName;
  } else {
  	hasReaction = false;
  }
}

void ParallelPlateADRSolver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{  
  m_faceSet = &(domain.m_feFaceManager.GetSet(m_faceSetName));
  m_nodeSet = &(domain.m_feNodeManager.GetSet(m_faceSetName));
  m_numFaces = m_faceSet->size();
  
  m_cumNodeCount = 0;
  
  // build face and node-dof map
  lSet::const_iterator fi=m_faceSet->begin();
  for(localIndex i =0; i < m_numFaces; ++i){
  	localIndex f = *(fi++);
  	m_faceDofMap[f] = i;
  	
  	const lArray1d& faceNodeMap = domain.m_feFaceManager.m_toNodesRelation[f];
    for( localIndex a=0 ; a<faceNodeMap.size(); ++a ){
      const localIndex nd = faceNodeMap[a];
  	  m_nodeDofMap[nd] = 0;
  	}
  	
  	m_cumNodeCount += faceNodeMap.size();
  }
  
  // edges to faces
  fi=m_faceSet->begin();
  for( localIndex i=0 ; i<m_numFaces ; ++i )
  {
    localIndex kf = *(fi++);
    localIndex numEdges = domain.m_feFaceManager.m_toEdgesRelation[kf].size();
    for(size_t a =0; a < numEdges; ++a){
      localIndex eg = domain.m_feFaceManager.m_toEdgesRelation[kf][a];
    	
      lSet& edgeFaces = domain.m_feEdgeManager.m_toFacesRelation[eg];
      lArray1d edgeList;
      
      for( lSet::iterator edgeFace=edgeFaces.begin() ; edgeFace!=edgeFaces.end() ; ++edgeFace ){
      	if(isMember(*edgeFace,m_faceDofMap)){
      	  edgeList.push_back(*edgeFace);
      	}
      }
      m_edgesToFaces[eg] = edgeList;
    }
  }
  
  std::map<localIndex,localIndex>::iterator itr;
  std::map<localIndex,localIndex>::iterator iEnd = m_nodeDofMap.end();
  localIndex count = 0;
  for(itr = m_nodeDofMap.begin(); itr != iEnd; ++itr){
  	itr->second = count;
  	count++;
  }
  m_numNodes = count;
  
}

void ParallelPlateADRSolver::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_TrilinosIndexStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(concentrationFieldName_);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_TrilinosIndexStr);
  partition.SetBufferSizes(syncedFields, CommRegistry::twoDARDSolver);
}


// w = Sum 1/|x_node - x_faceCenter|
void ParallelPlateADRSolver::GenerateFaceCenters(NodeManagerT& nodeManager,FaceManagerT& faceManager){
  
  localIndex setSize = m_faceSet->size();
  
  Array1dT<R1Tensor>& centers = faceManager.GetFieldData<R1Tensor>(FaceCenterStr);
  Array1dT<realT>& areas = faceManager.GetFieldData<realT>(FaceAreaStr);
  Array1dT<R1Tensor>& normals = faceManager.GetFieldData<R1Tensor>(FaceNormalStr);
    
  lSet::const_iterator fi=m_faceSet->begin();
  for(localIndex i =0; i < setSize; ++i){
    const localIndex f = *(fi++);
    
    // set face center
    faceManager.FaceCenter( nodeManager, f, centers[f]);
    
   // set face area and normals
  	R1Tensor eZ = faceManager.FaceNormal(nodeManager,f); 
  	areas[f] = faceManager.ProjectedArea(nodeManager,f,eZ);
  	normals[f] = eZ;
  }
  
}

/**
 * 
 * 
**/
double ParallelPlateADRSolver::TimeStep( const realT& time,
                          const realT& dt,
                          const int cycleNumber,
                          PhysicalDomainT& domain,
                          const sArray1d& namesOfSolverRegions ,
                          SpatialPartition& partition,
                          FractunatorBase* const fractunator )
{  

  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();

  // repeat solve for multiple species 
  // nb assumes same sparsity pattern - so may need to be modified if
  // new bc types are added
  for(unsigned i = 0; i < m_concentrationFieldNames.size(); ++i){
      SetConcentrationField(m_concentrationFieldNames[i],
                            m_reactionRateFieldNames[i],
                            rrDerivFieldName_);
      if(i==0) SetupSystem (domain,partition);
      Assemble  (domain,partition, time,dt);
      Solve     (domain,partition);
  }


  // limit maximum increase in timestep
  m_stabledt.m_maxdt *=  this->m_courant;
  if(m_stabledt.m_maxdt > 2*dt && dt > 0.0 && ~m_doSteadyState){
    m_stabledt.m_maxdt = 2*dt;
  }

  return dt;
}


void ParallelPlateADRSolver:: SetupSystem (PhysicalDomainT&  domain,
                                                SpatialPartition& partition)
{

  // calculate face centers and node weights
  GenerateFaceCenters(domain.m_feNodeManager,domain.m_feFaceManager);


  // Trilinos
  ///////////	
  iArray1d& face_trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& face_is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  iArray1d& node_trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& node_is_ghost       = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  
  
  // count ghost and local rows
  ///////////////////////////////
//  int n_hosted_rows = m_faceSet->size()+m_nodeSet->size();
//  int n_hosted_faces = m_faceSet->size();
  int n_ghost_rows = 0;
  int n_local_rows = 0;
  
  // faces
  for( lSet::const_iterator fItr=m_faceSet->begin() ; fItr!=m_faceSet->end() ; ++fItr )
  {
  	if( face_is_ghost[*fItr] < 0 ){
  		++n_local_rows;
  	} else {
  		++n_ghost_rows;
  	}
  }

  m_numLocalFaces = n_local_rows;

  // nodes
  for(lSet::const_iterator nItr = m_nodeSet->begin(); nItr != m_nodeSet->end(); ++nItr){
  	if( node_is_ghost[*nItr] < 0 ){
  		++n_local_rows;
  	} else {
  		++n_ghost_rows;
  	}
  }
  
  // determine the global/local degree of freedom distribution.
  ////////////////////////////////////////////////////////////

  std::vector<int> gather(n_mpi_processes);

  m_epetra_comm.GatherAll(&n_local_rows,
                        &gather.front(),
                        1);

  int first_local_row = 0;
  int n_global_rows = 0;

  for( int p=0; p<n_mpi_processes; ++p)
  {
    n_global_rows += gather[p];
    if(p<this_mpi_process)
      first_local_row += gather[p];
  }
  
   // create trilinos dof indexing
  //////////////////////////////////
  unsigned local_count = 0;
  // faces
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    if(face_is_ghost[*kf] < 0)
    {
      face_trilinos_index[*kf] = first_local_row+local_count;
      local_count++;
    }
    else
    {
      face_trilinos_index[*kf] = -INT_MAX;
    }
  }

  // nodes
  for( lSet::const_iterator kn=m_nodeSet->begin() ; kn!=m_nodeSet->end() ; ++kn )
  {
    if(node_is_ghost[*kn] < 0)
    {
      node_trilinos_index[*kn] = first_local_row+local_count;
      local_count++;
    }
    else
    {
      node_trilinos_index[*kn] = -INT_MAX;
    }
  }

  assert(local_count == static_cast<unsigned int>(n_local_rows) );

  partition.SynchronizeFields(syncedFields, CommRegistry::twoDARDSolver);

  // create epetra map
  ////////////////////

  row_map = Teuchos::rcp(new Epetra_Map(n_global_rows,n_local_rows,0,m_epetra_comm));

  // set up sparsity graph
  ////////////////////////

  sparsity = Teuchos::rcp(new Epetra_FECrsGraph(Copy,*row_map,0));
  
  iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

 
  // loop over edges
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {

    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0 ){   
      int numFaces = itr->second.size();
      if( numFaces >= 2){
        std::vector<int> dofIndex (numFaces); 
        std::vector<int> dofIndexB (2); 
  		
  	    eg = itr->first;
        for(int i =0; i < numFaces; ++i){
  	      localIndex kf = itr->second[i];
          dofIndex[i] = face_trilinos_index[kf];  
        }

        localIndex& nda = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
        localIndex& ndb = domain.m_feEdgeManager.m_toNodesRelation(eg,1);
      	dofIndexB[0] = node_trilinos_index[nda]; 
      	dofIndexB[1] = node_trilinos_index[ndb]; 

        // face -> face terms
        sparsity->InsertGlobalIndices(dofIndex.size(),
                                      &dofIndex.front(),
                                      dofIndex.size(),
                                      &dofIndex.front());

        // node (col) -> face (row) terms
        if( numFaces == 2){ // fixme currently only valid for 2 faces per edge
          sparsity->InsertGlobalIndices(dofIndex.size(),
                                        &dofIndex.front(),
                                        dofIndexB.size(),
                                        &dofIndexB.front());
        }

      }
    }	
  }
  
  // face (col) -> node (row) interactions
  // loop over faces
  for( lSet::const_iterator fi=m_faceSet->begin() ; fi!=m_faceSet->end() ; ++fi )
  {
    if(face_is_ghost[*fi] < 0)
    {      
      const lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[*fi];
      std::vector<int> dofIndex (1);
      std::vector<int> dofIndexB (2);
      dofIndexB[1] = face_trilinos_index[*fi];

      for( size_t a=0 ; a<nodeList.size(); ++a ){
        const localIndex nd = nodeList[a];
        dofIndex[0] = dofIndexB[0] = node_trilinos_index[nd];
      
        sparsity->InsertGlobalIndices(dofIndex.size(),
                                      &dofIndex.front(),
                                      dofIndexB.size(),
                                      &dofIndexB.front());
      }
    }   
       
  }

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();
    
}

/* Assemble */

void ParallelPlateADRSolver :: Assemble (PhysicalDomainT&  domain,
                                         SpatialPartition& partition ,
                                         const realT& time, const realT& dt)

{
	
    // (re-)init linear system

  matrix   = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*sparsity));
  solution = Teuchos::rcp(new Epetra_FEVector(*row_map));
  rhs      = Teuchos::rcp(new Epetra_FEVector(*row_map));

  // trilinos data

  iArray1d& face_trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& face_is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  iArray1d& node_trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);

  iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  

  // data
  Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  rArray1d& faceApertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );        
  rArray1d& faceAreas = domain.m_feFaceManager.GetFieldData<realT>(FaceAreaStr);
  Array1dT<R1Tensor>& faceNormals = domain.m_feFaceManager.GetFieldData<R1Tensor>(FaceNormalStr);
  
  rArray1d& fluidFluxes = domain.m_feEdgeManager.GetFieldData<realT>( VolumetricFluxStr );
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  rArray1d& faceConcentrations = domain.m_feFaceManager.GetFieldData<realT>(concentrationFieldName_);


  const Array1dT<realT>* reactionRateFieldPtr = NULL;  // reaction rate
  const Array1dT<realT>* dRRdCFieldPtr = NULL;   // reaction rate derivative wrt self dR/dC
  if(hasReaction){
    reactionRateFieldPtr = &(domain.m_feFaceManager.GetFieldData<realT>(reactionRateFieldName_));
    dRRdCFieldPtr = &(domain.m_feFaceManager.GetFieldData<realT>(rrDerivFieldName_));
  }  
  
    
  // loop over edges
  //////////////////
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0 ){   
  	  int numFaces = itr->second.size();

      Epetra_IntSerialDenseVector face_dof(numFaces);
      Epetra_IntSerialDenseVector node_dof(2);

      Epetra_SerialDenseMatrix face_K (numFaces,numFaces);
      Epetra_SerialDenseMatrix face_node_K (numFaces,2);

  	  if( numFaces == 2){

        localIndex& nda = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
        localIndex& ndb = domain.m_feEdgeManager.m_toNodesRelation(eg,1);
        node_dof[0] = node_trilinos_index[nda];
        node_dof[1] = node_trilinos_index[ndb];

  	    realT Ka = 0.0;
  	    realT Kb = 0.0;
  		
  	    localIndex kf = itr->second[0];
  	    localIndex kfb = itr->second[1];
        face_dof[0] = face_trilinos_index[kf];
        face_dof[1] = face_trilinos_index[kfb];
  	  
        realT edgeAperture = 0.5*(faceApertures[kf] + faceApertures[kfb] );
      
        R1Tensor eZ = faceNormals[kf]; // nb means dispersion must be defined in face a's coordinate system
      
        R1Tensor dL = faceCenters[kfb] - faceCenters[kf]; // branch vector (normalized eventually)
        R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,dT); // edge vector
      
        realT Lnorm = dL.L2_Norm();// branch vector length
        realT Tnorm = dT.L2_Norm();// edge length
      
        dL /= Lnorm;
        dT /= Tnorm;
      
        // edge norm
        R1Tensor dN; dN.Cross(dT,eZ); dN.Normalize();
            
        realT CosTheta = dN*dL; //abs(dN*dL);
      
        realT SinTheta = dT*dL;
            
        // diffusion terms  =  \int - D_{ij}n_{i}dC/dx_{j} dS
        //////////////////
      
        R1Tensor phiDijNi; phiDijNi.AijBi(m_D,dN);
        phiDijNi *= edgeAperture;
      
        realT phiDijNiNj = phiDijNi*dN;
        realT phiDijNiTj = phiDijNi*dT;
      
        // face diffusion term
        realT DT1 = Tnorm*phiDijNiNj/(abs(CosTheta)*Lnorm);
        Ka += DT1;
        Kb += DT1;

        //DEBUG(DT1)
      
        // node diffusion term (tangential contribution)
        realT DT2 = phiDijNiTj - phiDijNiNj*SinTheta/CosTheta;

      
        // upwind differencing for advection term
        //////////////////////////////////////////
        realT edgeFlux = fluidFluxes[eg];  // nb. flux = net volume flux between cells
        if(edgeFlux > 0){  // flowing a -> b
          	Ka += edgeFlux;
        } else { // flowing b -> a
          	Kb += -edgeFlux;
        }
      
        // Matrix entry for faces
        ////////////////////////////////
          
        // self terms
        //AddSparseMatrixTerm( dof,  dof, Ka, ki);
        //AddSparseMatrixTerm(dofB, dofB, Kb, ki);
        face_K(0,0) = Ka;
        face_K(1,1) = Kb;
 
        // off-axis terms
        //AddSparseMatrixTerm(dof, dofB, -Kb, ki);
        //AddSparseMatrixTerm(dofB, dof, -Ka, ki);
        face_K(0,1) = -Kb;
        face_K(1,0) = -Ka;

        matrix->SumIntoGlobalValues(face_dof, face_K);


        // Matrix entry for face/node terms
        ///////////////////////////////////

        //AddSparseMatrixTerm(dof,  dofNdA,  DT2, ki);
        //AddSparseMatrixTerm(dofB, dofNdA, -DT2, ki);
        face_node_K(0,0) =  DT2;
        face_node_K(1,0) = -DT2;
      	
        //AddSparseMatrixTerm(dof,  dofNdB, -DT2, ki);
        //AddSparseMatrixTerm(dofB, dofNdB,  DT2, ki);
        face_node_K(0,1) = -DT2;
        face_node_K(1,1) =  DT2;

        matrix->SumIntoGlobalValues(face_dof,
                                    node_dof,
                                    face_node_K);

      } else if(numFaces > 2){

          R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,dT); // edge vector
          realT Tnorm = dT.L2_Norm();// edge length
 
          //localIndex& nda = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
          //localIndex& ndb = domain.m_feEdgeManager.m_toNodesRelation(eg,1);

          R1Tensor edgeCenter; domain.m_feEdgeManager.EdgeCenter(domain.m_feNodeManager,eg,edgeCenter);
          lArray1d& faces = itr->second;

          // advection
          rArray1d qs(numFaces,0.0); // fluxes out of faces (into edge)
          rArray1d kappas(numFaces,0.0); // permeabilities
          rArray1d apertures(numFaces,0.0); // apertures
          std::vector<R1Tensor> dLs(numFaces); // dls
          rArray1d lNorms(numFaces,0.0); // dls
          rArray1d Ps(numFaces,0.0); // pressures
          //realT Cm = 0.0;// effective mixed concentration (for advection) Cm = \sum C^{u}_{i}Q^{u}_{i}/\sum Q^{u}_{i} (i.e. normalized by flux into edge)
          realT Pe = 0.0;// pressure at edge Pe = \sum K_{i}P_{i}/\sum K_{j}
          realT ksum = 0.0;
          realT qsum = 0.0;
          realT qsumb = 0.0;

          // gather data
          for(int i =0; i < numFaces; ++i){
              localIndex kf = faces[i];
              face_dof[i] = face_trilinos_index[kf];

              R1Tensor dL = edgeCenter - faceCenters[kf]; // branch vector (normalized eventually)
              apertures[i] = faceApertures[kf];
              lNorms[i] = dL.L2_Norm();// branch vector length
              dLs[i] = dL;
              Ps[i] = faceFluidPressure[kf];

              // edge pressure contribution
              kappas[i] = PPFS::CalculatePermeability( lNorms[i], apertures[i], Tnorm ,m_mu,m_SHP_FCT);
              ksum += kappas[i];
              Pe += kappas[i]*Ps[i];
          }
          // calculate edge pressure
          Pe /= ksum;

          // determine fluxes into faces
          for(int i =0; i < numFaces; ++i){
        	  qs[i] = kappas[i]*(Ps[i]-Pe);
        	  if(qs[i] < 0){ // downstream faces only
        	    qsum -= qs[i];
        	  } else {
        		qsumb += qs[i];

        	  }
          }

          // Sanity check
          if(true){
              std::cout << qsum << " "<< qsumb << std::endl; // flux in = flux out
              for(int i =0; i < numFaces; ++i){
                  localIndex kf = faces[i];
            	if(qs[i] <0){
            		std::cout <<"Downstream: edge " << eg <<" face " << kf ;// << std::cout;
            	}else {
            		std::cout <<"Upstream: edge " << eg  <<" face " <<  kf ; //<< std::cout;

            	}
                //lArray1d& edges = domain.m_feFaceManager.m_toEdgesRelation[kf];

                for (unsigned int ii = 0; ii < domain.m_feFaceManager.m_toEdgesRelation[kf].size(); ++ii){
              	  localIndex egb = domain.m_feFaceManager.m_toEdgesRelation[kf][ii];

              	  if(egb != eg) {
            		  std::cout << "egb: " << egb << " " << fluidFluxes[egb] << std::endl;
            	  }

                }
    		    std::cout << "eg: " << eg << " " << qs[0] << std::endl;
              }

          }


          // effective mixed concentration (for advection) Cm = \sum C^{u}_{i}Q^{u}_{i}/\sum Q^{u}_{i} (i.e. normalized by flux into edge)
          for(int i =0; i < numFaces; ++i){
        	  if(qs[i] > 0){ // downstream face, amount removed = concentration*flux
                  face_K(i,i) +=  qs[i];
        	  } else { //upstream face, amount injected*mixed concentration*flux in
        		  for(int ii =0; ii < numFaces; ++ii){
        		      if(qs[ii] > 0){
        		         face_K(i,ii) +=  qs[i]*qs[ii]/qsum;
        		      }
        		  }
        	  }
          }


/*
          // diffusion
          rArray1d ksD(numFaces,0.0);
          realT kSumD = 0.0;
          // calculate edge contribution
          for(int i =0; i < numFaces; ++i){
            localIndex kf = faces[i];

            realT app = faceApertures[kf];
            R1Tensor dL = edgeCenter - faceCenters[kf]; // branch vector (normalized eventually)
      
            realT Lnorm = dL.L2_Norm();// branch vector length
            dL /= Lnorm;

            R1Tensor phiDijLi; phiDijLi.AijBi(m_D,dL);
            phiDijLi *= app;
            realT phiDijLiLj = phiDijLi*dL;

            ksD[i] = phiDijLiLj/Lnorm;
            kSumD += ksD[i];
          }

          //realT invkSumADR = 1.0/(kSumADR+ TINY);
          realT invkSumD = 1.0/(kSumD + TINY);


          // add contributions to stiffness matrices
          for(int i =0; i < numFaces; ++i)
          {
              face_K(i,i) +=  -ksD[i]*(1.0 - ksD[i]*invkSumD);
              // cross terms
              for(int ii =0; ii < numFaces; ++ii)
              {
                if(i!=ii){
                    face_K(i,ii) +=  ksD[i]*ksD[ii]*invkSumD;
                }
              }
          }

        */
          // Matrix entry for faces
          ////////////////////////////////

          matrix->SumIntoGlobalValues(face_dof, face_K);



          // Matrix entry for face/node terms
          ///////////////////////////////////
          /*matrix->SumIntoGlobalValues(face_dof,
                                      node_dof,
                                      face_node_K);*/
  	  }

    }// ghost edges
  }
  
  // loop over faces
  //////////////////
  for( lSet::const_iterator fi=m_faceSet->begin() ; fi!=m_faceSet->end() ; ++fi )
  {
    if(face_is_ghost[*fi] < 0){
      localIndex fc = *fi;

      // Face self terms
      //////////////////

      R1Tensor center = faceCenters[fc];
      realT area = faceAreas[fc];
      realT vol = area*faceApertures[fc];

      Epetra_IntSerialDenseVector  face_dof (1);
      Epetra_SerialDenseVector  face_p (1);
      Epetra_SerialDenseMatrix face_K (1,1);
      
      face_dof[0] = face_trilinos_index[fc];
      
      if(~m_doSteadyState){
        face_K(0,0) = vol/(dt+TINY);
        face_p(0) = vol*faceConcentrations[fc]/(dt+TINY);

      }

      if(hasReaction){
        // V*[Cn - Co]/dt = A*Q(Cn);  Q(Cn) ~= Q(Co) + dQ(Co)/dC [Cn - Co]

        // limit maximum reaction amount to 0.5* concentration
        realT AQ = area*(*reactionRateFieldPtr)[fc];
        if(fabs(AQ) > 0.0 ){
          realT AQmax = 0.5*vol*faceConcentrations[fc]/(dt+TINY) ;

          if(AQ < -AQmax){
        	m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt,-dt*AQmax/AQ);
        	AQ = -AQmax;
        	//std::cout << "AQmax exceeded." <<std::endl;
          } else if(AQ > AQmax){
        	m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt,dt*AQmax/AQ);
          	AQ = AQmax;
        	//std::cout << "AQmax exceeded." << std::endl;
          };
        }

        //face_K(0,0) += vol/(dt+TINY) - area*(*dRRdCFieldPtr)[fc];
        //face_p(0) += vol*faceConcentrations[fc]/(dt+TINY) + AQ - area*(*dRRdCFieldPtr)[fc]*faceConcentrations[fc];

        face_K(0,0) += -area*(*dRRdCFieldPtr)[fc];
        face_p(0) += AQ - area*(*dRRdCFieldPtr)[fc]*faceConcentrations[fc];
      } 

      if(m_doSteadyState && face_K(0,0) == 0){
        face_K(0,0) += TINY; // avoid isolated DOF
      }

      // record matrix entry for face self terms
      matrix->SumIntoGlobalValues(face_dof, face_K);
      rhs->SumIntoGlobalValues(face_dof, face_p);
    

      // Face - node interactions
      /////////////////////////////
      Epetra_IntSerialDenseVector  node_dof (1);
      Epetra_IntSerialDenseVector  node_dof_face_dof (2);
      Epetra_SerialDenseMatrix node_face_matrix (1,2);
      node_dof_face_dof(1) = face_trilinos_index[fc];

      const lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[fc];
      for( size_t a=0 ; a<nodeList.size(); ++a ){
        const localIndex nd = nodeList[a];
        R1Tensor nodePos; domain.m_feNodeManager.GetPosition(nd,nodePos);
      
        realT denom = (nodePos -  center).L2_Norm() + TINY;
        realT Kn = 1.0/denom; 
        // record matrix entry
        node_dof(0) = node_trilinos_index[nd];
        node_dof_face_dof(0) = node_trilinos_index[nd];
        node_face_matrix(0,0) =  Kn+ TINY;// avoid isolated DOF
        node_face_matrix(0,1) = -Kn; 
        matrix->SumIntoGlobalValues(node_dof,
                                    node_dof_face_dof,
                                    node_face_matrix);
      
      }   
    }
       
  }
  
  
  // boundary conditions
  ApplyBoundaryCondition<realT>(this, &ParallelPlateADRSolver::DirichletBoundaryCondition, 
                                domain, domain.m_feEdgeManager, concentrationFieldName_,time);
  ApplyBoundaryCondition<realT>(this, &ParallelPlateADRSolver::OutflowBoundaryCondition, 
                                domain, domain.m_feEdgeManager, concentrationFieldName_,time);
   
  matrix->GlobalAssemble();
  rhs->GlobalAssemble();

  //rhs->Print(std::cout);
  //std::cout << matrix->NormInf() << std::endl;
  //EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*matrix);
  //exit(0);
}


/* Solve */

void ParallelPlateADRSolver:: Solve (PhysicalDomainT&  domain,
                                               SpatialPartition& partition)
{

  iArray1d& face_trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& face_is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d& face_concentration  = domain.m_feFaceManager.GetFieldData<realT>(concentrationFieldName_);

  iArray1d& node_trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& node_is_ghost       = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d& node_concentration  = domain.m_feNodeManager.GetFieldData<realT>(concentrationFieldName_);

  int dummy;
  double* local_solution = NULL;

  solution->ExtractView(&local_solution,&dummy);


  // copy previous solution from faces
  ///////////////////////////////////

  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    if(face_is_ghost[*kf] < 0)
    {
      int lid = row_map->LID(face_trilinos_index[*kf]);
      local_solution[lid] = face_concentration[*kf];
    }
  }

  // copy previous solution from nodes
  ////////////////////////////////////

  for( lSet::const_iterator kn=m_nodeSet->begin() ; kn!=m_nodeSet->end() ; ++kn )
  {
    if(node_is_ghost[*kn] < 0)
    {
      int lid = row_map->LID(node_trilinos_index[*kn]);
      local_solution[lid] = node_concentration[*kn];
    }
  }

 // EpetraExt::RowMatrixToMatlabFile("beforeScaling.dat",*matrix);

  // Scaling
  // Row scaling
  if(true)
  {
	//    EpetraExt::RowMatrixToMatlabFile("beforeScaling.dat",*matrix);
	//    EpetraExt::MultiVectorToMatrixMarketFile("rhsBeforeScaling.txt", *rhs);
    Epetra_Vector scaling(*row_map); //m_matrix->RowMap()  *row_map
    matrix->InvRowSums(scaling);
    matrix->LeftScale(scaling);

    Epetra_FEVector tmp (*rhs);
    rhs->Multiply(1.0,scaling,tmp,0.0);

   // EpetraExt::VectorToMatrixMarketFile("Scaling.dat", scaling);
    //EpetraExt::RowMatrixToMatlabFile("afterScaling.dat", *matrix );
    //EpetraExt::MultiVectorToMatrixMarketFile("rhsAfterScaling.txt", *rhs);
  //  EpetraExt::MultiVectorToMatrixMarketFile("scalingVector.txt",scaling);
  }


  // Krylov solver

  Epetra_LinearProblem problem(&(*matrix),
                               &(*solution),
                               &(*rhs));

  // ML preconditioner
  //////////////////////

  // create a parameter list for ML options
  Teuchos::ParameterList MLList;

  ML_Epetra::SetDefaults("SA",MLList);

  // create the preconditioning object.
  ML_Epetra::MultiLevelPreconditioner* MLPrec = NULL;
  if(numerics.useMLPrecond){
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*matrix, MLList);
  }


  //Solve

  AztecOO solver(problem);
  //        solver.SetAztecOption(AZ_solver,AZ_bicgstab);
 //         solver.SetAztecOption(AZ_precond, AZ_Jacobi);
          solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
          solver.SetAztecOption(AZ_solver,AZ_gmres);
          solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
  if(m_doSteadyState){
	  // steady state typically has very small RHS - use a different convergence criteria?
          solver.SetAztecOption(AZ_conv,AZ_rhs); // |r|_{\infty}/(|A|_{\infty}*|x| + |b|_{\infty}) Az_sol
  }else{
          solver.SetAztecOption(AZ_conv,AZ_rhs);// |r|_{2}/|b|_{2}
  }

  if(m_verboseFlag){
	  solver.SetAztecOption(AZ_output,AZ_all);
  }  else {
	  solver.SetAztecOption(AZ_output,AZ_none);
  }

  solver.Iterate(1000,numerics.krylov_tol);
  //solution->Print(std::cout);
  //EpetraExt::MultiVectorToMatrixMarketFile("solution.txt", *solution);


  // destroy the preconditioner
  if(numerics.useMLPrecond){
     delete MLPrec;
  }


  // copy solution to faces
  ////////////////////////
  
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    if(face_is_ghost[*kf] < 0)
    {
      int lid = row_map->LID(face_trilinos_index[*kf]);
      face_concentration[*kf] = local_solution[lid];
    }
  }

  // copy solution to nodes
  ////////////////////////
  
  for( lSet::const_iterator kn=m_nodeSet->begin() ; kn!=m_nodeSet->end() ; ++kn )
  {
    if(node_is_ghost[*kn] < 0)
    {
      int lid = row_map->LID(node_trilinos_index[*kn]);
      node_concentration[*kn] = local_solution[lid];
    }
  }


  // re-sync ghost nodes

  partition.SynchronizeFields(syncedFields, CommRegistry::twoDARDSolver);
}





// Add a Dirichlet boundary condition to a given set of objects (nodes faces etc).
void ParallelPlateADRSolver::DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bc, const lSet& set,realT time){

  
  EdgeManagerT *isEdgeBC = dynamic_cast<EdgeManagerT *>(&object);
  if(isEdgeBC && bc->GetComponent(time) != 2){ // hack to get outflow bc working fixme
  	
    iArray1d& face_trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
    rArray1d& faceApertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );  
    rArray1d& fluidFluxes = domain.m_feEdgeManager.GetFieldData<realT>( VolumetricFluxStr );

    Epetra_IntSerialDenseVector face_dof(1); 
    Epetra_SerialDenseMatrix face_K (1,1);
    Epetra_SerialDenseVector face_p(1); 
    
    iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
      
    for( lSet::const_iterator egItr=set.begin() ; egItr!=set.end() ; ++egItr )
    {
   
      localIndex eg = *egItr;
      if(edge_is_ghost[eg] < 0){
        lArray1d& edgeFaces = m_edgesToFaces[eg];
        R1Tensor edgeCenter; domain.m_feEdgeManager.EdgeCenter(domain.m_feNodeManager,eg,edgeCenter);
      
      
        // edge tangent vector
        R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,dT); 
        realT Tnorm = dT.L2_Norm();// edge length
        dT /= Tnorm; 
      
        realT volFlux = fluidFluxes[eg];
        realT value = bc->GetValue(domain.m_feEdgeManager,egItr,time);
      
        for(size_t ii =0; ii < edgeFaces.size(); ++ii){
      	
    	  localIndex fc = edgeFaces[ii];
    	
          realT aperture = faceApertures[fc]; // edge aperture = face aperture at boundary
          realT Kf = 0.0;
          realT Kc = 0.0;
        
    	  R1Tensor faceCenter; domain.m_feFaceManager.FaceCenter(domain.m_feNodeManager,fc,faceCenter);
    	      
          R1Tensor dL = edgeCenter-faceCenter; // branch vector (normalized eventually)
      
          realT Lnorm = dL.L2_Norm();// branch vector length
      
          dL /= Lnorm; 
      
          // edge norm
          R1Tensor dN; dN = dL - (dL*dT)*dT;  dN.Normalize();
            
         // realT SinTheta = dT*dL;
          realT CosTheta =dN*dL; // abs(dN*dL); 
            
          // diffusion terms  =  \int - D_{ij}n_{i}dC/dx_{j} dS 
          //////////////////
      
          R1Tensor phiDijNi; phiDijNi.AijBi(m_D,dN); 
          phiDijNi *= aperture;
      
          realT phiDijNiNj = phiDijNi*dN;
          //realT phiDijNiTj = phiDijNi*dT;
        
          // face diffusion term
          realT DT1 = Tnorm*phiDijNiNj/(CosTheta*Lnorm);
          Kc += DT1;
          Kf += DT1;
      
          // No tangential term - constant concentration along edge
      
          // upwind differencing for advection term
          //////////////////////////////////////////
        
      	  if(volFlux > 0){  // flowing out
      	    Kc += volFlux;
          } else { // flowing in 
        	  Kf += -volFlux;
          }
      
          // record matrix entry for face
          face_dof(0) = face_trilinos_index[fc];
          face_K(0,0) = Kc;
          matrix->SumIntoGlobalValues(face_dof, face_K);
        
          // record change to RHS
          face_p[0] = Kf*value;
          rhs->SumIntoGlobalValues(face_dof, face_p);
        }
      }
    }
    
  } 
}


// Add an outflow boundary condition to a given set of edges.
void ParallelPlateADRSolver::OutflowBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bc, const lSet& set, realT time){

  
  EdgeManagerT *isEdgeBC = dynamic_cast<EdgeManagerT *>(&object);
  if(isEdgeBC && bc->GetComponent(time) == 2){ // hack to get outflow bc working fixme
  
//    rArray1d& face_concentration  = domain.m_faceManager.GetFieldData<realT>(concentrationFieldName_);
  	 // set dC/dn = 0 at the boundary
    iArray1d& face_trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
//    rArray1d& faceApertures = domain.m_faceManager.GetFieldData<realT>( ApertureStr );
    rArray1d& fluidFluxes = domain.m_feEdgeManager.GetFieldData<realT>( VolumetricFluxStr );

    Epetra_IntSerialDenseVector face_dof(1); 
    Epetra_SerialDenseMatrix face_K (1,1);
    Epetra_SerialDenseVector face_p(1); 
    
    iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
      
    for( lSet::const_iterator egItr=set.begin() ; egItr!=set.end() ; ++egItr )
    {
   
      localIndex eg = *egItr;
      if(edge_is_ghost[eg] < 0){
        lArray1d& edgeFaces = m_edgesToFaces[eg];
      
        realT volFlux = fluidFluxes[eg];
      
        for(size_t ii =0; ii < edgeFaces.size(); ++ii){
      	
      	  localIndex fc = edgeFaces[ii];
        
          // upwind differencing for advection term
          //////////////////////////////////////////
          
      	  if(volFlux > 0){  // flowing out - ignored if flowing in
      	    //Kc = volFlux;
      
            // record matrix entry for face
            face_dof(0) = face_trilinos_index[fc];
            face_K(0,0) = volFlux;
            matrix->SumIntoGlobalValues(face_dof, face_K);
        
            // record change to RHS
            //face_p[0] = volFlux*face_concentration[fc];
            //rhs->SumIntoGlobalValues(face_dof, face_p);
          }
        }
      }
    }
    
  } 
}
  
//////////////////////////////////////////
/// Register solver in the solver factory
REGISTER_SOLVER( ParallelPlateADRSolver )

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////


/// Steady-State Advection diffusion reaction solver
SteadyStateParallelPlateADRSolver::SteadyStateParallelPlateADRSolver(  const std::string& name,
                                                                       ProblemManagerT* const pm ):
SolverBase(name,pm)
{
  m_solverMapPtr = &(pm->m_solvers);
}

SteadyStateParallelPlateADRSolver::~SteadyStateParallelPlateADRSolver()
{
  // TODO Auto-generated destructor stub
}

/*
 *  <SteadyStateParallelPlateADRSolver name="ssSolver"               * Name of the steady state solver 
 *                        adrSolver="adrSolver"         * Name of parallel plate advection diffusion reaction solver used to find the steady state
 *                        tol="0.01"                    * Convergence criteria. Will converge when (dConc)/Max(Conc) < tol for all species
 *                        faceset="Zmax"                * Faces to apply the solver to
 *                        speciesList="Ca CO2"          * List of solver species eg "Ca Fe CO2aq" etc
 *                        diffusiveLengthScale="1">     * Minimum lengthscale to resolve diffusion solution (used to set timstep for adrsolver).
 *    <ReactionRateFunctions>                           * Vector of reaction functions (one per species in species list)           
 *       <ReactionRateFunction                          
 *                    species="Ca"                      * Name of species     
 *                    function="CaReactionFunction"     * Name of function to calculate reaction rate
 *                    variables="Ca CO2"/>              * Reaction rate function variables (In order used by function)
 *       <ReactionRateFunction                          
 *                    species="CO2"                           
 *                    function="CO2ReactionFunction"     
 *                    variables="Ca CO2"/>
 *    </ReactionRateFunctions>
 * </SteadyStateParallelPlateADRSolver>
 * 
 */
void SteadyStateParallelPlateADRSolver::ReadXML(TICPP::HierarchicalDataNode* const hdn)
{
  SolverBase::ReadXML( hdn );

  m_tol = hdn->GetAttributeOrDefault<realT>("tol",0.01);
  m_dL = hdn->GetAttributeOrDefault<realT>("diffusiveLengthScale",1);
  
  std::string componentsStr = hdn->GetAttributeString("speciesList");
  m_ADRSolverName = hdn->GetAttributeString("adrSolver");
  std::cout << "   Parallel Plate ADRSolver " << m_ADRSolverName << std::endl;
  
  m_speciesList = Tokenize(componentsStr," ");
  
  // ReactionRateFunctions
  TICPP::HierarchicalDataNode* rrfsNode = hdn->GetChild("ReactionRateFunctions");
  if(rrfsNode)	{
  	
    for(TICPP::HierarchicalDataNode* fNode = rrfsNode->Next(true); fNode; fNode=rrfsNode->Next()){
  		
      rrfStruct rrf;
      rrf.species = fNode->GetAttributeString("species");
      rrf.function = fNode->GetAttributeString("function");
      std::string varsStr = fNode->GetAttributeString("variables");
      rrf.variables = Tokenize(varsStr," ");
  	  m_reactionRateFunctions[rrf.species] = rrf;
    }
  }
  
  // face set
  m_faceSetName = hdn->GetAttributeString("faceset");


}


void SteadyStateParallelPlateADRSolver::RegisterFields( PhysicalDomainT& domain )
{
 /* domain.m_nodeManager.AddKeyedDataField<FieldInfo::referencePosition>();
  domain.m_nodeManager.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_nodeManager.AddKeyedDataField<FieldInfo::incrementalDisplacement>();*/
  
  
  // register species and reaction rates
  for(sArray1d::size_type i =0; i < m_speciesList.size(); ++i){
    domain.m_feFaceManager.AddKeylessDataField<realT>( m_speciesList[i], true, true );
    domain.m_feFaceManager.AddKeylessDataField<realT>( m_speciesList[i]+"_ReactionRate", true, true );
    domain.m_feNodeManager.AddKeylessDataField<realT>( m_speciesList[i], true, true );
  }
  domain.m_feFaceManager.AddKeylessDataField<realT>( ReactionRateDerivStr, true, true );
  domain.m_feFaceManager.AddKeylessDataField<realT>( PreviousConcentrationStr, true, true );
}

void SteadyStateParallelPlateADRSolver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
	m_faceSet = &(domain.m_feFaceManager.GetSet(m_faceSetName));  
    m_numFaces = m_faceSet->size();
  
}

/**
 * 
 * 
**/
double SteadyStateParallelPlateADRSolver::TimeStep( const realT& time,
                                     const realT& dt ,
                                     const int cycleNumber,
                                     PhysicalDomainT& domain,
                                     const sArray1d& namesOfSolverRegions,
                                     SpatialPartition& partition,
                                     FractunatorBase* const fractunator )
{   
 
 FunctionManager& functionManager = FunctionManager::Instance();   
 localIndex numNodes = domain.m_feNodeManager.m_numNodes;
 
 std::map<std::string,SolverBase*>& solverMap = *m_solverMapPtr;
 SolverBase* solverPtr = solverMap[m_ADRSolverName]; 
 ParallelPlateADRSolver* adrSolverPtr = dynamic_cast<ParallelPlateADRSolver*>( solverPtr );
 if(adrSolverPtr == 0) throw GPException("SteadyStateParallelPlateADRSolver: " + m_ADRSolverName + " is not a ParallelPlateADRSolver.");
 
 // Solver timestep
 double D = adrSolverPtr->GetDiffusivity();
 const double dtau = 0.25*(m_dL*m_dL)/D; // internal timestep - fixme
 // DEBUG(dtau)
 
 // Internal timestep
 bool hasConverged = false;
 while(!hasConverged){
   std::cout << "Running steady state solver" << std::endl;
 	
   hasConverged = true;

   // Loop over species list
   ///////////////////////
   for( sArray1d::size_type i = 0; i < m_speciesList.size(); ++i){

     // Update reaction rate
     ///////////////////////
     const Array1dT<realT>& concFieldData = domain.m_feFaceManager.GetFieldData<realT>(m_speciesList[i]);
     Array1dT<realT>& prevConcFieldData = domain.m_feFaceManager.GetFieldData<realT>(PreviousConcentrationStr);
     
     if(m_reactionRateFunctions.find(m_speciesList[i]) != m_reactionRateFunctions.end() ){
       rrfStruct& rrf = m_reactionRateFunctions[m_speciesList[i]];
       Function& func = functionManager.GetFunction(rrf.function);
       
       Array1dT<realT>& rrFieldData = domain.m_feFaceManager.GetFieldData<realT>(m_speciesList[i]+"_ReactionRate");
       Array1dT<realT>& rrDerivFieldData = domain.m_feFaceManager.GetFieldData<realT>(ReactionRateDerivStr);
       
       // field pointers to reaction species
       int nVars = rrf.variables.size();
       std::vector<realT> x(nVars );
       std::vector<Array1dT<realT>*> fieldPtr(nVars );
       int c_indx = -1;  // index of this component -1 indicates reaction rate is not a function of component's concentration
       for(int j =0; j < nVars; ++j){
       	 fieldPtr[j] = &(domain.m_feFaceManager.GetFieldData<realT>( rrf.variables[j]) );
       	 if(rrf.variables[j] == m_speciesList[i]) c_indx = j;
       }
       
       for( lSet::const_iterator fi=m_faceSet->begin() ; fi!=m_faceSet->end() ; ++fi ){
       	  localIndex fc = *fi;
       	  for(int j =0; j < nVars; ++j){
       	  	x[j] = (*fieldPtr[j])[fc];
       	  }
          double rr = func(x[0]);
          rrFieldData[fc] = rr;
          if(c_indx != -1){ 
          	realT temp = x[c_indx];
          	realT dX = !isZero(temp)? fabs(temp) * SMALL : SMALL;
            x[c_indx] += dX;
            dX = x[c_indx] - temp; // reduce machine precision error
            double rrb = func(x[0]);
            rrDerivFieldData[fc] = (rrb - rr)/dX;
          } else { // reaction rate not a function of component concentration
          	rrDerivFieldData[fc] = 0.0;
          }
       }
       adrSolverPtr->SetConcentrationField(m_speciesList[i],m_speciesList[i]+"_ReactionRate",ReactionRateDerivStr);
     } else {
       adrSolverPtr->SetConcentrationField(m_speciesList[i]);
     }
     
     // Record concentration prior to advection 
     ////////////////////////////////////////

     for( lSet::const_iterator fi=m_faceSet->begin(); fi!=m_faceSet->end(); ++fi){
       localIndex fc = *fi;
       prevConcFieldData[fc] = concFieldData[fc];
     }
     
    
     // Solve advection diffusion reaction
     /////////////////////////////////////
     std::cout << "Running ADR solver: " << m_speciesList[i] << std::endl;
     adrSolverPtr->TimeStep(time, dtau, 0, domain, namesOfSolverRegions, partition, NULL );
    
     // Check convergence
     ////////////////////
     
     // err = |C_n -C_o|/max(C_n)
     bool calculateMaxErr = true;
     if(hasConverged || calculateMaxErr){
   
       double maxC = 0.0;
       double maxErr = 0.0;
       for( localIndex nd =0; nd < numNodes; ++nd){
       	 realT conc = concFieldData[nd];
       	 if(maxC < conc) maxC = conc;
       }
       std::cout << "Maximum concentration: " << maxC << std::endl;
       
       if(maxC > 0.0){
         for( lSet::const_iterator fi=m_faceSet->begin(); (fi!=m_faceSet->end())&& (hasConverged || calculateMaxErr); ++fi){
     
           localIndex fc = *fi;
           
           double err = fabs(concFieldData[fc] - prevConcFieldData[fc])/maxC;
                
       	   if(maxErr < err){
       	   	 maxErr = err;
             hasConverged = (err < m_tol);
       	   }
         }
       }
       
       std::cout << "Parallel-Plate Steady-State ADR Solver: Error - "<< m_speciesList[i] << " "<< maxErr << std::endl; 
       
     }
   } // component loop
   
   //hasConverged=true; // debug - if uncommented should revert to ADR solver
 } // Converged?
 
 std::cout << "Parallel-Plate Steady-State ADR Solver has converged.\n" << std::endl; 
 
 return dt;
}


/////////////////////////////////////////
/// Register solver in the solver factory
REGISTER_SOLVER( SteadyStateParallelPlateADRSolver )
