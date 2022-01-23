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

#include "TwoDADRSolver_Proppant.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"
#include "IO/ticpp/TinyXMLParser.h"

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
#include "ParallelPlateProppantFlowSolverImplicit.h"


using namespace BoundaryConditionFunctions;
using namespace PS_STR;
using namespace TDSSADR_STR;
namespace{
  std::string NodeWeightsStr =  "NodeWeights";
}


//#define DEBUG(x) std::cout << #x << " " << x << std::endl;

int ParallelPlateProppantSolver::m_instances = 0;

namespace{

 R2Tensor Eye3d(1.0, 0.0,  0.0,
                0.0, 1.0,  0.0,
                0.0, 0.0,  1.0);
 
  const realT TINY = 1e-64;
  const realT SMALL = 1e-8;
}


/// Least-Squares Finite-element advection diffusion reaction solver
ParallelPlateProppantSolver::ParallelPlateProppantSolver(const std::string& name,
                     ProblemManagerT* pm):
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
//m_proppantData(),
//m_slurryModel(),
hasReaction(),
m_D(),
concentrationFieldName_(),
reactionRateFieldName_(),
rrDerivFieldName_()
{
  m_solverMapPtr = &(pm->m_solvers);
  //ReadXML( hdn ); 
  m_TrilinosIndexStr = "TwoDPPADRS_" +  toString<int>(m_instances) + "_GlobalDof";  
  ++m_instances; 
  
}

ParallelPlateProppantSolver::~ParallelPlateProppantSolver()
{
  // TODO Auto-generated destructor stub
}

void ParallelPlateProppantSolver::ReadXML(TICPP::HierarchicalDataNode* hdn){
	


  m_FlowSolverName = hdn->GetAttributeString("flowSolver");
  std::cout << "   Proppant Flow Solver " << m_FlowSolverName << std::endl;

  realT D = hdn->GetAttributeOrDefault(DiffusivityStr,"1e-6m^2/s");
  m_D = Eye3d; m_D *= D;
  
  m_alpha_L = hdn->GetAttributeOrDefault("longitudialDispersivity","0");
  m_alpha_T = hdn->GetAttributeOrDefault("transverseDispersivity",0.1*m_alpha_L);

  m_verboseFlag = hdn->GetAttributeOrDefault<bool>(VerboseStr,true);

  m_doSteadyState = false; // = hdn->GetAttributeOrDefault<bool>("steadyState",false);
  
  m_doScreenOut = hdn->GetAttributeOrDefault<bool>("explicitPack",false);

  TICPP::HierarchicalDataNode* erosionNode = hdn->GetChild("ErosionModel");
  if(erosionNode)
  {
	  m_doErosion = true;
	  m_ErosionModel.ReadXML(erosionNode);
  } else {
	  m_doErosion = false;
	   std::cout << "No erosion model specified." <<std::endl;
  }

  m_useCollisionalSlipVelocity =  hdn->GetAttributeOrDefault<bool>("collisionalSlip",false);

  // face set
  m_faceSetName = hdn->GetAttributeString("faceset"); // if empty is based on "flow face" field (and recalculated every timestep)
  
  // Fieldnames
  //concentrationFieldName_ = hdn->GetAttributeStringOrDefault("species",ConcentrationStr);
  //reactionRateFieldName_ = hdn->GetAttributeStringOrDefault("reactionRateField","");

  concentrationFieldName_ = ProppantVolumeFractionStr;
  reactionRateFieldName_ = "";

  m_concentrationFieldNames.push_back(concentrationFieldName_);
  m_reactionRateFieldNames.push_back(reactionRateFieldName_);

  /**
  if(reactionRateFieldName_!=""){
  	hasReaction = true;
    rrDerivFieldName_ = hdn->GetAttributeStringOrDefault("reactionRateDerivField",ReactionRateDerivStr);
  } else {
  	hasReaction = false;
  }
  **/
  hasReaction =false;
  
  numerics.krylov_tol     = hdn->GetAttributeOrDefault("tol",1.0e-10);

  numerics.useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);
  
  this->m_courant  = hdn->GetAttributeOrDefault("courant",1.0);


  std::map<std::string,SolverBase*>& solverMap = *m_solverMapPtr;
  SolverBase* solverPtr = solverMap[m_FlowSolverName];
  m_flowSolverPtr = dynamic_cast<ParallelPlateProppantFlowSolverImplicit*>( solverPtr );
  if(m_flowSolverPtr == 0) throw GPException("Error ParallelPlateProppantSolver: No flow solver provided.");

  m_proppantDataPtr = &(m_flowSolverPtr->m_proppantData);
  m_slurryModelPtr = &(m_flowSolverPtr->m_proppantData.m_slurryModel);

  m_maxFluxCoefficient = 0.1; // limiter on max flux per timestep

  m_init_dt = hdn->GetAttributeOrDefault<realT>("initial_dt",0.001);
  
  m_upVector = R1Tensor(0,0,1);
  m_upVector.Normalize();
}


void ParallelPlateProppantSolver::RegisterFields( PhysicalDomainT& domain )
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
  domain.m_feEdgeManager.AddKeylessDataField<realT>( concentrationFieldName_, true, true );

  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>( FaceCenterStr, true, true );
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>( FaceNormalStr, true, true );
  domain.m_feFaceManager.AddKeylessDataField<realT>( ApertureStr, true, true );
  domain.m_feFaceManager.AddKeylessDataField<realT>( FaceAreaStr, true, true );

  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>( FluidVelocityStr, true, true );
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>( FluidVelocityStr+"_proppant", true, true );
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>( FluidVelocityStr+"_fluid", true, true );
  
  domain.m_feEdgeManager.AddKeylessDataField<realT>( VolumetricFluxStr, true, true );
  domain.m_feEdgeManager.AddKeylessDataField<realT>( VolumetricFluxStr+"_proppant", true, true );
  domain.m_feEdgeManager.AddKeylessDataField<realT>( VolumetricFluxStr+"_fluid", true, true );
  
  domain.m_feFaceManager.AddKeylessDataField<realT>(ProppantPackVolumeFractionStr ,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<realT>("ExcessProppantPackVolume",true,true);

  if(hasReaction){
    domain.m_feFaceManager.AddKeylessDataField<realT>( reactionRateFieldName_, true, true );
    domain.m_feFaceManager.AddKeylessDataField<realT>( rrDerivFieldName_, true, true );
  }
  
  domain.m_feFaceManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);
  domain.m_feNodeManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);

  if(m_doErosion){
    domain.m_feFaceManager.AddKeylessDataField<realT>("ErodedAperture", true, true);
    domain.m_feFaceManager.AddKeylessDataField<realT>("ErosionRate", true, true);
  }

     
}

void ParallelPlateProppantSolver::SetConcentrationField(std::string& concentrationFieldName,std::string reactionRateFieldName,std::string rrDerivFieldName){
  concentrationFieldName_ = concentrationFieldName;
  if(reactionRateFieldName != ""){
  	hasReaction = true;
  	reactionRateFieldName_ = reactionRateFieldName;
  	rrDerivFieldName_ = rrDerivFieldName;
  } else {
  	hasReaction = false;
  }
}

void ParallelPlateProppantSolver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{  
  if(m_faceSetName.empty()){
	    m_faceSet = &(domain.m_feFaceManager.m_Sets["flowFaceSet"]);
	    m_nodeSet = &(domain.m_feNodeManager.m_Sets["flowFaceNodeSet"]);

	    DefineFlowSets(domain);

  } else {
    m_faceSet = &(domain.m_feFaceManager.GetSet(m_faceSetName));
    m_nodeSet = &(domain.m_feNodeManager.GetSet(m_faceSetName));
  
    m_numFaces = m_faceSet->size();

    //m_cumNodeCount = 0;
  
    // build face and node-dof map
    lSet::const_iterator fi=m_faceSet->begin();
    for(localIndex i =0; i < m_numFaces; ++i){
    	localIndex f = *(fi++);
    	m_faceDofMap[f] = i;
  	
    	const lArray1d& faceNodeMap = domain.m_feFaceManager.m_toNodesRelation[f];
        for( localIndex a=0 ; a<faceNodeMap.size(); ++a ){
          //const localIndex nd = faceNodeMap[a];
          localIndex nd = domain.m_feNodeManager.GetParentIndex(faceNodeMap[a]);
  	      m_nodeDofMap[nd] = 0;
        }
  	
  	   // m_cumNodeCount += faceNodeMap.size();
    }

    std::map<localIndex,localIndex>::iterator itr;
    std::map<localIndex,localIndex>::iterator iEnd = m_nodeDofMap.end();
    localIndex count = 0;
    for(itr = m_nodeDofMap.begin(); itr != iEnd; ++itr){
  	  itr->second = count;
  	  count++;
    }
    m_numNodes = count;
  
    // edge parents to face parents
    m_edgesToFaces.clear();
    fi=m_faceSet->begin();
    for( localIndex i=0 ; i<m_numFaces ; ++i )
    {
      localIndex kf = *(fi++);
      localIndex numEdges = domain.m_feFaceManager.m_toEdgesRelation[kf].size();
      for(size_t a =0; a < numEdges; ++a){
        localIndex eg = domain.m_feFaceManager.m_toEdgesRelation[kf][a];
        localIndex egParent = domain.m_feEdgeManager.GetParentIndex(eg);
    	
        lSet& edgeFaces = domain.m_feEdgeManager.m_toFacesRelation[eg];
        lArray1d edgeList;
      
        for( lSet::iterator edgeFace=edgeFaces.begin() ; edgeFace!=edgeFaces.end() ; ++edgeFace ){
          localIndex faceParent = domain.m_feFaceManager.GetParentIndex(*edgeFace);
      	  if(isMember(faceParent,*m_faceSet) && !isMember( faceParent, m_edgesToFaces[egParent] ) ){
            m_edgesToFaces[egParent].push_back(faceParent);
      	  }
        }
      }
    }
  
  }
}

void ParallelPlateProppantSolver::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_TrilinosIndexStr);
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(concentrationFieldName_); // needed?
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_TrilinosIndexStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(concentrationFieldName_);
  partition.SetBufferSizes(syncedFields, CommRegistry::twoDARDSolver);
}


// w = Sum 1/|x_node - x_faceCenter|
void ParallelPlateProppantSolver::GenerateFaceCenters(NodeManagerT& nodeManager,FaceManagerT& faceManager){
  
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
double ParallelPlateProppantSolver::TimeStep( const realT& time,
                          const realT& dt,
                          const int cycleNumber,
                          PhysicalDomainT& domain,
                          const sArray1d& namesOfSolverRegions ,
                          SpatialPartition& partition,
                          FractunatorBase* const fractunator )
{  

  m_dt = dt;
  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();

  DefineFlowSets( domain );

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
      if(m_doScreenOut) CalculateScreenout(domain,partition,time,dt);
      if(m_doErosion) CalculateErosionRate(domain,partition,time,dt);
  }

  // limit maximum increase in timestep
  m_stabledt.m_maxdt *=  this->m_courant;
  if(m_stabledt.m_maxdt > 2*dt && dt > 0.0 && ~m_doSteadyState){
    m_stabledt.m_maxdt = 2*dt;
  }

  if(dt == 0.0){
	  m_stabledt.m_maxdt = m_init_dt;
  }

  return dt;
}

// this will only update every timestep if the face set is left blank
// may want to add a flag to update it or not instead...
// or update after mesh has been altered
// note that face set may be updated by another solver, but nodeset may not (nor may edges to faces)
void ParallelPlateProppantSolver::DefineFlowSets( PhysicalDomainT& domain )
{

  if( m_faceSetName.empty() )
  {
  FaceManagerT& faceManager = domain.m_feFaceManager;
  //EdgeManagerT& edgeManager = domain.m_feEdgeManager;

  const iArray1d& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");
 //  const iArray1d& flowEdgeType = edgeManager.GetFieldData<int>("flowEdgeType");
 // const Array1dT<lSet>& edgeToFlowFaces = edgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");


  m_faceSet->clear();
  m_nodeSet->clear();
  m_edgesToFaces.clear();

  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 0 && domain.m_feFaceManager.IsParent(kf)  )
    {
      m_faceSet->insert( kf );

      //m_nodeSet->insert(faceManager.m_toNodesRelation[kf].begin(),faceManager.m_toNodesRelation[kf].end());
      lArray1d::iterator ndend = faceManager.m_toNodesRelation[kf].end();
      for(lArray1d::iterator nd = faceManager.m_toNodesRelation[kf].begin(); nd != ndend; ++nd){
          m_nodeSet->insert(domain.m_feNodeManager.GetParentIndex(*nd) );
      }
    }
  }

  m_numFaces = m_faceSet->size();
/*
  m_faceDofMap.clear();
  lSet::const_iterator si=m_faceSet->begin();
  for(localIndex i =0; i < m_numFaces; ++i, ++si){
    localIndex f = *si;
    m_faceDofMap[f] = i;
  }
  */

  m_faceDofMap.clear();
  //m_cumNodeCount = 0;
  // build face and node-dof map
  lSet::const_iterator fi=m_faceSet->begin();
  for(localIndex i =0; i < m_numFaces; ++i){
	  localIndex f = *(fi++);
	  m_faceDofMap[f] = i;

	  const lArray1d& faceNodeMap = domain.m_feFaceManager.m_toNodesRelation[f];
	  for( localIndex a=0 ; a<faceNodeMap.size(); ++a ){
		  //const localIndex nd = faceNodeMap[a];
		  localIndex nd = domain.m_feNodeManager.GetParentIndex(faceNodeMap[a]);
		  m_nodeDofMap[nd] = 0;
	  }

	 // m_cumNodeCount += faceNodeMap.size();
  }

  // node dof
  std::map<localIndex,localIndex>::iterator itr;
  std::map<localIndex,localIndex>::iterator iEnd = m_nodeDofMap.end();
  localIndex count = 0;
  for(itr = m_nodeDofMap.begin(); itr != iEnd; ++itr){
	  itr->second = count;
  	  count++;
  }
  m_numNodes = count;

  //iArray1d& ffCount = edgeManager.GetFieldData<int>("FlowFaceCount"); // debug
/*
  for( localIndex ke=0 ; ke < edgeManager.DataLengths() ; ++ke )
  {
    if( flowEdgeType[ke] == 0 && !(edgeToFlowFaces[ke].empty()) )
    {
      m_edgesToFaces[ke].assign( edgeToFlowFaces[ke].begin(), edgeToFlowFaces[ke].end() ) ;
     // ffCount[ke] = m_edgesToFaces[ke].size();

    }
  }
  */

  // parent edge to parent face
  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 0 && domain.m_feFaceManager.IsParent(kf)  )
    {
      for( localIndex ke=0 ; ke<domain.m_feFaceManager.m_toEdgesRelation[kf].size() ; ++ke )
      {
        const localIndex edgeIndex = domain.m_feFaceManager.m_toEdgesRelation[kf][ke];
        //flowEdgeType[edgeIndex] = 0;

        const localIndex parentEdgeIndex = domain.m_feEdgeManager.GetParentIndex( edgeIndex );
        m_edgesToFaces[parentEdgeIndex].push_back( kf );
      }
    }
  }

  /*
  std::map<localIndex,localIndex>::iterator itr;
  std::map<localIndex,localIndex>::iterator iEnd = m_nodeDofMap.end();
  localIndex count = 0;
  for(itr = m_nodeDofMap.begin(); itr != iEnd; ++itr){
	  itr->second = count;
	  count++;
  }
  m_numNodes = count;
  */

  }

}


void ParallelPlateProppantSolver:: SetupSystem (PhysicalDomainT&  domain,
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
      if( numFaces == 2){ // fixme currently only valid for 2 faces per edge
        std::vector<int> dofIndex (numFaces); 
        std::vector<int> dofIndexB (2); 
  		
  	    eg = itr->first;
        for(int i =0; i < numFaces; ++i){
  	      localIndex kf = itr->second[i];
          dofIndex[i] = face_trilinos_index[kf];  
        }

        localIndex nda = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
        localIndex ndb = domain.m_feEdgeManager.m_toNodesRelation(eg,1);
        nda = domain.m_feNodeManager.GetParentIndex(nda);
        ndb = domain.m_feNodeManager.GetParentIndex(ndb);
      	dofIndexB[0] = node_trilinos_index[nda]; 
      	dofIndexB[1] = node_trilinos_index[ndb]; 

        // face -> face terms
        sparsity->InsertGlobalIndices(dofIndex.size(),
                                      &dofIndex.front(),
                                      dofIndex.size(),
                                      &dofIndex.front());

        // node (col) -> face (row) terms
        sparsity->InsertGlobalIndices(dofIndex.size(),
                                      &dofIndex.front(),
                                      dofIndexB.size(),
                                      &dofIndexB.front());
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
        localIndex nd = domain.m_feNodeManager.GetParentIndex(nodeList[a]);
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

void ParallelPlateProppantSolver :: Assemble (PhysicalDomainT&  domain,
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

  const Array1dT<R1Tensor>& mixedFluidVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  Array1dT<R1Tensor>& particlePhaseVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr+"_proppant");
  Array1dT<R1Tensor>& fluidPhaseVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr+"_fluid");

  rArray1d& faceConcentrations = domain.m_feFaceManager.GetFieldData<realT>(concentrationFieldName_);

  // needed for triple junctions
  rArray1d& edgeMus = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  // settled pack
  rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);


  // Calculate proppant volume fraction, slipVelocity/fluid velocity/particle velocity
  
  for( lSet::const_iterator fi=m_faceSet->begin() ; fi!=m_faceSet->end() ; ++fi )
  {

     localIndex kf = *fi;
     R1Tensor slipVelocity(0);
     realT vf = std::min(faceConcentrations[kf],m_proppantDataPtr->m_maxVf); // concentration suspended in fluid

     //realT mf = vf/(vf+ m_proppantDataPtr->m_rhoRatio*(1-vf)); // mass fraction

     //const realT& area = faceAreas[kf];
     //realT vol = area*faceApertures[kf]*(1-facePackVfs[kf]);

     // determine proppant volume fraction from proppant mass

     R1Tensor eN = faceNormals[kf];
     //realT vol = area*faceApertures[kf];
     CalculateSlipVelocity(vf, eN, mixedFluidVelocity[kf], slipVelocity);
     particlePhaseVelocity[kf] =  mixedFluidVelocity[kf] + (1.0-vf)*slipVelocity;
     fluidPhaseVelocity[kf] = particlePhaseVelocity[kf] -slipVelocity;

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

      	  localIndex nda = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
      	  localIndex ndb = domain.m_feEdgeManager.m_toNodesRelation(eg,1);
          nda = domain.m_feNodeManager.GetParentIndex(nda);
          ndb = domain.m_feNodeManager.GetParentIndex(ndb);

      	  node_dof[0] = node_trilinos_index[nda]; 
      	  node_dof[1] = node_trilinos_index[ndb]; 

  	  realT Ka = 0.0;  // matrix contributions
  	  realT Kb = 0.0;	
  		
          // face indicies
  	  localIndex kf = itr->second[0];  
  	  localIndex kfb = itr->second[1];
          face_dof[0] = face_trilinos_index[kf];
          face_dof[1] = face_trilinos_index[kfb];    

          realT facePackVfa = facePackVfs[kf];
          realT facePackVfb = facePackVfs[kfb];
  	  
          realT edgePackVf = std::max(facePackVfa,facePackVfb); // fixme should account for orientation

          realT edgeAperture = 0.5*(faceApertures[kf] + faceApertures[kfb] )*(1-edgePackVf);
      
          R1Tensor eZ = faceNormals[kf]; // nb eZ is used for reference plane which means dispersion must be defined in face a's coordinate system
      
          R1Tensor dL = faceCenters[kfb] - faceCenters[kf]; // branch vector (normalized eventually)
          R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,dT); // edge vector
      
          realT Lnorm = dL.L2_Norm();// branch vector length
          realT Tnorm = dT.L2_Norm();// edge length
      
          dL /= Lnorm; 
          dT /= Tnorm; 
      
          // edge norm
          R1Tensor dN; dN.Cross(dT,eZ); dN.Normalize(); 
            
          realT CosTheta = dN*dL; 
      
          realT SinTheta = dT*dL;

          // dispersion
          R2Tensor Dij(m_D);
          if(false){

            R1Tensor dv(mixedFluidVelocity[kf]);
            R1Tensor dvb(mixedFluidVelocity[kfb]);
            // rotate dvb into face a's plane
            R1Tensor eZb = faceNormals[kfb];
            eZb *= copysign(1.0,Dot(eZ,eZb) );  // pointing in same direction
            R1Tensor dNb; dNb.Cross(dT,eZb);dNb.Normalize();

            // find edge velocity (average velocity in face a's plane
            dv += Dot(dvb,dT)*dT + Dot(dvb,dNb)*dN;
            dv *= 0.5;

            realT vNorm = dv.L2_Norm();
            dv /= vNorm;
        	realT D_Diag = m_alpha_T*vNorm;
        	Dij(0,0) += D_Diag;
        	Dij(1,1) += D_Diag;
        	Dij(2,2) += D_Diag;
        	R2Tensor vivj; vivj.dyadic_aa(dv); vivj *= ((m_alpha_L-m_alpha_T)*vNorm );
        	Dij += vivj;
          }


          // diffusion/Dispersion terms  =  \int - D_{ij}n_{i}dC/dx_{j} dS
          //////////////////
      
          R1Tensor phiDijNi; phiDijNi.AijBi(Dij,dN);
          phiDijNi *= edgeAperture;
      
          realT phiDijNiNj = phiDijNi*dN;
          realT phiDijNiTj = phiDijNi*dT;
      
          // face diffusion term
          realT DT1 = Tnorm*phiDijNiNj/(fabs(CosTheta)*Lnorm);
      
          // node diffusion term (tangential contribution)
          realT DT2 = phiDijNiTj - phiDijNiNj*SinTheta/CosTheta;

      
          // upwind differencing for advection term
          //////////////////////////////////////////

          if(faceApertures[kf] >  m_proppantDataPtr->m_diameter  // no transport if aperture is too small
             && faceApertures[kfb] > m_proppantDataPtr->m_diameter
             && facePackVfa < 1.0  // no transport if occupied by proppant pack
             && facePackVfb < 1.0  ){
            realT edgeFlux = fluidFluxes[eg];  // nb. flux = net volume flux between cells

            R1Tensor slipVelocity(0);
            realT driftFlux(0.0);
            realT vf = std::min(faceConcentrations[kf],m_proppantDataPtr->m_maxVf);
      	    realT vfb = std::min(faceConcentrations[kfb],m_proppantDataPtr->m_maxVf);

      	    realT area = faceAreas[kf];
            realT vol = area*faceApertures[kf]*(1-facePackVfs[kf]);
      	    realT areab = faceAreas[kfb];
            realT volb = areab*faceApertures[kfb]*(1-facePackVfs[kfb]);

            // drift velocity
            if(edgeFlux > 0){  // flowing a -> b
                CalculateSlipVelocity(vf, eZ,mixedFluidVelocity[kf], slipVelocity);
                //realT mf = vf/(vf+ m_proppantDataPtr->m_rhoRatio*(1-vf)); // mass fraction
                driftFlux = edgeAperture*Tnorm*(1-vf)*Dot(slipVelocity,dN)*copysign(1.0,CosTheta);
            } else { // flowing b -> a
                R1Tensor eZb = faceNormals[kfb];
                R1Tensor dNb; dNb.Cross(dT,eZb); dNb.Normalize();
                CalculateSlipVelocity(vfb, eZb,mixedFluidVelocity[kfb], slipVelocity);
                //realT mfb = vfb/(vfb+ m_proppantDataPtr->m_rhoRatio*(1-vfb)); // mass fraction
                realT CosThetab = dNb*dL;
                driftFlux =  edgeAperture*Tnorm*(1-vfb)*Dot(slipVelocity,dNb)*copysign(1.0,CosThetab);
            }
            edgeFlux += driftFlux;

            if(edgeFlux > 0){  // flowing a -> b
              realT maxFlux = std::min( vol , volb*(m_proppantDataPtr->m_maxVf- vfb)/(vf+TINY) )/(dt+TINY);
              if(maxFlux < 0) maxFlux = 0;
              maxFlux *= m_maxFluxCoefficient;
              if(maxFlux > edgeFlux){
                Ka += edgeFlux;
              } else {
                Ka += maxFlux;

                // cancel diffusion
                DT1 = 0;
                DT2 = 0;
              }
            } else { // flowing b -> a
              realT maxFlux = std::min( volb , vol*(m_proppantDataPtr->m_maxVf- vf)/(vfb+TINY) )/(dt+TINY);
              maxFlux *= m_maxFluxCoefficient;
              if(maxFlux < 0) maxFlux = 0;
              if(maxFlux > -edgeFlux){
                Kb += -edgeFlux;
              } else {
                Kb += maxFlux;

                // cancel diffusion
                DT1 = 0;
                DT2 = 0;
              }
        	  // adjust timestep
              realT minVol = std::min(area*faceApertures[kf],areab*faceApertures[kfb]);
              m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt, fabs(  0.5*m_maxFluxCoefficient*minVol/(edgeFlux + TINY)  ) );
            }
          }

          // Self diffusion terms
          Ka += DT1;
          Kb += DT1;
      
          // Matrix entry for faces
          ////////////////////////////////
          
          // self terms
          face_K(0,0) = Ka;
          face_K(1,1) = Kb;
 
          // off-axis terms
          face_K(0,1) = -Kb;
          face_K(1,0) = -Ka;

          matrix->SumIntoGlobalValues(face_dof, face_K);


          // Matrix entry for face/node terms
          ///////////////////////////////////

          face_node_K(0,0) =  DT2;
          face_node_K(1,0) = -DT2;
      	
          face_node_K(0,1) = -DT2;
          face_node_K(1,1) =  DT2;

          matrix->SumIntoGlobalValues(face_dof, 
                                      node_dof,
                                      face_node_K);

      
        }   else if(numFaces > 2){

            realT mu = edgeMus[eg];
            realT SHP_FCT =  1.0;

            R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,dT); // edge vector
            realT Tnorm = dT.L2_Norm();// edge length

            //localIndex& nda = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
            //localIndex& ndb = domain.m_feEdgeManager.m_toNodesRelation(eg,1);
           // nda = domain.m_feNodeManager.getParentIndex(nda);
           // ndb = domain.m_feNodeManager.getParentIndex(ndb);

            R1Tensor edgeCenter; domain.m_feEdgeManager.EdgeCenter(domain.m_feNodeManager,eg,edgeCenter);
            lArray1d& faces = itr->second;

            // advection
            rArray1d qs(numFaces,0.0); // fluxes out of faces (into edge)
            rArray1d kappas(numFaces,0.0); // permeabilities
            rArray1d apertures(numFaces,0.0); // apertures
            std::vector<R1Tensor> dLs(numFaces); // branch vectors
            rArray1d lNorms(numFaces,0.0); // branch vector norms
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
                kappas[i] = PPFS::CalculatePermeability( lNorms[i], apertures[i], Tnorm , mu,SHP_FCT);
                ksum += kappas[i];
                Pe += kappas[i]*Ps[i];
            }
            // calculate edge pressure
            Pe /= (ksum+TINY);

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
            /*
            if(true){
                std::cout << qsum << " "<< qsumb << std::endl; // flux in = flux out
                for(int i =0; i < numFaces; ++i){
                    localIndex kf = faces[i];
                  if(qs[i] <0){
              		std::cout <<"Downstream: edge " << eg <<" face " << kf  << std::cout;
              	  }else {
              		std::cout <<"Upstream: edge " << eg  <<" face " <<  kf << std::cout;

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

            }*/


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
        
     } // ghost edge
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

      face_K(0,0) = vol/(dt+TINY);
      face_p(0) = vol*faceConcentrations[fc]/(dt+TINY);

      /*
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

      */
      


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
        localIndex nd = domain.m_feNodeManager.GetParentIndex(nodeList[a]);

        R1Tensor nodePos; domain.m_feNodeManager.GetPosition(nd,nodePos);
      
        realT denom = (nodePos -  center).L2_Norm() + TINY;
        realT Kn = 1.0/denom; 
        // record matrix entry
        node_dof(0) = node_trilinos_index[nd];
        node_dof_face_dof(0) = node_trilinos_index[nd];
        node_face_matrix(0,0) =  Kn;
        node_face_matrix(0,1) = -Kn; 
        matrix->SumIntoGlobalValues(node_dof,
                                    node_dof_face_dof,
                                    node_face_matrix);
      
      }   

    }
       
  }
  
  
  // boundary conditions
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantSolver::DirichletBoundaryCondition,
                                domain, domain.m_feEdgeManager, concentrationFieldName_,time);
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantSolver::OutflowBoundaryCondition,
                                domain, domain.m_feEdgeManager, concentrationFieldName_,time);
   
  matrix->GlobalAssemble();
  rhs->GlobalAssemble();

  //rhs->Print(std::cout);
  //std::cout << matrix->NormInf() << std::endl;
  //EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*matrix);
  //exit(0);
}


/* Solve */

void ParallelPlateProppantSolver:: Solve (PhysicalDomainT&  domain,
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
          solver.SetAztecOption(AZ_solver,AZ_bicgstab);
          solver.SetAztecOption(AZ_precond, AZ_Jacobi);
 //         solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
 //         solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
  if(m_doSteadyState){
	  // steady state typically has very small RHS - use a different convergence criteria?
          solver.SetAztecOption(AZ_conv,AZ_rhs); // |r|_{\infty}/(|A|_{\infty}*|x| + |b|_{\infty}) Az_sol
  }else{
          solver.SetAztecOption(AZ_conv,AZ_rhs);// |r|_{2}/|b|_{2}
  }

  if(~m_verboseFlag) solver.SetAztecOption(AZ_output,AZ_none);

  solver.Iterate(1000,numerics.krylov_tol);
  //solution->Print(std::cout);


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
void ParallelPlateProppantSolver::DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bc, const lSet& set,realT time){

  
  EdgeManagerT *isEdgeBC = dynamic_cast<EdgeManagerT *>(&object);



  ::OutflowBoundaryCondition* ofbc  = dynamic_cast< ::OutflowBoundaryCondition*> ( bc);
  bool isOutflowBC = (isEdgeBC && bc->GetComponent(time) == 2) || (ofbc != 0); // old way, new way

  //if(isEdgeBC && bc->GetComponent(time) != 2){ // hack to get outflow bc working fixme
  if(isEdgeBC && ~isOutflowBC){
  	
    iArray1d& face_trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
    rArray1d& faceApertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
    rArray1d& faceAreas = domain.m_feFaceManager.GetFieldData<realT>(FaceAreaStr);
    rArray1d& fluidFluxes = domain.m_feEdgeManager.GetFieldData<realT>( VolumetricFluxStr );
    Array1dT<R1Tensor>& faceNormals = domain.m_feFaceManager.GetFieldData<R1Tensor>(FaceNormalStr);

    rArray1d& faceConcentrations  = domain.m_feFaceManager.GetFieldData<realT>(concentrationFieldName_);
    rArray1d& edgeConcentrations  = domain.m_feEdgeManager.GetFieldData<realT>(concentrationFieldName_);

    const Array1dT<R1Tensor>& mixedFluidVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);

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
        realT boundaryVf = bc->GetValue(domain.m_feEdgeManager,egItr,time);
      
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
          //realT phiDijNiTj = phiDijNi;
        
          // face diffusion term
          realT DT1 = Tnorm*phiDijNiNj/(fabs(CosTheta)*Lnorm);
          Kc += DT1;
          Kf += DT1;
      
          // No tangential term - constant concentration along edge
      
          // upwind differencing for advection term
          //////////////////////////////////////////
          if(m_proppantDataPtr->m_diameter < aperture ){
            realT vf = std::min(faceConcentrations[fc],m_proppantDataPtr->m_maxVf);
        	if(volFlux > 0){  // flowing out

                R1Tensor eZ = faceNormals[fc]; // nb means dispersion must be defined in face a's coordinate system
        		R1Tensor slipVelocity;
                realT driftFlux(0.0);
                CalculateSlipVelocity(vf, eZ,mixedFluidVelocity[fc], slipVelocity);
                //realT mf = vf/(vf+ m_proppantDataPtr->m_rhoRatio*(1-vf)); // mass fraction
                driftFlux = aperture*Tnorm*(1-vf)*Dot(slipVelocity,dN);

        	    Kc += volFlux+driftFlux;
            } else { // flowing in - ignore drift flux

              realT vol = faceAreas[fc]*aperture;
              realT maxFlux = vol*(m_proppantDataPtr->m_maxVf- vf)/((boundaryVf+TINY)*(m_dt+TINY));
             // if(maxFlux < 0) maxFlux = 0;
              maxFlux *= m_maxFluxCoefficient;
              if(maxFlux > -volFlux){
        	    Kf += -volFlux;
              } else {
          	    Kf += maxFlux;
               // m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt, fabs( 0.5*m_dt*maxFlux/(volFlux+TINY) ) );
                m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt, fabs(  m_maxFluxCoefficient*vol/(volFlux + TINY)  ) );
              }
            }
          }

  	      edgeConcentrations[eg] = boundaryVf;
      
          // record matrix entry for face
          face_dof(0) = face_trilinos_index[fc];
          face_K(0,0) = Kc;
          matrix->SumIntoGlobalValues(face_dof, face_K);
        
          // record change to RHS
          face_p[0] = Kf*boundaryVf;
          rhs->SumIntoGlobalValues(face_dof, face_p);
        }
      }
    }
    
  } 
}


// Add an outflow boundary condition to a given set of edges.
void ParallelPlateProppantSolver::OutflowBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bc, const lSet& set, realT time){

  rArray1d& faceApertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );

  //rArray1d& faceAreas = domain.m_feFaceManager.GetFieldData<realT>(FaceAreaStr);
  Array1dT<R1Tensor>& faceNormals = domain.m_feFaceManager.GetFieldData<R1Tensor>(FaceNormalStr);
  const Array1dT<R1Tensor>& mixedFluidVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  
  EdgeManagerT *isEdgeBC = dynamic_cast<EdgeManagerT *>(&object);

  ::OutflowBoundaryCondition* ofbc  = dynamic_cast< ::OutflowBoundaryCondition*> ( bc);
  bool isOutflowBC = (isEdgeBC && bc->GetComponent(time) == 2) || (ofbc != 0); // old way, new way

  if(isOutflowBC){
  
    rArray1d& faceConcentrations  = domain.m_feFaceManager.GetFieldData<realT>(concentrationFieldName_);
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
    	R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,dT);
        realT Tnorm = dT.L2_Norm();// edge length
        dT/=Tnorm;

        for(size_t ii =0; ii < edgeFaces.size(); ++ii){
      	
      	  localIndex fc = edgeFaces[ii];

      	  // upwind differencing for advection term
          //////////////////////////////////////////
          
      	  if(volFlux > 0 && m_proppantDataPtr->m_diameter < faceApertures[fc]){  // flowing out - ignored if flowing in
      	    //Kc = volFlux;

            realT vf = faceConcentrations[fc];
            R1Tensor eZ = faceNormals[fc];
            R1Tensor dN; dN.Cross(dT,eZ); dN.Normalize();
            R1Tensor slipVelocity;
            CalculateSlipVelocity(vf, eZ,mixedFluidVelocity[fc], slipVelocity);
            // realT mf = vf/(vf+ m_proppantDataPtr->m_rhoRatio*(1-vf)); // mass fraction
            realT driftFlux = faceApertures[fc]*Tnorm*(1-vf)*Dot(slipVelocity,dN);
      
            // record matrix entry for face
            face_dof(0) = face_trilinos_index[fc];
            face_K(0,0) = volFlux+driftFlux;
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


// Slip velocity = V_Particle-V_fluid
// vf = proppant volume fraction
void ParallelPlateProppantSolver::CalculateSlipVelocity(const realT vf, const R1Tensor& normal, const R1Tensor& mixtureVelocity, R1Tensor& slipVelocity){

    slipVelocity = 0.0;

    if(vf < m_proppantDataPtr->m_maxVf && vf > 0.0){
      // settling velocity
      realT sz = m_slurryModelPtr->HinderedSettlingSpeed(vf);
      //slipVelocity[2] = -sz*(1.0-normal[2]);

      slipVelocity = -sz*m_upVector;
      slipVelocity -= Dot(slipVelocity,normal)*normal;  // remove component into face


      // wall drag/friction
      if(m_useCollisionalSlipVelocity){
    	  realT lambda = SuperadvectiveTransportFactor(vf);
    	  realT scaledVf = vf/( std::max( m_proppantDataPtr->m_maxVf - vf,TINY ));
    	  lambda /= (1+scaledVf*scaledVf);
    	  slipVelocity += ( (1-lambda*vf)/(1-vf) )*mixtureVelocity;
      }
    } else if (vf >= m_proppantDataPtr->m_maxVf){
  	  slipVelocity += ( 1.0/(1-vf) )*mixtureVelocity;
    }

}



void ParallelPlateProppantSolver:: CalculateScreenout(PhysicalDomainT&  domain,
                                                      SpatialPartition& partition, const realT& time,const realT& dt){
    rArray1d& faceConcentrations  = domain.m_feFaceManager.GetFieldData<realT>(concentrationFieldName_);
    const rArray1d&  faceAreas = domain.m_feFaceManager.GetFieldData<realT>( PS_STR::FaceAreaStr );
    rArray1d& faceApertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
    rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);
    Array1dT<R1Tensor>& faceNormals = domain.m_feFaceManager.GetFieldData<R1Tensor>(FaceNormalStr);
    Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

    rArray1d& edgeExcessPackVolume = domain.m_feEdgeManager.GetFieldData<realT>("ExcessProppantPackVolume");
    edgeExcessPackVolume = 0.0; // reset edgeExcess

    const Array1dT<R1Tensor>& nodePosition = *domain.m_feNodeManager.m_refposition;

    Array1dT<R1Tensor>& particlePhaseVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr+"_proppant");

	// search for neighbors of "screened out" faces - make sure vol frac is > 0
	//////////////////////////////////////////////////////////////////////////
	for( lSet::const_iterator fi=m_faceSet->begin() ; fi!=m_faceSet->end() ; ++fi )
	{
	  localIndex fc = *fi;
	  realT facePackVf = facePackVfs[fc];
      if(facePackVf >= 1.0 ){
          lArray1d& faceEdges = domain.m_feFaceManager.m_toEdgesRelation[fc];
          localIndex numEdges = faceEdges.size();
          for(unsigned  i =1; i < numEdges; ++i){
          	localIndex eg = domain.m_feEdgeManager.GetParentIndex(faceEdges[i]);
          	lArray1d& edgeFaces = m_edgesToFaces[eg];
          	unsigned numFaces =  edgeFaces.size();
          	for(unsigned ii =0; ii < numFaces; ++ii){
          	  localIndex nbrfc = edgeFaces[ii];
          	  if( facePackVfs[nbrfc] <= 0 ) facePackVfs[nbrfc] = TINY;
          	}
          }
      }
	}

	// loop over faces
	//////////////////
	for( lSet::const_iterator fi=m_faceSet->begin() ; fi!=m_faceSet->end() ; ++fi )
	{
	  localIndex fc = *fi;
      realT facePackVf = facePackVfs[fc];
      bool doScreenout = (facePackVf > 0.0 && facePackVf < 1.0);
	  R1Tensor center = faceCenters[fc];

      lArray1d& faceEdges = domain.m_feFaceManager.m_toEdgesRelation[fc];
      localIndex numEdges = faceEdges.size();

      if(!doScreenout){
    	// check if face is at fracture base (lowest edge is non-flow)
    	localIndex eg = domain.m_feEdgeManager.GetParentIndex(faceEdges[0]);
        localIndex lowestEg = eg;
        R1Tensor edgeCenter; domain.m_feEdgeManager.EdgeCenter(domain.m_feNodeManager,eg,edgeCenter);
        realT minZ = Dot(edgeCenter,m_upVector);
        for(unsigned  i =1; i < numEdges; ++i){
          eg = domain.m_feEdgeManager.GetParentIndex(faceEdges[i]);
          domain.m_feEdgeManager.EdgeCenter(domain.m_feNodeManager,eg,edgeCenter);
          realT z = Dot(edgeCenter,m_upVector);
          if(z < minZ){
        	minZ = z;
        	lowestEg = eg;
          }
        }
        if(m_edgesToFaces[lowestEg].size() ==  1 ){
          doScreenout = true;
        }

      }

      if(doScreenout){
        const lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[fc];

	    realT area = faceAreas[fc];
	    realT vol = area*faceApertures[fc];

        R1Tensor eN = faceNormals[fc];

        // min max z coordinates (relative to up vector)
        realT minZ = 1e64; // unlikely value
        realT maxZ = -1e64; // unlikely value
        for( size_t a=0 ; a<nodeList.size(); ++a ){
            localIndex nd = domain.m_feNodeManager.GetParentIndex(nodeList[a]);
        	realT z = Dot(nodePosition[nd],m_upVector);
        	if (z < minZ) minZ =  z;
        	if (z > maxZ) maxZ = z;
        }
        realT lz = maxZ-minZ+TINY;

        // fill rate is based on average cross section of element - may want to change this
        //  dx/dt * (phi_settled - phi_suspended) = vz*phi_suspended  i.e.  settling rate = flux in
        //  dvolFrac/dt = average cross section  * dh/dt/area (where h is the fill height)
        //  average cross section = area*|zxn|/lz
        // dvolFrac = (dt*dh/dt)*|zxn|/lz
        realT deltaPhi = (faceConcentrations[fc] < m_proppantDataPtr->m_maxVf)? std::max(m_proppantDataPtr->m_maxVf - faceConcentrations[fc],1e-3) : 1e-3;
        realT dxdt = -Dot(particlePhaseVelocity[fc],m_upVector)*faceConcentrations[fc]/deltaPhi;
        realT cosTheta = Dot(eN,m_upVector);
        realT sinTheta = sqrt(1-cosTheta*cosTheta);
        realT deltaVf = dt*dxdt*sinTheta/lz;  // change in proppant volume fraction

        realT newFacePackVf = facePackVfs[fc];
        newFacePackVf += deltaVf;        // pack volume frac increased by this much
        if (newFacePackVf > 1.0){
          realT excessPackVol = vol*(newFacePackVf - 1.0); // excess
          facePackVfs[fc] = 1.0;
          faceConcentrations[fc] = m_proppantDataPtr->m_maxVf;
          R1Tensor& faceCenter = faceCenters[fc];
          // push excess to edges in fill direction
          /////////////////////////////////////////////
          rArray1d weights(numEdges,0.0);
          realT wTotal = 0.0;
          for(unsigned  i =0; i < numEdges; ++i){
        	  localIndex eg = domain.m_feEdgeManager.GetParentIndex(faceEdges[i]);
              R1Tensor edgeCenter; domain.m_feEdgeManager.EdgeCenter(domain.m_feNodeManager,eg,edgeCenter);
        	  R1Tensor l = edgeCenter-faceCenter;
        	  if( Dot(l,m_upVector) > 0 && m_edgesToFaces[eg].size() > 1 ){
            	  R1Tensor edgeVector;
            	  domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,edgeVector);
        		  realT w = (Cross(edgeVector,m_upVector)).L2_Norm();
        		  weights[i] = w;
        		  wTotal += w;
        	  }
          }
          // distribute excess
          for(unsigned i =0; i < numEdges; ++i){
        	  localIndex eg = domain.m_feEdgeManager.GetParentIndex(faceEdges[i]);
        	  edgeExcessPackVolume[eg] += excessPackVol*weights[i]/wTotal;
          }
        } else if(newFacePackVf < 0.0){
          facePackVfs[fc] = 0.0;

        } else {
          //  std::cout << "Got Here! Doing screenout " <<  particlePhaseVelocity[fc] << std::endl;

        /*
         // don't think we need to do this - as particles settle out they should be replaced by particles above

          realT avSuspendedProppantConc = faceConcentrations[fc]*(1-newFacePackVf);  // vf of suspended proppant averaged over entire cell

          avSuspendedProppantConc -= deltaVf*m_proppantDataPtr->m_maxVf; // this much dropped out of suspension

          if(avSuspendedProppantConc < 0 ){ // not enough in cell - may need to adjust timestep
        	newFacePackVf -= deltaVf;
        	newFacePackVf += faceConcentrations[fc]*(1-newFacePackVf)/m_proppantDataPtr->m_maxVf; // everything drops out
            faceConcentrations[fc] = 0;
          } else {
            faceConcentrations[fc] =  std::min( avSuspendedProppantConc/(1-newFacePackVf), m_proppantDataPtr->m_maxVf);
          }
          */
          facePackVfs[fc] = newFacePackVf;
        }


      }
	}

	// redistribute excess pack volume on edges
	///////////////////////////////////////////

	// loop over edges
	std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
	for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr ){
		localIndex eg = itr->first;
		realT excess = edgeExcessPackVolume[eg];
		if(excess > 0.0){
			lArray1d& edgeFaces = itr->second;
			unsigned numFaces =  edgeFaces.size();
			// order edge faces by z height
			std::set< std::pair<realT,localIndex> > centerHeights;
			for(unsigned ii =0; ii < numFaces; ++ii){
				localIndex fc = edgeFaces[ii];
				realT h = Dot(faceCenters[fc],m_upVector);
				centerHeights.insert(std::pair<realT,localIndex>(h,fc) );
			}

			// distribute excess across faces in order of z height
			for(std::set< std::pair<realT,localIndex> >::iterator hItr = centerHeights.begin();
					excess > 0.0 && hItr != centerHeights.end() ; ++hItr){
				localIndex fc = hItr->second;

				realT faceVolume = faceAreas[fc]*faceApertures[fc];
				realT facePackVf = facePackVfs[fc];
				//realT suspendedMaterialEquivalentVf = (faceConcentrations[fc]/m_proppantDataPtr->m_maxVf)*(1-facePackVf); // vol occupied by suspended material
				realT availVol = (1-facePackVf)*faceVolume; // -suspendedMaterialEquivalentVf)*faceVolume;
				if(availVol < excess){
					excess -= availVol;
					facePackVfs[fc] = 1.0;
					faceConcentrations[fc] = m_proppantDataPtr->m_maxVf; // this doesn't destroy mass as suspended material is automatically added to vf.
				} else {
					realT newFacePackVf = facePackVf + excess/faceVolume;
					excess = 1e-64;
					//realT newFaceConc = faceConcentrations[fc] * (1-facePackVf)/(1 - newFacePackVf + TINY);
					//faceConcentrations[fc] = std::min( newFaceConc, m_proppantDataPtr->m_maxVf); // std::min should not be required
					facePackVfs[fc] = newFacePackVf;
				}
			}
			edgeExcessPackVolume[eg] = excess; // record any remaining excess pack volume
		}
	}
}




// calculates the erosion rate and updates the eroded aperture and the actual aperture
void ParallelPlateProppantSolver:: CalculateErosionRate(PhysicalDomainT&  domain,
                                                      SpatialPartition& partition, const realT& time,const realT& dt){
    const rArray1d& faceConcentrations  = domain.m_feFaceManager.GetFieldData<realT>(concentrationFieldName_);
    const Array1dT<R1Tensor>& particlePhaseVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr+"_proppant");

    rArray1d& faceApertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
    rArray1d& erodedApertures = domain.m_feFaceManager.GetFieldData<realT>( "ErodedAperture" );
    rArray1d& erosionRates = domain.m_feFaceManager.GetFieldData<realT>( "ErosionRate" );

    const rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);


    const Array1dT<R1Tensor>& faceVelocities = domain.m_feFaceManager.GetFieldData<R1Tensor>( FluidVelocityStr );

	// loop over faces
	//////////////////
	for( lSet::const_iterator fi=m_faceSet->begin() ; fi!=m_faceSet->end() ; ++fi )
	{
		  localIndex fc = *fi;
		  realT newErosionRate = 0;
		  realT h = faceApertures[fc];

		  // estimate pressure gradient from
		  realT velMag = faceVelocities[fc].L2_Norm();
          realT qMag = velMag*h;
		  realT K = m_flowSolverPtr->CalculateFacePermeability(h,qMag,faceConcentrations[fc],facePackVfs[fc]);


		  // shear stress on fracture wall
		  // deltaPMag = qMag/K;  // estimate pressure gradient from
		  // tau = 0.5*deltaPMag*h; // force on fracture faces balanced by pressure gradient
		  realT tau = 0.5*qMag*h/(K+TINY);


		  newErosionRate = 2*m_ErosionModel.ErosionRate(tau, faceConcentrations[fc]); //

		  /*
		  if(velMag > 80){
			  std::cout << "K "<< K << std::endl;
			  std::cout << "qMag "<< qMag << std::endl;
			  std::cout << "h "<< h << std::endl;
			  std::cout << "tau "<< tau << std::endl;
			  std::cout << "newErosionRate "<< newErosionRate << std::endl;
		  }
          */

		  realT deltaAperture =  0.5*dt*(erosionRates[fc]+newErosionRate);


          // supress spurious erosion rate spikes near max volume fraction
		  if(faceConcentrations[fc] > 0.95*m_proppantDataPtr->m_maxVf){
			  newErosionRate = 0.0;
			  deltaAperture = 0.0;
		  }

		  // limit max aperture change per timestep
		  deltaAperture = std::min(deltaAperture,0.1*faceApertures[fc]);



		  //faceApertures[fc] += deltaAperture;
		  erodedApertures[fc] += deltaAperture;
		  erosionRates[fc] = newErosionRate;
	}
}


realT ParallelPlateProppantSolver::SuperadvectiveTransportFactor(const realT vf){
    // empirical model that captures entrainment in center of fracture
	// see Barree & Conway 1994
    return 1.27- pow(fabs(vf-0.1),1.5);
}
  
//////////////////////////////////////////
/// Register solver in the solver factory
REGISTER_SOLVER( ParallelPlateProppantSolver )

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

