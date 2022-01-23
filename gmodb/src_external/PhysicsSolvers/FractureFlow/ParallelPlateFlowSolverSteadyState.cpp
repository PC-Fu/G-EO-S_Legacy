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
 * @file ParallelPlateFlowSolverSteadyState.cpp
 * @author walsh24
 * @date June 1, 2011
 */

#include "ParallelPlateFlowSolverSteadyState.h"
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
using namespace PPFS;

namespace{
	realT TINY = 1e-64;	
}


ParallelPlateFlowSolverSteadyState::ParallelPlateFlowSolverSteadyState( const std::string& name, ProblemManagerT* const pm):
ParallelPlateFlowSolverBase(name,pm),
m_faceSet(NULL),
m_numFaces(0),
m_faceDofMap(),
m_edgesToFaces(),
m_useMLPrecond(false),
m_flux_bc_counter(0),
m_flux_bc_trilinos_host(),
m_flux_bc_trilinos_index(),
m_flux_bc_pressure(),
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
m_numerics()
{
  //ReadXML(hdn);
  ++m_instances; 
  m_TrilinosIndexStr = "PPFS_SS_" +  toString<int>(m_instances) + "_GlobalDof";
}

ParallelPlateFlowSolverSteadyState::~ParallelPlateFlowSolverSteadyState()
{
  // TODO Auto-generated destructor stub
}

void ParallelPlateFlowSolverSteadyState::ReadXML( TICPP::HierarchicalDataNode* hdn )
{ 
  ParallelPlateFlowSolverBase::ReadXML(hdn);
  // flags
  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("updateAperture",false);

  // Linear Solver
  m_numerics.m_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-10);
  m_numerics.m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);
  m_numerics.m_verboseFlag = hdn->GetAttributeOrDefault<bool>(VerboseStr,true);
  
  m_useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);

}


void ParallelPlateFlowSolverSteadyState::RegisterFields( PhysicalDomainT& domain )
{
	
  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(VolumetricFluxStr,true,true);
    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::pressure>();    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::density>(); 
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::mass>(); 
  
  domain.m_feEdgeManager.AddKeylessDataField<realT>(PermeabilityStr,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(VolumetricFluxStr,true,true);
  

  //domain.m_edgeManager.AddKeylessDataField<int>("FlowFaceCount",true,true);// debug
  

  domain.m_feFaceManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);
  
}

void ParallelPlateFlowSolverSteadyState::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  EdgeManagerT& edgeManager = domain.m_feEdgeManager;
  
  m_faceSet = &(faceManager.GetSet(m_flowFaceSetName));
  m_numFaces = m_faceSet->size();
  
  // build face-dof map

  lSet::const_iterator si=m_faceSet->begin();
  for(localIndex i =0; i < m_numFaces; ++i, ++si){
  	localIndex f = *si;
  	m_faceDofMap[f] = i;
  }
  
 //  iArray1d& ffCount = domain.m_edgeManager.GetFieldData<int>("FlowFaceCount"); // debug
  

  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    localIndex numEdges = faceManager.m_toEdgesRelation[*kf].size();
    for(localIndex a =0; a < numEdges; ++a){
      localIndex eg = faceManager.m_toEdgesRelation[*kf][a];
    	
      lSet& edgeFaces = edgeManager.m_toFacesRelation[eg];
      lArray1d edgeList;
      
      for( lSet::iterator edgeFace=edgeFaces.begin() ; edgeFace!=edgeFaces.end() ; ++edgeFace ){
        if(isMember(*edgeFace,m_faceDofMap)){
          edgeList.push_back(*edgeFace);
        }
      }
      m_edgesToFaces[eg] = edgeList;
     // ffCount[eg] = edgeList.size();
    }
  }


  GenerateParallelPlateGeometricQuantities( domain, 0,0 );


  using namespace BoundaryConditionFunctions;
  ApplyMultiSetBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::SinglePartitionBC_NeighborUpdate,
                                           domain, domain.m_feEdgeManager,
                                           SinglePartitionPeriodicBoundaryCondition::BoundaryConditionName(), 0.0 );
    
}


/**
 * @author walsh24
 * @brief  Update neighbors for single partition periodic boundary condition.
 *
 */
void ParallelPlateFlowSolverSteadyState::
SinglePartitionBC_NeighborUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                                BoundaryConditionBase* bc, realT time  ){

  SinglePartitionPeriodicBoundaryCondition* spbc = dynamic_cast<SinglePartitionPeriodicBoundaryCondition*> (bc);
  if(spbc){
   spbc->SetNeighborMaps(domain);
   std::cout << "Neighbor update" << std::endl;
  }

}

/**
 * @author walsh24
 * @brief  Update sparsity pattern for single partition periodic boundary condition.
 *
 */
void ParallelPlateFlowSolverSteadyState::
SinglePartitionBC_Sparsity(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                                BoundaryConditionBase* bc, realT time  ){

   std::cout << "SparsityA" << std::endl;
  SinglePartitionPeriodicBoundaryCondition* spbc = dynamic_cast<SinglePartitionPeriodicBoundaryCondition*> (bc);
  if(spbc){
   std::cout << "Sparsity" << std::endl;
    iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
    std::map<localIndex, localIndex>& edgeSetMap = spbc->m_edgeNeighborMapA;

    std::map<localIndex, localIndex>::iterator itr = edgeSetMap.begin();
    std::map<localIndex, localIndex>::iterator iend = edgeSetMap.end();

    std::vector<int> dofIndex(2);

    // loop over edges
    for(;itr!=iend; ++itr){
      localIndex edgeA = itr->first;
      localIndex edgeB = itr->second;

      // fixme assumes that there is only one face per periodic boundary edge
      // otherwise would have to loop over all possible combinations of faceA and faceB
      // not hard to implement here - but painful when calculating relative permeabilities later
      // and likely redundant once proper periodic bcs come in
      std::map<localIndex,lArray1d>::iterator itrA = m_edgesToFaces.find(edgeA);
      lArray1d& facesA = itrA->second;
      const localIndex kfA = facesA[0];
      dofIndex[0] = trilinos_index[kfA];

      std::map<localIndex,lArray1d>::iterator itrB = m_edgesToFaces.find(edgeB);
      lArray1d& facesB = itrB->second;
      const localIndex kfB = facesB[0];
      dofIndex[1] = trilinos_index[kfB];

      sparsity->InsertGlobalIndices(2,
                                    &dofIndex.front(),
                                    2,
                                    &dofIndex.front());

    }
  }

}

/**
 * @author walsh24
 * @brief  Single partition periodic boundary condition.
 *
 */
void ParallelPlateFlowSolverSteadyState::
SinglePartitionBC(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                                BoundaryConditionBase* bc, realT time  ){

   std::cout << "SinglePartitionBCA" << std::endl;
  SinglePartitionPeriodicBoundaryCondition* spbc = dynamic_cast<SinglePartitionPeriodicBoundaryCondition*> (bc);
  if(spbc){
   std::cout << "SinglePartitionBC" << std::endl;

    iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
    rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);

    const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
    const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );

    std::map<localIndex, localIndex>& edgeSetMap = spbc->m_edgeNeighborMapA;

    std::map<localIndex, localIndex>::iterator itr = edgeSetMap.begin();
    std::map<localIndex, localIndex>::iterator iend = edgeSetMap.end();

    Epetra_IntSerialDenseVector  edgeDofIndex (2);
    Epetra_SerialDenseVector     edge_rhs     (2);
    Epetra_SerialDenseMatrix     edge_matrix  (2,2);

    // loop over edges
    for(;itr!=iend; ++itr){
      localIndex eg = itr->first;
      localIndex egb = itr->second;

      // fixme assumes that there is only one face per periodic boundary edge
      // otherwise would have to loop over all possible combinations of faceA and faceB
      // and remove contributions from previous perm calculations
      std::map<localIndex,lArray1d>::iterator itrA = m_edgesToFaces.find(eg);
      lArray1d& facesA = itrA->second;
      const localIndex kf = facesA[0];

      std::map<localIndex,lArray1d>::iterator itrB = m_edgesToFaces.find(egb);
      lArray1d& facesB = itrB->second;
      const localIndex kfb = facesB[0];

      // calculate edge permeability
      R1Tensor edgeCenter,edgeCenterB, la, lb;
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter );
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, egb , edgeCenterB );

      la = edgeCenter;
      la -= faceCenters[kf];

      lb = edgeCenterB;
      lb -= faceCenters[kfb];

      realT w = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, eg);

      realT app = (apertures[kf] <  m_max_aperture)? apertures[kf]: m_max_aperture;
      realT appb = (apertures[kfb] <  m_max_aperture)? apertures[kfb]: m_max_aperture;
      if (app <  m_min_aperture) app = m_min_aperture;
      if (appb <  m_min_aperture) appb = m_min_aperture;

      realT kappa = CalculatePermeability( la.L2_Norm() , lb.L2_Norm(),
                                           app,appb,
                                           w,m_mu,m_SHP_FCT);

      edgePermeabilities[eg] = kappa;
      edgePermeabilities[egb] = kappa;

      // build stiffness matrix
      if( kappa > 0.0 ){

        // diagonal
        edge_matrix(0,0) = kappa;
        edge_matrix(1,1) = kappa;

        // cross terms
        edge_matrix(0,1) = -kappa;
        edge_matrix(1,0) = -kappa;

        edgeDofIndex[0] = trilinos_index[kf];
        edgeDofIndex[1] = trilinos_index[kfb];

        matrix->SumIntoGlobalValues(edgeDofIndex, edge_matrix);
      }

    }
  }


}

/**
 * @author walsh24
 * @brief  Update velocity pattern for single partition periodic boundary condition.
 *
 */
void ParallelPlateFlowSolverSteadyState::
SinglePartitionBC_VelocityUpdate(PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
                                 BoundaryConditionBase* bc, realT time  ){

  SinglePartitionPeriodicBoundaryCondition* spbc = dynamic_cast<SinglePartitionPeriodicBoundaryCondition*> (bc);
  if(spbc){
   std::cout << "Velocity update" << std::endl;

      rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);

      const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

      Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);

      Array1dT<R1Tensor>& faceVolFlux = domain.m_feFaceManager.GetFieldData<R1Tensor>(VolumetricFluxStr);
      const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

      iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

      std::map<localIndex, localIndex>& edgeSetMap = spbc->m_edgeNeighborMapA;

      std::map<localIndex, localIndex>::iterator itr = edgeSetMap.begin();
      std::map<localIndex, localIndex>::iterator iend = edgeSetMap.end();

      for( ; itr!=iend ; ++itr ) {
        localIndex eg = itr->first;
        localIndex egb = itr->second;

        std::map<localIndex,lArray1d>::iterator itrA = m_edgesToFaces.find(eg);
        std::map<localIndex,lArray1d>::iterator itrB = m_edgesToFaces.find(egb);
        if(itrA != m_edgesToFaces.end() && itrB != m_edgesToFaces.end()  ){

          lArray1d& facesA = itrA->second;
          lArray1d& facesB = itrB->second;

          R1Tensor edgeCenter, edgeCenterB;
          domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter );
          domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, egb , edgeCenterB );

          // nb only uses first face to calculate volume flux
          // flux will be incorrect for other edges if at junction.
          {
            localIndex fc = facesA[0];
            localIndex fcb = facesB[0];

            const realT Pa = pressures[fc];
            const realT Pb = pressures[fcb];

            realT volumeFlux =  edgePermeabilities[eg]*(Pa-Pb); // flux out of a into b
            volFlux[eg] = volumeFlux;
            volFlux[egb] = -volumeFlux;

            // vecFlux a

            R1Tensor vecFlux;

            vecFlux = edgeCenter;
            vecFlux -= faceCenters[fc];

            vecFlux.Normalize();
            vecFlux *= volumeFlux;

            if(is_ghost[fc] < 0) faceVolFlux[fc] += vecFlux;

            // vecFlux b
            vecFlux = edgeCenterB;
            vecFlux -= faceCenters[fcb];

            vecFlux.Normalize();
            vecFlux *= -volumeFlux;

            if(is_ghost[fcb] < 0) faceVolFlux[fcb] += vecFlux;


          }
        }
      }

  }

}



void ParallelPlateFlowSolverSteadyState:: SetupSystem (PhysicalDomainT&  domain,
                                                SpatialPartition& partition, const realT& time)
{
	
  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  
  // count local dof
  ///////////////////////////////
  
  // local rows
  int n_local_rows = 0;
  
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
  	if( is_ghost[*kf] < 0 ){
  		++n_local_rows;
  	} 
  }
  
  // local flux dof
  FluxControlledBoundaryCondition_Setup(domain, n_local_rows,time);
  
  // determine the global/local degree of freedom distribution.
  ////////////////////////////////////////////////////////////

  std::vector<int> gather(n_mpi_processes);
  std::vector<int> cum_global_rows(n_mpi_processes);

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
    cum_global_rows[p] = n_global_rows; 
  }
  
  // create trilinos dof indexing
  //////////////////////////////////
  unsigned local_count = 0;
  // faces
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      trilinos_index[*kf] = first_local_row+local_count;
      local_count++;
    }
    else
    {
      trilinos_index[*kf] = -INT_MAX;
    }
  }
  
  // flux bc's    
  localIndex numBCs = m_flux_bc_trilinos_host.size();
  iArray1d offsets(n_mpi_processes);  
  offsets = -1;             
  for(localIndex i =0; i < numBCs; ++i){
    int proc = m_flux_bc_trilinos_host(i);
    m_flux_bc_trilinos_index(i) = cum_global_rows[proc]+ offsets(proc);
    offsets(proc) -=  1;
    if(proc == this_mpi_process){
      local_count++; // nb local_count no longer increases with trilinos index 
    }
  }

  assert(static_cast<int>(local_count) == n_local_rows);

  partition.SynchronizeFields(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);

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
  	  unsigned int numFaces = itr->second.size();
      if( numFaces > 1){
        std::vector<int> dofIndex (numFaces);
        for(unsigned i=0; i<numFaces; ++i)
        {
          localIndex kf = itr->second[i];
          dofIndex[i] = trilinos_index[kf];
        }

        sparsity->InsertGlobalIndices(dofIndex.size(),
                                      &dofIndex.front(),
                                      dofIndex.size(),
                                      &dofIndex.front());
      }
    }
  }
  
  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition_Sparsity,
                                domain, domain.m_feEdgeManager, "NetFlux", time );

  using namespace BoundaryConditionFunctions;
  ApplyMultiSetBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::SinglePartitionBC_Sparsity,
                                           domain, domain.m_feEdgeManager,
                                           SinglePartitionPeriodicBoundaryCondition::BoundaryConditionName(), time );
     
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();
  
}

/* Assemble */

void ParallelPlateFlowSolverSteadyState :: Assemble (PhysicalDomainT&  domain,
                                                     SpatialPartition& partition ,
                                                     const realT& time)

{
	
    // (re-)init linear system

  matrix   = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*sparsity));
  solution = Teuchos::rcp(new Epetra_FEVector(*row_map));
  rhs      = Teuchos::rcp(new Epetra_FEVector(*row_map));

  // basic face data ( = dof data for our problem)

  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

  iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  
  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  
    
  // loop over edges
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0 ){
      int numFaces = itr->second.size();
   	
      Epetra_IntSerialDenseVector  edgeDofIndex (numFaces);
      Epetra_SerialDenseVector     edge_rhs     (numFaces);
      Epetra_SerialDenseMatrix     edge_matrix  (numFaces,numFaces);
   	
      if( numFaces == 2){


  		
        eg = itr->first;
        localIndex kf = itr->second[0];
        localIndex kfb = itr->second[1]; 	
  	
        // calculate edge permeability
        R1Tensor edgeCenter, la, lb;
        domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter );

        la = edgeCenter;
        la -= faceCenters[kf];

        lb = edgeCenter;
        lb -= faceCenters[kfb];

        realT w = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, eg);
        
        realT app = (apertures[kf] <  m_max_aperture)? apertures[kf]: m_max_aperture;
        realT appb = (apertures[kfb] <  m_max_aperture)? apertures[kfb]: m_max_aperture;
        if (app <  m_min_aperture) app = m_min_aperture;
        if (appb <  m_min_aperture) appb = m_min_aperture;

        realT kappa = CalculatePermeability( la.L2_Norm() , lb.L2_Norm(),
                                             app,appb,
                                             w,m_mu,m_SHP_FCT);
                                           
        edgePermeabilities[eg] = kappa;
      
  	    // build stiffness matrix       
        if( kappa > 0.0 ){
      	  
  	      // diagonal
  	      edge_matrix(0,0) = kappa;
  	      edge_matrix(1,1) = kappa;
  	    	 
  	      // cross terms
  	      edge_matrix(0,1) = -kappa;
  	      edge_matrix(1,0) = -kappa;
  	  
  	      edgeDofIndex[0] = trilinos_index[kf];
  	      edgeDofIndex[1] = trilinos_index[kfb];
  	  
  	      matrix->SumIntoGlobalValues(edgeDofIndex, edge_matrix);
        }

      
      } else if(numFaces > 2){

  	
  	    eg = itr->first;
        lArray1d& faces = itr->second; 
  	    R1Tensor edgeCenter;
        domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter );
        realT w = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, eg);
      
        rArray1d kappas(numFaces,0.0);
        realT kappaSum =0.0;
        for(int i =0; i < numFaces; ++i){
  	      localIndex kf = faces[i];
  	      edgeDofIndex[i] = trilinos_index[kf];
  	      R1Tensor la = edgeCenter;
  	      la -= faceCenters[kf];
  	      //          realT app = (apertures[kf] <  m_max_aperture)? apertures[kf]: m_max_aperture;
  	      kappas[i] = CalculatePermeability( la.L2_Norm(), apertures[kf], w,m_mu,m_SHP_FCT);
  	      kappaSum += kappas[i];
        }
        realT invKappaSum = 1.0/(kappaSum + TINY);

        // add contributions to stiffness matrices
        for(int i =0; i < numFaces; ++i){
          realT selfTerm = kappas[i]*(1.0 - kappas[i]*invKappaSum);
          // diagonal
          edge_matrix(i,i) =  selfTerm;

          // cross terms
          for(int ii =0; ii < numFaces; ++ii){
            if(i!=ii){
              realT crossTerm = -kappas[i]*kappas[ii]*invKappaSum;
              edge_matrix(i,ii) =  crossTerm;
            }
          }
        }
  	  
        matrix->SumIntoGlobalValues(edgeDofIndex, edge_matrix);
      } // numFaces > 2 ?
    } // ghost edge
   
  } // edge loop
  
  
  // boundary conditions
  using namespace BoundaryConditionFunctions;
  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::PressureBoundaryCondition,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time );

  ApplyMultiSetBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::SinglePartitionBC,
                                           domain, domain.m_feEdgeManager,
                                           SinglePartitionPeriodicBoundaryCondition::BoundaryConditionName(), time );
                              
  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition,
                                domain, domain.m_feEdgeManager, "NetFlux",time );
  	  

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  //rhs->Print(std::cout);
  //std::cout << matrix->NormInf() << std::endl;
  //EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*matrix);
  //exit(0);
}


/// Apply a pressure boundary condition to a given set of edges
void ParallelPlateFlowSolverSteadyState::PressureBoundaryCondition(PhysicalDomainT& domain,
                                                                    ObjectDataStructureBaseT& object ,
                                                                    BoundaryConditionBase* bc, const lSet& set, realT time){
 

  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  //iArray1d& face_is_ghost       = domain.m_faceManager.GetFieldData<FieldInfo::ghostRank>();
  iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
 
  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
 
  	  	
//  int numBCs = set.size();
  
  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);
  Epetra_SerialDenseMatrix     face_matrix(1,1);
  	
  // loop over edges, find permeabilities and apply boundary conditions.

  lSet::const_iterator eg=set.begin() ;

  for( localIndex i =0; i < set.size() ; ++i, ++eg ){

    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){
      lArray1d& faces = itr->second; 
      
      R1Tensor edgeCenter;
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *eg , edgeCenter );
      realT value = bc->GetValue(domain.m_feEdgeManager,eg,time);
      for(size_t ii =0; ii < faces.size(); ++ii){
      	
      	
      	const localIndex kf = faces[ii];

        // calculate edge permeability for face
      	/*old - use vector from face center to edge center
        R1Tensor la = edgeCenter; la -= faceCenters[kf];
        realT w = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, *eg);
        */

      	// new project face-edge vector normal to face
        R1Tensor la = edgeCenter; la -= faceCenters[kf];
      	R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,*eg,dT); // edge vector
      	realT w = dT.Normalize();
      	la = la- (la*dT)*dT;

        const realT la_norm = la.L2_Norm();
        const realT app = apertures[kf];
        const realT kappa = CalculatePermeability( la_norm , app, w,m_mu,m_SHP_FCT);
        
        if(ii==0) edgePermeabilities[*eg] = kappa;  // only assign first permeability (at junction)
        
        face_dof(0) = trilinos_index[kf]; 
        face_matrix(0,0) = kappa;
        matrix->SumIntoGlobalValues(face_dof, face_matrix);
    
        face_rhs[0] =  value*kappa;
        rhs->SumIntoGlobalValues(face_dof, face_rhs);
    
      }
    }
  }

}




/* Solve */

void ParallelPlateFlowSolverSteadyState:: Solve (PhysicalDomainT&  domain,
                                               SpatialPartition& partition)
{

  // face fields
  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d& pressure       = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  // set initial guess
  int dummy;
  double* local_solution = NULL;

  solution->ExtractView(&local_solution,&dummy);
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      long long int rowTmp = static_cast<long long int>(trilinos_index[*kf]);
      int lid = row_map->LID(rowTmp);
      local_solution[lid] = pressure[*kf];
    }
  }
  
  unsigned int numBCs = m_flux_bc_trilinos_host.size();
  if(numBCs > 0){
    for( unsigned int i =0; i < numBCs;++i){
      if(m_flux_bc_trilinos_host[i] == this_mpi_process){
        long long int rowTmp = static_cast<long long int>(m_flux_bc_trilinos_index[i]);
        int lid = row_map->LID(rowTmp);
  	    local_solution[lid] = m_flux_bc_pressure[i];
  	  }
    }
  }

  // krylov solver

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
  if(m_useMLPrecond){
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*matrix, MLList);
  }


  // Solver
  AztecOO solver(problem);

          if(m_useMLPrecond) solver.SetPrecOperator(MLPrec);

          solver.SetAztecOption(AZ_solver,AZ_bicgstab);
          solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
          solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
          solver.SetAztecOption(AZ_conv,AZ_rhs);
          //solver.SetAztecOption(AZ_output,AZ_none);
          
          /*
          // Symmetric preconditioner
          // - domain decomposition preconditioner
          // - ICC factorization on each subdomain
          solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
          solver.SetAztecOption(AZ_overlap,0);
          solver.SetAztecOption(AZ_subdomain_solve, AZ_icc);
          
          
          solver.SetAztecOption(AZ_conv,AZ_rhs);
          solver.SetAztecOption(AZ_solver,AZ_cg);
         // solver.SetAztecOption(AZ_precond, AZ_Jacobi);*/

          if(m_numerics.m_verboseFlag){
        	  solver.SetAztecOption(AZ_output,AZ_all);
          }  else {
        	  solver.SetAztecOption(AZ_output,AZ_none);
          }
          

          solver.Iterate(m_numerics.m_maxIters,m_numerics.m_tol);
  //solution->Print(std::cout);

  // destroy the preconditioner
  if(m_useMLPrecond){
     delete MLPrec;
  }

  // copy solution to faces
  ////////////////////////
  
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      long long int rowTmp = static_cast<long long int>(trilinos_index[*kf]);
      int lid = row_map->LID(rowTmp);
      pressure[*kf] = local_solution[lid];
    }
  }
  
  // distribute flux BC pressures
  ///////////////////////////
  if(numBCs > 0){
    std::vector<double> aBuffer(numBCs,0.0);
    for( unsigned int i =0; i < numBCs;++i){
      if(m_flux_bc_trilinos_host[i] == this_mpi_process){
        long long int rowTmp = static_cast<long long int>(m_flux_bc_trilinos_index[i]);
        int lid = row_map->LID(rowTmp);
        aBuffer[i] = local_solution[lid];
      }
    }
    m_epetra_comm.SumAll(&aBuffer[0],&m_flux_bc_pressure[0],numBCs); // fixme not very efficient
  }

  // re-sync ghost nodes
  partition.SynchronizeFields(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);
}


void ParallelPlateFlowSolverSteadyState::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  //syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(FaceCenterStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::pressure>::Name());
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(ApertureStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_TrilinosIndexStr);

  //syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(VolumetricFluxStr);
  syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(PermeabilityStr);

  syncedFieldsB.clear();
  syncedFieldsB[PhysicalDomainT::FiniteElementFaceManager].push_back(VolumetricFluxStr);

  partition.SetBufferSizes(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);
  partition.SetBufferSizes(syncedFieldsB, CommRegistry::steadyStateParallelPlateFlowSolverB);

}


double ParallelPlateFlowSolverSteadyState::TimeStep( const realT& time,
                                                        const realT& dt,
                                                        const int cycleNumber,
                                                        PhysicalDomainT& domain,
                                                        const sArray1d& namesOfSolverRegions ,
                                                        SpatialPartition& partition,
                                                        FractunatorBase* const fractunator )
{

  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();

  GenerateParallelPlateGeometricQuantities( domain, time,dt );

  CalculateSteadyStatePressureDistribution( domain, partition, time );

  UpdateEOS( dt, domain );
  UpdateFlux( time, dt, domain, partition);

  return dt;

}





void ParallelPlateFlowSolverSteadyState::GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,realT time,
        realT dt  )
{

  Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);
  
  // get apertures of flow path (loop over faces)
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf ) {

    domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *kf, faceCenter[*kf]);

    if(m_doApertureUpdate){
      R1Tensor gap;
      R1Tensor N;
      N = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *kf );

      gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, *kf );
      apertures[*kf] = Dot(gap,N) ;

      if( apertures[*kf]<m_min_aperture )
        apertures[*kf] = m_min_aperture;
      else if( apertures[*kf] > m_max_aperture )
        apertures[*kf] = m_max_aperture;
    }
     
  }
}



void ParallelPlateFlowSolverSteadyState::CalculateSteadyStatePressureDistribution( PhysicalDomainT& domain,
                                                                                       SpatialPartition& partition, realT time )
{
  
  SetupSystem (domain,partition, time);
  Assemble    (domain,partition, time);
  Solve       (domain,partition);
  
}



void ParallelPlateFlowSolverSteadyState::UpdateEOS( const realT dt , PhysicalDomainT& domain)
{


  const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  rArray1d& density = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& mass    = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& pressure    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();


  for( lSet::const_iterator fc=m_faceSet->begin() ; fc!=m_faceSet->end() ; ++fc ) {
    density[*fc] = rho_EOS(pressure[*fc],m_bulk_modulus,m_rho_o);
    realT area = domain.m_feFaceManager.SurfaceArea(domain.m_feNodeManager, *fc);
    mass[*fc] = density[*fc] *area*apertures[*fc];

    // propagate pressure to children
    lArray1d& childFaces = domain.m_feFaceManager.m_childIndices[*fc];
    for(unsigned i =0; i < childFaces.size(); ++i){
      pressure[childFaces[i]] = pressure[*fc];
    }
  }
}


void ParallelPlateFlowSolverSteadyState::UpdateFlux( const realT time, const realT dt , PhysicalDomainT& domain,
                                                     SpatialPartition& partition
){
	
  const rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  const rArray1d& face_aperture = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);
  volFlux = 0.0;
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<R1Tensor>& faceVolFlux = domain.m_feFaceManager.GetFieldData<R1Tensor>( VolumetricFluxStr );
  faceVolFlux = 0.0;
  
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  
  //iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  // update edge fluxes
  {
    std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
    for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
    {
      const int numFaces = itr->second.size();
      if( numFaces > 1) {
        localIndex eg = itr->first;
        localIndex kfa = itr->second[0];
        localIndex kfb = itr->second[1];

        realT Pa = pressures[kfa];
        realT Pb = pressures[kfb];

        // will be incorrect at junctions of 3 or more faces
        volFlux[eg] =  edgePermeabilities[eg]*(Pa-Pb); // flux A -> B

        R1Tensor edgeCenter; domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter );
        R1Tensor vecFluxA,vecFluxB;

        vecFluxA = edgeCenter - faceCenters[kfa];
        vecFluxB = faceCenters[kfb] - edgeCenter;

        vecFluxA.Normalize();
        vecFluxB.Normalize();
        vecFluxA *= volFlux[eg];
        vecFluxB *= volFlux[eg];

        //if(is_ghost[kfa] < 0) faceVolFlux[kfa] += vecFluxA;
        //if(is_ghost[kfb] < 0) faceVolFlux[kfb] += vecFluxB;


      }
    }
  }
  
  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::PressureBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time );
                                
  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, "NetFlux", time );


  ApplyMultiSetBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::SinglePartitionBC_VelocityUpdate,
                                           domain, domain.m_feEdgeManager,
                                           SinglePartitionPeriodicBoundaryCondition::BoundaryConditionName(), time );

   /*
  // synchronize ghost edge fluxes// this is problematic because edge-face orderings are not preserved for ghost edges
  partition.SynchronizeFields(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);

  // update face fluxes

  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    const int numFaces = itr->second.size();

    localIndex eg = itr->first;
    R1Tensor edgeCenter;
    domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter );
    for(int i = 0; i < numFaces; ++i){
      localIndex kf = itr->second[i];
      R1Tensor vecFlux;

      vecFlux = edgeCenter;
      vecFlux -= faceCenters[kf];

      vecFlux.Normalize();
      vecFlux *= volFlux[eg];  // this is problematic because edge-face orderings are not preserved for ghost edges
      if(i!=0)  vecFlux *= -1;

      faceVolFlux[kf] += vecFlux;
    }
  }
  faceVolFlux *= 0.5;
  */

  // synchronize ghost face fluxes
  partition.SynchronizeFields(syncedFieldsB, CommRegistry::steadyStateParallelPlateFlowSolverB);


  // convert flux volume per element  on edges to flux per cross-sectional area on faces
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
     {

 	  R1Tensor faceNorm = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *kf );

 	  //realT weight = 0.0;
 	  localIndex numEdges = domain.m_feFaceManager.m_toEdgesRelation[*kf].size();

 	  // new flux calculation
 	  R1TensorT<2> NQ,q;
 	  R2SymTensorT<2> NN;
 	  
 	  // choose plane for calculation (xi,yi - plane for calculating least squares minimization, zi - out of plane coordinate)
      int xi = 0;// fixme need to project into coordinates on plane of face - currently assumes vertical orientation
      int yi = 1;
      int zi = 2;

      if(abs(faceNorm[zi]) < abs(faceNorm[xi]) ){
    	  int temp = xi;
    	  xi = zi; xi = temp;
      }
      if(abs(faceNorm[zi]) < abs(faceNorm[yi]) ){
    	  int temp = yi;
    	  yi = zi; zi = temp;
      }

  	 //std::cout << " face calc " << std::endl << std::endl;

 	  for(size_t a =0; a < numEdges; ++a){
 	          localIndex eg = domain.m_feFaceManager.m_toEdgesRelation[*kf][a];

 	          R1Tensor edgeCenter;
 	          domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter );

 	          R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,dT); // edge vector
 	          realT l = dT.L2_Norm();

 	          R1Tensor dN;

              const bool useNormalizedNetFlux = true; // much better than other method
 	          if(useNormalizedNetFlux){
 	            // calculate flux using net flux divided by face area
 	 	        // Solve q_{j} n^{a}_{j} = Q^{a}/l

                realT Q = volFlux[eg]/l;


                localIndex numfaces = m_edgesToFaces[eg].size();
                if(numfaces == 2){  // fixme 2 faces only
                  localIndex fcb = m_edgesToFaces[eg][0];
                  if( fcb != *kf){
                	  Q = -Q; // neighbor correct, flux is not
                  } else {
                    fcb = m_edgesToFaces[eg][1]; // flux is correct but neighbor is not
                  }
          	      dN = faceCenters[fcb] - faceCenters[*kf];
                } else if(numfaces == 1){
                  dN = edgeCenter - faceCenters[*kf];
                  dN = dN- ((dN*dT)/(l*l))*dT; // project normal to face
                } else {
                	dN = edgeCenter - faceCenters[*kf];
                }
        	    dN.Normalize();

 	            NN(0,0) += dN[xi]*dN[xi]+TINY;
 	            NN(0,1) += dN[xi]*dN[yi];
 	            NN(1,1) += dN[yi]*dN[yi]+TINY; // TINY -> prevent singular inverse



                NQ(0) += Q*dN[xi];
                NQ(1) += Q*dN[yi];

               // DEBUG(dN)
               // DEBUG(NQ)
               // DEBUG(Q)

 	          } else {

 	            // calculate average flux weighting edge contributions by length
 	        	// Solve lq_{j} n^{a}_{j} = Q^{a}
 	            realT Q = volFlux[eg];
                localIndex numfaces = m_edgesToFaces[eg].size();
 	            if(numfaces == 2){  // fixme 2 faces only
 	              localIndex fcb = m_edgesToFaces[eg][0];
 	              if( fcb != *kf){
 	                Q = -Q; // neighbor correct, flux is not
 	              } else {
 	                fcb = m_edgesToFaces[eg][1]; // flux is correct but neighbor is not
 	              }
 	              dN = faceCenters[fcb] - faceCenters[*kf];
 	            } else if(numfaces == 1){
 	              // ultimately may be that best method is to enforce bc constraints directly - recalculate for pressure/outflow bcs
 	              dN = edgeCenter - faceCenters[*kf];
 	              dN = dN- ((dN*dT)/(l*l))*dT; // project normal to face
 	            } else {
 	              dN = edgeCenter - faceCenters[*kf];
 	            }
 	            dN.Normalize();
 	            dN*l;

 	          	NN(0,0) += dN[xi]*dN[xi]+TINY;
 	          	NN(0,1) += dN[xi]*dN[yi];
 	          	NN(1,1) += dN[yi]*dN[yi]+TINY; // prevent 0 inverse

                NQ(0) += Q*dN[xi];
 	            NQ(1) += Q*dN[yi];
 	          }
 	  }

 	  // invert NN
 	  NN.Inverse();

 	  // calculate flux
 	  q.AijBj(NN,NQ);

 	  // out of plane component
 	  realT qz = (-q[0]*faceNorm[xi] - q[1]*faceNorm[yi])/ faceNorm[zi];

 	  realT app = face_aperture[*kf];

 	  // Project back to xyz from face coordinates
 	  faceVolFlux[*kf][xi] = q[0]/app;
 	  faceVolFlux[*kf][yi] = q[1]/app;
 	  faceVolFlux[*kf][zi] = qz/app;

  }
    
}



/// Update the velocity fluxes on the boundary faces
void ParallelPlateFlowSolverSteadyState::
       PressureBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                 ObjectDataStructureBaseT& object ,
                                                 BoundaryConditionBase* bc, const lSet& set, realT time){
  
  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  
  Array1dT<R1Tensor>& faceVolFlux = domain.m_feFaceManager.GetFieldData<R1Tensor>(VolumetricFluxStr);
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
    
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  for( lSet::const_iterator eg=set.begin() ; eg!=set.end() ; ++eg ) {
     
    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if(itr != m_edgesToFaces.end() ){
    	
      lArray1d& faces = itr->second; 
      
      R1Tensor edgeCenter;
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *eg , edgeCenter );
      
      // nb only uses first face to calculate volume flux
      // flux will be incorrect for other edges if at junction.
      {            
        localIndex fc = faces[0];
                
        const realT Pa = pressures[fc];
        const realT Pb = bc->GetValue(domain.m_feEdgeManager,eg,time);
    	
        volFlux[*eg] =  edgePermeabilities[*eg]*(Pa-Pb); // flux out of a into b


        R1Tensor vecFlux;

        vecFlux = edgeCenter;
        vecFlux -= faceCenters[fc];

        vecFlux.Normalize();
        vecFlux *= volFlux[*eg];

        if(is_ghost[fc] < 0) faceVolFlux[fc] += vecFlux;


      }
    }	
  }

}



// this becomes more expensive as number of processors increases - may just want to do it once.
void ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition_Setup(PhysicalDomainT& domain,
                                                                                    int& n_local_rows, const realT& time){
                                                                          
  m_flux_bc_counter = 0;
  m_flux_bc_trilinos_host.resize(0);
  m_flux_bc_trilinos_index.resize(0);
  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition_Count,
                                domain, domain.m_feEdgeManager, "NetFlux",time );
  m_flux_bc_counter = 0; // reset after count
  
  localIndex numBCs = m_flux_bc_trilinos_host.size();
  m_flux_bc_pressure.resize(numBCs);
  m_flux_bc_trilinos_index.resize(numBCs);
  if(numBCs > 0){
    std::vector<int> gather(n_mpi_processes*numBCs);

    m_epetra_comm.GatherAll(&m_flux_bc_trilinos_host.front(),
                            &gather.front(),
                            numBCs);
    
    // determine which processor has the dummy dof
    for(localIndex i = 0; i < numBCs; ++i){
      int maxIndexJ = 0;
      localIndex maxVal = gather[i*n_mpi_processes];
      for( int j = 1; j < n_mpi_processes; ++j){
        localIndex val = gather[i*n_mpi_processes+j];
    	if(val > maxVal){
          maxIndexJ = j;
          maxVal = val;
    	}
      }   	
      m_flux_bc_trilinos_host[i] = maxIndexJ; // record winning processor
      if(maxIndexJ == this_mpi_process) n_local_rows++;
      
    }
    
  }                                                    
                                	
}

/// Count number of flux controlled boundary conditions
/// and store the global index of the first element in the set
void ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition_Count( PhysicalDomainT& domain ,
                                                                                     ObjectDataStructureBaseT& object,
                                                                                     BoundaryConditionBase* bc ,
                                                                                     const lSet& set,realT time  )
{
                                     	
  m_flux_bc_trilinos_host.resize(m_flux_bc_counter+1);
  
  // store global index of 1st element in set 
  // The host is determined later - given to the processor with the largest non-ghost element in the set
  iArray1d& is_ghost       = object.GetFieldData<FieldInfo::ghostRank>();
  m_flux_bc_trilinos_host[m_flux_bc_counter] = -1; // empty 
  
  if(set.size() > 0){
    lSet::const_iterator itr=set.begin(); 
    while(itr != set.end() && m_flux_bc_trilinos_host[m_flux_bc_counter] < 0 ){
      
      if(is_ghost[*itr] < 0){
      	 localIndex indx = object.m_localToGlobalMap[*itr];
      	 m_flux_bc_trilinos_host[m_flux_bc_counter] = indx; // store global index
      }
      ++itr;
    }
  }     
  
  m_flux_bc_counter++;                 	
}

void ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition_Sparsity(
                                     PhysicalDomainT& domain,
                                     ObjectDataStructureBaseT& object ,
                                     BoundaryConditionBase* bc ,
                                     const lSet& set,realT time ){
                                     	 
  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  //iArray1d& face_is_ghost       = domain.m_faceManager.GetFieldData<FieldInfo::ghostRank>();
  iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
 
 
  std::vector<int> dofIndex(2);
  dofIndex[0] = m_flux_bc_trilinos_index[m_flux_bc_counter]; 
  	
  // loop over edges
  lSet::const_iterator eg=set.begin() ;

  for( localIndex i =0; i < set.size() ; ++i, ++eg ){

    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){
      lArray1d& faces = itr->second; 
      
      for(size_t ii =0; ii < faces.size(); ++ii){
      	const localIndex kf = faces[ii];
        dofIndex[1] = trilinos_index[kf];
        
        sparsity->InsertGlobalIndices(2,
                                      &dofIndex.front(),
                                      2,
                                      &dofIndex.front());
      }
    }
  }
              
  m_flux_bc_counter++;      
  if( m_flux_bc_counter == m_flux_bc_trilinos_index.size() ) m_flux_bc_counter=0;                 	
}

void ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition( PhysicalDomainT& domain,
                                                                               ObjectDataStructureBaseT& object ,
                                                                               BoundaryConditionBase* bc,
                                                                               const lSet& set,
                                                                               realT time){
                                     	 
  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  //iArray1d& face_is_ghost       = domain.m_faceManager.GetFieldData<FieldInfo::ghostRank>();
  iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
 
  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
 
  	  	
//  int numBCs = set.size();
  
  Epetra_IntSerialDenseVector  bc_dof_face_dof(2);
  Epetra_SerialDenseMatrix     bc_face_matrix(2,2);
  
  Epetra_IntSerialDenseVector  bc_dof(1);
  Epetra_SerialDenseVector     bc_rhs(1);
  
  bc_dof(0) = m_flux_bc_trilinos_index[m_flux_bc_counter]; 
  bc_dof_face_dof(0) = m_flux_bc_trilinos_index[m_flux_bc_counter]; 
  	
  // loop over edges, find permeabilities and apply boundary conditions.

  lSet::const_iterator eg=set.begin() ;

  for( localIndex i =0; i < set.size() ; ++i, ++eg ){

    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){
      lArray1d& faces = itr->second; 
      
      R1Tensor edgeCenter;
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *eg , edgeCenter );
      for(size_t ii =0; ii < faces.size(); ++ii){
      	
      	
      	const localIndex kf = faces[ii];

        // calculate edge permeability for face
        R1Tensor la = edgeCenter; la -= faceCenters[kf];
        realT w = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, *eg);

        const realT la_norm = la.L2_Norm();
        const realT app = apertures[kf];
        const realT kappa = CalculatePermeability( la_norm , app, w,m_mu,m_SHP_FCT);
        
        if(ii==0) edgePermeabilities[*eg] = kappa;  // only assign first permeability (at junction)
        
        bc_dof_face_dof(1) = trilinos_index[kf]; 
        
  	    bc_face_matrix(0,0) = kappa;
  	    bc_face_matrix(1,1) = kappa;
  	    bc_face_matrix(0,1) = -kappa;
  	    bc_face_matrix(1,0) = -kappa;
  	    
  	    matrix->SumIntoGlobalValues(bc_dof_face_dof, bc_face_matrix);
    
      }
    }
  }
    
  if(m_flux_bc_trilinos_host(m_flux_bc_counter) == this_mpi_process){
    bc_rhs(0) =  bc->GetValue(domain.m_feEdgeManager,eg,time);
    rhs->SumIntoGlobalValues(bc_dof, bc_rhs);
  }       
     
  m_flux_bc_counter++;    
  if( m_flux_bc_counter == m_flux_bc_trilinos_index.size() ) m_flux_bc_counter=0;                     	
}
                                  
void ParallelPlateFlowSolverSteadyState::FluxControlledBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                                                              ObjectDataStructureBaseT& object ,
                                                                                              BoundaryConditionBase* bc ,
                                                                                              const lSet& set,
                                                                                              realT time ){

 rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  
  Array1dT<R1Tensor>& faceVolFlux = domain.m_feFaceManager.GetFieldData<R1Tensor>(VolumetricFluxStr);
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
    
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  for( lSet::const_iterator eg=set.begin() ; eg!=set.end() ; ++eg ) {
     
    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if(itr != m_edgesToFaces.end() ){
    	
      lArray1d& faces = itr->second; 
      
      R1Tensor edgeCenter;
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *eg , edgeCenter );
      
      // nb only uses first face to calculate volume flux
      // flux will be incorrect for other edges if at junction.
      {            
        localIndex fc = faces[0];
                
        const realT Pa = pressures[fc];
        const realT Pb = m_flux_bc_pressure[m_flux_bc_counter];
    	
        volFlux[*eg] =  edgePermeabilities[*eg]*(Pa-Pb); // flux out of a into b


        R1Tensor vecFlux;

        vecFlux = edgeCenter;
        vecFlux -= faceCenters[fc];

        vecFlux.Normalize();
        vecFlux *= volFlux[*eg];

        if(is_ghost[fc] < 0) faceVolFlux[fc] += vecFlux;


      }
    }	
  }
  
  m_flux_bc_counter++;    
  if( m_flux_bc_counter == m_flux_bc_trilinos_index.size() ) m_flux_bc_counter=0;         

}



/// Register solver in the solver factory
REGISTER_SOLVER( ParallelPlateFlowSolverSteadyState )
