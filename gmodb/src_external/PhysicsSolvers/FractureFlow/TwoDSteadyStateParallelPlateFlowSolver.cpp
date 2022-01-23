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
 * @file TwoDSteadyStateParallelPlateFlowSolver.cpp
 * @author walsh24
 * @date June 1, 2011
 */

#include "TwoDSteadyStateParallelPlateFlowSolver.h"
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

namespace{
	realT TINY = 1e-64;	
}


static realT CalculatePermeability( const realT la, const realT lb, const realT apa,
                                    const realT apb, const realT w, const realT mu, const realT SHP_FAC) ;
                                    
static realT CalculatePermeability( const realT la, const realT h, const realT w, const realT mu, const realT SHP_FAC) ;

SteadyStateParallelPlateFlowSolver_TwoD::SteadyStateParallelPlateFlowSolver_TwoD( const std::string& name,
                                                                                  ProblemManagerT* const pm ):
SolverBase(name,pm),
m_faceSet(NULL),
m_faceSetName(),
m_numFaces(0),
m_faceDofMap(),
m_edgesToFaces(),
m_tol(0.0),
m_maxSolverIters(0),
doDataWrite(true),
verboseFlag(false),
m_SHP_FCT(0.0),
m_mu(0.0),
m_bulk_modulus(0.0),
m_min_aperture(0.0),
m_max_aperture(0.0),
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
numerics()
{
  ++m_instances; 
  m_TrilinosIndexStr = "TwoDSSPPFS_" +  toString<int>(m_instances) + "_GlobalDof"; 
}

SteadyStateParallelPlateFlowSolver_TwoD::~SteadyStateParallelPlateFlowSolver_TwoD()
{
}

void SteadyStateParallelPlateFlowSolver_TwoD::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML( hdn );

  m_SHP_FCT = hdn->GetAttributeOrDefault<realT>("shapefactor",1.0);
  m_mu = hdn->GetAttributeOrDefault("mu","1.0e-3 N.s/m^2");
  m_bulk_modulus = hdn->GetAttributeOrDefault(BulkModulusStr,"2.0e9 Pa");
  m_min_aperture = hdn->GetAttributeOrDefault<realT>(MinimumApertureStr,0.000);
  m_max_aperture = hdn->GetAttributeOrDefault<realT>(MaximumApertureStr,1.000e64);
  
  m_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-10);
  m_maxSolverIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);
  doDataWrite = hdn->GetAttributeOrDefault<bool>("doDataWrite",false);
  verboseFlag = hdn->GetAttributeOrDefault<bool>(VerboseStr,false);
  
  updateApertureFlag = hdn->GetAttributeOrDefault<bool>("updateAperture",false);
  
  m_faceSetName = hdn->GetAttributeString("faceset");
  
  
  numerics.krylov_tol     = m_tol;
  
}


void SteadyStateParallelPlateFlowSolver_TwoD::RegisterFields( PhysicalDomainT& domain )
{
	
  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FluidVelocityStr,true,true);
    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::pressure>();    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::density>(); 
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::mass>(); 
  
  domain.m_feEdgeManager.AddKeylessDataField<realT>(PermeabilityStr,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(VolumetricFluxStr,true,true);
  
  
  domain.m_feEdgeManager.AddKeylessDataField<int>("FlowFaceCount",true,true);// debug
  

  domain.m_feFaceManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);
  
}

void SteadyStateParallelPlateFlowSolver_TwoD::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  EdgeManagerT& edgeManager = domain.m_feEdgeManager;
  
  m_faceSet = &(faceManager.GetSet(m_faceSetName));
  m_numFaces = m_faceSet->size();
  
  // build face-dof map

  lSet::const_iterator si=m_faceSet->begin();
  for(localIndex i =0; i < m_numFaces; ++i, ++si){
  	localIndex f = *si;
  	m_faceDofMap[f] = i;
  }
  
   iArray1d& ffCount = domain.m_feEdgeManager.GetFieldData<int>("FlowFaceCount"); // debug
  

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
      ffCount[eg] = edgeList.size();
    }
  }
    
}


/**
 * Transfer aperture data from the external face manager to the face manager. 
 * 
 * 
 */
void SteadyStateParallelPlateFlowSolver_TwoD::UpdateAperture(PhysicalDomainT& domain)
{

  const iArray1d& isExternal = domain.m_feFaceManager.m_isExternal;
  const lArray1d& externalFaceIndex = domain.m_feFaceManager.GetFieldData<localIndex>(
      "externalFaceIndex");
  const rArray1d& normal_approach = domain.m_externalFaces.GetFieldData<realT>("normalApproach");

  rArray1d& face_aperture = domain.m_feFaceManager.GetFieldData<realT>("Aperture");

  for (lSet::const_iterator si = m_faceSet->begin(); si != m_faceSet->end(); ++si)
  {
    localIndex kf = *si;
    if (isExternal[kf])
    {
      face_aperture[kf] = std::max(-normal_approach[externalFaceIndex[kf]], 0.0);
    }
  }
}

void SteadyStateParallelPlateFlowSolver_TwoD:: SetupSystem (PhysicalDomainT&  domain,
                                                SpatialPartition& partition, const realT& time)
{
	
  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  if(updateApertureFlag){
  	UpdateAperture(domain);
  }
  
  
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
  
  ApplyBoundaryCondition<realT>(this, &SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition_Sparsity, 
                                domain, domain.m_feEdgeManager, "FixedFlux", time );
     
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();
  
}

/* Assemble */

void SteadyStateParallelPlateFlowSolver_TwoD :: Assemble (PhysicalDomainT&  domain,
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
  ApplyBoundaryCondition<realT>(this, &SteadyStateParallelPlateFlowSolver_TwoD::PressureBoundaryCondition, 
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time );
                              
  ApplyBoundaryCondition<realT>(this, &SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition, 
                                domain, domain.m_feEdgeManager, "FixedFlux",time );
  	  

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  //rhs->Print(std::cout);
  //std::cout << matrix->NormInf() << std::endl;
  //EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*matrix);
  //exit(0);
}


/// Apply a pressure boundary condition to a given set of edges
void SteadyStateParallelPlateFlowSolver_TwoD::PressureBoundaryCondition(PhysicalDomainT& domain,
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
        R1Tensor la = edgeCenter; la -= faceCenters[kf];
        realT w = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, *eg);

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

void SteadyStateParallelPlateFlowSolver_TwoD:: Solve (PhysicalDomainT&  domain,
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

  // ml preconditioner

  //Teuchos::ParameterList amg_params;
  //int *options    = new int[AZ_OPTIONS_SIZE];
  //double *params  = new double[AZ_PARAMS_SIZE];
  //ML_Epetra::SetDefaults("SA",amg_params,options,params);

  AztecOO solver(problem);
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
          
          
  solver.Iterate(1000,numerics.krylov_tol);
  //solution->Print(std::cout);


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


void SteadyStateParallelPlateFlowSolver_TwoD::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  //syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(FaceCenterStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::pressure>::Name());
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(ApertureStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_TrilinosIndexStr);
  //syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(VolumetricFluxStr);

  partition.SetBufferSizes(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);

}


double SteadyStateParallelPlateFlowSolver_TwoD::TimeStep( const realT& time,
                                                        const realT& dt,
                                                        const int cycleNumber,
                                                        PhysicalDomainT& domain,
                                                        const sArray1d& namesOfSolverRegions ,
                                                        SpatialPartition& partition,
                                                        FractunatorBase* const fractunator )
{
  GenerateParallelPlateGeometricQuantities( domain );

  CalculateSteadyStatePressureDistribution( domain, partition, time );

  UpdateEOS( dt, domain );
  UpdateFlux( time, dt, domain);

  return dt;
}





void SteadyStateParallelPlateFlowSolver_TwoD::GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain )
{

//  rArray1d& apertures = domain.m_faceManager.GetFieldData<realT>( ApertureStr );
  Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  
  // get apertures of flow path (loop over faces)
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf ) {

    domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *kf, faceCenter[*kf]);
     
  }
}



void SteadyStateParallelPlateFlowSolver_TwoD::CalculateSteadyStatePressureDistribution( PhysicalDomainT& domain,
                                                                                       SpatialPartition& partition, realT time )
{
  
  SetupSystem (domain,partition, time);
  Assemble    (domain,partition, time);
  Solve       (domain,partition);
  
}



void SteadyStateParallelPlateFlowSolver_TwoD::UpdateEOS( const realT dt , PhysicalDomainT& domain)
{


  const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  rArray1d& density = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& mass    = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& pressure    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();


  for( lSet::const_iterator fc=m_faceSet->begin() ; fc!=m_faceSet->end() ; ++fc ) {
  	
    density[*fc] = 1.0 + pressure[*fc]/m_bulk_modulus; // pressure is +ve in compression
    realT area = domain.m_feFaceManager.SurfaceArea(domain.m_feNodeManager, *fc);
    mass[*fc] = density[*fc] *area*apertures[*fc];
  }
}


void SteadyStateParallelPlateFlowSolver_TwoD::UpdateFlux( const realT time, const realT dt , PhysicalDomainT& domain
){
	
  const rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  volFlux = 0.0;
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>( FluidVelocityStr );
  flowVelocity = 0.0;
  
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  // loop over edges
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
      
      
      R1Tensor vecFlux;

      vecFlux = faceCenters[kfb];
      vecFlux -= faceCenters[kfa];

      vecFlux.Normalize();
      vecFlux *= volFlux[eg];

      if(is_ghost[kfa] < 0) flowVelocity[kfa] += vecFlux;
      if(is_ghost[kfb] < 0) flowVelocity[kfb] += vecFlux;
    	
    }  
  }
  
  ApplyBoundaryCondition<realT>(this, &SteadyStateParallelPlateFlowSolver_TwoD::PressureBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time );
                                
  ApplyBoundaryCondition<realT>(this, &SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, "FixedFlux", time );
                                
  flowVelocity *= 0.5;
    
}



/// Update the velocity fluxes on the boundary faces
void SteadyStateParallelPlateFlowSolver_TwoD::
       PressureBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                 ObjectDataStructureBaseT& object ,
                                                 BoundaryConditionBase* bc, const lSet& set, realT time){
  
  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  
  Array1dT<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
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

        if(is_ghost[fc] < 0) flowVelocity[fc] += vecFlux;


      }
    }	
  }

}



// this becomes more expensive as number of processors increases - may just want to do it once.
void SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition_Setup(PhysicalDomainT& domain, 
                                                                                    int& n_local_rows, const realT& time){
                                                                          
  m_flux_bc_counter = 0;
  m_flux_bc_trilinos_host.resize(0);
  m_flux_bc_trilinos_index.resize(0);
  ApplyBoundaryCondition<realT>(this, &SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition_Count,
                                domain, domain.m_feEdgeManager, "FixedFlux",time );                                        
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
void SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition_Count( PhysicalDomainT& domain ,
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

void SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition_Sparsity(
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

void SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition( PhysicalDomainT& domain,
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
                                  
void SteadyStateParallelPlateFlowSolver_TwoD::FluxControlledBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                                                              ObjectDataStructureBaseT& object ,
                                                                                              BoundaryConditionBase* bc ,
                                                                                              const lSet& set,
                                                                                              realT time ){

 rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  
  Array1dT<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
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

        if(is_ghost[fc] < 0) flowVelocity[fc] += vecFlux;


      }
    }	
  }
  
  m_flux_bc_counter++;    
  if( m_flux_bc_counter == m_flux_bc_trilinos_index.size() ) m_flux_bc_counter=0;         

}

/**
 *
 * @param la length of the path between centers of two flow elements projected on normal attached to element a
 * @param lb length of the path between centers of two flow elements projected on normal attached to element b
 * @param apa aperture of the first flow element
 * @param apb aperture of the second flow element
 * @param w width of the flow path
 *
 * The "Permeability" returned is scaled by the cell face area
 * and divided by the product of the distance between the cells and kinematic viscosity. 
 * 
 *   k' = k*(h*w)/(mu*L) 
 * 
 * so that the mass flux from cell a to cell b is given by
 * 
 *   q_Mass = k' * rho * (Pa-Pb)
 *
 */
static realT CalculatePermeability(const realT la,
                                   const realT lb,
                                   const realT apa,
                                   const realT apb,
                                   const realT w,
                                   const realT mu,
                                   const realT SHP_FCT)
{
  if (apa <= 0 || apb <= 0)
    return 0.0;

  // realT h = 0.5*(apa+apb);
  // realT permeability = h*h*h *w / ( 12.0 * mu * ( la + lb ) );

  const realT ka = apa*apa*apa;
  const realT kb = apb*apb*apb;
  
  realT permeability = ka * kb * w/(12.0 * mu *(ka*lb + kb*la)); // harmonic mean

  permeability *= SHP_FCT;
  return permeability;
}

/**
 * One sided permeability calculation
 * 
 * @author walsh24
 * 
 * @param la length of the path between centers of two flow elements projected on normal attached to element a\
 * @param h aperture of the flow element
 * @param w width of the flow path
 * 
 * NB formulation differs from that in Johnson Morris 2009 
 * which uses distance between adjacent elements (2*la) rather than element and boundary (la). 
 * 
 * The "Permeability" returned is scaled by the cell face area
 * and divided by the product of the distance between the cells and kinematic viscosity. 
 * 
 *   k' = k*(h*w)/(mu*L) 
 * 
 * so that the mass flux from cell a to cell b is given by
 * 
 *   q_Mass = k' * rho * (Pa-Pb)
 *
 */
static realT CalculatePermeability(const realT la,
                                   const realT h,
                                   const realT w,
                                   const realT mu,
                                   const realT SHP_FCT)
{
  if (h <= 0) return 0.0;

  realT permeability = h*h*h *w / ( 12.0 * mu * la );

  permeability *= SHP_FCT;
  return permeability;
}

/**
 * @author walsh24
 * 
 * @brief Method to calculate the fluxes at joints with more than two flow elements.
 * 
 * We introduce a joint pressure such that flux into the joint is balanced by flux out:
 * 
 * \sum_i K_i(Pi-Po) = 0;
 * 
 * Then calculate the individual fluxes from
 * Q_i = K(Pi-Po);
 * 
 */
/*
static void CalculateEdgeFluxes(const rArray1d& ls,
                                const rArray1d& hs,
                                const rArray1d& ws,
                                const rArray1d& Ps,
                                const realT mu,
                                const realT SHP_FCT,
                                rArray1d& qs)
{
	// get one-sided permeabilities for each flow element
	realT TINY = 1e-64;
	const int numElems = ls.size();
	rArray1d ks(numElems);
	qs.resize(numElems);
	realT ksum = TINY;
	realT kPsum = 0.0;
	for(int i =0; i < numElems; ++i){
		ks[i] = CalculatePermeability(ls[i],hs[i],ws[i],mu,SHP_FCT);
		ksum +=ks[i];
		kPsum += ks[i] + Ps[i];	
	}
	
	// calculate joint pressure
	realT Po = kPsum/ksum;
	
	// find fluxes out of i into joint
	for(int i =0; i < numElems; ++i) qs[i] = ks[i]*(Ps[i]-Po);
	
}
*/

/// Register solver in the solver factory
REGISTER_SOLVER( SteadyStateParallelPlateFlowSolver_TwoD )
