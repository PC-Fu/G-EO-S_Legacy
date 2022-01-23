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
 * @file TwoDSteadyStateParallelPlateFlowSolver.cpp
 * @author walsh24
 * @date February 21, 2012
 */

#include "ParallelPlateFlowSolverFEM.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"


// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"









ParallelPlateFlowSolverFEM::ParallelPlateFlowSolverFEM(  const std::string& name,
                                                         ProblemManagerT* const pm ):
ParallelPlateFlowSolverBase(name,pm),
m_faceSet(),
m_phi(1.0),
this_mpi_process(pm->m_epetraComm.MyPID()),
n_mpi_processes(pm->m_epetraComm.NumProc()),
m_epetra_comm(pm->m_epetraComm),
row_map(),
sparsity(),
m_TrilinosIndexStr()
{
  m_TrilinosIndexStr = "ParallelPlateFlowSolverFEM_GlobalDof";
}

ParallelPlateFlowSolverFEM::~ParallelPlateFlowSolverFEM()
{
  // TODO Auto-generated destructor stub
}

void ParallelPlateFlowSolverFEM::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  ParallelPlateFlowSolverBase::ReadXML(hdn);

  // Mixed difference parameter
  m_phi = hdn->GetAttributeOrDefault<realT>("phi",1.0); // backward difference by default.
  
  // Linear Solver
  m_linearSolverParams.m_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-10);
  m_linearSolverParams.m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);

  // Flags
  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("UpdateAperture",false);
}


void ParallelPlateFlowSolverFEM::RegisterFields( PhysicalDomainT& domain )
{
  
  ParallelPlateFlowSolverBase::RegisterFields( domain.m_feFaceManager, domain.m_feEdgeManager );


    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::pressure>();    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::density>(); 
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::mass>(); 
  

  domain.m_feFaceManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);


  
}

void ParallelPlateFlowSolverFEM::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");

  flowFaceType = -1;
  flowEdgeType = -1;

  if( !(m_flowFaceSetName.empty()) )
  {
    m_faceSet = faceManager.GetSet(m_flowFaceSetName);
  }
  else
  {
    DefineFlowSets(domain);
  }
}


void ParallelPlateFlowSolverFEM:: SetupSystem (PhysicalDomainT&  domain,
                                                SpatialPartition& partition, const realT& time)
{
  
  const iArray1d& faceGhostRank  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d& nodeGhostRank  = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);
  
  //if(m_doApertureUpdate) UpdateAperture(domain);
  
  
  // count local dof
  ///////////////////////////////
  
  // local rows
  int n_local_rows = 0;
  
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if( faceGhostRank[*kf] < 0 )
    {
      for( lArray1d::const_iterator a=domain.m_feFaceManager.m_toNodesRelation[*kf].begin() ;
                                    a!=domain.m_feFaceManager.m_toNodesRelation[*kf].end() ; ++ a )
      {
        if( nodeGhostRank[*a] < 0 )
          ++n_local_rows;
      }
    } 
  }

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
  trilinos_index = -INT_MAX;
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(faceGhostRank[*kf] < 0)
    {
      for( lArray1d::const_iterator a=domain.m_feFaceManager.m_toNodesRelation[*kf].begin() ;
                                    a!=domain.m_feFaceManager.m_toNodesRelation[*kf].end() ; ++ a )
      {
        if( nodeGhostRank[*a] < 0 )
        {
          trilinos_index[*kf] = first_local_row+local_count;
          local_count++;
        }
      }
    }
  }

  assert(static_cast<int>(local_count) == n_local_rows);

  partition.SynchronizeFields(m_syncedFields, CommRegistry::parallelPlateFlowSolver);

  // create epetra map
  ////////////////////

  row_map = Teuchos::rcp(new Epetra_Map(n_global_rows,n_local_rows,0,m_epetra_comm));

  // set up sparsity graph
  ////////////////////////

  sparsity = Teuchos::rcp(new Epetra_FECrsGraph(Copy,*row_map,0));
  

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(faceGhostRank[*kf] < 0)
    {

      const lArray1d& nodelist = domain.m_feFaceManager.m_toNodesRelation[*kf];
      iArray1d dofIndex (nodelist.size());
      for( lArray1d::size_type a=0 ; a<=nodelist.size() ; ++ a )
      {
        localIndex b = nodelist[b];
        dofIndex[a] = trilinos_index[b] ;
      }
      sparsity->InsertGlobalIndices(dofIndex.size(),
                                    &dofIndex.front(),
                                    dofIndex.size(),
                                    &dofIndex.front());

    }
  }

  
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // (re-)init linear system
  matrix   = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*sparsity));
  solution = Teuchos::rcp(new Epetra_FEVector(*row_map));
  rhs      = Teuchos::rcp(new Epetra_FEVector(*row_map));

}

/* Assemble */

void ParallelPlateFlowSolverFEM :: Assemble (PhysicalDomainT&  domain,
                                                  SpatialPartition& partition ,
                                                  const realT& time,
                                                  const realT& dt)
{
  const iArray1d& faceGhostRank  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);


  matrix->Scale(0.0);
  solution->Scale(0.0);
  rhs->Scale(0.0);
/*
  realT Phi[4][4];
  {
  realT Eta[4] = { };
  for( int q=0 ; q<4 ; ++q )
  {
    for( int a=0 ; a<4 ; ++a )
    {
      Phi[q][a] =
    }
  }
  }
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(faceGhostRank[*kf] < 0)
    {
      const int fluid_dof_per_face = domain.m_feFaceManager.m_toNodesRelation[*kf].size();
      Epetra_IntSerialDenseVector  faceLocalDofIndex   (fluid_dof_per_face);
      Epetra_SerialDenseVector     face_rhs     (fluid_dof_per_face);
      Epetra_SerialDenseMatrix     face_matrix  (fluid_dof_per_face,fluid_dof_per_face);

      const lArray1d& nodelist = domain.m_feFaceManager.m_toNodesRelation[*kf];

      for( lArray1d::size_type a=0 ; a<=nodelist.size() ; ++ a )
      {
        faceLocalDofIndex[a] = trilinos_index[ nodelist[a] ];
      }

      for( int q=0 ; q<4 ; ++q )
      {
        for( lArray1d::size_type a=0 ; a<=nodelist.size() ; ++ a )
        {
          for( lArray1d::size_type b=0 ; b<=nodelist.size() ; ++ b )
          {
            face_rhs[a] += Phi[a]*
          }
        }
      }

    }
  }
*/

    
  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();


}



/* Solve */

void ParallelPlateFlowSolverFEM:: Solve (PhysicalDomainT&  domain,
                                               SpatialPartition& partition)
{

  // face fields
  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  rArray1d& massRate = domain.m_feFaceManager.GetFieldData<realT>("massRate");

  // set initial guess
  int dummy;
  double* local_solution = NULL;

  solution->ExtractView(&local_solution,&dummy);
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      int lid = row_map->LID(trilinos_index[*kf]);
      local_solution[lid] = massRate[*kf];
    }
  }
  


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
          //solver.SetAztecOption(AZ_conv,AZ_rhs);
          //solver.SetAztecOption(AZ_output,AZ_none);
          
          
  solver.Iterate(m_linearSolverParams.m_maxIters, m_linearSolverParams.m_tol);
  //solution->Print(std::cout);


  // copy solution to faces
  ////////////////////////
  
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      int lid = row_map->LID(trilinos_index[*kf]);
      massRate[*kf] = local_solution[lid];
    }
  }
  

  // re-sync ghost nodes
  partition.SynchronizeFields(m_syncedFields, CommRegistry::parallelPlateFlowSolver);
}


void ParallelPlateFlowSolverFEM::InitializeCommunications( PartitionBase& partition )
{
  m_syncedFields.clear();
  m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_TrilinosIndexStr);

  partition.SetBufferSizes(m_syncedFields, CommRegistry::parallelPlateFlowSolver);
}


void ParallelPlateFlowSolverFEM::TimeStep( const realT& time,
                                                        const realT& dt,
                                                        PhysicalDomainT& domain,
                                                        const sArray1d& namesOfSolverRegions ,
                                                        SpatialPartition& partition,
                                                        FractunatorBase* const fractunator )
{

}

void ParallelPlateFlowSolverFEM::DefineFlowSets( PhysicalDomainT& domain )
{
  if( m_flowFaceSetName.empty() )
  {
    FaceManagerT& faceManager = domain.m_feFaceManager;
    EdgeManagerT& edgeManager = domain.m_feEdgeManager;

    const iArray1d& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");

    m_faceSet.clear();

    for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
    {
      if( flowFaceType[kf] == 0 )
      {
        m_faceSet.insert( kf );
      }
    }
  }
}



/// Register solver in the solver factory
REGISTER_SOLVER( ParallelPlateFlowSolverFEM )
