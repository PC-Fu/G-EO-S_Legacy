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

#include "ParallelPlateFlowSolverFV.h"
#include "../SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"


// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"



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


using namespace BoundaryConditionFunctions;
using namespace PS_STR;
using namespace PPFS;


int ParallelPlateFlowSolverFV::m_instances = 0;






realT  clear_row ( Epetra_FECrsMatrix* const matrix,
                   const unsigned int row,
                   const realT factor );

/**
 * @author settgast
 * @param nodeIndex the local node index relative to the face
 * @param ndofPerNode number of degrees of freedom associated with a node
 *
 * This function gives the dof index for a node-pair associated with an original
 * parent node on a parent face. This function is used for consistency when
 * forming and applying contributions to the system.
 */
static void faceNodePairIndexing( const int nodeIndex,
                                    const int ndofPerNode,
                                    int index[2])
{
  index[0] = ndofPerNode*(nodeIndex*2);
  index[1] = ndofPerNode*(nodeIndex*2+1);
}

ParallelPlateFlowSolverFV::ParallelPlateFlowSolverFV(  const std::string& name,
                                                       ProblemManagerT* const pm ):
ParallelPlateFlowSolverBase(name,pm),
m_faceSet(),
m_numFaces(0),
m_faceDofMap(),
m_phi(1.0),
this_mpi_process(pm->m_epetraComm.MyPID()),
n_mpi_processes(pm->m_epetraComm.NumProc()),
syncedFields()
{
  ++m_instances; 
  m_TrilinosIndexStr = "TwoDIMPPFS_" +  toString<int>(m_instances) + "_GlobalDof";
}

ParallelPlateFlowSolverFV::~ParallelPlateFlowSolverFV()
{
  // TODO Auto-generated destructor stub
}

void ParallelPlateFlowSolverFV::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  ParallelPlateFlowSolverBase::ReadXML(hdn);

  // Mixed difference parameter
  m_phi = hdn->GetAttributeOrDefault<realT>("phi",1.0); // backward difference by default.
  
  // Linear Solver
  m_numerics.krylov_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-10);
  m_numerics.m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);

  // Flags
  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("UpdateAperture",false);
}


void ParallelPlateFlowSolverFV::RegisterFields( PhysicalDomainT& domain )
{
  
  ParallelPlateFlowSolverBase::RegisterFields( domain.m_feFaceManager, domain.m_feEdgeManager );

  const bool plotOldValues = false; // no need to plot old values unless debugging - overwritten at end of timestep.

  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr,true,true);

  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FluidVelocityStr,true,true);
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::pressure>();    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::density>(); 
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::mass>(); 
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::volume>();
  domain.m_feFaceManager.AddKeylessDataField<realT>("BoundedAperture",true,true);

  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("faceNormal0",true,true);


  domain.m_feFaceManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<int>("FlowFaceCount",true,true);// debug
  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>(EdgeCenterStr,true,plotOldValues);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(EdgeLengthStr,true,true);

}

void ParallelPlateFlowSolverFV::RegisterTemporaryFields( PhysicalDomainT& domain )
{
  domain.m_feFaceManager.AddKeylessDataField<realT>("FluidMass_n",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("FluidVolume_n",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("FluidDensity_n",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("Pressure_n",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("dPdM",false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr+std::string("_n"),false,false);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr+std::string("_n"),false,false);
  domain.m_feFaceManager.AddKeylessDataField<realT>("BoundedAperture_n",false,false);

  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>(EdgeCenterStr+std::string("_n"),false,false);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(EdgeLengthStr+std::string("_n"),false,false);


}


void ParallelPlateFlowSolverFV::DeregisterTemporaryFields( PhysicalDomainT& domain )
{

  domain.m_feFaceManager.RemoveDataField<realT>("FluidMass_n");
  domain.m_feFaceManager.RemoveDataField<realT>("FluidVolume_n");
  domain.m_feFaceManager.RemoveDataField<realT>("FluidDensity_n");
  domain.m_feFaceManager.RemoveDataField<realT>("Pressure_n");
  domain.m_feFaceManager.RemoveDataField<realT>("dPdM");
  domain.m_feFaceManager.RemoveDataField<realT>(ApertureStr+std::string("_n"));
  domain.m_feFaceManager.RemoveDataField<R1Tensor>(FaceCenterStr+std::string("_n"));
  domain.m_feFaceManager.RemoveDataField<R1Tensor>("BoundedAperture_n");

  domain.m_feEdgeManager.RemoveDataField<R1Tensor>("EdgeCenter_n");
  domain.m_feEdgeManager.RemoveDataField<realT>(EdgeLengthStr+std::string("_n"));
}


void ParallelPlateFlowSolverFV::FillTemporaryFields( PhysicalDomainT& domain )
{
  const rArray1d& mass_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& mass_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");
  mass_n = mass_np1;

  const rArray1d& volume_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& volume_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidVolume_n");
  volume_n = volume_np1;

  const rArray1d& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& density_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");
  density_n = density_np1;

  const rArray1d& pressure_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& pressure_n   = domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");
  pressure_n = pressure_np1;

  const rArray1d& apertures_np1 = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  rArray1d& apertures_n   = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
  apertures_n = apertures_np1;

  const rArray1d& boundedApertures_np1 = domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture" );
  rArray1d& boundedApertures_n   = domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture_n" );
  boundedApertures_n = boundedApertures_np1;

  const Array1dT<R1Tensor>& faceCenters_np1 = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  Array1dT<R1Tensor>& faceCenters_n   = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+std::string("_n") );
  faceCenters_n = faceCenters_np1;

  const Array1dT<R1Tensor>& edgeCenters_np1 = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  Array1dT<R1Tensor>& edgeCenters_n   = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr+std::string("_n"));
  edgeCenters_n = edgeCenters_np1;

  const rArray1d& edgeLengths_np1 = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  rArray1d& edgeLengths_n   = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr+std::string("_n") );
  edgeLengths_n = edgeLengths_np1;

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if( isZero(apertures_n[*kf]) )
    {
//      apertures_n[*kf] = this->m_min_aperture;
    }

    if( isZero(volume_n[*kf]) )
      volume_n[*kf] = (m_zeroApertureVolume + boundedApertures_np1[*kf]) * domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, *kf );

    if( isZero(mass_n[*kf]) )
      mass_n[*kf] = volume_n[*kf] * this->m_rho_o;

    if( isZero(density_n[*kf]) )
      density_n[*kf] = m_rho_o;

  }


}


void ParallelPlateFlowSolverFV::OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain )
{
  rArray1d& mass_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const rArray1d& mass_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");
  mass_np1 = mass_n;

  rArray1d& volume_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  const rArray1d& volume_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidVolume_n");
  volume_np1 = volume_n;

  rArray1d& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const rArray1d& density_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");
  density_np1 = density_n;

  rArray1d& pressure_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& pressure_n   = domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");
  pressure_np1 = pressure_n;

  rArray1d& apertures_np1 = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const rArray1d& apertures_n   = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
  apertures_np1 = apertures_n;

  rArray1d& boundedApertures_np1 = domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture" );
  const rArray1d& boundedApertures_n   = domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture_n" );
  boundedApertures_np1 = boundedApertures_n;

  Array1dT<R1Tensor>& faceCenters_np1 = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const Array1dT<R1Tensor>& faceCenters_n   = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+std::string("_n") );
  faceCenters_np1 = faceCenters_n;

  Array1dT<R1Tensor>& edgeCenters_np1 = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  const Array1dT<R1Tensor>& edgeCenters_n   = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr+std::string("_n"));
  edgeCenters_np1 = edgeCenters_n;

  rArray1d& edgeLengths_np1 = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  const rArray1d& edgeLengths_n   = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr+std::string("_n") );
  edgeLengths_np1 = edgeLengths_n;

}


void ParallelPlateFlowSolverFV::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  EdgeManagerT& edgeManager = domain.m_feEdgeManager;
  
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");

  flowFaceType = -1;
  flowEdgeType = -1;

  if( !(m_flowFaceSetName.empty()) )
  {
    m_faceSet = faceManager.GetSet(m_flowFaceSetName);
    m_numFaces = m_faceSet.size();

    // build face-dof map
  
    lSet::const_iterator si=m_faceSet.begin();
    for(localIndex i =0; i < m_numFaces; ++i, ++si){
      localIndex f = *si;
      m_faceDofMap[f] = i;
    }

    iArray1d& ffCount = domain.m_feEdgeManager.GetFieldData<int>("FlowFaceCount"); // debug

  
    for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
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
  else
  {
    DefineFlowSets(domain);
  }

  GenerateParallelPlateGeometricQuantities( domain,0,0 );
//  FillTemporaryFields( domain );

  InitializeDensity( domain);
}




void ParallelPlateFlowSolverFV::SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                                              SpatialPartition& partition,
                                                              int& numLocalRows,
                                                              int& numGlobalRows )
{

  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  //if(m_doApertureUpdate) UpdateAperture(domain);
  
  
  // count local dof
  ///////////////////////////////
  
  // local rows
  numLocalRows = 0;
  
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if( is_ghost[*kf] < 0 )
    {
      ++numLocalRows;
    } 
  }

  // determine the global/local degree of freedom distribution.
  ////////////////////////////////////////////////////////////

  std::vector<int> gather(n_mpi_processes);
  std::vector<int> cum_global_rows(n_mpi_processes);

  epetra_comm->GatherAll(&numLocalRows,
                        &gather.front(),
                        1);

  int first_local_row = 0;
  numGlobalRows = 0;

  for( int p=0; p<n_mpi_processes; ++p)
  {
    numGlobalRows += gather[p];
    if(p<this_mpi_process)
      first_local_row += gather[p];
    cum_global_rows[p] = numGlobalRows;
  }
  
  // create trilinos dof indexing
  //////////////////////////////////
  unsigned local_count = 0;
  // faces
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
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

  assert(static_cast<int>(local_count) == numLocalRows);

  partition.SynchronizeFields(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);

}

void ParallelPlateFlowSolverFV:: SetupSystem (PhysicalDomainT&  domain,
                                              SpatialPartition& partition,
                                              const realT& time)
{
  int n_local_rows;
  int n_global_rows;
  SetNumRowsAndTrilinosIndices( domain, partition,
                                n_local_rows,
                                n_global_rows );

  const iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);


  #if USECPP11==1
  // create epetra map
  ////////////////////

  m_rowMap = std::make_shared<Epetra_Map>(n_global_rows,n_local_rows,0,*epetra_comm);

  // set up sparsity graph
  ////////////////////////

  m_sparsity = std::make_shared<Epetra_FECrsGraph>(Copy,*m_rowMap,0);
#else
  m_rowMap = new Epetra_Map(n_global_rows,n_local_rows,0,*epetra_comm);
  m_sparsity = new Epetra_FECrsGraph(Copy,*m_rowMap,0);
  
#endif
  iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

  // loop over edges
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {    
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0 )
    {
      unsigned int numFaces = itr->second.size();
      if( numFaces > 1)
      {
        std::vector<int> dofIndex (numFaces);
        for(unsigned i=0; i<numFaces; ++i)
        {
          localIndex kf = itr->second[i];
          dofIndex[i] = trilinos_index[kf];
        }

        m_sparsity->InsertGlobalIndices(dofIndex.size(),
                                      &dofIndex.front(),
                                      dofIndex.size(),
                                      &dofIndex.front());
      }
    }
  }
  
  m_sparsity->GlobalAssemble();
  m_sparsity->OptimizeStorage();
  
}

/* Assemble */

realT ParallelPlateFlowSolverFV :: Assemble ( PhysicalDomainT&  domain,
                                              Epetra_System& epetraSystem,
                                              const realT& time,
                                              const realT& dt )
{
  realT maxMassScale = 0.0;
  // (re-)init linear system

  // basic face data ( = dof data for our problem)

  iArray1d& trilinosIndexFace = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

  iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  iArray1d& face_is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  const Array1dT<R1Tensor>& faceCenters_np1 = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const Array1dT<R1Tensor>& faceCenters_n   = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+std::string("_n") );


  const Array1dT<R1Tensor>& edgeCenters_np1 = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  const Array1dT<R1Tensor>& edgeCenters_n   = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr+std::string("_n"));

  const rArray1d& edgeLengths_np1 = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  const rArray1d& edgeLengths_n   = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr+std::string("_n") );

  const rArray1d& apertures_np1 = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const rArray1d& apertures_n   = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
  
  const rArray1d& mass_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const rArray1d& mass_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");

  const rArray1d& volume_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
//  const rArray1d& volume_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidVolume_n");

  const rArray1d& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const rArray1d& density_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");

  const rArray1d& pressure_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& pressure_n   = domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");

  const rArray1d& dPdM = domain.m_feFaceManager.GetFieldData<realT>("dPdM");

  rArray1d& edgePermeability = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);

  const int dim = 3;

  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");
//  const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );


//  Epetra_FECrsMatrix& dRdDisp = epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement );


  // Calculate the self(i.e. diagonal) contributions of the mass Residual, and dRdM. Also find the maximum mass in the
  // problem for scaling of the residual norm.
  {
    localIndex numFaces = domain.m_feFaceManager.DataLengths();
    Epetra_IntSerialDenseVector  faceDofIndex (numFaces);
    Epetra_SerialDenseVector     flowRHS     (numFaces);
    Epetra_SerialDenseMatrix     flowMatrix  (numFaces,numFaces);

    for( localIndex r=0 ; r<numFaces ; ++r )
    {
      faceDofIndex[r] = trilinosIndexFace[r];
      if( flowFaceType[r] == 0 && face_is_ghost[r] < 0 )
      {
        // self contribution to the residual
        flowRHS[r] -= -( mass_np1[r] - mass_n[r] );

        // derivative of the residual wrt mass_np1[r]
        flowMatrix(r,r) -= 1.0;

        // find maximum fluid mass in problem
        maxMassScale = std::max( maxMassScale, fabs( mass_np1[r] ) );
        maxMassScale = std::max( maxMassScale, fabs( mass_n[r] ) );
      }
    }

    m_matrix->SumIntoGlobalValues(faceDofIndex, flowMatrix);
    m_rhs->SumIntoGlobalValues(faceDofIndex,flowRHS);

  }

    CalculateApertureDerivatives( domain.m_feFaceManager,
                                  domain.m_feNodeManager );


  // loop over edges and process all flow through edges between faces attached to this edge
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0  && flowEdgeType[eg] == 0 )
    {

      // number of faces attached to edge.
      const lArray1d& facelist = itr->second;
      const int numFaces = facelist.size();

      // allocate space for dof map, residual, and dRdM.
      Epetra_IntSerialDenseVector  faceDofIndex (numFaces);
      Epetra_SerialDenseVector     flowRHS     (numFaces);
      Epetra_SerialDenseMatrix     flowMatrix  (numFaces,numFaces);
     
      flowRHS.Scale(0.0);
      flowMatrix.Scale(0.0);

      Epetra_IntSerialDenseVector  displacementDofIndex;
      Epetra_SerialDenseMatrix     matrix_10;


      rArray1d dVr_du;
      rArray1d dVs_du;

      rArray1d dKappa_du_rterm;
      rArray1d dKappa_du_sterm;

      for( unsigned int kr=0 ; kr<facelist.size() ; ++kr )
      {
        const unsigned int r = facelist[kr];
        if( flowFaceType[r] == 0 )
        {

          faceDofIndex[kr] = trilinosIndexFace[r];

          realT area_r = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, r, false );

          for( unsigned int ks=kr+1 ; ks<facelist.size() ; ++ks )
          {
            const unsigned int s = facelist[ks];
            if( flowFaceType[s] == 0 )
            {


              realT kappa_n  = TwoFacePermeability( edgeCenters_n, edgeLengths_n,
                                                          faceCenters_n, apertures_n,
                                                          eg, r, s,
                                                          NULL, NULL, NULL );

              realT kappa_np1  = TwoFacePermeability( edgeCenters_np1, edgeLengths_np1,
                                                      faceCenters_np1, apertures_np1,
                                                      eg, r, s, &m_dwdu, &dKappa_du_rterm, &dKappa_du_sterm );

//              kappa_n = 0;
//              kappa_np1 = kappa_n;
//              dKappa_du_rterm = 0.0;
//              dKappa_du_sterm = 0.0;

              edgePermeability[eg] = kappa_np1;






              // set variables for upwinding

              const localIndex ku_n   =   pressure_n[r] > pressure_n[s]   ? kr : ks ;
              const localIndex ku_np1 = pressure_np1[r] > pressure_np1[s] ? kr : ks ;

              const localIndex u_n   = facelist[ku_n];
              const localIndex u_np1 = facelist[ku_np1];

              const realT rho_n   = density_n[u_n] ;
              const realT rho_np1 = density_np1[u_np1] ;

              const realT drho_dm = 1.0 / volume_np1[u_np1] ;

              realT area_s = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, s, false );


              const realT dp_n   = pressure_n[s] - pressure_n[r];
              const realT dp_np1 = pressure_np1[s] - pressure_np1[r];


              // while adding contributions, we are doing this for a face pair (s,r) on the residual for r.s
              // add the contributions to the residual
              realT R_sr = ( ( 1.0 - m_phi ) * kappa_n * rho_n * dp_n + m_phi * kappa_np1 * rho_np1 * dp_np1 ) * dt;
              flowRHS[kr] -=  R_sr;
              flowRHS[ks] -= -R_sr;


              const realT dRdm_drho = ( m_phi * kappa_np1 * ( dp_np1 * drho_dm ) ) * dt;
              flowMatrix(kr,ku_np1) += dRdm_drho;

              // the s-face term is just the negative of the r-face term
              flowMatrix(ks,ku_np1) -= dRdm_drho;


              // when q=r
              const realT dRdm_dpr = m_phi * kappa_np1 * ( - rho_np1 * dPdM[r] ) * dt;
              // when q=s
              const realT dRdm_dps = m_phi * kappa_np1 * (   rho_np1 * dPdM[s] ) * dt;

              flowMatrix(kr,kr) +=  dRdm_dpr;
              flowMatrix(kr,ks) +=  dRdm_dps;

              // The s-face terms are the negative of the r-face term, and they
              // are also flipped.
              flowMatrix(ks,ks) += -dRdm_dps;
              flowMatrix(ks,kr) += -dRdm_dpr;





              if( epetraSystem.GetBlockID(EpetraBlock::displacement) != -1 )
              {
                // form coupled matrix for flow DOF's wrt displacement DOF's
                if( epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement ) )
                {
                  // calculate displacement DOF for these two faces
                  const localIndex rnodes = domain.m_feFaceManager.m_toNodesRelation[r].size();
                  const localIndex snodes = domain.m_feFaceManager.m_toNodesRelation[s].size();
                  rArray1d drhoRdu = m_dwdu(r);
                  drhoRdu *= -density_np1[r]/volume_np1[r]*area_r * m_dwdw(r);
                  rArray1d drhoSdu = m_dwdu(s);
                  drhoSdu *= -density_np1[s]/volume_np1[s]*area_s * m_dwdw(s);


                  const realT dPdrho_r = dPdM[r] * volume_np1[r] ;
                  const realT dPdrho_s = dPdM[s] * volume_np1[s] ;

                  Epetra_IntSerialDenseVector  dispSubFaceDofIndex(2);
                  dispSubFaceDofIndex(0) = trilinosIndexFace[r];
                  dispSubFaceDofIndex(1) = trilinosIndexFace[s];

                  displacementDofIndex.Resize( dim*2*(rnodes + snodes ) );
                  matrix_10.Reshape( 2, dim*2*(rnodes + snodes ) );
                  matrix_10.Scale(0.0);



                  // loop over the nodes on the r-face
                  for( localIndex a=0 ; a<rnodes ; ++a )
                  {
                    int rnodeDofIndex[2] ;
                    faceNodePairIndexing( a, dim, rnodeDofIndex);
                    for( int i=0 ; i<dim ; ++i )
                    {
                      for( int side=0 ; side<2 ; ++side )
                      {
                        const int ai=rnodeDofIndex[side]+i;
                        displacementDofIndex(ai) = m_dwdu_dof(r)(ai);

                        // r-face residual derivative
                        realT dRrdu = m_phi * ( rho_np1*dp_np1*dKappa_du_rterm(ai)
                                              + kappa_np1*rho_np1*(-dPdrho_r*drhoRdu(ai) ) ) * dt;

                        // upwind density derivative only applied if r is the upwind face.
                        if( ku_np1==kr )
                        {
                          dRrdu += m_phi * ( kappa_np1*drhoRdu(ai) * dp_np1 ) * dt;
                        }

                        // s-face residual derivative
                        realT dRsdu = - dRrdu;

                        matrix_10(0,ai) += dRrdu;
                        matrix_10(1,ai) += dRsdu;

                      }
                    }
                  }

                  // loop over the nodes on the s-face
                  for( localIndex a=0 ; a<snodes ; ++a )
                  {
                    int snodeDofIndex[2] ;
                    faceNodePairIndexing( a, dim, snodeDofIndex);
                    for( int i=0 ; i<dim ; ++i )
                    {
                      for( int side=0 ; side<2 ; ++side )
                      {
                        const int ai=snodeDofIndex[side]+i;
                        const int aiOffset = ai + rnodes*dim*2;
                        displacementDofIndex(aiOffset) = m_dwdu_dof(s)(ai);

                        // r-face residual derivative
                        realT dRrdu = m_phi * ( rho_np1*dp_np1*dKappa_du_sterm(ai)
                                              + kappa_np1*rho_np1*(dPdrho_s*drhoSdu(ai) ) ) * dt;

                        // upwind density derivative only applied if r is the upwind face.
                        if( ku_np1==ks )
                        {
                          dRrdu += m_phi * ( kappa_np1*drhoSdu(ai) * dp_np1 ) * dt;
                        }

                        // s-face residual derivative
                        realT dRsdu = - dRrdu;

                        matrix_10(0,aiOffset) += dRrdu;
                        matrix_10(1,aiOffset) += dRsdu;

                      }
                    }
                  }



                  epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement )->SumIntoGlobalValues(dispSubFaceDofIndex, displacementDofIndex, matrix_10);
                }
              }

            }
          }
        }
      }

      
      m_matrix->SumIntoGlobalValues(faceDofIndex, flowMatrix);
      m_rhs->SumIntoGlobalValues(faceDofIndex,flowRHS);

    } // ghost edge
  } // edge loop
  
  
  // boundary conditions
//  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverFV::PressureBoundaryCondition,
//                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()),
//                                time, dt, flowDofOffset );

  ApplyBoundaryCondition<R1Tensor>(this, &ParallelPlateFlowSolverFV::MassRateBC,
                                   domain, domain.m_feFaceManager,
                                   "MassRate", time + 0.5*dt, dt );

  ApplyBoundaryCondition<R1Tensor>(this, &ParallelPlateFlowSolverFV::MassBC,
                                   domain, domain.m_feFaceManager,
                                   Field<FieldInfo::mass>::Name(), time + dt );

  m_matrix->GlobalAssemble(true);
  m_rhs->GlobalAssemble();


  if( epetraSystem.GetBlockID(EpetraBlock::displacement) != -1 )
  if( epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement ) )
  {
//    epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement )->Print(std::cout);

#if USECPP11==1
    const std::shared_ptr<Epetra_Map> dispRowMap = epetraSystem.GetRowMap(EpetraBlock::displacement);
    const std::shared_ptr<Epetra_Map> flowRowMap = epetraSystem.GetRowMap(EpetraBlock::fluidMass);
#else
    const Epetra_Map* const dispRowMap = epetraSystem.GetRowMap(EpetraBlock::displacement);
    const Epetra_Map* const flowRowMap = epetraSystem.GetRowMap(EpetraBlock::fluidMass);
#endif
//    epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement )->Print(std::cout);
    epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement )->GlobalAssemble(*dispRowMap,*flowRowMap);
//    epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement )->Print(std::cout);
  }


//  m_rhs->Print(std::cout);
//  m_matrix->Print( std::cout);
//  std::cout << m_matrix->NormInf() << std::endl;
//  EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*(m_matrix.get()));
  //exit(0);

  {
//    std::shared_ptr<Epetra_FECrsMatrix> scratch = epetraSystem.GetScratch( EpetraBlock::fluidMass, EpetraBlock::fluidMass );
//    scratch = std::make_shared<Epetra_FECrsMatrix>( epetraSystem.GetMatrix( EpetraBlock::fluidMass, EpetraBlock::fluidMass ) );

  }



  return maxMassScale;
}


void ParallelPlateFlowSolverFV::MassBC( PhysicalDomainT& domain,
                                        ObjectDataStructureBaseT& object,
                                        BoundaryConditionBase* bc,
                                        const lSet& set,
                                        realT time)
{

  realT LARGE ;

  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& faceGhostRank  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const rArray1d& mass_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();


  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);

  for( lSet::const_iterator faceID=set.begin() ; faceID!=set.end() ; ++faceID )
  {
    face_dof(0) = trilinos_index[*faceID];

    if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]==0 )
    {
#if USECPP11==1
      LARGE = clear_row( m_matrix.get(), face_dof(0), 1.0 );
      clear_row( this->m_system->GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement ).get(), face_dof(0),1.0);
#else
      LARGE = clear_row( m_matrix, face_dof(0), 1.0 );
      clear_row( this->m_system->GetMatrix( EpetraBlock::fluidMass, EpetraBlock::displacement ), face_dof(0),1.0);
#endif

      face_rhs(0) = LARGE*( bc->GetValue(domain.m_feFaceManager,faceID,time) - mass_np1[*faceID] );
    }
    else
    {
//      LARGE = clear_row( matrix, node_dof(0), 0.0 );
//      node_rhs(0) = 0.0;
    }
    m_rhs->ReplaceGlobalValues(face_dof, face_rhs);
  }
}



void ParallelPlateFlowSolverFV::MassRateBC( PhysicalDomainT& domain,
                                            ObjectDataStructureBaseT& object,
                                            BoundaryConditionBase* bc,
                                            const lSet& set,
                                            realT time,
                                            realT dt )
{


  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& faceGhostRank  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");


  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);


  for( lSet::const_iterator faceID=set.begin() ; faceID!=set.end() ; ++faceID )
  {
    face_dof(0) = trilinos_index[*faceID];

    if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]==0 )
    {

      const realT massRate = bc->GetValue(domain.m_feFaceManager,faceID,time) ;

      if( faceGhostRank[*faceID] < 0 && flowFaceType[*faceID]==0 )
      {
        face_rhs(0) -= massRate * dt;
      }
      else
      {
  //      LARGE = clear_row( matrix, node_dof(0), 0.0 );
  //      node_rhs(0) = 0.0;
      }
      m_rhs->SumIntoGlobalValues(face_dof, face_rhs);
    }
  }
}


realT ParallelPlateFlowSolverFV::TwoFacePermeability(const Array1dT<R1Tensor>& edgeCenters,
                                                           const rArray1d& edgeLengths,
                                                           const Array1dT<R1Tensor>& faceCenters,
                                                           const rArray1d& apertures,
                                                           const localIndex eg,
                                                           const localIndex r,
                                                           const localIndex s,
                                                           const Array1dT<rArray1d>* const dwdu,
                                                           rArray1d* const dkdu_r,
                                                           rArray1d* const dkdu_s)
{

  R1Tensor edgeCenter = edgeCenters[eg];
  R1Tensor lr, ls;

  lr = edgeCenter;
  lr -= faceCenters[r];

  ls = edgeCenter;
  ls -= faceCenters[s];

  realT l_edge = edgeLengths[eg];

  realT wr = BoundedAperture(apertures[r]);
  realT ws = BoundedAperture(apertures[s]);

  realT kappa = CalculatePermeability(lr.L2_Norm(), ls.L2_Norm(), wr, ws, l_edge, m_mu, m_SHP_FCT);


  if( dkdu_r!=NULL && dkdu_s!=NULL && dwdu!=NULL )
  {
    if( dwdu->size() > std::max(r,s) )
    {
      const realT denom = lr.L2_Norm()*ws*ws*ws+ls.L2_Norm()*wr*wr*wr;

      *dkdu_r  = (*dwdu)[r] ;
      *dkdu_r *= lr.L2_Norm()*ws*ws*ws/wr;

      *dkdu_s  = (*dwdu)[s] ;
      *dkdu_s *= ls.L2_Norm()*wr*wr*wr/ws;

      *dkdu_r *= 3 * kappa / denom  * BoundedApertureDerivative(apertures[r]);
      *dkdu_s *= 3 * kappa / denom  * BoundedApertureDerivative(apertures[s]);
    }
  }

  return kappa;
}


/// Apply a pressure boundary condition to a given set of edges
void ParallelPlateFlowSolverFV::PressureBoundaryCondition( PhysicalDomainT& domain,
                                                           ObjectDataStructureBaseT& object ,
                                                           BoundaryConditionBase* const bc,
                                                           const lSet& set,
                                                           const realT time,
                                                           const realT dt,
                                                           const int dofOffset )
{
 

  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

  iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  
  const Array1dT<R1Tensor>& faceCenters_n  = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+std::string("_n") );
  const Array1dT<R1Tensor>& faceCenters_np1  = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

//  const rArray1d& faceFluidVolume_n  = domain.m_feFaceManager.GetFieldData<realT>( "FluidVolume_n" );
  const rArray1d& faceFluidVolume_np1  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  const Array1dT<R1Tensor>& edgeCenters_n  = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr+std::string("_n") );
  const Array1dT<R1Tensor>& edgeCenters_np1  = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );

  const rArray1d& edgeLengths_n  = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr+std::string("_n") );
  const rArray1d& edgeLengths_np1  = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  const rArray1d& apertures_n  = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr+std::string("_n") );
  const rArray1d& apertures_np1  = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );

//  const rArray1d& boundedApertures_n  = domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture_n" );
//  const rArray1d& boundedApertures_np1  = domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture" );

//  const rArray1d& faceFluidMass_np1  = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();


  const rArray1d& faceFluidPressure_n   = domain.m_feFaceManager.GetFieldData<realT>("Pressure_n");
  const rArray1d& faceFluidPressure_np1  = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  const rArray1d& faceDensity_n   = domain.m_feFaceManager.GetFieldData<realT>("FluidDensity_n");
  const rArray1d& faceDensity_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();

  const rArray1d& dPdM = domain.m_feFaceManager.GetFieldData<realT>("dPdM");

  rArray1d& edgePermeability = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);

//  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
 

  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     edge_rhs(1);
  Epetra_SerialDenseMatrix     edge_matrix(1,1);
    
  Epetra_IntSerialDenseVector  displacementDofIndex;
  Epetra_IntSerialDenseVector  dispSubFaceDofIndex(1);
  Epetra_SerialDenseMatrix     edge_matrix_dispDOF;

  // loop over edges, find permeabilities and apply boundary conditions.

  lSet::const_iterator eg=set.begin() ;

  for( localIndex iset=0; iset < set.size() ; ++iset, ++eg ){

    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) )
    {

      lArray1d& faces = itr->second; 
      
      realT bc_pressure_n = bc->GetValue(domain.m_feEdgeManager,eg,time);
      realT bc_pressure_np1 = bc->GetValue(domain.m_feEdgeManager,eg,time+dt);

      // rhoP = rho * Pressure at BC
      realT rhoP_n = rho_EOS(bc_pressure_n,m_bulk_modulus,m_rho_o);
      realT rhoP_np1 = rho_EOS(bc_pressure_np1,m_bulk_modulus,m_rho_o);

      for(size_t kf = 0; kf < faces.size(); ++kf)
      {
        
        const localIndex r = faces[kf];


        const int dim=3;



        // calculate displacement DOF for these two faces
        const localIndex rnodes = domain.m_feFaceManager.m_toNodesRelation[r].size();
        {
          dispSubFaceDofIndex(0) = trilinos_index[r] + dofOffset;
          displacementDofIndex.Resize( dim*2*rnodes );
          edge_matrix_dispDOF.Reshape( 1, dim*2*rnodes );
          edge_matrix_dispDOF.Scale(0.0);
          for( localIndex a=0 ; a<rnodes ; ++a )
          {
            for( int i=0 ; i<dim ; ++i )
            {
              displacementDofIndex(dim*a*2+i) = m_dwdu_dof(r)(dim*a*2+i);
              displacementDofIndex(dim*a*2+dim+i) = m_dwdu_dof(r)(dim*a*2+dim+i);

            }
          }
        }


        // calculate edge permeability for face
        R1Tensor la_old = edgeCenters_n[*eg]; la_old -= faceCenters_n[r];
        R1Tensor la_new = edgeCenters_np1[*eg]; la_new -= faceCenters_np1[r];

        realT kappa_n = CalculatePermeability( la_old.L2_Norm(), BoundedAperture(apertures_n[r]), edgeLengths_n[*eg], m_mu, m_SHP_FCT);
        realT kappa_np1 = CalculatePermeability( la_new.L2_Norm(), BoundedAperture(apertures_np1[r]), edgeLengths_np1[*eg], m_mu, m_SHP_FCT);
        rArray1d dKappa_du_rterm = m_dwdu(r);
        dKappa_du_rterm *= 3 * kappa_np1 / BoundedAperture(apertures_np1[r]) * BoundedApertureDerivative(apertures_np1[r]);

        edgePermeability[*eg] = kappa_np1;

//        kappa_n = 0;
//        kappa_np1 = kappa_n;
//        dKappa_du_rterm *= 0.0;


        const realT dPdM_r = dPdM[r];


        // matrix
        const realT dp_sr = bc_pressure_np1 - faceFluidPressure_np1[r];

        realT rho_sr = faceDensity_np1[r] + rhoP_np1;

        if( bc_pressure_np1 > faceFluidPressure_np1[r] )
        {

        }
        else
        {

        }


        edge_matrix(0,0) = 0.5 * m_phi * kappa_np1 * ( dp_sr / faceFluidVolume_np1[r] - rho_sr * dPdM_r );

        // rhs
        edge_rhs(0) = -(  0.5 * ( 1.0 - m_phi ) * kappa_n * ( rhoP_n + faceDensity_n[r] ) * ( bc_pressure_n - faceFluidPressure_n[r] )
                        + 0.5 * m_phi * kappa_np1 * rho_sr * dp_sr );


        face_dof(0) =  dofOffset + trilinos_index[r];
        m_matrix->SumIntoGlobalValues(face_dof, edge_matrix);
    
        m_rhs->SumIntoGlobalValues(face_dof, edge_rhs);


        realT area_r = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, r, false );
        const realT dPdV_r = - dPdM_r * faceDensity_np1[r] ;
        for( unsigned int i=0 ; i<dim*2*rnodes ; ++i )
        {
          edge_matrix_dispDOF(0,i) = 0.5 * m_phi * ( rho_sr*dp_sr*dKappa_du_rterm(i)
                                                    - kappa_np1*( faceDensity_np1[r]/faceFluidVolume_np1[r]*area_r * dp_sr
                                                                  + rho_sr*dPdV_r*area_r ) * m_dwdu(r)(i) * m_dwdw(r) );
        }

        m_matrix->SumIntoGlobalValues(dispSubFaceDofIndex, displacementDofIndex, edge_matrix_dispDOF);


    
      }
    }
  }
}



realT ParallelPlateFlowSolverFV::CheckSolution( const realT* const local_solution,
                                                const PhysicalDomainT& domain,
                                                const localIndex dofOffset )
{
  realT rval = 1.0;


  // face fields
  const iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  const iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  const rArray1d& mass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      const int intDofOffset = dofOffset;
      int lid = m_rowMap->LID(intDofOffset+trilinos_index[*kf]);
//      if( mass[*kf] + local_solution[lid] <= 0.0 )
//      {
//        rval =
//      }
      realT ratio = fabs(local_solution[lid]) / mass[*kf];
      if( ratio > 0.9 )
      {
        realT temp = 0.9 / ratio;
        if( temp < rval )
        {
          rval = temp;
//          std::cout<<"dmass = "<<local_solution[lid]<<std::endl;
//          std::cout<<"mass  = "<<mass[*kf]<<std::endl;
//          std::cout<<"scalingFactor = "<<rval<<std::endl;
        }
      }
    }
  }


  return rval;
}



void ParallelPlateFlowSolverFV::PropagateSolution( const realT* const local_solution,
                                                   const realT scalingFactor,
                                                   PhysicalDomainT& domain,
                                                   const localIndex dofOffset  )
{
  // face fields
  const iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  const iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d& mass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();

  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(is_ghost[*kf] < 0)
    {
      const int intDofOffset = dofOffset;
      int lid = m_rowMap->LID(intDofOffset+trilinos_index[*kf]);
      mass[*kf] += scalingFactor*local_solution[lid];
//      if( mass[*kf] < 0.0 )
//        mass[*kf] = 0.0;
    }
  }

}

void ParallelPlateFlowSolverFV::PostSyncConsistency( PhysicalDomainT& domain,
                                                     SpatialPartition& partition )
{
  // re-sync ghost nodes
  partition.SynchronizeFields(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);

}




void ParallelPlateFlowSolverFV::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(FaceCenterStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("massRate");
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(ApertureStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_TrilinosIndexStr);
  syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(VolumetricFluxStr);

  partition.SetBufferSizes(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);
}


double ParallelPlateFlowSolverFV::TimeStep( const realT& time,
                                                        const realT& dt,
                                                        const int cycleNumber,
                                                        PhysicalDomainT& domain,
                                                        const sArray1d& namesOfSolverRegions ,
                                                        SpatialPartition& partition,
                                                        FractunatorBase* const fractunator )
{
  realT dt_return = dt;

  DefineFlowSets( domain );

  RegisterTemporaryFields( domain );

  GenerateParallelPlateGeometricQuantities( domain, time,dt );

  FillTemporaryFields( domain );



  SetupSystem (domain,partition, time);


#if USECPP11==1
  m_matrix   = std::make_shared<Epetra_FECrsMatrix>(Copy,*m_sparsity);
  m_solution = std::make_shared<Epetra_FEVector>(*m_rowMap);
  m_rhs      = std::make_shared<Epetra_FEVector>(*m_rowMap);
#else
  m_matrix   = new Epetra_FECrsMatrix(Copy,*m_sparsity);
  m_solution = new Epetra_FEVector(*m_rowMap);
  m_rhs      = new Epetra_FEVector(*m_rowMap);

#endif

  int dummy;
  double* res = NULL;
  m_rhs->ExtractView(&res,&dummy);

  for( int i=0 ; i<this->m_numerics.m_maxIterNewton ; ++i )
  {
    m_rhs->Scale(0.0);
    m_matrix->Scale(0.0);

    GenerateParallelPlateGeometricQuantities( domain, time, dt );
    UpdateEOS( time, dt, domain );

    Epetra_System junk;
    const realT massScale = Assemble(domain, junk, time, dt );

    realT residual = 0;
    m_rhs->Norm2( &residual);

    if( residual / massScale < this->m_numerics.m_tolNewton )
      break;

//    m_matrix->Print(std::cout);
//    m_rhs->Print(std::cout);

    Solve(domain,partition, time, dt );

//    m_solution->Print(std::cout);

  }

  UpdateEOS( time, dt, domain );

  UpdateFlux( time, dt, domain,partition );

  DeregisterTemporaryFields( domain );
  return dt_return;
}

void ParallelPlateFlowSolverFV::DefineFlowSets( PhysicalDomainT& domain )
{
  if( m_flowFaceSetName.empty() )
  {
  FaceManagerT& faceManager = domain.m_feFaceManager;
  EdgeManagerT& edgeManager = domain.m_feEdgeManager;

  const iArray1d& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");
  const iArray1d& flowEdgeType = edgeManager.GetFieldData<int>("flowEdgeType");
  const Array1dT<lSet>& edgeToFlowFaces = edgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");





  m_faceSet.clear();
  m_edgesToFaces.clear();

  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 0 )
    {
      m_faceSet.insert( kf );
    }
  }

  m_numFaces = m_faceSet.size();

  m_faceDofMap.clear();
  lSet::const_iterator si=m_faceSet.begin();
  for(localIndex i =0; i < m_numFaces; ++i, ++si){
    localIndex f = *si;
    m_faceDofMap[f] = i;
  }

  iArray1d& ffCount = edgeManager.GetFieldData<int>("FlowFaceCount"); // debug

  for( localIndex ke=0 ; ke < edgeManager.DataLengths() ; ++ke )
  {
    if( flowEdgeType[ke] == 0 && !(edgeToFlowFaces[ke].empty()) )
    {
      m_edgesToFaces[ke].assign( edgeToFlowFaces[ke].begin(), edgeToFlowFaces[ke].end() ) ;
      ffCount[ke] = m_edgesToFaces[ke].size();

    }
  }
  }

}


void ParallelPlateFlowSolverFV::GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,realT time,realT dt )
{

  Array1dT<R1Tensor>& faceCenter  = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr );
  rArray1d& mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& boundedAperture = domain.m_feFaceManager.GetFieldData<realT>( "BoundedAperture" );



  // update face quantities
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf ) {

    domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *kf, faceCenter[*kf]);

    if(m_doApertureUpdate){
      R1Tensor gap;
      R1Tensor N;
      N = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *kf );

      gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, *kf );
      aperture[*kf] = Dot(gap,N) + m_zeroApertureVolume;

      if( m_boundPhysicalAperture )
      {
        boundedAperture[*kf] = BoundedAperture(aperture[*kf]);
      }
      else
      {
        boundedAperture[*kf] = aperture[*kf];
      }

    }


    fluidVolume[*kf] = ( boundedAperture[*kf]) * domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, *kf );

    if( fluidVolume[*kf]<=0.0 )
      throw GPException("you have a negative volume finite volume element!!");

    if( isZero( mass[*kf] ) )
    {
      mass[*kf] = fluidVolume[*kf] * this->m_rho_o;
    }

  }


  // update edge properties
  iArray1d& edge_is_ghost            = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  Array1dT<R1Tensor>& edgeCenter_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  Array1dT<realT>& edgeLength_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0 )
    {
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter_new[eg] );
      edgeLength_new[eg] = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, eg);
    }
  }
}


void ParallelPlateFlowSolverFV::CalculateMassRate( PhysicalDomainT& domain,
                                                   SpatialPartition& partition,
                                                   realT time, realT dt )
{

}

// set initial fluid density based on known pressure;
void ParallelPlateFlowSolverFV::InitializeDensity( PhysicalDomainT& domain)
{
  rArray1d& mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& density  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const rArray1d& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  for( lSet::const_iterator fc=m_faceSet.begin() ; fc!=m_faceSet.end() ; ++fc ) {
    density[*fc] = rho_EOS(pressure[*fc],m_bulk_modulus,m_rho_o );
    mass[*fc] = density[*fc]*fluidVolume[*fc];
  }
}


void ParallelPlateFlowSolverFV::UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain )
{
  rArray1d& mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& density  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  rArray1d& dPdM = domain.m_feFaceManager.GetFieldData<realT>("dPdM");

  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  for( lSet::const_iterator fc=m_faceSet.begin() ; fc!=m_faceSet.end() ; ++fc )
  {
    density[*fc] = mass[*fc] / fluidVolume[*fc];
    realT P = P_EOS(density[*fc],m_bulk_modulus,m_rho_o );

    dPdM[*fc] = dPdRho_EOS(m_bulk_modulus,m_rho_o)/fluidVolume[*fc];

    const realT Pmin = 1.0e5;
    if( false&&P<0 )
    {

      const realT pi = 3.14159265358979323846;
      const realT a = 2 * Pmin / pi;
//      dPdM[*fc] *= 1 / ( 1 + pow( (P-Pmin)/a, 2 ) ) ;
//      P = a * atan( (P-Pmin) / a ) + Pmin;

      dPdM[*fc] *= 1 / ( 1 + pow( P / a, 2 ) ) ;
      P = a * atan( P / a );


    }
    pressure[*fc] = P;
    // propagate pressure to children
    lArray1d& childFaces = domain.m_feFaceManager.m_childIndices[*fc];
    for(unsigned i =0; i < childFaces.size(); ++i){
      pressure[childFaces[i]] = P;
    }
  }
}


void ParallelPlateFlowSolverFV::UpdateFlux( const realT time,
                                            const realT dt,
                                            PhysicalDomainT& domain,
                                            SpatialPartition& partition  )
{
  
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
  
  ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolverFV::PressureBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time, dt );

  flowVelocity *= 0.5;
    
}




/// Update the velocity fluxes on the boundary faces
void ParallelPlateFlowSolverFV::
       PressureBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                 ObjectDataStructureBaseT& object ,
                                                 BoundaryConditionBase* const bc,
                                                 const lSet& set,
                                                 const realT time,
                                                 const realT dt )
{
  
  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  
  Array1dT<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
    
  iArray1d& is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
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


void ParallelPlateFlowSolverFV::CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain )
{

  const rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  rArray1d& edgePermeability = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& flowRate = domain.m_feFaceManager.GetFieldData<realT>("flowRate");
//  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");


  iArray1d& tracer = domain.m_feFaceManager.GetFieldData<int>( "tracer" );


  rArray1d& massRate = domain.m_feEdgeManager.GetFieldData<realT>("massRate");


  const rArray1d& edgeLength = domain.m_feEdgeManager.GetFieldData<realT>("length");
  const Array1dT<R1Tensor>& edgeCenter = domain.m_feEdgeManager.GetFieldData<R1Tensor>("center");

  //  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");


  const Array1dT<lSet>& edgeToFlowFaces = domain.m_feEdgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");

  R1Tensor la, lb;

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();



  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    flowRate[kf] = 0.0;
  }


  // loop over all edges
  for( localIndex ke=0 ; ke<domain.m_feEdgeManager.DataLengths() ; ++ke )
  {

    // only do work on edges that are active flow edges
    if( flowEdgeType[ke]== 0 )
    {
      realT w = edgeLength[ke];


      // we have two different ways of calculating the flow. One for edges attached to only 2
      // faces, and one for more than two faces.
      const unsigned int numFlowFaces = edgeToFlowFaces[ke].size();

      if( numFlowFaces == 2 )
      {
        lSet::const_iterator edgeToFlowFace = edgeToFlowFaces[ke].begin();
        const localIndex faceIndex0 = *edgeToFlowFace;
        const localIndex faceIndex1 = *(++edgeToFlowFace);


        la = edgeCenter[ke];
        la -= faceCenter[faceIndex0];

        lb = edgeCenter[ke];
        lb -= faceCenter[faceIndex1];

        realT norm_la = la.L2_Norm();
        realT norm_lb = lb.L2_Norm();

        realT permeabilityAperture0 = std::min(aperture[faceIndex0], m_max_aperture);
        realT permeabilityAperture1 = std::min(aperture[faceIndex1], m_max_aperture);

        edgePermeability[ke] = PPFS::CalculatePermeability( norm_la , norm_lb,
                                                      permeabilityAperture0, permeabilityAperture1,
                                                      w,m_mu,m_SHP_FCT);
        // determine the mass flux across edge at t_n+1/2
        massRate[ke] = edgePermeability[ke] * ( -( faceFluidDensity[faceIndex0] * faceFluidPressure[faceIndex0])
            + ( faceFluidDensity[faceIndex1] * faceFluidPressure[faceIndex1] ) );

        if (domain.m_feFaceManager.m_toEdgesRelation[faceIndex0].size() == 2 )
        {//2D
          if (domain.m_feFaceManager.m_toEdgesRelation[faceIndex0][0] == ke)
          {
            flowRate[faceIndex0] += 0.5*massRate[ke] / m_rho_o;
          }
          else
          {
            flowRate[faceIndex0] -= 0.5*massRate[ke] / m_rho_o;
          }

          if (domain.m_feFaceManager.m_toEdgesRelation[faceIndex1][0] == ke)
          {
            flowRate[faceIndex1] -= 0.5*massRate[ke] / m_rho_o;
          }
          else
          {
            flowRate[faceIndex1] += 0.5*massRate[ke] / m_rho_o;
          }
        }
        else
        {//3D
          flowRate[faceIndex0] += 0.5 * fabs(massRate[ke]) / m_rho_o;
          flowRate[faceIndex1] += 0.5 * fabs(massRate[ke]) / m_rho_o;
        }


        faceFluidMass[faceIndex0] +=  massRate[ke] * dt;
        faceFluidMass[faceIndex1] -=  massRate[ke] * dt;

        if( tracer[faceIndex0] == 1 )
        {
          tracer[faceIndex1] = 1;
        }
        if( tracer[faceIndex1] == 1 )
        {
          tracer[faceIndex0] = 1;
        }


        realT thisdt = 6 * m_mu *( norm_la + norm_lb ) * ( norm_la + norm_lb )
                     / ( m_bulk_modulus * 0.25 * (permeabilityAperture0+permeabilityAperture1)* (permeabilityAperture0+permeabilityAperture1) ) ;

        if( thisdt < m_stabledt.m_maxdt )
          m_stabledt.m_maxdt = thisdt;



      }
      else if( numFlowFaces > 2 )
      {
        realT rhoP = 0.0;
        realT sumK = 0.0;
        Array1dT<R1Tensor> length(numFlowFaces);
        rArray1d k(numFlowFaces);
        rArray1d q(numFlowFaces);
        rArray1d kRhoP(numFlowFaces);

        lSet::const_iterator faceIndex=edgeToFlowFaces[ke].begin();

        for( localIndex kf=0 ; kf<numFlowFaces ; ++kf, ++faceIndex)
        {

          length[kf] = edgeCenter[ke];
          length[kf] -= faceCenter[*faceIndex];

          realT permeabilityAperture = std::min(aperture[*faceIndex], m_max_aperture);

          k[kf] = PPFS::CalculatePermeability( length[kf].L2_Norm(),
                                         permeabilityAperture,
                                         w, m_mu, m_SHP_FCT );

          sumK += k[kf];

          kRhoP[kf] = k[kf] * faceFluidDensity[*faceIndex] * faceFluidPressure[*faceIndex];
          rhoP += kRhoP[kf];



          realT thisdt = 6 * m_mu * 4 * Dot(length[kf],length[kf])
          / ( m_bulk_modulus * permeabilityAperture * permeabilityAperture ) ;

          if( thisdt < m_stabledt.m_maxdt )
            m_stabledt.m_maxdt = thisdt;
        }


        rhoP /= sumK ;
        faceIndex=edgeToFlowFaces[ke].begin();

        int tracerActive = 0;
        for( localIndex kf=0 ; kf<numFlowFaces ; ++kf, ++faceIndex )
        {
          q[kf] = k[kf] * ( faceFluidDensity[*faceIndex] * faceFluidPressure[*faceIndex] - rhoP );

          faceFluidMass[*faceIndex] -= q[kf] * dt;

          if (domain.m_feFaceManager.m_toEdgesRelation[*faceIndex].size() == 2 )
          {//2D
            if (domain.m_feFaceManager.m_toEdgesRelation[*faceIndex][0] == ke)
            {
              flowRate[*faceIndex] -= 0.5*q[kf] / m_rho_o;
            }
            else
            {
              flowRate[*faceIndex] += 0.5*q[kf] / m_rho_o;
            }
          }
          else //3D
          {
            flowRate[*faceIndex] += 0.5* fabs(q[kf]) / m_rho_o;
          }




          if( tracer[*faceIndex] == 1 )
          {
            tracerActive = 1;
          }

        }

        for( faceIndex=edgeToFlowFaces[ke].begin() ; faceIndex!=edgeToFlowFaces[ke].end() ; ++faceIndex )
        {
          if( tracerActive == 1 && tracer[*faceIndex] == 0 )
          {
            tracer[*faceIndex] = 1;
          }
        }

      }
    }
  }



  if (m_stabledt.m_maxdt != std::numeric_limits<double>::max())
  {
    m_stabledt.m_maxdt *= this->m_courant;
  }
}

void ParallelPlateFlowSolverFV::CalculateCarterLeakOff( const realT time,
                                                      const realT dt,
                                                      PhysicalDomainT& domain )
{
  rArray1d* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");
  rArray1d* totalLeakedVolume = domain.m_feFaceManager.GetFieldDataPointer<realT>("totalLeakedVolume");
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();



  if (totalLeakedVolume != NULL)
  {
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
    {
      if (flowFaceType[kf] == 0 && faceFluidPressure[kf] > 0.0 && (*initialSaturatedTime)[kf] != std::numeric_limits<realT>::max())
      {
        localIndex face0, face1;
        if (domain.m_feFaceManager.m_toEdgesRelation[kf].size() == 2 ) //2D
        {
          face0 = domain.m_feFaceManager.m_childIndices[kf][0];
          face1 = domain.m_feFaceManager.m_childIndices[kf][1];
        }
        else
        {
          face0 = kf;
          face1 = domain.m_feFaceManager.m_childIndices[kf][0];
        }
        realT leakOffVelocity = m_leakoffCoef / sqrt( time + dt * 0.5 - (*initialSaturatedTime)[kf]);
        if (m_pressureDependentLeakoff == 1)
        {
          leakOffVelocity *= faceFluidPressure[kf] - m_farFieldPorePressure;
        }
        realT leakOffMassInc = leakOffVelocity * dt *  domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, kf ) * m_rho_o;
        leakOffMassInc *= 2;  //Counting for both sides.

        if (leakOffMassInc < faceFluidMass[kf])
        {
          faceFluidMass[kf] -= leakOffMassInc;
          (*totalLeakedVolume)[face0] += leakOffVelocity * dt ;
          (*totalLeakedVolume)[face1] += leakOffVelocity * dt ;
        }
        else
        {
          (*totalLeakedVolume)[face0] += faceFluidMass[kf] * 0.5;
          (*totalLeakedVolume)[face1] += faceFluidMass[kf] * 0.5;
          faceFluidMass[kf] = 0.0;
        }


      }
    }
  }


}


// We know the combined flow rate and have to distribute the given flow rate among the faces in the set.
void ParallelPlateFlowSolverFV::ApplyFluxBoundaryCondition( const realT time,
                                                          const realT dt,
                                                          const int cycleNumber,
                                                          const int rank,
                                                          PhysicalDomainT& domain )
{
  m_dT = dt;

  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  //const rArray1d& edgeLength = domain.m_feEdgeManager.GetFieldData<realT>("length");
  const Array1dT<R1Tensor>& edgeCenter = domain.m_feEdgeManager.GetFieldData<R1Tensor>("center");
  const Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();


  for( Array1dT<BoundaryConditionBase*>::const_iterator bcItr=domain.m_feFaceManager.m_bcData.begin() ; bcItr!=domain.m_feFaceManager.m_bcData.end() ; ++ bcItr )
  {
    // check to see if the requested field has a boundary condition applied to it.
    BoundaryConditionBase* bc = *bcItr;
    if( streq( bc->GetFieldName(time), "combinedFlowRate") )
    {
      for(localIndex i =0; i < bc->m_setNames.size(); ++i)
      {
        std::map< std::string, lSet >::iterator setMap = domain.m_feFaceManager.m_Sets.find( bc->m_setNames[i] );
        if( setMap != domain.m_feFaceManager.m_Sets.end() )
        {
          lSet& set = setMap->second;

          lSet::const_iterator b=set.begin();
          realT qTotal = bc->GetValue(domain.m_feFaceManager,b,time); // The iterator b in this function is not doing anything.
          int nInlets = set.size();

          realT sumKP = 0.0;
          realT sumK = 0.0;
          R1Tensor length;
          rArray1d k(nInlets);
          realT kP;

          nInlets = 0;
          localIndex kf=0;
          for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a, ++kf)
          {
            if (flowFaceType[*a] == 0 && isGhost[*a]<0)
            {
              nInlets++;
              localIndex ke = domain.m_feFaceManager.m_toEdgesRelation[*a][0];
              length = edgeCenter[ke];
              length -= faceCenter[*a];
              realT w = 1.0;
              realT permeabilityAperture = std::min(aperture[*a], m_max_aperture);
              k[kf] = PPFS::CalculatePermeability( length.L2_Norm(),
                                             permeabilityAperture,
                                             w, m_mu, m_SHP_FCT );
              sumK += k[kf];

              kP = k[kf] * faceFluidPressure[*a];
              sumKP += kP;
            }
          }

          int myNInlets = nInlets;
          realT mySumK = sumK;
          realT mySumKP = sumKP;
          int myRank0 = nInlets * rank;
          int rank0;

          MPI_Allreduce(&myNInlets, &nInlets, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(&mySumK, &sumK, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(&mySumKP, &sumKP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

          if (cycleNumber%100 == 0) MPI_Allreduce(&myRank0, &rank0, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

          if (myNInlets >=1)
          {
            realT pInj;
            sumKP += qTotal;
            pInj = sumKP / sumK;

            kf = 0;
            for( lSet::const_iterator a=set.begin() ; a!=set.end() ; ++a, ++kf)
            {
              if (flowFaceType[*a] == 0)
              {
                realT q = k[kf] * (pInj - faceFluidPressure[*a]);
                faceFluidMass[*a] += q * dt * std::max(faceFluidDensity[*a], m_rho_o);
              }
            }
            if (cycleNumber%100 == 0 && myRank0 == rank0)
            {
              std::cout <<"Flow rate BC applied to set: " << bc->m_setNames[i] << ", t=" << time << ", P_inj=" << pInj <<" rank: " <<rank << std::endl;
            }

          }
        }
      }
    }
  }



  // Edge fluxes
//  BoundaryConditionFunctions::ApplyBoundaryCondition<realT>(this, &ParallelPlateFlowSolver::FlowControlledBoundaryCondition,
//                                domain, domain.m_feEdgeManager, "FixedFlowRate", time );

}

// constant flow into faces
void ParallelPlateFlowSolverFV::FlowControlledBoundaryCondition( PhysicalDomainT& domain,
                                                               ObjectDataStructureBaseT& object ,
                                                               BoundaryConditionBase* bc ,
                                                               const lSet& set,
                                                               realT time ){


  //iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();


  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  const rArray1d& edgeLengths = domain.m_feEdgeManager.GetFieldData<realT>("length");
  //const Array1dT<R1Tensor>& edgeCenter = domain.m_feEdgeManager.GetFieldData<R1Tensor>("center");
  //const Array1dT<R1Tensor>& faceCenter = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );
  const rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  //const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  const Array1dT<lSet>& edgeToFlowFaces = domain.m_feEdgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");

  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);
  Epetra_SerialDenseMatrix     face_matrix(1,1);

  // loop over edges, apply flux to faces

  for( lSet::const_iterator eg=set.begin(); eg != set.end(); ++eg  ){

    realT qdt = bc->GetValue(domain.m_feEdgeManager,eg,time)*m_dT;

    for( lSet::const_iterator fc=edgeToFlowFaces[*eg].begin() ; fc!=edgeToFlowFaces[*eg].end() ; ++fc )
    {
        if(flowFaceType[*fc] == 0){
          const realT area = edgeLengths[*eg]*BoundedAperture(apertures[*fc]);

          const realT massFlux = qdt*area*std::max(faceFluidDensity[*fc], m_rho_o);
          faceFluidMass[*fc] += massFlux;
        }
    }
  }
}


void ParallelPlateFlowSolverFV::CalculateApertureDerivatives( const FaceManagerT& faceManager,
                                                              const NodeManagerT& nodeManager )
{
  const iArray1d* const trilinosIndexNode = nodeManager.GetFieldDataPointer<int>("IMS_0_GlobalDof");

  if( trilinosIndexNode!=NULL )
  {
    m_dwdu.resize(faceManager.DataLengths());
    m_dwdu_dof.resize( faceManager.DataLengths() );
    m_dwdw.resize( faceManager.DataLengths() );
    m_dwdw = 1.0;

    const int dim=3;

    const iArray1d& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");
    const rArray1d& apertures_np1 = faceManager.GetFieldData<realT>( ApertureStr );
    const OrderedVariableOneToManyRelation& childFaceIndex = faceManager.GetVariableOneToManyMap( "childIndices" );


    const Array1dT<R1Tensor>& faceNormal = faceManager.GetFieldData<R1Tensor>("faceNormal0");

    // set aperture derivatives
    for( localIndex r=0 ; r<faceManager.DataLengths() ; ++r )
    {
      if( /*face_is_ghost[r] < 0 &&*/ flowFaceType[r] == 0 )
      {
        m_dwdu(r).resize( faceManager.m_toNodesRelation[r].size() * dim * 2 );
        m_dwdu_dof(r).resize( faceManager.m_toNodesRelation[r].size() * dim * 2 );

        if( m_boundPhysicalAperture )
        {
          m_dwdw(r) = BoundedApertureDerivative( apertures_np1[r] );
        }
        const localIndex numNodes = faceManager.m_toNodesRelation[r].size();
        for( localIndex a=0 ; a<numNodes ; ++a )
        {
          const localIndex faceIndex[2] = { r, childFaceIndex[r][0] };

          const R1Tensor N[2] = { faceManager.FaceNormal( nodeManager, faceIndex[0] ),
                                  faceManager.FaceNormal( nodeManager, faceIndex[1] )};

//          const R1Tensor N[2] = { faceNormal[faceIndex[0]],
 //                                 faceNormal[faceIndex[1]] };

          R1Tensor Nbar = N[0];
          Nbar -= N[1];
          Nbar.Normalize();

          const localIndex aa = a == 0 ? a : numNodes - a;
          const localIndex node0 = faceManager.m_toNodesRelation[faceIndex[0]][a];
          const localIndex node1 = faceManager.m_toNodesRelation[faceIndex[1]][aa];


          int nodeDofIndex[2] ;
          faceNodePairIndexing( a, dim, nodeDofIndex);

          for( int i=0 ; i<dim ; ++i )
          {
            m_dwdu(r)(nodeDofIndex[0]+i) = -Nbar[i]/faceManager.m_toNodesRelation[r].size() ;
            m_dwdu(r)(nodeDofIndex[1]+i) =  Nbar[i]/faceManager.m_toNodesRelation[r].size() ;

            m_dwdu_dof(r)(nodeDofIndex[0]+i) = dim*(*trilinosIndexNode)[node0]+i;
            m_dwdu_dof(r)(nodeDofIndex[1]+i) = dim*(*trilinosIndexNode)[node1]+i;
          }
        }
      }
    }
  }
}

/// Register solver in the solver factory
REGISTER_SOLVER( ParallelPlateFlowSolverFV )
