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
 * @file ParallelPlateProppantFlowSolverImplicit.cpp
 * @author walsh24
 * @date February 21, 2012
 */

#include "ParallelPlateProppantFlowSolverImplicit.h"
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
  realT LARGE = 1e64;
  const std::string oldStr = "_old";
}

int ParallelPlateProppantFlowSolverImplicit::m_instances = 0;


ParallelPlateProppantFlowSolverImplicit::ParallelPlateProppantFlowSolverImplicit( const std::string& name, ProblemManagerT* const pm):
ParallelPlateFlowSolverBase(name,pm),
m_faceSet(),
m_numFaces(0),
m_faceDofMap(),
m_edgesToFaces(),
m_phi(1.0),
this_mpi_process(pm->m_epetraComm.MyPID()),
n_mpi_processes(pm->m_epetraComm.NumProc()),
m_epetra_comm(pm->m_epetraComm),
row_map(),
sparsity(),
matrix(),
solution(),
rhs(),
syncedFields(),
syncedFieldsB(),
m_TrilinosIndexStr(),
m_numerics(),
m_cycleNumber(0)
{
 // ReadXML(hdn);
  ++m_instances; 
  m_TrilinosIndexStr = "TwoDIMPPFS_" +  toString<int>(m_instances) + "_GlobalDof";
}

ParallelPlateProppantFlowSolverImplicit::~ParallelPlateProppantFlowSolverImplicit()
{
  // TODO Auto-generated destructor stub
}

void ParallelPlateProppantFlowSolverImplicit::ReadXML( TICPP::HierarchicalDataNode* hdn )
{
  ParallelPlateFlowSolverBase::ReadXML(hdn);

  // Mixed difference parameter
  m_phi = hdn->GetAttributeOrDefault<realT>("phi",1.0); // backward difference by default.
  
  // Linear Solver
  m_numerics.m_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-10);
  m_numerics.m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);

  // Flags
  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("UpdateAperture",false);
  m_verboseFlag = hdn->GetAttributeOrDefault<bool>("Verbose",false);
  m_doErosion = hdn->GetAttributeOrDefault<bool>("UseErosionModel",false);

  //m_relPermFunctionStr = hdn->GetAttributeOrDefault("RelativePermeabilityFunction","");

  m_apertureRelaxationTime = hdn->GetAttributeOrDefault("apertureRelaxationTime","10s");

  // Proppant

  TICPP::HierarchicalDataNode* proppantNode = hdn->GetChild("Proppant");
   if(!proppantNode) throw GPException("Error ParallelPlateProppantFlowSolverImplicit: No proppant data provided.");
   m_proppantData.ReadXML(proppantNode);


   TICPP::HierarchicalDataNode* leakoffNode = hdn->GetChild("LeakoffModel");

   if(leakoffNode)
   {
	   m_leakoffModel.ReadXML(leakoffNode);
   } else {
	   std::cout << "No leakoff model specified." <<std::endl;
   }

   // Barton Joint Parameters
   std::string temp = hdn->GetAttributeString("BartonJointParameters"); // aperture at zero effective stress; reference stress; aperture at ref stress
     if( !temp.empty() )
     {
       R1Tensor tempArray;
       tempArray.StrVal( temp );
       m_wZeroStress = tempArray[0];
       realT stressRef = tempArray[1];
       realT wRef = tempArray[2];

       if (wRef < 0.99 * m_wZeroStress)
       {
         m_bBarton = (m_wZeroStress - wRef) / stressRef / wRef;
         m_aBarton = m_wZeroStress * (m_wZeroStress - wRef) / stressRef / wRef;
       }
       else
       {
         m_bBarton = 0;
         m_aBarton = 0;
         m_min_aperture = m_wZeroStress;
       }
     }
     else
     {
       m_bBarton = 0;
       m_aBarton = 0;
     }


}


void ParallelPlateProppantFlowSolverImplicit::RegisterFields( PhysicalDomainT& domain )
{
  
  ParallelPlateFlowSolverBase::RegisterFields( domain.m_feFaceManager, domain.m_feEdgeManager );

  const bool plotOldValues = false; // no need to plot old values unless debugging - overwritten at end of timestep.

  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr+oldStr,true,plotOldValues);
  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr+oldStr,true,plotOldValues);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr,true,true);

  domain.m_feFaceManager.AddKeylessDataField<realT>(PS_STR::FaceAreaStr ,true,true);

  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FluidVelocityStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R2Tensor>("ShearStrainRate",true,true);
    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::pressure>();    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::density>(); 
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::mass>(); 
  
  domain.m_feFaceManager.AddKeylessDataField<realT>("massRate",true,true); // face mass rate
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::volume>();
  domain.m_feFaceManager.AddKeylessDataField<realT>("Volume_old",true,plotOldValues);

  domain.m_feFaceManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);


  domain.m_feEdgeManager.AddKeylessDataField<realT>(PermeabilityStr,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(VolumetricFluxStr,true,true);
  
  domain.m_feEdgeManager.AddKeylessDataField<int>("FlowFaceCount",true,true);// debug
  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>(EdgeCenterStr+oldStr,true,plotOldValues);
  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>(EdgeCenterStr,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(EdgeLengthStr+oldStr,true,plotOldValues);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(EdgeLengthStr,true,true);



  domain.m_feEdgeManager.AddKeylessDataField<realT>(ProppantVolumeFractionStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>(ProppantVolumeFractionStr,true,true);

  domain.m_feFaceManager.AddKeylessDataField<realT>(ProppantPackVolumeFractionStr ,true,true); // the volume fraction occupied by the proppant pack (fluid and particles)

  domain.m_feEdgeManager.AddKeylessDataField<realT>(ViscosityStr,true,true);

  domain.m_feEdgeManager.AddKeylessDataField<realT>("DeltaP",true,true);

  domain.m_feFaceManager.AddKeylessDataField<realT>("initialSaturatedTime", true, true);

  // for visualization
  domain.m_feNodeManager.AddKeylessDataField<realT>("nodalVolumeFraction", true, true);
  domain.m_feNodeManager.AddKeylessDataField<realT>("nodalAperture", true, true);

  // erosion

  domain.m_feFaceManager.AddKeylessDataField<realT>("ErodedAperture", true, true);


}

void ParallelPlateProppantFlowSolverImplicit::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
  iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");

  flowFaceType = -1;
  flowEdgeType = -1;

  if( !(m_flowFaceSetName.empty()) )
  {
    m_faceSet = faceManager.GetSet(m_flowFaceSetName);
    m_numFaces = m_faceSet.size();

    // build face-dof map
   /*
    lSet::const_iterator si=m_faceSet.begin();
    for(localIndex i =0; i < m_numFaces; ++i, ++si){
      localIndex f = *si;
      m_faceDofMap[f] = i;
    }
*/
    //iArray1d& ffCount = domain.m_feEdgeManager.GetFieldData<int>("FlowFaceCount"); // debug

  /*
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
    */
    // edge parents to face parents
      m_edgesToFaces.clear();
      lSet::const_iterator fi=m_faceSet.begin();
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
          	if(isMember(faceParent,m_faceSet) && !isMember( faceParent, m_edgesToFaces[egParent] ) ){
              m_edgesToFaces[egParent].push_back(faceParent);
          	}
          }
        }
      }
  }
  else
  {
	//m_faceSet = &(domain.m_feFaceManager.m_Sets["flowFaceSet"]);
    DefineFlowSets(domain);
  }

  rArray1d& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");
  initialSaturatedTime = std::numeric_limits<realT>::max();

  GenerateParallelPlateGeometricQuantities( domain,0,0 );
  OverwriteOldGeometricQuantities(domain);
  GenerateSlurryParameters( domain );

  InitializeDensity( domain);

}


/**
 * Transfer aperture data from the external face manager to the face manager. 
 * 
 * 
 */
void ParallelPlateProppantFlowSolverImplicit:: UpdateAperture(PhysicalDomainT&  domain){
  
  iArray1d& isExternal = domain.m_feFaceManager.GetFieldData<FieldInfo::isExternal>();
  lArray1d& externalFaceIndex = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");
  //rArray1d& external_aperture = domain.m_externalFaces.GetFieldData<realT>("aperture");
  rArray1d& normal_approach = domain.m_externalFaces.GetFieldData<realT>("normalApproach");
  rArray1d& face_aperture = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);
  
  for( lSet::const_iterator si=m_faceSet.begin() ; si!=m_faceSet.end() ; ++si ){
    localIndex kf = *si;
    if(isExternal[kf])
    {
      face_aperture[kf] = std::max(-normal_approach[externalFaceIndex[kf]],0.0);
    }
  }
}

void ParallelPlateProppantFlowSolverImplicit:: SetupSystem (PhysicalDomainT&  domain,
                                                SpatialPartition& partition, const realT& time)
{
  
  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  //if(m_doApertureUpdate) UpdateAperture(domain);
  
  
  // count local dof
  ///////////////////////////////
  
  // local rows
  int n_local_rows = 0;
  
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if( is_ghost[*kf] < 0 )
    {
      ++n_local_rows;
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

        sparsity->InsertGlobalIndices(dofIndex.size(),
                                      &dofIndex.front(),
                                      dofIndex.size(),
                                      &dofIndex.front());
      }
    }
  }
  
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();
  
}

/* Assemble */

void ParallelPlateProppantFlowSolverImplicit :: Assemble (PhysicalDomainT&  domain,
                                                  SpatialPartition& partition ,
                                                  const realT& time,
                                                  const realT& dt)
{
  
  // (re-)init linear system
  matrix   = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*sparsity));
  solution = Teuchos::rcp(new Epetra_FEVector(*row_map));
  rhs      = Teuchos::rcp(new Epetra_FEVector(*row_map));

  // basic face data ( = dof data for our problem)

  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

  iArray1d& edge_is_ghost       = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  iArray1d& face_is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  const Array1dT<R1Tensor>& faceCenters_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr + oldStr );
  const Array1dT<R1Tensor>& faceCenters_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  const rArray1d& faceFluidVolume_old  = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
  const rArray1d& faceFluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  const rArray1d&  faceArea = domain.m_feFaceManager.GetFieldData<realT>( PS_STR::FaceAreaStr );

  const Array1dT<R1Tensor>& edgeCenters_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr + oldStr );
  const Array1dT<R1Tensor>& edgeCenters_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );

  const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
  const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr  );
  
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  //rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);

  rArray1d& edgeMus = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);
  rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);
  rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");

  const Array1dT<R1Tensor>& faceVelocities = domain.m_feFaceManager.GetFieldData<R1Tensor>( FluidVelocityStr );

  rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);

  // leakoff
  rArray1d& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");

  // loop over faces and create an identity matrix
  Epetra_IntSerialDenseVector  faceDofIndex (1);
  Epetra_SerialDenseMatrix     face_matrix  (1,1);
  Epetra_SerialDenseVector     face_rhs (1);
  face_matrix(0,0) = 1;
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {
    if(face_is_ghost[*kf] < 0)
    {
      faceDofIndex(0) = trilinos_index[*kf];
      matrix->SumIntoGlobalValues(faceDofIndex, face_matrix);

      // calculate leakoff
      realT P = faceFluidPressure[*kf];
      // adjust for weight of water - no need
      // realT dz = faceCenters_new[*kf][2];
      // P -= m_proppantData.BuoyancyDeltaP(0.0,dz); // pore fluid pressure
      realT q = m_leakoffModel.CalculateFlux(P,time,initialSaturatedTime[*kf]); // flux out
      realT A = faceArea[*kf];
      realT m = faceFluidMass[*kf];
      //realT rho = std::min( m/(faceFluidVolume_old[*kf]+TINY), 1.5*m_rho_o);  // cap max leakoff density
      realT rho = m/(faceFluidVolume_old[*kf]+TINY);  
      face_rhs(0) = std::max(-2*q*A*rho,-m*0.1); // limit max leakoff to 10% of mass per timestep, factor of 2 accounts for both sides of fracture
      rhs->SumIntoGlobalValues(faceDofIndex,face_rhs);

    }
  }

    
  // loop over edges
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {

	localIndex eg = itr->first;

    int numFaces = itr->second.size();
    m_mu = edgeMus[eg];

    if( edge_is_ghost[eg] < 0 ) // need edge deltaP for ghost edges
    {
     
      Epetra_IntSerialDenseVector  edgeDofIndex (numFaces);
      Epetra_SerialDenseVector     edge_rhs     (numFaces);
      Epetra_SerialDenseMatrix     edge_matrix  (numFaces,numFaces);
     
      if( numFaces == 2)
      {
        localIndex kf = itr->second[0];
        localIndex kfb = itr->second[1];   

        // calculate edge permeability
        realT kappa_old,kappa_new;

        if(m_usePowerlawFluid){
            kappa_old = TwoFacePermeability_PowerLaw(edgeCenters_old,edgeLengths_old,
                                                  faceCenters_old, apertures_old,
                                                  faceVelocities,facePackVfs,
                                                  eg,kf,kfb);
            kappa_new = TwoFacePermeability_PowerLaw(edgeCenters_new,edgeLengths_new,
                                                     faceCenters_new, apertures_new,
                                                     faceVelocities,facePackVfs,
                                                      eg,kf,kfb);
        } else {
            kappa_old = TwoFacePermeability(edgeCenters_old,edgeLengths_old,
                                               faceCenters_old, apertures_old,facePackVfs,
                                               eg,kf,kfb);
            kappa_new = TwoFacePermeability(edgeCenters_new,edgeLengths_new,
                                               faceCenters_new, apertures_new,facePackVfs,
                                               eg,kf,kfb);
        }


        edgePermeabilities[eg] = kappa_new;

      
        // build stiffness matrix and rhs
        if( kappa_new > 0.0 ){

            // stiffness matrix contribution
            /////////////////////////////////
                 /**
            realT volume_old[2];
            volume_old[0] = faceFluidVolume_old[kf];
            volume_old[1] = faceFluidVolume_old[kfb];

            realT volume_new[2];
            volume_new[0] = faceFluidVolume_new[kf];
            volume_new[1] = faceFluidVolume_new[kfb];

            realT mass[2];
            mass[0] = faceFluidMass[kf];
            mass[1] = faceFluidMass[kfb];

            realT rhoAv_old = 0.5*(mass[0]/volume_old[0] + mass[1]/volume_old[1]);
            realT rhoAv_s = 0.5*(mass[0]/volume_new[0] + mass[1]/volume_new[1]);

            realT pressure[2];
            realT pressure_s[2];
            realT dPdM_s[2];

            pressure[0] = faceFluidPressure[kf];
            pressure[1] = faceFluidPressure[kfb];

            pressure_s[0] = P_EOS(mass[0]/volume_new[0],m_bulk_modulus,m_rho_o);
            pressure_s[1] = P_EOS(mass[1]/volume_new[1],m_bulk_modulus,m_rho_o);

            dPdM_s[0] = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_new[0]); // dP/dm = dP/dRho * dRho/dm
            dPdM_s[1] = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_new[1]);

            // Body force
            realT dz = faceCenters_new[kfb][2] - faceCenters_new[kf][2];
            realT DeltaP = m_proppantData.BuoyancyDeltaP(edgeVfs[eg],dz);// body force contribution   -g*dz*((rho_p - rho_l)*phi+rho_l)
            // DeltaP +=  m_proppantData.m_rholG *dz; // weight of fluid
            edgeDeltaP[eg] = DeltaP;

            // Assemble
            ///////////

            edgeDofIndex[0] = trilinos_index[kf];
            edgeDofIndex[1] = trilinos_index[kfb];
            **/
            /**
            // self terms
            edge_matrix(0,0) = m_phi*kappa_new*dt
                *(rhoAv_s*dPdM_s[0] - 0.5*(pressure_s[1]-pressure_s[0])/volume_new[0]);

            edge_matrix(1,1) = m_phi*kappa_new*dt
                *(rhoAv_s*dPdM_s[1] - 0.5*(pressure_s[0]-pressure_s[1])/volume_new[1]);

            // cross terms
            edge_matrix(0,1) = -m_phi*kappa_new*dt*(rhoAv_s*dPdM_s[1] + 0.5*(pressure_s[1]-pressure_s[0])/volume_new[1]);
            edge_matrix(1,0) = -m_phi*kappa_new*dt*(rhoAv_s*dPdM_s[0] + 0.5*(pressure_s[0]-pressure_s[1])/volume_new[0]);

            // rhs
            edge_rhs(0) = (1.0-m_phi) * kappa_old*rhoAv_old*(pressure[1] - pressure[0])
                                  + m_phi * kappa_new*rhoAv_s*(pressure_s[1]- pressure_s[0] )
                                  - kappa_old*DeltaP*rhoAv_old;
            edge_rhs(1) = -edge_rhs(0);

            **/

        	realT volume_old[2];
        	          volume_old[0] = faceFluidVolume_old[kf];
        	          volume_old[1] = faceFluidVolume_old[kfb];

        	          realT volume_new[2];
        	          volume_new[0] = faceFluidVolume_new[kf];
        	          volume_new[1] = faceFluidVolume_new[kfb];

        	          realT deltaV[2];
        	          deltaV[0] = volume_new[0] - volume_old[0];
        	          deltaV[1] = volume_new[1] - volume_old[1];

        	          realT mass[2];
        	          mass[0] = faceFluidMass[kf];
        	          mass[1] = faceFluidMass[kfb];

        	          realT pressure[2];
        	          realT dPdM[2];
        	          realT dPdV[2];
        	          pressure[0] = faceFluidPressure[kf];
        	          pressure[1] = faceFluidPressure[kfb];

        	          /*
        	           * No pressure cap
        	           *

        	          dPdM[0] = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_old[0]); // dP/dm = dP/dRho * dRho/dm
        	          dPdM[1] = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_old[1]);

        	          dPdV[0] = -(mass[0]/volume_old[0])* dPdM[0]; // dP/dV = -(m/V^2) * dP/dRho = -(m/V)*dP/dm
        	          dPdV[1] = -(mass[1]/volume_old[1])* dPdM[1];
                      */


        	         // Pressure cap
        	          realT rho0 = mass[0]/volume_old[0];
        	          realT rho1 = mass[1]/volume_old[1];

        	          //dPdM[0] = dPdRho_EOS(rho0,m_bulk_modulus,m_rho_o,m_pressureCap)*(1.0/volume_old[0]); // dP/dm = dP/dRho * dRho/dm
        	          //dPdM[1] = dPdRho_EOS(rho1,m_bulk_modulus,m_rho_o,m_pressureCap)*(1.0/volume_old[1]);
        	          dPdM[0] = m_fluidEOS->dPdRho(rho0)*(1.0/volume_old[0]); // dP/dm = dP/dRho * dRho/dm
        	          dPdM[1] = m_fluidEOS->dPdRho(rho1)*(1.0/volume_old[1]);


        	          dPdV[0] = -rho0* dPdM[0]; // dP/dV = -(m/V^2) * dP/dRho = -(m/V)*dP/dm
        	          dPdV[1] = -rho1* dPdM[1];


        	          // Body force
        	          realT dz = faceCenters_new[kfb][2] - faceCenters_new[kf][2];
        	          realT DeltaP = m_proppantData.BuoyancyDeltaP(edgeVfs[eg],dz);// body force contribution   -g*dz*((rho_p - rho_l)*phi+rho_l)
        	          // DeltaP +=  m_proppantData.m_rholG *dz; // weight of fluid
        	          edgeDeltaP[eg] = DeltaP;
        	          //realT rhoAv_old = 0.5*( mass[1]/volume_old[1] + mass[0]/volume_old[0]);
        	          pressure[0] += DeltaP/2;
                	  pressure[1] -= DeltaP/2;

        	          // self terms
        	          realT self_0 = m_phi*(kappa_new/volume_new[0])*(pressure[0] + dPdM[0]*mass[0])*dt;
        	          realT self_1 = m_phi*(kappa_new/volume_new[1])*(pressure[1] + dPdM[1]*mass[1])*dt;
        	          edge_matrix(0,0) = self_0;
        	          edge_matrix(1,1) = self_1;

        	          // cross terms
        	          edge_matrix(0,1) = -self_1;
        	          edge_matrix(1,0) = -self_0;

        	          edgeDofIndex[0] = trilinos_index[kf];
        	          edgeDofIndex[1] = trilinos_index[kfb];

        	          // rhs
        	       //    edge_rhs(0) = (1.0-m_phi) * kappa_old*( pressure[1] * mass[1]/volume_old[1] - pressure[0] * mass[0]/volume_old[0])
        	       //                      + m_phi * kappa_new*( (pressure[1] + deltaV[1]*dPdV[1]) * mass[1]/volume_new[1]
        	       //                                          - (pressure[0] + deltaV[0]*dPdV[0]) * mass[0]/volume_new[0] )
        	       //                                          - kappa_old*DeltaP*rhoAv_old;
        	          edge_rhs(0) = (1.0-m_phi) * kappa_old*( pressure[1] * mass[1]/volume_old[1] - pressure[0] * mass[0]/volume_old[0])
        	                  	                            + m_phi * kappa_new*( (pressure[1] + deltaV[1]*dPdV[1]) * mass[1]/volume_new[1]
        	                  	                                                - (pressure[0] + deltaV[0]*dPdV[0]) * mass[0]/volume_new[0] );
        	           edge_rhs(1) = -edge_rhs(0);

          // assemble
          matrix->SumIntoGlobalValues(edgeDofIndex, edge_matrix);
          rhs->SumIntoGlobalValues(edgeDofIndex,edge_rhs);
        }
  
      } else if(numFaces > 2){
    
        lArray1d& faces = itr->second; 
        R1Tensor edgeCenter_old = edgeCenters_old[eg];
        R1Tensor edgeCenter_new = edgeCenters_new[eg];
        //realT w_old = edgeLengths_old[eg];
        //realT w_new = edgeLengths_new[eg];
      
        rArray1d kappas_old(numFaces,0.0);
        rArray1d kappas_new(numFaces,0.0);

        rArray1d volume_old(numFaces,0.0);
        rArray1d volume_new(numFaces,0.0);
        rArray1d deltaV(numFaces,0.0);


        rArray1d pressure(numFaces,0.0);
        rArray1d dPdM(numFaces,0.0);
        rArray1d dPdV(numFaces,0.0);
        rArray1d mass(numFaces,0.0);

        realT kappaSum_old =0.0;
        realT kappaSum_new =0.0;
        for(int i =0; i < numFaces; ++i){
          localIndex kf = faces[i];
          edgeDofIndex[i] = trilinos_index[kf];

          R1Tensor la_old = edgeCenters_old[eg];
          la_old -= faceCenters_old[kf];

          // kappas_old[i] = CalculatePermeability( la_old.L2_Norm(), BoundedAperture(apertures_old[kf]), w_old,m_mu,m_SHP_FCT);
          kappas_old[i] = OneFacePermeability(edgeCenters_old, edgeLengths_old,faceCenters_old, apertures_old,facePackVfs,eg,kf);
          kappaSum_old += kappas_old[i];

          R1Tensor la_new = edgeCenters_new[eg];
          la_new -= faceCenters_new[kf];
          //kappas_new[i] = CalculatePermeability( la_new.L2_Norm(), BoundedAperture(apertures_new[kf]), w_new,m_mu,m_SHP_FCT);
          kappas_new[i] = OneFacePermeability(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new,facePackVfs,eg,kf);
          kappaSum_new += kappas_new[i];

          // volumes
          volume_old[i] =  faceFluidVolume_old[kf];
          volume_new[i] =  faceFluidVolume_new[kf]; 
          deltaV[i] = volume_new[i] - volume_old[i];

          mass[i] = faceFluidMass[kf];
          pressure[i] = faceFluidPressure[kf];

          // Body force
          realT dz = edgeCenter_new[2] - faceCenters_new[kf][2];
          realT DeltaP = m_proppantData.BuoyancyDeltaP(edgeVfs[eg],dz);// body force contribution   -g*dz*((rho_p - rho_l)*phi+rho_l)
          pressure[i] += DeltaP;


          /*
           * No pressure cap
           *

          dPdM[i] = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_old[i]);
          dPdV[i] = -dPdM[i]*(mass[i]/volume_old[i]);
           */
          // Pressure cap
          realT rho = mass[i]/volume_old[i];
          // dPdM[i] = dPdRho_EOS(rho,m_bulk_modulus,m_rho_o,m_pressureCap)*(1.0/volume_old[i]);
          dPdM[i] = m_fluidEOS->dPdRho(rho)*(1.0/volume_old[i]);
          dPdV[i] = -dPdM[i]*rho;
        }

        realT invKappaSum_old = 1.0/(kappaSum_old + TINY);
        realT invKappaSum_new = 1.0/(kappaSum_new + TINY);
      
        // add contributions to stiffness matrices
        for(int i =0; i < numFaces; ++i)
        {
          realT selfCoeff_old = -kappas_old[i]*(1.0 - kappas_old[i]*invKappaSum_old);
          realT selfCoeff_new = -kappas_new[i]*(1.0 - kappas_new[i]*invKappaSum_new);
          // diagonal
          edge_matrix(i,i) =  m_phi*(selfCoeff_new/volume_new[i])
                                   *(pressure[i] + dPdM[i]*mass[i])*dt;
          realT eRHS = (1.0-m_phi) * selfCoeff_old * pressure[i] * mass[i]/volume_old[i]
                           + m_phi * selfCoeff_new * (pressure[i] + deltaV[i]*dPdV[i]) * mass[i]/volume_new[i];
          
          // cross terms
          for(int ii =0; ii < numFaces; ++ii)
          {
            if(i!=ii){
              realT crossCoeff_old = kappas_old[i]*kappas_old[ii]*invKappaSum_old;
              realT crossCoeff_new = kappas_new[i]*kappas_new[ii]*invKappaSum_new;

              edge_matrix(i,ii) =  m_phi*(crossCoeff_new/volume_new[ii])
                                        *(pressure[ii] + dPdM[ii]*mass[ii])*dt;

              eRHS +=  (1.0-m_phi) * crossCoeff_old * pressure[ii]*mass[ii]/volume_old[ii]
                           + m_phi * crossCoeff_new * (pressure[ii] + deltaV[ii]*dPdV[ii])*mass[ii]/volume_new[ii];
            }
          }

          edge_rhs(i) += eRHS;
        }
      
        matrix->SumIntoGlobalValues(edgeDofIndex, edge_matrix);
        rhs->SumIntoGlobalValues(edgeDofIndex,edge_rhs);

      } // numFaces > 2 ?
    } else { // ghost edge

        if( numFaces == 2)
        {
            localIndex kf = itr->second[0];
            localIndex kfb = itr->second[1];

	        realT dz = faceCenters_new[kfb][2] - faceCenters_new[kf][2];
	        realT DeltaP = m_proppantData.BuoyancyDeltaP(edgeVfs[eg],dz);// body force contribution   -g*dz*((rho_p - rho_l)*phi+rho_l)
	        // DeltaP +=  m_proppantData.m_rholG *dz; // weight of fluid
	        edgeDeltaP[eg] = DeltaP;

        }

    }// ghost edge
   
  } // edge loop
  
  
  // boundary conditions
  ///////////////////////

  // Pressure
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantFlowSolverImplicit::PressureBoundaryCondition,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time, dt );

  // Constant flux
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantFlowSolverImplicit::FixedFluxBoundaryCondition,
                                domain, domain.m_feEdgeManager, "FixedFlux", time, dt );
  // Constant net flux
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantFlowSolverImplicit::FixedNetFluxBoundaryCondition,
                                domain, domain.m_feEdgeManager, "FixedNetFlux", time, dt );

  // Calculate fixed pressure field
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantFlowSolverImplicit::FixedFacePressureBoundaryCondition,
                                domain, domain.m_feFaceManager, std::string(Field<FieldInfo::pressure>::Name()), time, dt );


  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  //rhs->Print(std::cout);
  //std::cout << matrix->NormInf() << std::endl;
  //EpetraExt::RowMatrixToMatlabFile("system-matrix.dat",*matrix);
  //exit(0);
}

realT ParallelPlateProppantFlowSolverImplicit::TwoFacePermeability(const Array1dT<R1Tensor>& edgeCenters,
                                                           const rArray1d& edgeLengths,
                                                           const Array1dT<R1Tensor>& faceCenters,
                                                           const rArray1d& apertures,
                                                           const rArray1d& packVfs,
                                                           localIndex eg,localIndex kf, localIndex kfb)
{

  R1Tensor edgeCenter = edgeCenters[eg];
  R1Tensor la, lb;

  la = edgeCenter;
  la -= faceCenters[kf];

  lb = edgeCenter;
  lb -= faceCenters[kfb];

  realT w = edgeLengths[eg];

  realT app = (apertures[kf] < m_max_aperture) ? apertures[kf] : m_max_aperture;
  realT appb = (apertures[kfb] < m_max_aperture) ? apertures[kfb] : m_max_aperture;

  if (app < m_min_aperture)
    app = m_min_aperture;
  if (appb < m_min_aperture)
    appb = m_min_aperture;
  realT lla = la.L2_Norm();
  realT llb = lb.L2_Norm();

  realT kappa = CalculatePermeability(lla,llb, app, appb, w, m_mu, m_SHP_FCT);

  // get proppant pack volume fraction
  realT packVf = std::max(packVfs[kf],packVfs[kfb]);

  if(packVf > 0.0){
	realT kappa_pack = ProppantPackPermeability(0.5*(app+appb), w,lla+llb, m_proppantData.m_maxVf,m_proppantData.m_fluidViscosity);
	kappa = PartiallyFilledEdgePermeability(packVf,kappa,kappa_pack);
  }

  return kappa;
}

realT ParallelPlateProppantFlowSolverImplicit::TwoFacePermeability_PowerLaw(const Array1dT<R1Tensor>& edgeCenters,
                                                           const rArray1d& edgeLengths,
                                                           const Array1dT<R1Tensor>& faceCenters,
                                                           const rArray1d& apertures,
                                                           const Array1dT<R1Tensor>& fluidVelocity,
                                                           const rArray1d& packVfs,
                                                           localIndex eg,localIndex kf, localIndex kfb)
{
	  R1Tensor edgeCenter = edgeCenters[eg];
	  R1Tensor la, lb;

	  la = edgeCenter;
	  la -= faceCenters[kf];

	  lb = edgeCenter;
	  lb -= faceCenters[kfb];

	  realT w = edgeLengths[eg];

	  realT app = (apertures[kf] < m_max_aperture) ? apertures[kf] : m_max_aperture;
	  realT appb = (apertures[kfb] < m_max_aperture) ? apertures[kfb] : m_max_aperture;

	  if (app < m_min_aperture)
	    app = m_min_aperture;
	  if (appb < m_min_aperture)
	    appb = m_min_aperture;
	  realT lla = la.L2_Norm();
	  realT llb = lb.L2_Norm();

	  realT avApp = 0.5*(app+appb);
	  realT qMag = 0.5*( fluidVelocity[kf].L2_Norm() + fluidVelocity[kfb].L2_Norm() + TINY )*avApp;
	  realT kappa = m_fluidModelPtr->CalculatePermeability(lla, llb, app, appb, w, qMag, m_SHP_FCT);

	  realT packVf = std::max(packVfs[kf],packVfs[kfb]);
	  if(packVf > 0.0){
		realT kappa_pack = ProppantPackPermeability(avApp, w,lla+llb, m_proppantData.m_maxVf,m_proppantData.m_fluidViscosity);
		kappa = PartiallyFilledEdgePermeability(packVf,kappa,kappa_pack);
	  }


   return kappa;
}

/*

realT ParallelPlateProppantFlowSolverImplicit::TwoFacePermeability_HerschelBulkley(const Array1dT<R1Tensor>& edgeCenters,
                                                           const rArray1d& edgeLengths,
                                                           const Array1dT<R1Tensor>& faceCenters,
                                                           const rArray1d& apertures,
                                                           const Array1dT<R1Tensor>& fluidVelocity,
                                                           const rArray1d& packVfs,
                                                           localIndex eg,localIndex kf, localIndex kfb)
{
	  R1Tensor edgeCenter = edgeCenters[eg];
	  R1Tensor la, lb;

	  la = edgeCenter;
	  la -= faceCenters[kf];

	  lb = edgeCenter;
	  lb -= faceCenters[kfb];

	  realT w = edgeLengths[eg];

	  realT app = (apertures[kf] < m_max_aperture) ? apertures[kf] : m_max_aperture;
	  realT appb = (apertures[kfb] < m_max_aperture) ? apertures[kfb] : m_max_aperture;

	  if (app < m_min_aperture)
	    app = m_min_aperture;
	  if (appb < m_min_aperture)
	    appb = m_min_aperture;
	  realT lla = la.L2_Norm();
	  realT llb = lb.L2_Norm();

	  realT avApp = 0.5*(app+appb);
	  realT qMag = 0.5*( fluidVelocity[kf].L2_Norm() + fluidVelocity[kfb].L2_Norm() + TINY )*avApp;
	  realT kappa = m_fluidModelPtr->CalculatePermeability(lla, llb, app, appb, w, qMag, m_SHP_FCT);

	  realT packVf = std::max(packVfs[kf],packVfs[kfb]);
	  if(packVf > 0.0){
		realT kappa_pack = ProppantPackPermeability(avApp, w,lla+llb, m_proppantData.m_maxVf,m_proppantData.m_fluidViscosity);
		kappa = PartiallyFilledEdgePermeability(packVf,kappa,kappa_pack);
	  }


   return kappa;
}
*/

realT ParallelPlateProppantFlowSolverImplicit::OneFacePermeability(const Array1dT<R1Tensor>& edgeCenters,
                                                                   const rArray1d& edgeLengths,
                                                                   const Array1dT<R1Tensor>& faceCenters,
                                                                   const rArray1d& apertures,
                                                                   const rArray1d& packVfs,
                                                                   localIndex eg,localIndex kf)
{

  R1Tensor edgeCenter = edgeCenters[eg];
  R1Tensor la, lb;

  la = edgeCenter;
  la -= faceCenters[kf];

  realT w = edgeLengths[eg];

  realT app = (apertures[kf] < m_max_aperture) ? apertures[kf] : m_max_aperture;

  if (app < m_min_aperture)
    app = m_min_aperture;

  realT lla = la.L2_Norm();
  realT kappa = CalculatePermeability(lla, app, w, m_mu, m_SHP_FCT);

  // get proppant pack volume fraction
  realT packVf = packVfs[kf];

  if(packVf > 0.0){
	realT kappa_pack = ProppantPackPermeability(app, w,lla, m_proppantData.m_maxVf,m_proppantData.m_fluidViscosity);
	kappa = PartiallyFilledEdgePermeability(packVf,kappa,kappa_pack);
  }

  return kappa;
}

realT ParallelPlateProppantFlowSolverImplicit::OneFacePermeability_PowerLaw(const Array1dT<R1Tensor>& edgeCenters,
                                                           const rArray1d& edgeLengths,
                                                           const Array1dT<R1Tensor>& faceCenters,
                                                           const rArray1d& apertures,
                                                           const Array1dT<R1Tensor>& fluidVelocity,
                                                           const rArray1d& packVfs,
                                                           localIndex eg,localIndex kf)
{
	  R1Tensor edgeCenter = edgeCenters[eg];
	  R1Tensor la, lb;

	  la = edgeCenter;
	  la -= faceCenters[kf];

	  realT w = edgeLengths[eg];

	  realT app = (apertures[kf] < m_max_aperture) ? apertures[kf] : m_max_aperture;

	  if (app < m_min_aperture)
	    app = m_min_aperture;

	  realT qMag = (fluidVelocity[kf].L2_Norm()+TINY)*app;
	  realT lla = la.L2_Norm();
	  realT kappa = m_fluidModelPtr->CalculatePermeability(lla, app, w, qMag, m_SHP_FCT);

	  realT packVf = packVfs[kf];
	  if(packVf > 0.0){
		realT kappa_pack = ProppantPackPermeability(app, w,lla, m_proppantData.m_maxVf,m_proppantData.m_fluidViscosity);
		kappa = PartiallyFilledEdgePermeability(packVf,kappa,kappa_pack);
	  }


   return kappa;
}

/*
realT ParallelPlateProppantFlowSolverImplicit::OneFacePermeability_HerschelBulkley(const Array1dT<R1Tensor>& edgeCenters,
                                                           const rArray1d& edgeLengths,
                                                           const Array1dT<R1Tensor>& faceCenters,
                                                           const rArray1d& apertures,
                                                           const Array1dT<R1Tensor>& fluidVelocity,
                                                           const rArray1d& packVfs,
                                                           localIndex eg,localIndex kf)
{
	  R1Tensor edgeCenter = edgeCenters[eg];
	  R1Tensor la, lb;

	  la = edgeCenter;
	  la -= faceCenters[kf];

	  realT w = edgeLengths[eg];

	  realT app = (apertures[kf] < m_max_aperture) ? apertures[kf] : m_max_aperture;

	  if (app < m_min_aperture)
	    app = m_min_aperture;

	  realT qMag = (fluidVelocity[kf].L2_Norm()+TINY)*app;
	  realT lla = la.L2_Norm();
	  realT kappa = m_fluidModelPtr->CalculatePermeability(lla, app, w, qMag, m_SHP_FCT);

	  realT packVf = packVfs[kf];
	  if(packVf > 0.0){
		realT kappa_pack = ProppantPackPermeability(app, w,lla, m_proppantData.m_maxVf,m_proppantData.m_fluidViscosity);
		kappa = PartiallyFilledEdgePermeability(packVf,kappa,kappa_pack);
	  }


   return kappa;
}
*/

realT ParallelPlateProppantFlowSolverImplicit:: CalculateFacePermeability(realT app,realT qMag,realT concentration,realT packVf){
	  realT const w = 1.0;
	  realT const l = 1.0;

      m_mu =  m_proppantData.m_slurryModel.GetViscosity(concentration);

      /* fixme - need to set the effective viscosity as a function of concentration for powerlaw fluids
       *  could treat it as a multiplier
      m_fluidModelPtr->SetViscosity(m_mu);
	  realT kappa = m_fluidModelPtr->CalculatePermeability(l, app, w, qMag, m_SHP_FCT);
      */

      realT kappa = CalculatePermeability(l, app, w, m_mu, m_SHP_FCT);  // fixme only implemented for newtonian

	  if(packVf > 0.0){
		realT kappa_pack = ProppantPackPermeability(app, w,l, m_proppantData.m_maxVf,m_proppantData.m_fluidViscosity);
		kappa = PartiallyFilledEdgePermeability(packVf,kappa,kappa_pack);
	  }
	  return kappa;

}

/// Apply a pressure boundary condition to a given set of edges
void ParallelPlateProppantFlowSolverImplicit::PressureBoundaryCondition( PhysicalDomainT& domain,
                                                                 ObjectDataStructureBaseT& object,
                                                                 BoundaryConditionBase* const bc,
                                                                 const lSet& set,
                                                                 const realT time,
                                                                 const realT dt )
{

  ::InflowBoundaryCondition* ifbc  = dynamic_cast< ::InflowBoundaryCondition*> ( bc);
  bool isInflowBC = (ifbc !=0);

  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

  iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  
  const Array1dT<R1Tensor>& faceCenters_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr + oldStr );
  const Array1dT<R1Tensor>& faceCenters_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  const rArray1d& faceFluidVolume_old  = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
  const rArray1d& faceFluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  const Array1dT<R1Tensor>& edgeCenters_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr + oldStr );
  const Array1dT<R1Tensor>& edgeCenters_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );

  const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
  const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );

  const rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  const Array1dT<R1Tensor>& fluidVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  const rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);

  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  //rArray1d& faceVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantVolumeFractionStr);
  rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);
  rArray1d& edgeMus = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);

  rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");
 

  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);
  Epetra_SerialDenseMatrix     face_matrix(1,1);
    
  // loop over edges, find permeabilities and apply boundary conditions.

  lSet::const_iterator eg=set.begin() ;

  for( localIndex i=0; i < set.size() ; ++i, ++eg ){
	realT edgeVf = edgeVfs[*eg];
	m_mu =  m_proppantData.m_slurryModel.GetViscosity(edgeVf);
	edgeMus[*eg] = m_mu;
    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){

      lArray1d& faces = itr->second; 
      
      realT z = edgeCenters_old[*eg][2];

      realT bc_pressure_old = bc->GetValue(domain.m_feEdgeManager,eg,time) + m_proppantData.BuoyancyDeltaP(edgeVfs[*eg],z);
      realT bc_pressure_new = bc->GetValue(domain.m_feEdgeManager,eg,time+dt)+ m_proppantData.BuoyancyDeltaP(edgeVfs[*eg],z);


      // rhoP = rho * Pressure at BC
      //realT rhoP_old = rho_EOS(bc_pressure_old,m_bulk_modulus,m_rho_o,m_pressureCap)*bc_pressure_old;
      //realT rhoP_new = rho_EOS(bc_pressure_new,m_bulk_modulus,m_rho_o,m_pressureCap)*bc_pressure_new;
      realT rhoP_old = m_fluidEOS->density(bc_pressure_old)*bc_pressure_old;
      realT rhoP_new = m_fluidEOS->density(bc_pressure_new)*bc_pressure_new;

      for(size_t ii = 0; ii < faces.size(); ++ii){
    	  /**

          const localIndex kf = faces[ii];
        
    	  // calculate edge permeability for face
    	  R1Tensor la_old = edgeCenters_old[*eg]; la_old -= faceCenters_old[kf];
    	  R1Tensor la_new = edgeCenters_new[*eg]; la_new -= faceCenters_new[kf];

    	  const realT kappa_old = CalculatePermeability( la_old.L2_Norm(), apertures_old[kf], edgeLengths_old[*eg], m_mu, m_SHP_FCT);
    	  const realT kappa_new = CalculatePermeability( la_new.L2_Norm(), apertures_new[kf], edgeLengths_new[*eg], m_mu, m_SHP_FCT);

    	  const realT volume_old = faceFluidVolume_old[kf];
    	  const realT volume_new = faceFluidVolume_new[kf];
    	          //const realT deltaV = volume_new-volume_old;

    	  const realT mass = faceFluidMass[kf];
    	          //const realT pressure = faceFluidPressure[kf];

    	  const realT fc_pressure = faceFluidPressure[kf];
    	  const realT fc_pressure_s = P_EOS(mass/volume_new,m_bulk_modulus,m_rho_o);

    	  const realT rhoAv_old = bc_rho_old;//0.5*(bc_rho_old+mass/volume_old);
    	  const realT rhoAv_s = bc_rho_new; //0.5*(bc_rho_new+mass/volume_new);

    	  //        realT dPdM = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_old);
    	  realT dPdM_s = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_new); // dP/dm = dP/dRho * dRho/dm
    	  //        realT dPdV = -(mass/volume_old)*dPdM;
    	  // Body force
    	  realT dz = la_old[2];
    	  realT DeltaP = m_proppantData.BuoyancyDeltaP(edgeVfs[*eg],dz);// body force contribution   -g*dz*((rho_p - rho_l)*phi+rho_l)
    	  // DeltaP +=  m_proppantData.m_rholG *dz; // weight of fluid

    	  if(ii==0){
    	            edgeDeltaP[*eg] = DeltaP; // will be incorrect at boundary edges with 2 or more faces
    	            edgePermeabilities[*eg] = kappa_new;// only assign first permeability (at junction)
    	  }
    	  **/


    	          /*
    	          // matrix
    	          face_matrix(0,0) = m_phi*(kappa_new/volume_new)*(pressure + dPdM*mass)*m_dt;

    	          // rhs
    	          face_rhs(0) = (1.0-m_phi) * kappa_old*( rhoP_old - pressure * mass/volume_old)
    	                            + m_phi * kappa_new*( rhoP_new - (pressure + deltaV*dPdV) * mass/volume_new)
    	                                    - kappa_old* deltaP*rhoAv_old;
    	          */

    	  /** old - incorrect?
    	    // matrix
    	    face_matrix(0,0) = m_phi*kappa_new*m_dt
    	                             *(rhoAv_s*dPdM_s- 0.5*(bc_pressure_new-fc_pressure_s)/volume_new);

    	    // rhs
    	    face_rhs(0) = (1.0-m_phi) * kappa_old*rhoAv_old*(bc_pressure_old - fc_pressure)
    	                                      + m_phi * kappa_new*rhoAv_s*(bc_pressure_new - fc_pressure_s )
    	                                      - kappa_old*DeltaP*rhoAv_old;
            **/
    	  /**
  	    // matrix
  	    face_matrix(0,0) = m_phi*kappa_new*m_dt
  	                             *(rhoAv_s*dPdM_s- 0.5*(bc_pressure_new-fc_pressure_s)/volume_new);

  	    // rhs
  	    face_rhs(0) = (1.0-m_phi) * kappa_old*rhoAv_old*(bc_pressure_old - fc_pressure)
  	                                      + m_phi * kappa_new*rhoAv_s*(bc_pressure_new - fc_pressure_s )
  	                                      - kappa_old*DeltaP*rhoAv_old;
                **/
    	  ////
    	  const localIndex kf = faces[ii];
    	  R1Tensor la_old = edgeCenters_old[*eg]; la_old -= faceCenters_old[kf];
    	  //R1Tensor la_new = edgeCenters_new[*eg]; la_new -= faceCenters_new[kf];

    	  // calculate edge permeability for face


    	  realT kappa_new,kappa_old;
          if(m_usePowerlawFluid){

            kappa_old =  OneFacePermeability_PowerLaw(edgeCenters_old, edgeLengths_old,faceCenters_old, apertures_old, fluidVelocity,facePackVfs,*eg,kf);
            kappa_new =  OneFacePermeability_PowerLaw(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new, fluidVelocity,facePackVfs,*eg,kf);

          } /*else if(m_useHBPowerlawFluid){
              kappa_old =  OneFacePermeability_HerschelBulkley(edgeCenters_old, edgeLengths_old,faceCenters_old, apertures_old, fluidVelocity,facePackVfs,*eg,kf);
              kappa_new =  OneFacePermeability_HerschelBulkley(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new, fluidVelocity,facePackVfs,*eg,kf);

          }*/
          else {
        	  /*
    	    kappa_old = CalculatePermeability( la_old.L2_Norm(), BoundedAperture(apertures_old[kf]), edgeLengths_old[*eg], m_mu, m_SHP_FCT);
    	    kappa_new = CalculatePermeability( la_new.L2_Norm(), BoundedAperture(apertures_new[kf]), edgeLengths_new[*eg], m_mu, m_SHP_FCT);
    	    */
            kappa_old =  OneFacePermeability(edgeCenters_old, edgeLengths_old,faceCenters_old, apertures_old,facePackVfs,*eg,kf);
            kappa_new =  OneFacePermeability(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new,facePackVfs,*eg,kf);

          }

    	  const realT volume_old = faceFluidVolume_old[kf];
    	  const realT volume_new = faceFluidVolume_new[kf];
    	  const realT deltaV = volume_new-volume_old;

    	  const realT mass = faceFluidMass[kf];
    	  realT pressure = faceFluidPressure[kf];
    	  /*
    	   * No pressure cap


    	  realT dPdM = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_old);
    	  realT dPdV = -( mass/volume_old)*dPdM;

    	   */

          realT rho = mass/volume_old;
    	  //realT dPdM = dPdRho_EOS(rho,m_bulk_modulus,m_rho_o,m_pressureCap)*(1.0/volume_old);
          realT dPdM = m_fluidEOS->dPdRho(rho)*(1.0/volume_old);
    	  realT dPdV = -rho*dPdM;

    	  //if(ii==0) edgePermeabilities[*eg] = kappa_new;  // only assign first permeability (at junction)

    	  // body force (across edge)
    	  realT dz = la_old[2];
    	  realT DeltaP = m_proppantData.BuoyancyDeltaP(edgeVfs[*eg],dz);// body force contribution (across edge)   -g*dz*((rho_p - rho_l)*phi+rho_l)
    	  
    	  pressure += DeltaP;

    	  bool isActiveBC = true;
    	  if (isInflowBC  && pressure + deltaV*dPdV > bc_pressure_new){ // i.e. wants to flow out
    		  isActiveBC = false;
    	  }

    	  if(isActiveBC ) {
    	    if(ii==0){
    	       edgeDeltaP[*eg] = DeltaP; // will be incorrect at boundary edges with 2 or more faces
    	       edgePermeabilities[*eg] = kappa_new;// only assign first permeability (at junction)
    	    }
    	    // matrix
    	    face_matrix(0,0) = m_phi*(kappa_new/volume_new)*(pressure + dPdM*mass)*dt;

    	    // rhs
    	    face_rhs(0) = (1.0-m_phi) * kappa_old*( rhoP_old - pressure * mass/volume_old)
    	                              + m_phi * kappa_new*( rhoP_new - (pressure + deltaV*dPdV) * mass/volume_new);
    	                           //  - kappa_old*DeltaP*rhoP_old;

    	    face_dof(0) = trilinos_index[kf];
    	    matrix->SumIntoGlobalValues(face_dof, face_matrix);
 	        rhs->SumIntoGlobalValues(face_dof, face_rhs);
    	  } else {
    		edgePermeabilities[*eg] = 0;
    	  }
    
      }
    }
  }
}


/// Apply a fixed aperture to a given set of faces
void ParallelPlateProppantFlowSolverImplicit::FixedApertureBoundaryCondition( PhysicalDomainT& domain,
                                                                 ObjectDataStructureBaseT& object,
                                                                 BoundaryConditionBase* const bc,
                                                                 const lSet& set,
                                                                 const realT time,
                                                                 const realT dt )
{


  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

  iArray1d& face_is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  const Array1dT<R1Tensor>& faceCenters_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr + oldStr );
  const Array1dT<R1Tensor>& faceCenters_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  rArray1d& faceFluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& faceFluidAperture  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>( PS_STR::FaceAreaStr  );

  rArray1d& faceProppantVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantVolumeFractionStr);
  rArray1d& faceProppantPackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);

  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");




  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);
  Epetra_SerialDenseMatrix     face_matrix(1,1);

  // loop over faces, and fix aperture.

  lSet::const_iterator fc=set.begin() ;

  for( localIndex i=0; i < set.size() ; ++i, ++fc ){
	localIndex kf = *fc;
    if(  ( face_is_ghost[kf] < 0 ) && flowFaceType[kf] ==0 ){

      realT originalAperture = faceFluidAperture[kf];
      realT newAperture = bc->GetValue(domain.m_feFaceManager,fc,time);


      // need to update proppant volume fractions
      realT proppantPackVf = faceProppantPackVfs[kf];
      realT proppantVf = faceProppantVfs[kf];
	  realT volFRatio = originalAperture/(newAperture+TINY);

  	  if(newAperture < originalAperture){
    	// prevent compression beyond max allowable for proppant pack + proppant concentration
    	realT minInvVFRatio =  proppantPackVf + proppantVf/m_proppantData.m_maxVf + TINY;
    	// don't allow aperture to halve in one timestep.
    	//minInvVFRatio =  std::max(0.9,minInvVFRatio); <- bad idea - seems to mess with pressure distribution under one timestep
    	if( volFRatio  > 1.0/minInvVFRatio ){
    		// minInvVFRatio = 0.5 + 0.5*minInvVFRatio;  // damp out changes in vol - doesn't seem to help much
    		volFRatio = 1.0/minInvVFRatio;
    		newAperture = originalAperture*minInvVFRatio;
    		// apply force ... ?FIXME - should do this when pressure is applied to fracture
    	}
      }

      // update aperture, fluid volume, proppant volume fractions
      if(newAperture != originalAperture){
        faceFluidAperture[kf] = newAperture;
        faceFluidVolume[kf] = newAperture*faceArea[kf];

        faceProppantPackVfs[kf] = proppantPackVf*volFRatio;
        faceProppantVfs[kf] = proppantVf*volFRatio;
      }

    }
  }
}


/// Apply a pressure boundary condition to a given set of faces
void ParallelPlateProppantFlowSolverImplicit::FixedFacePressureBoundaryCondition( PhysicalDomainT& domain,
                                                                 ObjectDataStructureBaseT& object,
                                                                 BoundaryConditionBase* const bc,
                                                                 const lSet& set,
                                                                 const realT time,
                                                                 const realT dt )
{


  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

  iArray1d& face_is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  const Array1dT<R1Tensor>& faceCenters_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr + oldStr );
  const Array1dT<R1Tensor>& faceCenters_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  const rArray1d& faceFluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  const rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();


  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");


  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);
  Epetra_SerialDenseMatrix     face_matrix(1,1);

  // loop over faces, and fix pressure.

  lSet::const_iterator fc=set.begin() ;

  for( localIndex i=0; i < set.size() ; ++i, ++fc ){
	localIndex kf = *fc;
    if(  ( face_is_ghost[kf] < 0 ) && flowFaceType[kf] ==0 ){

      realT bc_pressure_new = bc->GetValue(domain.m_feFaceManager,fc,time);
      //realT bc_rho_new = rho_EOS(bc_pressure_new,m_bulk_modulus,m_rho_o, m_pressureCap ); //rho_EOS(bc_pressure_new,m_bulk_modulus,m_rho_o);
      realT bc_rho_new = m_fluidEOS->density(bc_pressure_new);

      const realT volume_new = faceFluidVolume_new[kf];

      const realT mass = faceFluidMass[kf];
      const realT newMass = bc_rho_new*volume_new;

      realT dmdt = (dt > 0.0)? (newMass-mass)/dt : 0.0;

      // matrix
      face_matrix(0,0) = LARGE;

      // rhs
      face_rhs(0) = LARGE*dmdt;

      face_dof(0) = trilinos_index[kf];

      if(face_dof(0) == 0) {
    	  std::cout << "dmdt: " << dmdt << std::endl;
    	  std::cout << "bc_pressure_new: " << bc_pressure_new << std::endl;
    	  std::cout << "volume_new: " << volume_new << std::endl;

      }
      matrix->SumIntoGlobalValues(face_dof, face_matrix);
 	  rhs->SumIntoGlobalValues(face_dof, face_rhs);

    }
  }
}

/// Apply a FixedFlux boundary condition to a given set of edges
/// BC value sets the velocity of the fluid entering - need to change the name
void ParallelPlateProppantFlowSolverImplicit::FixedFluxBoundaryCondition( PhysicalDomainT& domain,
                                                                 ObjectDataStructureBaseT& object,
                                                                 BoundaryConditionBase* const bc,
                                                                 const lSet& set,
                                                                 const realT time,
                                                                 const realT dt ){

	  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

	  iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

	  //const Array1dT<R1Tensor>& faceCenters_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr + oldStr );
	  const Array1dT<R1Tensor>& faceCenters_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

	  const rArray1d& faceFluidVolume_old  = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
	  //const rArray1d& faceFluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

	  //const Array1dT<R1Tensor>& edgeCenters_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr + oldStr );
	  const Array1dT<R1Tensor>& edgeCenters_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );

	  //const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
	  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

	  //const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
	  const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );

	  const rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
	  //const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

	  const Array1dT<R1Tensor>& fluidVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
	  const rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);

	  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
	  //rArray1d& faceVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantVolumeFractionStr);
	  rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);
	  rArray1d& edgeMus = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);

	  //rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");


	  Epetra_IntSerialDenseVector  face_dof(1);
	  Epetra_SerialDenseVector     face_rhs(1);
	  Epetra_SerialDenseMatrix     face_matrix(1,1);

	  // loop over edges, find permeabilities and apply boundary conditions.

	  lSet::const_iterator eg=set.begin() ;

	  for( localIndex i=0; i < set.size() ; ++i, ++eg ){
		realT edgeVf = edgeVfs[*eg];
		m_mu =  m_proppantData.m_slurryModel.GetViscosity(edgeVf);
		edgeMus[*eg] = m_mu;
	    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
	    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){

	      lArray1d& faces = itr->second;

	      //realT z = edgeCenters_old[*eg][2];

	      realT flux = bc->GetValue(domain.m_feEdgeManager,eg,time);



	      for(size_t ii = 0; ii < faces.size(); ++ii){

	    	  ////
	    	  const localIndex kf = faces[ii];

	    	  // calculate edge permeability for face
	    	  //R1Tensor la_old = edgeCenters_old[*eg]; la_old -= faceCenters_old[kf];
	    	  R1Tensor la_new = edgeCenters_new[*eg]; la_new -= faceCenters_new[kf];

	    	  realT kappa_new;
              if(m_usePowerlawFluid){
    	        kappa_new = OneFacePermeability_PowerLaw(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new, fluidVelocity,facePackVfs,*eg,kf);
              } /*else if(m_useHBPowerlawFluid){
            	  kappa_new = OneFacePermeability_HerschelBulkley(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new, fluidVelocity,facePackVfs,*eg,kf);
              } */else {
    	          kappa_new = OneFacePermeability(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new,facePackVfs,*eg,kf);
              }

	    	  const realT volume_old = faceFluidVolume_old[kf];

	    	  const realT mass = faceFluidMass[kf];
	    	  //realT pressure = faceFluidPressure[kf];

	    	  const realT massFlux = flux*edgeLengths_new[*eg]*BoundedAperture(apertures_new[kf])*mass/volume_old;

   	          if(ii==0) edgePermeabilities[*eg] = kappa_new;  // only assign first permeability (at junction)


	    	  // rhs
	    	  face_rhs(0) = massFlux;//fixme *dt?
	    	  face_dof(0) = trilinos_index[kf];
	    	  //matrix->SumIntoGlobalValues(face_dof, face_matrix);
	 	      rhs->SumIntoGlobalValues(face_dof, face_rhs);


	      }
	    }
	  }

}

/*
 *  Not working - need to adjust velocity update also
 */
/// Apply a FixedNetFlux boundary condition to a given set of edges (constant flux per unit length of edge)
void ParallelPlateProppantFlowSolverImplicit::FixedNetFluxBoundaryConditionB( PhysicalDomainT& domain,
                                                                 ObjectDataStructureBaseT& object,
                                                                 BoundaryConditionBase* const bc,
                                                                 const lSet& set,
                                                                 const realT time,
                                                                 const realT dt ){

	  m_cycleNumber++;// NQR

	  int rank, size ;
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &size);

	  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

	  iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

	  //const Array1dT<R1Tensor>& faceCenters_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr + oldStr );
	  const Array1dT<R1Tensor>& faceCenters_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

	  const rArray1d& faceFluidVolume_old  = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
	  //const rArray1d& faceFluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

	  //const Array1dT<R1Tensor>& edgeCenters_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr + oldStr );
	  const Array1dT<R1Tensor>& edgeCenters_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );

	  //const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
	  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

	  //const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
	  const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );

	  const rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
	  const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

	  const Array1dT<R1Tensor>& fluidVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
	  const rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);

	  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
	  //rArray1d& faceVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantVolumeFractionStr);
	  rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);
	  rArray1d& edgeMus = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);

	  rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");


	  Epetra_IntSerialDenseVector  face_dof(1);
	  Epetra_SerialDenseVector     face_rhs(1);
	  Epetra_SerialDenseMatrix     face_matrix(1,1);

	  // loop over edges, sum permeabilities/ pressures to find injection pressure.

	  lSet::const_iterator eg=set.begin() ;
      realT qTotal = bc->GetValue(domain.m_feEdgeManager,eg,time);
      int nInlets = set.size();
      realT sumKP = 0.0;
      realT sumK = 0.0;

      nInlets = 0;
      for(localIndex i=0; i < set.size() ; ++i, ++eg){
    	realT edgeVf = edgeVfs[*eg];
		m_mu =  m_proppantData.m_slurryModel.GetViscosity(edgeVf);
		edgeMus[*eg] = m_mu;
	    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
	    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){

	      lArray1d& faces = itr->second;

	      for(size_t ii = 0; ii < faces.size(); ++ii){

	    	  ////
	    	  const localIndex kf = faces[ii];

	    	  // calculate edge permeability for face
	    	  //R1Tensor la_old = edgeCenters_old[*eg]; la_old -= faceCenters_old[kf];
	    	  R1Tensor la_new = edgeCenters_new[*eg]; la_new -= faceCenters_new[kf];

	    	  realT kappa_new; //,kappa_old;
              if(m_usePowerlawFluid){
  	            kappa_new = OneFacePermeability_PowerLaw(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new, fluidVelocity,facePackVfs,*eg,kf);
              } else {
  	            kappa_new = OneFacePermeability(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new,facePackVfs,*eg,kf);
              }

	    	  realT pressure = faceFluidPressure[kf];

	    	  sumK +=kappa_new;
	    	  sumKP += kappa_new*pressure;
	    	  nInlets++;
	      }
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

      if (m_cycleNumber%100 == 0) MPI_Allreduce(&myRank0, &rank0, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD);

      sumKP += qTotal;
      m_fixedFlux_pInj = sumKP / (sumK+1e-64);

      if (myNInlets >=1)
      {

        // No pressure cap
        //realT rhoInj =  rho_EOS(pInj,m_bulk_modulus,m_rho_o );

        //realT rhoInj = m_rho_o; // rho_EOS(pInj,m_bulk_modulus,m_rho_o, m_pressureCap ); using real rho seems to lead to instabilities
    	realT rhoInj = m_fluidEOS->rho_o();
  	    // loop over edges, find permeabilities and apply boundary conditions.
	    eg=set.begin() ;

	    for( localIndex i=0; i < set.size() ; ++i, ++eg ){
 		  realT edgeVf = edgeVfs[*eg];
		  m_mu =  m_proppantData.m_slurryModel.GetViscosity(edgeVf);
		  edgeMus[*eg] = m_mu;
	      std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
	      if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){

	        lArray1d& faces = itr->second;

	        for(size_t ii = 0; ii < faces.size(); ++ii){

	    	  ////
	    	  const localIndex kf = faces[ii];

	    	  // calculate edge permeability for face
	    	  //R1Tensor la_old = edgeCenters_old[*eg]; la_old -= faceCenters_old[kf];
	    	  R1Tensor la_new = edgeCenters_new[*eg]; la_new -= faceCenters_new[kf];

	    	  realT kappa_new; //,kappa_old;
              if(m_usePowerlawFluid){
    	        kappa_new = OneFacePermeability_PowerLaw(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new, fluidVelocity,facePackVfs,*eg,kf);
              } else {
            	kappa_new = OneFacePermeability(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new,facePackVfs,*eg,kf);
              }

	    	  const realT volume_old = faceFluidVolume_old[kf];

	    	  const realT mass = faceFluidMass[kf];
	    	  realT pressure = faceFluidPressure[kf];

	    	  const realT volFlux = kappa_new*(m_fixedFlux_pInj-pressure);
	    	  const realT massFlux = volFlux*rhoInj; // mass/volume_old;

   	          if(ii==0){
   	        	  edgePermeabilities[*eg] = kappa_new;  // only assign first permeability (at junction)
   	          }

	    	  // rhs
	    	  face_rhs(0) = massFlux;
	    	  face_dof(0) = trilinos_index[kf];
	    	  //matrix->SumIntoGlobalValues(face_dof, face_matrix);
	 	      rhs->SumIntoGlobalValues(face_dof, face_rhs);

	        }// loop over faces
	      } //not ghost & edge in flow set
	    }//loop over edges

        if (m_cycleNumber%100 == 0 && myRank0 == rank0)
        {
          std::cout <<"Fixed net flux rate BC applied: t=" << time << ", P_inj=" << m_fixedFlux_pInj <<" rank: " <<rank << std::endl;
        }
	  }// ninlets > 0

}

/***/


/*
 *  Not net flux boundary condition - quick hack to try to get working like fixed flux BC - need to implement a distributed flux (see above) but not clear why that is
 *  not working - poss because velocity update not implemented correctly
 *
 *   actually fixed flux per unit length (regardless of aperture
 */

void ParallelPlateProppantFlowSolverImplicit::FixedNetFluxBoundaryCondition( PhysicalDomainT& domain,
                                                                 ObjectDataStructureBaseT& object,
                                                                 BoundaryConditionBase* const bc,
                                                                 const lSet& set,
                                                                 const realT time,
                                                                 const realT dt ){

	  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);

	  iArray1d& edge_is_ghost = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

	  //const Array1dT<R1Tensor>& faceCenters_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr + oldStr );
	  const Array1dT<R1Tensor>& faceCenters_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

	  const rArray1d& faceFluidVolume_old  = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
	  //const rArray1d& faceFluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

	  //const Array1dT<R1Tensor>& edgeCenters_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr + oldStr );
	  const Array1dT<R1Tensor>& edgeCenters_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );

	  //const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
	  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

	  //const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
	  // const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );

	  const rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
	  //const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

	  //const Array1dT<R1Tensor>& fluidVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
	  //const rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);

	  //rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
	  //rArray1d& faceVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantVolumeFractionStr);
	  rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);
	  rArray1d& edgeMus = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);

	  //rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");


	  Epetra_IntSerialDenseVector  face_dof(1);
	  Epetra_SerialDenseVector     face_rhs(1);
	  Epetra_SerialDenseMatrix     face_matrix(1,1);

	  realT rho_o = m_fluidEOS->rho_o();

	  // loop over edges, find permeabilities and apply boundary conditions.

	  lSet::const_iterator eg=set.begin() ;

	  for( localIndex i=0; i < set.size() ; ++i, ++eg ){
		realT edgeVf = edgeVfs[*eg];
		m_mu =  m_proppantData.m_slurryModel.GetViscosity(edgeVf);
		edgeMus[*eg] = m_mu;
	    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
	    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){

	      lArray1d& faces = itr->second;

	      //realT z = edgeCenters_old[*eg][2];

	      realT flux = bc->GetValue(domain.m_feEdgeManager,eg,time); // fixed flux per length



	      for(size_t ii = 0; ii < faces.size(); ++ii){

	    	  ////
	    	  const localIndex kf = faces[ii];

	    	  // calculate edge permeability for face
	    	  //R1Tensor la_old = edgeCenters_old[*eg]; la_old -= faceCenters_old[kf];
	    	  R1Tensor la_new = edgeCenters_new[*eg]; la_new -= faceCenters_new[kf];

	    	  /*
	    	  realT kappa_new;
              if(m_usePowerlawFluid){
    	        kappa_new = OneFacePermeability_PowerLaw(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new, fluidVelocity,facePackVfs,*eg,kf);
              } else {
    	            //kappa_new = CalculatePermeability( la_new.L2_Norm(), BoundedAperture(apertures_new[kf]), edgeLengths_new[*eg], m_mu, m_SHP_FCT);
            	kappa_new = OneFacePermeability(edgeCenters_new, edgeLengths_new,faceCenters_new, apertures_new,facePackVfs,*eg,kf);
              }

   	          if(ii==0) edgePermeabilities[*eg] = kappa_new;  // only assign first permeability (at junction)
              */

	    	  const realT volume_old = faceFluidVolume_old[kf];

	    	  const realT mass = faceFluidMass[kf];
	    	  //realT pressure = faceFluidPressure[kf];

	    	  //There is a bit of a problem here - if we use the actual density
	    	  realT massFlux;
	    	  const bool useActualDensity = false;
	    	  if(useActualDensity){
	    	    massFlux = flux*edgeLengths_new[*eg]*mass/volume_old;
	    	  } else{
	    		//massFlux = flux*edgeLengths_new[*eg]*m_rho_o;
	    		massFlux = flux*edgeLengths_new[*eg]*rho_o;
	    	  }

	    	  // rhs
	    	  face_rhs(0) = massFlux;
	    	  face_dof(0) = trilinos_index[kf];
	    	  //matrix->SumIntoGlobalValues(face_dof, face_matrix);
	 	      rhs->SumIntoGlobalValues(face_dof, face_rhs);


	      }
	    }
	  }

}




/* Solve */

void ParallelPlateProppantFlowSolverImplicit:: Solve (PhysicalDomainT&  domain,
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
          solver.SetAztecOption(AZ_conv,AZ_rhs);
          if( !m_verboseFlag ) solver.SetAztecOption(AZ_output,AZ_none);
          
          
  solver.Iterate(m_numerics.m_maxIters,m_numerics.m_tol);
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
  partition.SynchronizeFields(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);
}


void ParallelPlateProppantFlowSolverImplicit::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(FaceCenterStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("massRate"); // needed?
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(ApertureStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::pressure>::Name());// needed?
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_TrilinosIndexStr);
  syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(VolumetricFluxStr);
  syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(PermeabilityStr);
  //syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back("DeltaP");

  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("initialSaturatedTime");


  syncedFieldsB.clear();
  syncedFieldsB[PhysicalDomainT::FiniteElementFaceManager].push_back(FluidVelocityStr);

  partition.SetBufferSizes(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);
  partition.SetBufferSizes(syncedFieldsB, CommRegistry::steadyStateParallelPlateFlowSolverB);
}


double ParallelPlateProppantFlowSolverImplicit::TimeStep( const realT& time,
                                                        const realT& dt,
                                                        const int cycleNumber,
                                                        PhysicalDomainT& domain,
                                                        const sArray1d& namesOfSolverRegions ,
                                                        SpatialPartition& partition,
                                                        FractunatorBase* const fractunator )
{

  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();

  DefineFlowSets( domain );

  GenerateParallelPlateGeometricQuantities( domain,time,dt );
  GenerateSlurryParameters( domain );

  CalculateMassRate( domain, partition, time,dt );

  UpdateEOS( domain,time,dt );

  UpdateFlux( time, dt, domain, partition);

  OverwriteOldGeometricQuantities(domain);

  return dt;
}

// this will only update every timestep if the face set is left blank
// may want to add a flag to update it or not instead...
// or update after mesh has been altered
// note that face set may be updated by another solver, but nodeset may not (nor may edges to faces)
void ParallelPlateProppantFlowSolverImplicit::DefineFlowSets( PhysicalDomainT& domain )
{

  if( m_flowFaceSetName.empty() )
  {
  FaceManagerT& faceManager = domain.m_feFaceManager;
//  EdgeManagerT& edgeManager = domain.m_feEdgeManager;

  const iArray1d& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");
 // const iArray1d& flowEdgeType = edgeManager.GetFieldData<int>("flowEdgeType");
 //  const Array1dT<lSet>& edgeToFlowFaces = edgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");


  m_faceSet.clear();
  //m_nodeSet->clear();
  m_edgesToFaces.clear();

  for( localIndex kf=0 ; kf<faceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 0 && domain.m_feFaceManager.IsParent(kf) )
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

  }

}


void ParallelPlateProppantFlowSolverImplicit::GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,
		                                                                        realT time,
                                                                                realT dt )
{

  Array1dT<R1Tensor>& faceCenter_old  = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr + oldStr );
  rArray1d&           fluidVolume_old = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
  rArray1d&           aperture_old    = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
  rArray1d&           erodedApertures = domain.m_feFaceManager.GetFieldData<realT>( "ErodedAperture" );
  rArray1d&           mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();


  Array1dT<R1Tensor>& faceCenter_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>( PS_STR::FaceAreaStr  );
  rArray1d& fluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& aperture_new = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);

  rArray1d& faceProppantVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantVolumeFractionStr);
  rArray1d& faceProppantPackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);

  realT originalAperture,newAperture; // used to update the aperture ("aperture_new") within this timestep

  const rArray1d* effectiveStressN;
  if(m_doApertureUpdate && m_bBarton != 0.0){
	  effectiveStressN = &domain.m_feFaceManager.GetFieldData<realT>("effectiveStressN");
  }

  // update face quantities
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf ) {

    domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *kf, faceCenter_new[*kf]);
    if( isZero(faceCenter_old[*kf].L2_Norm()) )
    {
      faceCenter_old[*kf] = faceCenter_new[*kf];
    }


    if(m_doApertureUpdate){

    	originalAperture = aperture_new[*kf];
    	newAperture = originalAperture;

    	//domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *kf, faceCenter[kf] );
    	R1Tensor gap;
    	R1Tensor N;

    	localIndex numChildren = domain.m_feFaceManager.m_childIndices[*kf].size();
    	if (numChildren <= 1)
    	{
    		N = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *kf );
    	}
    	else
    	{
    		N = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, domain.m_feFaceManager.m_childIndices[*kf][0] );
    	}

    	gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, *kf );

    	if( numChildren > 0){ // don't set aperture of unsplit faces
    		newAperture = Dot(gap,N) ;
    	}

    	if(m_doErosion) newAperture += erodedApertures[*kf];

    	if (newAperture < 0.0)
    	{
    		if (m_bBarton != 0.0)
    		{
    			newAperture = m_wZeroStress - m_aBarton * (*effectiveStressN)[*kf] / (1 + m_bBarton * (*effectiveStressN)[*kf]);
    		}
    		else
    		{
    			newAperture = m_min_aperture;
    		}
    	}
    	else
    	{
    		if( numChildren > 0){ // don't set aperture of unsplit faces
    			newAperture += m_wZeroStress;
    		}
    	}

    }



    if(m_doApertureUpdate){

      // damp aperture change
      newAperture = originalAperture + (newAperture - originalAperture)*dt/m_apertureRelaxationTime;


      // need to update proppant volume fractions
      realT proppantPackVf = faceProppantPackVfs[*kf];
      realT proppantVf = faceProppantVfs[*kf];
	  realT volFRatio = originalAperture/(newAperture+TINY);

  	  if(newAperture < originalAperture){
    	// prevent compression beyond max allowable for proppant pack + proppant concentration
    	realT minInvVFRatio =  proppantPackVf + proppantVf/m_proppantData.m_maxVf + TINY;
    	// don't allow aperture to halve in one timestep.
    	//minInvVFRatio =  std::max(0.9,minInvVFRatio); <- bad idea - seems to mess with pressure distribution under one timestep
    	if( volFRatio  > 1.0/minInvVFRatio ){
    		// minInvVFRatio = 0.5 + 0.5*minInvVFRatio;  // damp out changes in vol - doesn't seem to help much
    		volFRatio = 1.0/minInvVFRatio;
    		newAperture = originalAperture*minInvVFRatio;
    		// apply force ... ?FIXME - should do this when pressure is applied to fracture
    	}
      }

      faceProppantPackVfs[*kf] = proppantPackVf*volFRatio;
      faceProppantVfs[*kf] = proppantVf*volFRatio;

      // update aperture
      if(newAperture != originalAperture) aperture_new[*kf] = newAperture;
    }


    realT A = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, *kf );
    faceArea[*kf] = A;
    fluidVolume_new[*kf] = BoundedAperture(newAperture) * A;


    if( isZero(aperture_old[*kf]) )
    {
      aperture_old[*kf] = newAperture;
    }

    if( isZero(fluidVolume_old[*kf]) )
    {
      fluidVolume_old[*kf] = fluidVolume_new[*kf];
    }

    if( isZero( mass[*kf] ) )
    {
      //mass[*kf] = fluidVolume_new[*kf] * this->m_rho_o;
      mass[*kf] = fluidVolume_new[*kf] * m_fluidEOS->rho_o();
    }

  }

  // apply fixed aperture boundary condition
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantFlowSolverImplicit::FixedApertureBoundaryCondition,
                                domain, domain.m_feFaceManager, "Aperture", time, dt );


  // update edge properties
  iArray1d& edge_is_ghost            = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  Array1dT<R1Tensor>& edgeCenter_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  Array1dT<realT>& edgeLength_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  Array1dT<R1Tensor>& edgeCenter_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr + oldStr );
  Array1dT<realT>& edgeLength_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );

  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0 )
    {
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter_new[eg] );
      edgeLength_new[eg] = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, eg);

      if( isZero(edgeCenter_old[eg].L2_Norm()) )
      {
        edgeCenter_old[eg] = edgeCenter_new[eg];
      }

      if( isZero(edgeLength_old[eg]) )
      {
        edgeLength_old[eg] = edgeLength_new[eg];
      }

    }
  }
}


void ParallelPlateProppantFlowSolverImplicit::GenerateSlurryParameters( PhysicalDomainT& domain ){

  //const iArray1d& flowEdgeType = domain.m_feEdgeManager.GetFieldData<int>("flowEdgeType");
  const rArray1d& faceVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantVolumeFractionStr);
  rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);
  rArray1d& edgeMus = domain.m_feEdgeManager.GetFieldData<realT>(ViscosityStr);

  //const Array1dT<lSet>& edgeToFlowFaces = domain.m_feEdgeManager.GetUnorderedVariableOneToManyMap("edgeToFlowFaces");

  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
   for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
   {

	  localIndex ke = itr->first;
      const unsigned int numFlowFaces = itr->second.size();

      lArray1d& faces = itr->second;


      if(numFlowFaces > 1){ // i.e. not a boundary edge
          realT edgeVf(0.0);
          realT maxVf(1e-10);

          for( localIndex kf=0 ; kf<numFlowFaces ; ++kf)
          {
        	 realT faceVf = faceVfs[faces[kf]];
             edgeVf += faceVf;
             if(maxVf < faceVf) maxVf = faceVf;
          }
           edgeVf /= numFlowFaces;
           edgeVf =  std::min(edgeVf,m_proppantData.m_maxVf);

           edgeVfs[ke] = edgeVf;
           edgeMus[ke] =  m_proppantData.m_slurryModel.GetViscosity(maxVf);
      }
  }
}


void ParallelPlateProppantFlowSolverImplicit::OverwriteOldGeometricQuantities( PhysicalDomainT& domain){

  Array1dT<R1Tensor>& faceCenter_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  rArray1d& fluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& aperture_new = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);

  Array1dT<R1Tensor>& faceCenter_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+oldStr );
  rArray1d& fluidVolume_old  = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
  rArray1d& aperture_old = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr+oldStr);

  // update face quantities
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf ) {
    faceCenter_old[*kf] = faceCenter_new[*kf];
    aperture_old[*kf] = aperture_new[*kf];
    fluidVolume_old[*kf] = fluidVolume_new[*kf];
  }

  // update edge properties
  iArray1d& edge_is_ghost            = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  Array1dT<R1Tensor>& edgeCenter_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  Array1dT<realT>& edgeLength_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  Array1dT<R1Tensor>& edgeCenter_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr +oldStr );
  Array1dT<realT>& edgeLength_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr +oldStr);

  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    if( edge_is_ghost[eg] < 0 )
    {
      edgeCenter_old[eg] = edgeCenter_new[eg];
      edgeLength_old[eg] = edgeLength_new[eg];
    }
  }
}

 

void ParallelPlateProppantFlowSolverImplicit::CalculateMassRate( PhysicalDomainT& domain,
                                                         SpatialPartition& partition, realT time, realT dt )
{

  SetupSystem (domain,partition, time);
  Assemble    (domain,partition, time, dt);
  Solve       (domain,partition);
  
}

// set initial fluid density based on known pressure;
void ParallelPlateProppantFlowSolverImplicit::InitializeDensity( PhysicalDomainT& domain)
{
  rArray1d& mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& density  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const rArray1d& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  for( lSet::const_iterator fc=m_faceSet.begin() ; fc!=m_faceSet.end() ; ++fc ) {

	/* No pressure cap
    density[*fc] = rho_EOS(pressure[*fc],m_bulk_modulus,m_rho_o );
    */

    //density[*fc] = rho_EOS(pressure[*fc],m_bulk_modulus,m_rho_o, m_pressureCap );
	density[*fc] = m_fluidEOS->density(pressure[*fc] );
    mass[*fc] = density[*fc]*fluidVolume[*fc];
  }
}


void ParallelPlateProppantFlowSolverImplicit::UpdateEOS( PhysicalDomainT& domain, const realT time,const realT dt )
{
  rArray1d& mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& density  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  const rArray1d& massRate = domain.m_feFaceManager.GetFieldData<realT>("massRate");

  rArray1d* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");

  const bool updateMass = true;

    for( lSet::const_iterator fc=m_faceSet.begin() ; fc!=m_faceSet.end() ; ++fc )
    {
      if( updateMass )
      {
        mass[*fc] += massRate[*fc] *dt;
      }
      density[*fc] = mass[*fc] / (fluidVolume[*fc]+1e-64);

      /*
       * No Pressure cap
        realT P = P_EOS(density[*fc],m_bulk_modulus,m_rho_o );
       */

      //realT P = P_EOS(density[*fc],m_bulk_modulus,m_rho_o,m_pressureCap);

      realT P = m_fluidEOS->pressure(density[*fc]);

      if( P<0.0 )
        P = 0.0;
      pressure[*fc] = P;
      // propagate pressure to children
      lArray1d& childFaces = domain.m_feFaceManager.m_childIndices[*fc];
      for(unsigned i =0; i < childFaces.size(); ++i){
        pressure[childFaces[i]] = P;
      }


      //Simple Carter's leakoff model; initializing face
      if (pressure[*fc] > 0.0 && (*initialSaturatedTime)[*fc] == std::numeric_limits<realT>::max())
      {
            (*initialSaturatedTime)[*fc] = time - dt/2;
      }

    }

}


void ParallelPlateProppantFlowSolverImplicit::UpdateFlux( const realT time, const realT dt, PhysicalDomainT& domain,
                                                          SpatialPartition& partition)
{
  
  const rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  Array1dT<R1Tensor>& faceVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  const rArray1d& face_aperture = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);
  volFlux = 0.0;
  faceVelocity = 0.0;
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& faceVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  Array1dT<R2Tensor>& shearStrainRate = domain.m_feFaceManager.GetFieldData<R2Tensor>("ShearStrainRate");
  
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  

  rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");

  R1Tensor la, lb;

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

      realT DeltaP = edgeDeltaP[eg]; // body force contribution   -g*dz*((rho_p - rho_l)*phi+rho_l)
  
      // will be incorrect at junctions of 3 or more faces
      volFlux[eg] =  edgePermeabilities[eg]*(Pa-Pb+DeltaP); // flux A -> B
      
      R1Tensor vecFlux;

      vecFlux = faceCenters[kfb];
      vecFlux -= faceCenters[kfa];

      vecFlux.Normalize();
      vecFlux *= volFlux[eg];

      
    }  
  }
  
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantFlowSolverImplicit::PressureBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time );

  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantFlowSolverImplicit::FixedFluxBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, "FixedFlux", time );

  // use same velocity update for fixed net flux
  ApplyBoundaryCondition<realT>(this, &ParallelPlateProppantFlowSolverImplicit::FixedNetFluxBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, "FixedNetFlux", time );

  partition.SynchronizeFields(syncedFieldsB, CommRegistry::steadyStateParallelPlateFlowSolverB);


  // convert flux volume per element on edges to flux per cross-sectional area (velocity) on faces
  for( lSet::const_iterator kf=m_faceSet.begin() ; kf!=m_faceSet.end() ; ++kf )
  {

	  realT app = face_aperture[*kf];
	  if (app > 0.0){

		  R1Tensor faceNorm = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *kf );

		  //realT weight = 0.0;
		  localIndex numEdges = domain.m_feFaceManager.m_toEdgesRelation[*kf].size();

		  // new flux calculation
		  R1TensorT<2> NQ,q;
		  R2SymTensorT<2> NN;
		  R2SymTensorT<3> eij; // strain rate estimate

		  // choose plane for calculation (xi,yi - plane for calculating least squares minimization, zi - out of plane coordinate)
		  int xi = 0;
		  int yi = 1;
		  int zi = 2;

		  if(fabs(faceNorm[zi]) < fabs(faceNorm[xi]) ){
			  int temp = xi;
			  xi = zi; xi = temp;
		  }
		  if(fabs(faceNorm[zi]) < fabs(faceNorm[yi]) ){
			  int temp = yi;
			  yi = zi; zi = temp;
		  }

		  //std::cout << " face calc " << std::endl << std::endl;

		  for(size_t a =0; a < numEdges; ++a){
			  localIndex eg = domain.m_feFaceManager.m_toEdgesRelation[*kf][a];
			  eg = domain.m_feEdgeManager.GetParentIndex(eg);

			  R1Tensor edgeCenter;
			  domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter );

			  R1Tensor dT; domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager,eg,dT); // edge vector
			  dT[zi] = 0; // project dT onto xi,yi plane
			  realT l = dT.L2_Norm();

			  R1Tensor dN;


			  // Calculate flux using net flux divided by face area:
			  // Find least squares solution to q_{j} n^{a}_{j} = Q^{a}/l for q_{j}
			  // \sum n^{a}_{i}n^{a}_{j}  q_{j} = \sum Q^{a}n^{a}_{i}/l
			  // q_{j} = [\sum n^{a}_{i}n^{a}_{j} ]^{-1} * \sum Q^{a}n^{a}_{i}/l

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

			  // use dN to calculate in-plane shear strain contribution
			  dN.Normalize();
        eij += (l*Q)*DyadicProduct(dN);

			  // return to velocity calculation
			  dN[zi] = 0; // project dN onto xi,yi plane
			  dN.Normalize();

			  NN(0,0) += dN[xi]*dN[xi]+TINY;
			  NN(0,1) += dN[xi]*dN[yi];
			  NN(1,1) += dN[yi]*dN[yi]+TINY; // TINY -> prevent singular inverse

			  NQ(0) += Q*dN[xi];
			  NQ(1) += Q*dN[yi];

		  }

		  // invert NN
		  NN.Inverse();

		  // calculate flux
		  q.AijBj(NN,NQ);

		  // out of plane component (flux must lie in plane)
		  realT qz = (-q[0]*faceNorm[xi] - q[1]*faceNorm[yi])/ faceNorm[zi];

		  // Project back to xyz from face coordinates
		  faceVelocity[*kf][xi] = q[0]/app;
		  faceVelocity[*kf][yi] = q[1]/app;
		  faceVelocity[*kf][zi] = qz/app;

		  // Shear strain rate tensor
		  if(faceVolume[*kf] > 0.0){
			  realT costheta = sqrt(1-faceNorm[zi]*faceNorm[zi]);
			  eij *= costheta/faceVolume[*kf]; // divide by (projected) area to convert to average gradient, aperture to convert flux/unit length to velocity
			  realT eii = (1.0/3.0)*( eij(0,0)+eij(1,1)+eij(2,2) ); // not clear if this is correct or 2d version should be used instead
			  eij(0,0) -= eii;
			  eij(1,1) -= eii;
			  eij(2,2) -= eii;

			  shearStrainRate[*kf] = eij;

		  } else {
			  shearStrainRate[*kf] = 0.0;
		  }

	  } else {
		  // aperture <= 0
		  shearStrainRate[*kf] = 0.0;
		  faceVelocity[*kf] = 0.0;

	  }

  }

  partition.SynchronizeFields(syncedFieldsB, CommRegistry::steadyStateParallelPlateFlowSolverB);
    
}


/// Update the velocity fluxes on the boundary faces
void ParallelPlateProppantFlowSolverImplicit::
       PressureBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                 ObjectDataStructureBaseT& object ,
                                                 BoundaryConditionBase* bc, const lSet& set, realT time){

  //::InflowBoundaryCondition* ifbc  = dynamic_cast< ::InflowBoundaryCondition*> ( bc);
  //bool isInflowBC = (ifbc !=0);
  
  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  Array1dT<R1Tensor>& faceVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
    
  iArray1d& is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");
  
  for( lSet::const_iterator eg=set.begin() ; eg!=set.end() ; ++eg ) {
     
    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if(itr != m_edgesToFaces.end() ){
      
      lArray1d& faces = itr->second; 
      
      R1Tensor edgeCenter;
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *eg , edgeCenter );
      const realT z = edgeCenter[2];
      
      // nb only uses first face to calculate volume flux
      // flux will be incorrect for other edges if at junction.
      {            
        localIndex fc = faces[0];
                
        const realT Pa = pressures[fc];
        const realT Pb = bc->GetValue(domain.m_feEdgeManager,eg,time) + m_proppantData.BuoyancyDeltaP(edgeVfs[*eg],z);
        const realT DeltaP = edgeDeltaP[*eg];
      
        volFlux[*eg] =  edgePermeabilities[*eg]*(Pa-Pb+DeltaP); // flux out of a into b

        R1Tensor vecFlux;

        vecFlux = edgeCenter;
        vecFlux -= faceCenters[fc];

        vecFlux.Normalize();
        vecFlux *= volFlux[*eg];

        if(is_ghost[fc] < 0) faceVelocity[fc] += vecFlux;

      }
    }  
  }
}


/// Update the velocity fluxes on the boundary faces
void ParallelPlateProppantFlowSolverImplicit::
    FixedFluxBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                 ObjectDataStructureBaseT& object ,
                                                 BoundaryConditionBase* bc, const lSet& set, realT time){


  //rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  //rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);

  //const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  Array1dT<R1Tensor>& faceVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);

  //Array1dT<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  iArray1d& is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  //const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  //const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
  const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr  );

  for( lSet::const_iterator eg=set.begin() ; eg!=set.end() ; ++eg ) {

    realT flux = bc->GetValue(domain.m_feEdgeManager,eg,time);

    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if(itr != m_edgesToFaces.end() ){

      lArray1d& faces = itr->second;

      R1Tensor edgeCenter;
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *eg , edgeCenter );
      //const realT z = edgeCenter[2];

      // nb only uses first face to calculate volume flux
      // flux will be incorrect for other edges if at junction.
      for(size_t ii = 0; ii < faces.size(); ++ii){
        localIndex fc = faces[ii];

        // volFlux[*eg] = -flux*edgeLengths_new[*eg]*BoundedAperture(apertures_new[fc]); // flux out of a into b (hence -ve)
        realT volumeFlux = -flux*edgeLengths_new[*eg]*BoundedAperture(apertures_new[fc]);
        if(ii==0)volFlux[*eg] = volumeFlux; // FIXME just flux into first face

        R1Tensor vecFlux;

        vecFlux = edgeCenter;
        vecFlux -= faceCenters[fc];

        vecFlux.Normalize();
        vecFlux *= volumeFlux;

        if(is_ghost[fc] < 0) faceVelocity[fc] += vecFlux;

      }
    }
  }
}
/*
 * Fix me - still need to update to distribute fluxes across elements
 */
/// Update the velocity fluxes on the boundary faces
void ParallelPlateProppantFlowSolverImplicit::
    FixedNetFluxBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                  ObjectDataStructureBaseT& object ,
                                                  BoundaryConditionBase* bc, const lSet& set, realT time){


  //rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  //rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);

  //const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  Array1dT<R1Tensor>& faceVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);

  //Array1dT<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );

  iArray1d& is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  //const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  //const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
  //const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr  );

  for( lSet::const_iterator eg=set.begin() ; eg!=set.end() ; ++eg ) {

    realT flux = bc->GetValue(domain.m_feEdgeManager,eg,time); //FIXME

    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if(itr != m_edgesToFaces.end() ){

      lArray1d& faces = itr->second;

      R1Tensor edgeCenter;
      domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *eg , edgeCenter );
      //const realT z = edgeCenter[2];


      for(size_t ii = 0; ii < faces.size(); ++ii){
        //localIndex fc = faces[0];
        localIndex fc = faces[ii];

        //volFlux[*eg] = -flux*edgeLengths_new[*eg]; // FIXME

  	    realT volumeFlux = -flux*edgeLengths_new[*eg]; // flux out
        if(ii==0)volFlux[*eg] = volumeFlux; // FIXME just flux into first face

        R1Tensor vecFlux;

        vecFlux = edgeCenter;
        vecFlux -= faceCenters[fc];

        vecFlux.Normalize();
        vecFlux *= volumeFlux;

        if(is_ghost[fc] < 0) faceVelocity[fc] += vecFlux;

      }
    }
  }
}

void ParallelPlateProppantFlowSolverImplicit::
    FixedNetFluxBoundaryConditionB_VelocityUpdate( PhysicalDomainT& domain,
                                                  ObjectDataStructureBaseT& object ,
                                                  BoundaryConditionBase* bc, const lSet& set, realT time){


  //rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  //rArray1d& edgeVfs = domain.m_feEdgeManager.GetFieldData<realT>(ProppantVolumeFractionStr);

  const rArray1d& faceFluidPressure    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  Array1dT<R1Tensor>& faceVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);

  //Array1dT<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>(FluidVelocityStr);
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  const Array1dT<R1Tensor>& faceVelocities = domain.m_feFaceManager.GetFieldData<R1Tensor>( FluidVelocityStr );

  iArray1d& is_ghost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  //const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  const Array1dT<R1Tensor>& edgeCenters_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );


  //const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
  const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr  );

  const rArray1d& facePackVfs = domain.m_feFaceManager.GetFieldData<realT>(ProppantPackVolumeFractionStr);


  for( lSet::const_iterator eg=set.begin() ; eg!=set.end() ; ++eg ) {

    realT flux = bc->GetValue(domain.m_feEdgeManager,eg,time); //FIXME

    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if(itr != m_edgesToFaces.end() ){

      lArray1d& faces = itr->second;

      R1Tensor edgeCenter = edgeCenters_new[*eg];
      //domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, *eg , edgeCenter );
      //const realT z = edgeCenter[2];


      for(size_t ii = 0; ii < faces.size(); ++ii){

    	const localIndex kf = faces[ii];

    	//R1Tensor la_new = edgeCenter; la_new -= faceCenters[kf];

    	realT kappa_new;
    	if(m_usePowerlawFluid){
    	      	        kappa_new = OneFacePermeability_PowerLaw(edgeCenters_new, edgeLengths_new,faceCenters, apertures_new, faceVelocities,facePackVfs,*eg,kf);
    	} else {
    	              	kappa_new = OneFacePermeability(edgeCenters_new, edgeLengths_new,faceCenters, apertures_new,facePackVfs,*eg,kf);
    	}


    	realT pressure = faceFluidPressure[kf];
    	realT volumeFlux = kappa_new*(pressure-m_fixedFlux_pInj); // flux out


        if(ii==0)volFlux[*eg] = volumeFlux; // FIXME just flux into first face

        R1Tensor vecFlux;

        vecFlux = edgeCenter;
        vecFlux -= faceCenters[kf]; // (pointing towards edge)

        vecFlux.Normalize();
        vecFlux *= volumeFlux;

        if(is_ghost[kf] < 0) faceVelocity[kf] += vecFlux;

      }
    }
  }
}



// packVf = volume fraction filled by pack (solids+fluid) = vf_stationary_proppant/max_proppant_vf (i.e. 0<packVf<1)
// max_proppant_vf = max packing fraction of proppant (packing fraction in pack)
// fluid_mu = viscosity of fluid (w/o proppant)
// returns value assuming elements are entirely filled with proppant
realT ParallelPlateProppantFlowSolverImplicit::ProppantPackPermeability(realT h, realT w,realT l,const realT& max_proppant_vf, const realT& fluid_mu)
{
  realT k = m_proppantData.PackPermeability(max_proppant_vf,fluid_mu);
  return k*h*w/(l+1e-64);

}

// packVf = volume fraction filled by pack (solids+fluid) = vf_stationary_proppant/max_proppant_vf (i.e. 0<packVf<1)
realT ParallelPlateProppantFlowSolverImplicit::PartiallyFilledEdgePermeability(realT packVf, realT kappaOpen, realT kappaPack){
  return (1-packVf)*kappaOpen + packVf*kappaPack;
}



void ParallelPlateProppantFlowSolverImplicit::PostProcess (PhysicalDomainT& domain,
                  SpatialPartition& partition,
                  const sArray1d& namesOfSolverRegions){

	  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

	  rArray1d& aperture = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);

	  rArray1d& na = domain.m_feNodeManager.GetFieldData<realT>("nodalAperture");

	  //Calculate nodal aperture and proppant concentration
	  iArray1d nFace(domain.m_feNodeManager.DataLengths());
	  nFace = 0;
	  na = 0.0;
	  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
	  {
	    if (flowFaceType[kf] == 0)
	    {
	      for( lArray1d::iterator j = domain.m_feFaceManager.m_toNodesRelation[kf].begin() ;
	          j!=domain.m_feFaceManager.m_toNodesRelation[kf].end() ; ++j )
	      {
	        na[*j] += aperture[kf];
	        nFace[*j]++;
	      }
	    }
	  }

	  for (localIndex i = 0; i<domain.m_feNodeManager.DataLengths(); ++i)
	  {
	    if (nFace[i]>0)
	    {
	        na[i] /= nFace[i];
	    }
	  }

	  // Get nodal values of child nodes from their parents.
	  for (localIndex i = 0; i<domain.m_feNodeManager.DataLengths(); ++i)
	  {
	    localIndex ancestor = i;
	    while (domain.m_feNodeManager.m_parentIndex[ancestor] < domain.m_feNodeManager.DataLengths())
	    {
	      ancestor = domain.m_feNodeManager.m_parentIndex[ancestor];
	    }
	    if (ancestor != i)
	    {
	      na[i] = na[ancestor];
	    }
	  }


}


/// Register solver in the solver factory
REGISTER_SOLVER( ParallelPlateProppantFlowSolverImplicit )
