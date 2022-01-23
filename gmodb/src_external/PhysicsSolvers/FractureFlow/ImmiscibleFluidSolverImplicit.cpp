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
 * @date February 21, 2012
 */

#include "ImmiscibleFluidSolverImplicit.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"


//#include "../Materials/MaterialBaseStateDataT.h"

// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"





using namespace BoundaryConditionFunctions;
using namespace PS_STR;
using namespace PPFS;

namespace{
//  realT TINY = 1e-64;
  const std::string oldStr = "_old";
  const std::string InPlaneCurvatureStr = "InPlaneCurvature";
}

ImmiscibleFluidSolverImplicit::ImmiscibleFluidSolverImplicit( const std::string& name,
                                                              ProblemManagerT* const pm ):
ParallelPlateFlowSolverBase(name,pm),
m_faceSet(NULL),
m_numFaces(0),
m_faceDofMap(),
m_edgesToFaces(),
m_verboseFlag(false),
m_useMLPrecond(true),
m_phi(1.0),
m_dt(0.0),
this_mpi_process(pm->m_epetraComm.MyPID()),
n_mpi_processes(pm->m_epetraComm.NumProc()),
m_epetra_comm(pm->m_epetraComm),
row_map(),
sparsity(),
matrix(),
solution(),
rhs(),
syncedFields(),
syncedImmiscibleFields(),
m_TrilinosIndexStr(),
m_numerics()
{
  ++m_instances; 
  m_TrilinosIndexStr = "TwoDIMPPFS_" +  toString<int>(m_instances) + "_GlobalDof";
}

ImmiscibleFluidSolverImplicit::~ImmiscibleFluidSolverImplicit()
{
  // TODO Auto-generated destructor stub
}

void ImmiscibleFluidSolverImplicit::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  ParallelPlateFlowSolverBase::ReadXML(hdn);

  // Mixed difference parameter
  m_phi = hdn->GetAttributeOrDefault<realT>("phi",1.0); // backward difference by default.
  
  // Linear Solver
  m_numerics.m_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-10);
  m_numerics.m_maxIters = hdn->GetAttributeOrDefault<int>("max_solver_iterations",1000);

  // Flags
  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("UpdateAperture",false);
  m_verboseFlag = hdn->GetAttributeOrDefault<bool>("verbose",false);
  m_useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPrecoditioner",false);

  // Faceset
   m_flowFaceSetName = hdn->GetAttributeString("flowFaceSet");
   if(m_flowFaceSetName.empty()) m_flowFaceSetName = hdn->GetAttributeString("faceset");

   m_colorFunctionData.ReadXML(hdn);
   m_surfaceTensionData.ReadXML(hdn);
   m_permeabilityData.ReadXML(hdn);


   TICPP::HierarchicalDataNode* rfd_hdn = hdn->GetChild("RedFluid");
   TICPP::HierarchicalDataNode* bfd_hdn = hdn->GetChild("BlueFluid");

   if (!rfd_hdn)
     throw GPException("ImmiscibleFluidSolverColorFunction: Must have RedFluid defined in the input file.");
   if (!bfd_hdn)
     throw GPException("ImmiscibleFluidSolverColorFunction: Must have BlueFluid defined in the input file.");

   m_redFluidData.ReadXML(rfd_hdn);
   m_blueFluidData.ReadXML(bfd_hdn);


   m_min_dt = hdn->GetAttributeOrDefault<realT>("initial_dt",0.001);
   m_max_dt = hdn->GetAttributeOrDefault<realT>("max_dt",1);

}


void ImmiscibleFluidSolverImplicit::RegisterFields( PhysicalDomainT& domain )
{
  
  const bool plotOldValues = false; // no need to plot old values unless debugging - overwritten at end of timestep.

  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr+oldStr,true,plotOldValues);
  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr,true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr+oldStr,true,plotOldValues);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr,true,true);

  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>(FluidVelocityStr,true,true);

  // debug
  //domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("PressureVelocity",true,true);
  //domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("deltaPVelocity",true,true);
    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::pressure>();    
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::density>(); 
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::mass>(); 
  
  domain.m_feFaceManager.AddKeylessDataField<realT>("massRate",true,true); // face mass rate
  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::volume>();
  domain.m_feFaceManager.AddKeylessDataField<realT>("Volume_old",true,plotOldValues);

  domain.m_feFaceManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);


  domain.m_feEdgeManager.AddKeylessDataField<realT>(PermeabilityStr,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(VolumetricFluxStr,true,true);
  
 // domain.m_edgeManager.AddKeylessDataField<int>("FlowFaceCount",true,true);// debug
  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>(EdgeCenterStr+oldStr,true,plotOldValues);
  domain.m_feEdgeManager.AddKeylessDataField<R1Tensor>(EdgeCenterStr,true,true);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(EdgeLengthStr+oldStr,true,plotOldValues);
  domain.m_feEdgeManager.AddKeylessDataField<realT>(EdgeLengthStr,true,true);



  // color function
  //domain.m_feNodeManager.AddKeylessDataField<realT>("ColorFunction",true,true);
  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("ColorFunctionGradient",true,true);
  //domain.m_feNodeManager.AddKeylessDataField<realT>("NetFaceWeight",true,true);

  domain.m_feFaceManager.AddKeylessDataField<realT>("ColorFunction",true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("ColorFunctionIncrement",true,true);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("ColorFunctionGradient",true,true);


  domain.m_feFaceManager.AddKeylessDataField<realT>(InPlaneCurvatureStr,true,true);

  domain.m_feEdgeManager.AddKeylessDataField<realT>("DeltaP",true,true);
}


void ImmiscibleFluidSolverImplicit::Initialize( PhysicalDomainT& domain,SpatialPartition& partition )
{
  FaceManagerT& faceManager = domain.m_feFaceManager;
  EdgeManagerT& edgeManager = domain.m_feEdgeManager;
  
  m_faceSet = &(faceManager.GetSet(m_flowFaceSetName));
  m_nodeSet = &(domain.m_feNodeManager.GetSet(m_flowFaceSetName));
  m_numFaces = m_faceSet->size();
  
  // build face-dof map

  lSet::const_iterator si=m_faceSet->begin();
  for(localIndex i =0; i < m_numFaces; ++i, ++si){
    localIndex f = *si;
    m_faceDofMap[f] = i;
  }
  
 // iArray1d& ffCount = domain.m_edgeManager.GetFieldData<int>("FlowFaceCount"); // debug
  

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

  partition.SynchronizeFields(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);

  GenerateParallelPlateGeometricQuantities( domain,&partition,0 );
  OverwriteOldGeometricQuantities(domain);

  InitializeDensity( domain);

  partition.SynchronizeFields(syncedImmiscibleFields, CommRegistry::immiscibleFluidSolver);

}


/**
 * 
 * 
 */
void ImmiscibleFluidSolverImplicit:: UpdateAperture(PhysicalDomainT&  ){
}

void ImmiscibleFluidSolverImplicit:: SetupSystem (PhysicalDomainT&  domain,
                                                SpatialPartition& partition, const realT& )
{
  
  iArray1d& trilinos_index = domain.m_feFaceManager.GetFieldData<int>(m_TrilinosIndexStr);
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  
  // count local dof
  ///////////////////////////////
  
  // local rows
  int n_local_rows = 0;
  
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
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

void ImmiscibleFluidSolverImplicit :: Assemble (PhysicalDomainT&  domain,
                                                  SpatialPartition& partition __attribute__((unused)),
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

  const Array1dT<R1Tensor>& edgeCenters_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr + oldStr );
  const Array1dT<R1Tensor>& edgeCenters_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );

  const rArray1d& edgeLengths_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr + oldStr );
  const rArray1d& edgeLengths_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  const rArray1d& apertures_old = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr + oldStr );
  const rArray1d& apertures_new = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr  );
  
  rArray1d& faceFluidMass = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
//  rArray1d& faceFluidDensity = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();

  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);

//  const Array1dT<R1Tensor>& gradPhis = domain.m_feFaceManager.GetFieldData<R1Tensor>("ColorFunctionGradient");
  const rArray1d& colorFunction = domain.m_feFaceManager.GetFieldData<realT>("ColorFunction");

  rArray1d& curvatures = domain.m_feFaceManager.GetFieldData<realT>(InPlaneCurvatureStr);

  rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");


  // loop over faces and create an identity matrix
  Epetra_IntSerialDenseVector  faceDofIndex (1);
  Epetra_SerialDenseMatrix     face_matrix  (1,1);
  face_matrix(0,0) = 1;
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    if(face_is_ghost[*kf] < 0)
    {
      faceDofIndex(0) = trilinos_index[*kf];
      matrix->SumIntoGlobalValues(faceDofIndex, face_matrix);
    }
  }

    
  // loop over edges
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    //  if( edge_is_ghost[eg] < 0 )
    //  {
    int numFaces = itr->second.size();

    if( numFaces == 2)
    {
      localIndex kf = itr->second[0];
      localIndex kfb = itr->second[1];

      // calculate edge permeability
      realT kappa_old = TwoFacePermeability(edgeCenters_old,edgeLengths_old,
                                            faceCenters_old, apertures_old,
                                            eg,kf,kfb);
      realT kappa_new = TwoFacePermeability(edgeCenters_new,edgeLengths_new,
                                            faceCenters_new, apertures_new,
                                            eg,kf,kfb);

      edgePermeabilities[eg] = kappa_new;


      // build stiffness matrix and rhs
      if( kappa_new > 0.0 ){

        Epetra_IntSerialDenseVector  edgeDofIndex (numFaces);
        Epetra_SerialDenseVector     edge_rhs     (numFaces);
        Epetra_SerialDenseMatrix     edge_matrix  (numFaces,numFaces);

        // Color function
        /////////////////

        realT phiA = colorFunction[kf];
        realT phiB = colorFunction[kfb];
        realT phiAv = 0.5*(phiA+phiB);


        // FIXME - use new apertures/face centers or old? - assuming constant for time being
        R1Tensor branchVector = faceCenters_new[kfb]-faceCenters_new[kf];
        realT app = (apertures_new[kf] < m_permeabilityData.m_max_aperture) ? apertures_new[kf] : m_permeabilityData.m_max_aperture;
        realT appb = (apertures_new[kfb] < m_permeabilityData.m_max_aperture) ? apertures_new[kfb] : m_permeabilityData.m_max_aperture;


        realT convAngle = 0.0;
        if(m_surfaceTensionData.m_useConvergenceAngle){
          //FIXME - need to calculate angle of convergence based on local color gradient
          throw GPException("ImmiscibleFluidSolverColorFunction: ConvergenceAngle calculation is not yet implemented.");
        }
        realT h = 0.5*(app + appb);
        realT r1 = std::max(m_surfaceTensionData.CalculateNormalRadiusOfCurvature(h, convAngle),1e-64);

        realT curvA =1.0/r1;
        realT curvB = curvA;
        if(m_surfaceTensionData.m_calculateInPlaneCurvature){
          curvA += curvatures[kf];
          curvB += curvatures[kfb];
          // std::cout << "curvatures " <<  curv << " " << curvatures[kf] << " " << curvatures[kfb] << std::endl;
        }

        // surface tension forcing term
        //realT DeltaP = 0.5*m_surfaceTensionData.m_sigma*(curvA*gradPhis[kf]+curvB*gradPhis[kfb])*branchVector;
        realT DeltaP = 0.5*m_surfaceTensionData.m_sigma*(curvA+curvB)*(phiB-phiA);

        //realT DeltaP = 0.5*m_surfaceTensionData.m_sigma*curv*(gradPhis[kf]+gradPhis[kfb])*branchVector;

        // body force term
        DeltaP += (m_redFluidData.m_bodyForce*branchVector)*phiAv
            +(m_blueFluidData.m_bodyForce*branchVector)*(1.0-phiAv);

        edgeDeltaP[eg] = DeltaP;


        if( edge_is_ghost[eg] < 0 )
        {

          // stiffness matrix contribution
          /////////////////////////////////

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






          // Assemble
          ///////////

          edgeDofIndex[0] = trilinos_index[kf];
          edgeDofIndex[1] = trilinos_index[kfb];

          // self terms
          edge_matrix(0,0) = m_phi*kappa_new*dt
              *(rhoAv_s*dPdM_s[0] - 0.5*(pressure_s[1]-pressure_s[0])/volume_new[0]);

          edge_matrix(1,1) = m_phi*kappa_new*dt
              *(rhoAv_s*dPdM_s[1] - 0.5*(pressure_s[0]-pressure_s[1])/volume_new[1]);

          // cross terms
          edge_matrix(0,1) = -m_phi*kappa_new*dt
              *(rhoAv_s*dPdM_s[1] + 0.5*(pressure_s[1]-pressure_s[0])/volume_new[1]);
          edge_matrix(1,0) = -m_phi*kappa_new*dt
              *(rhoAv_s*dPdM_s[0] + 0.5*(pressure_s[0]-pressure_s[1])/volume_new[0]);

          // rhs
          edge_rhs(0) = (1.0-m_phi) * kappa_old*rhoAv_old*(pressure[1] - pressure[0])
                                + m_phi * kappa_new*rhoAv_s*(pressure_s[1]- pressure_s[0] )
                                - kappa_old*DeltaP*rhoAv_old;
          edge_rhs(1) = -edge_rhs(0);

          // assemble
          matrix->SumIntoGlobalValues(edgeDofIndex, edge_matrix);
          rhs->SumIntoGlobalValues(edgeDofIndex,edge_rhs);
        }
      } // ghost edge

    } else if(numFaces > 2){
      // fixme - not apparent how best to implement color function forcing on more than 2 faces.

    } // numFaces > 2 ?
    //} // ghost edge

  } // edge loop


  // boundary conditions
  ApplyBoundaryCondition<realT>(this, &ImmiscibleFluidSolverImplicit::PressureBoundaryCondition,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time );
  ApplyBoundaryCondition<realT>(this, &ImmiscibleFluidSolverImplicit::PressureBoundaryCondition,
                                domain, domain.m_feEdgeManager, "OutflowBoundary", time );

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  //rhs->Print(std::cout);
  //std::cout << matrix->NormInf() << std::endl;
  //std::string fileStr = "system-matrix_" + toString(time) + ".dat";
  //EpetraExt::RowMatrixToMatlabFile(fileStr.c_str(),*matrix);
  //exit(0);
}




realT ImmiscibleFluidSolverImplicit::TwoFacePermeability(const Array1dT<R1Tensor>& edgeCenters,
                                                           const rArray1d& edgeLengths,
                                                           const Array1dT<R1Tensor>& faceCenters,
                                                           const rArray1d& apertures,
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

  realT kappa = CalculatePermeability(la.L2_Norm(), lb.L2_Norm(), app, appb, w, m_mu, m_SHP_FCT);

  return kappa;
}


/// Apply a pressure boundary condition to a given set of edges
void ImmiscibleFluidSolverImplicit::PressureBoundaryCondition(PhysicalDomainT& domain,
                                                                ObjectDataStructureBaseT& object __attribute__((unused)),
                                                                BoundaryConditionBase* bc, const lSet& set, realT time){
 

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

  const rArray1d& colorFunction = domain.m_feFaceManager.GetFieldData<realT>("ColorFunction");


  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");
 

  Epetra_IntSerialDenseVector  face_dof(1);
  Epetra_SerialDenseVector     face_rhs(1);
  Epetra_SerialDenseMatrix     face_matrix(1,1);
    
  // loop over edges, find permeabilities and apply boundary conditions.

  lSet::const_iterator eg=set.begin() ;

  for( localIndex i=0; i < set.size() ; ++i, ++eg ){

    std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
    if( ( itr != m_edgesToFaces.end() ) && ( edge_is_ghost[*eg] < 0 ) ){

      lArray1d& faces = itr->second; 
      
      realT bc_pressure_old = bc->GetValue(domain.m_feEdgeManager,eg,time);
      realT bc_pressure_new = bc->GetValue(domain.m_feEdgeManager,eg,time+m_dt);


      realT bc_rho_old = rho_EOS(bc_pressure_old,m_bulk_modulus,m_rho_o);
      realT bc_rho_new = rho_EOS(bc_pressure_new,m_bulk_modulus,m_rho_o);

      // rhoP = rho * Pressure at BC
      //realT rhoP_old = rho_old*bc_pressure_old;
      //realT rhoP_new = rho_EOS(bc_pressure_new,m_bulk_modulus,m_rho_o)*bc_pressure_new;

      for(size_t ii = 0; ii < faces.size(); ++ii){
        
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

        const realT rhoAv_old = 0.5*(bc_rho_old+mass/volume_old);
        const realT rhoAv_s = 0.5*(bc_rho_new+mass/volume_new);

//        realT dPdM = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_old);
        realT dPdM_s = dPdRho_EOS(m_bulk_modulus,m_rho_o)*(1.0/volume_new); // dP/dm = dP/dRho * dRho/dm
//        realT dPdV = -(mass/volume_old)*dPdM;

        // body force
        realT phi = colorFunction[kf];

        realT deltaP =  (m_redFluidData.m_bodyForce*la_new)*phi
                       +(m_blueFluidData.m_bodyForce*la_new)*(1.0-phi);

        if(ii==0){
          edgeDeltaP[*eg] = deltaP; // will be incorrect at boundary edges with 2 or more faces
          edgePermeabilities[*eg] = kappa_new;// only assign first permeability (at junction)
        }


        /*
        // matrix
        face_matrix(0,0) = m_phi*(kappa_new/volume_new)*(pressure + dPdM*mass)*m_dt;

        // rhs
        face_rhs(0) = (1.0-m_phi) * kappa_old*( rhoP_old - pressure * mass/volume_old)
                          + m_phi * kappa_new*( rhoP_new - (pressure + deltaV*dPdV) * mass/volume_new)
                                  - kappa_old* deltaP*rhoAv_old;
        */

        // matrix
        face_matrix(0,0) = m_phi*kappa_new*m_dt
                           *(rhoAv_s*dPdM_s- 0.5*(bc_pressure_new-fc_pressure_s)/volume_new);

        // rhs
        face_rhs(0) = (1.0-m_phi) * kappa_old*rhoAv_old*(bc_pressure_old - fc_pressure)
                                    + m_phi * kappa_new*rhoAv_s*(bc_pressure_new - fc_pressure_s )
                                    - kappa_old*deltaP*rhoAv_old;



        face_dof(0) = trilinos_index[kf]; 
        matrix->SumIntoGlobalValues(face_dof, face_matrix);
    
        rhs->SumIntoGlobalValues(face_dof, face_rhs);
    
      }
    }
  }
}




/* Solve */

void ImmiscibleFluidSolverImplicit:: Solve (PhysicalDomainT&  domain,
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
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
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


  // ML preconditioner
  //////////////////////

  // create a parameter list for ML options
  Teuchos::ParameterList MLList;

  ML_Epetra::SetDefaults("SA",MLList);

  // create the preconditioning object.

/*
  ML_Epetra::MultiLevelPreconditioner* MLPrec(NULL);
  if(m_useMLPrecond){
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*matrix, MLList);
  }
*/


  Teuchos::RCP <ML_Epetra::MultiLevelPreconditioner> MLPrec;
  if(m_useMLPrecond){
    MLPrec =  Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*matrix, MLList));
  }



  // Aztec solver
  /////////////////
  AztecOO solver(problem);

          // ML preconditioner
          if(m_useMLPrecond) solver.SetPrecOperator( &(*MLPrec) );

          solver.SetAztecOption(AZ_solver,AZ_bicgstab);
          // solver.SetAztecOption(AZ_solver,AZ_gmres);
          solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
          solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
          solver.SetAztecOption(AZ_conv,AZ_rhs);
          if(!m_verboseFlag) solver.SetAztecOption(AZ_output,AZ_none);

          
          
  solver.Iterate(m_numerics.m_maxIters,m_numerics.m_tol);
  //solution->Print(std::cout);


  // destroy the preconditioner
/*
  if(m_useMLPrecond){
    delete MLPrec;
  }
*/


  // copy solution to faces
  ////////////////////////
  
  // refresh the view - not sure if needed (trying to eliminate a bug)
  solution->ExtractView(&local_solution,&dummy);

  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
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


void ImmiscibleFluidSolverImplicit::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(FaceCenterStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("massRate");
  //syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("Mass"); //needed?
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(PressureStr); // needed?
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(ApertureStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(ApertureStr+oldStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_TrilinosIndexStr);
//  syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(VolumetricFluxStr); // can't use - not oriented
//  syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back("DeltaP");// can't use - not oriented

  syncedImmiscibleFields.clear();
  syncedImmiscibleFields[PhysicalDomainT::FiniteElementFaceManager].push_back("ColorFunction");
  syncedImmiscibleFields[PhysicalDomainT::FiniteElementFaceManager].push_back("ColorFunctionGradient");
  syncedImmiscibleFields[PhysicalDomainT::FiniteElementNodeManager].push_back("ColorFunctionGradient");

  partition.SetBufferSizes(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);
  partition.SetBufferSizes(syncedImmiscibleFields, CommRegistry::immiscibleFluidSolver);
}


double ImmiscibleFluidSolverImplicit::TimeStep( const realT& time,
                                                        const realT& dt,
                                                        const int cycleNumber,
                                                        PhysicalDomainT& domain,
                                                        const sArray1d& namesOfSolverRegions ,
                                                        SpatialPartition& partition,
                                                        FractunatorBase* const fractunator )
{


  m_stabledt.m_maxdt = std::numeric_limits<double>::max();

  m_dt = dt; // Needed for pressure BC's

  GenerateParallelPlateGeometricQuantities( domain, &partition,dt );

  CalculateMassRate( domain, partition, time,dt );


  UpdateEOS( domain,dt );


  UpdateFlux( time, partition,domain);


  UpdateColorFunction( domain,partition,time, dt );


  OverwriteOldGeometricQuantities(domain);

  if(m_stabledt.m_maxdt < m_min_dt || m_stabledt.m_maxdt >= std::numeric_limits<double>::max()){
    m_stabledt.m_maxdt = m_min_dt;
  }

  // limit maximum increase in timestep
  if(m_stabledt.m_maxdt > 2*dt && dt > 0.0){
    m_stabledt.m_maxdt = 2*dt;
  }

  if(m_stabledt.m_maxdt > m_max_dt){
    m_stabledt.m_maxdt = m_max_dt;
  }

  return dt;
}


void ImmiscibleFluidSolverImplicit::GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain, SpatialPartition* partition, realT dt )
{

  // synchronize fields
  if(partition)
    partition->SynchronizeFields(syncedImmiscibleFields, CommRegistry::immiscibleFluidSolver);

  rArray1d& phis = domain.m_feFaceManager.GetFieldData<realT>("ColorFunction");
  rArray1d& deltaPhis = domain.m_feFaceManager.GetFieldData<realT>("ColorFunctionIncrement");
  Array1dT<R1Tensor>& gradPhis = domain.m_feFaceManager.GetFieldData<R1Tensor>("ColorFunctionGradient");
  rArray1d& curvatures = domain.m_feFaceManager.GetFieldData<realT>(InPlaneCurvatureStr);

  //rArray1d& node_phis = domain.m_feNodeManager.GetFieldData<realT>("ColorFunction");
  Array1dT<R1Tensor>& node_gradPhis = domain.m_feNodeManager.GetFieldData<R1Tensor>("ColorFunctionGradient");
 // rArray1d& node_weight = domain.m_feNodeManager.GetFieldData<realT>("NetFaceWeight");

  Array1dT<R1Tensor>& faceCenter_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  rArray1d& fluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& aperture_new = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);

  // reset fields
 // node_weight = 0;
 // node_phis = 0;
  deltaPhis = 0;
  node_gradPhis = 0;

  // debug
//  iArray1d& faceGhostRank            = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
//  iArray1d& nodeGhostRank            = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();

  // update face quantities
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf ) {


    R1Tensor& center = faceCenter_new[*kf];
    domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *kf, center);

    if(m_doApertureUpdate){
      R1Tensor gap;
      R1Tensor N;
      N = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *kf );

      gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, *kf );
      aperture_new[*kf] = Dot(gap,N) ;

      if( aperture_new[*kf]<m_min_aperture )
          aperture_new[*kf] = m_min_aperture;
      else if( aperture_new[*kf] > m_max_aperture )
          aperture_new[*kf] = m_max_aperture;
    }

    fluidVolume_new[*kf] = aperture_new[*kf] * domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, *kf );

 // loop over face nodes - add contribution to node color function and node color function weight
    // FIXME - maybe better done at the end of the timestep - but will need additional routine to initialize node color function
    const lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[*kf];
    for( size_t a=0 ; a<nodeList.size(); ++a ){
      const localIndex nd = nodeList[a];
      R1Tensor nodePos; domain.m_feNodeManager.GetPosition(nd, nodePos);


      R1Tensor l = center - nodePos;

     // realT ll = std::max(l.L2_Norm(),1e-64);
     // realT w = 1.0/ll;

     // node_weight[nd] += w;
      //node_phis[nd] += w*phis[*kf];

      // fixme assumes regular grid aligned in x,y
      l[0] = 0.25/l[0];
      l[1] = 0.25/l[1];
      node_gradPhis[nd] += phis[*kf]*l;

    }

  }


  // loop over nodes - adjust node color function
/*  for( lSet::const_iterator kn=m_nodeSet->begin(); kn!=m_nodeSet->end() ; ++kn ){
    const localIndex nd = *kn;
    node_phis[nd] /= node_weight[nd];
  }*/

  // update edge properties
  iArray1d& edge_is_ghost            = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();
  Array1dT<R1Tensor>& edgeCenter_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  Array1dT<realT>& edgeLength_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );

  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    domain.m_feEdgeManager.EdgeCenter( domain.m_feNodeManager, eg , edgeCenter_new[eg] );
    edgeLength_new[eg] = domain.m_feEdgeManager.EdgeLength( domain.m_feNodeManager, eg);
    if( true || edge_is_ghost[eg] < 0 )
    {

      // correct color gradients on boundary nodes
      if( itr -> second.size() == 1)
      {
        localIndex fc = itr->second[0];
        R1Tensor edgeNorm = (edgeCenter_new[eg] - faceCenter_new[fc]).UnitVector();

        localIndex ndA = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
        localIndex ndB = domain.m_feEdgeManager.m_toNodesRelation(eg,1);

        node_gradPhis[ndA] -= (node_gradPhis[ndA]*edgeNorm)*edgeNorm;
        node_gradPhis[ndB] -= (node_gradPhis[ndB]*edgeNorm)*edgeNorm;

        node_gradPhis[ndA] *= 2.0; // fixme hacked in for regular grid

      }
    }
  }

  // gradient at boundary edges - part  B - correct gradient along edge norm
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;

    if( edge_is_ghost[eg] < 0 )
    {
      // correct color gradients on boundary nodes
      if( itr -> second.size() == 1)
      {
        localIndex fc = itr->second[0];
        R1Tensor edgeNorm = (edgeCenter_new[eg] - faceCenter_new[fc]).UnitVector();

        localIndex ndA = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
        localIndex ndB = domain.m_feEdgeManager.m_toNodesRelation(eg,1);

        int count = 0;
        R1Tensor faceGrad;
        for( lArray1d::size_type i = 0; i< domain.m_feFaceManager.m_toNodesRelation[fc].size(); ++i){
          localIndex nd = domain.m_feFaceManager.m_toNodesRelation[fc][i];
          if(nd != ndA && nd != ndB){
             ++count;
             faceGrad += node_gradPhis[nd];
          }
        }
        faceGrad = (1.0/count)*faceGrad;

        node_gradPhis[ndA] += 0.5*(faceGrad*edgeNorm)*edgeNorm;
        node_gradPhis[ndB] += 0.5*(faceGrad*edgeNorm)*edgeNorm;

      }
    }
  }


  // synchronize fields
  if(partition)
    partition->SynchronizeFields(syncedImmiscibleFields, CommRegistry::immiscibleFluidSolver);



  // curvature calculation
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    const localIndex fc = *kf;


    if(phis[fc] < 0.0 || phis[fc] >  1.0){

      curvatures[fc] = 0.0;

    } else {

      R1Tensor& center = faceCenter_new[fc];
      R1Tensor gradPhi;
      //realT straightCurve = 0.0;
      //    realT dLNdX = 0.0;
      //    realT dLNdY = 0.0;
      //    realT divN = 0.0;
    //  const lArray1d& nodeList = domain.m_feFaceManager.m_FaceToNodeMap[fc];
    const lArray1d& nodeList = domain.m_feFaceManager.m_toNodesRelation[fc];

      // curvature calculation - Brackbill
      ////////////////////////////////////

      /*

    for( size_t a=0 ; a<nodeList.size(); ++a ){
      const localIndex nd = nodeList[a];
      R1Tensor nodePos; domain.m_feNodeManager.GetPosition(nd, nodePos);
      R1Tensor& node_gradPhi = node_gradPhis[nd];

      gradPhi += node_gradPhi;

      // FIXME assumes regular grid
      R1Tensor l = nodePos-center;
      realT ln = node_gradPhi.L2_Norm();
      dLNdX += ln/l[0];
      dLNdY += ln/l[1];
      divN += node_gradPhi[0]/l[0] + node_gradPhi[1]/l[1];

      // debug - straight curvature calculation
    //R1Tensor norm_node_gradPhi = node_gradPhi;
    // norm_node_gradPhi /=  ln+1e-64;
    //  straightCurve -= norm_node_gradPhi[0]/l[0] + norm_node_gradPhi[1]/l[1];
    }
    // FIXME assumes 4 nodes per face
    dLNdX /= 4.0;
    dLNdY /= 4.0;
    divN  /= 4.0;
    gradPhi/= 4.0;




    realT lnn = gradPhi.L2_Norm();

    gradPhis[fc] = gradPhi;



   // if(phis[fc] > 1e-8 && phis[fc] <  1.0-1e-8){
    if(phis[fc] > 0.0 && phis[fc] <  1.0){
      realT termA = (dLNdX*gradPhi[0] + dLNdY*gradPhi[1])/lnn;

      curvatures[fc] = (termA - divN)/lnn;

     // curvatures[fc] = straightCurve;
    } else {
      curvatures[fc] = 0.0;
    }
       */

      // straightCurve/= 4.0;
      /* DEBUG */
      /*
     if(gradPhi.L2_Norm() > 0.1){
       for( size_t a=0 ; a<nodeList.size(); ++a ){
         const localIndex nd = nodeList[a];
         std::cout << a << std::endl;
         R1Tensor nodePos; domain.m_feNodeManager.GetPosition(nd, nodePos);
         R1Tensor& node_gradPhi = node_gradPhis[nd];
         R1Tensor l = nodePos-center;
         realT ln = node_gradPhi.L2_Norm();
         std::cout << "  Position " << l << std::endl;
         std::cout << "  ln " << ln << std::endl;
         std::cout << "  node_gradPhi " << node_gradPhi << std::endl;

       }
       std::cout << "dLNdX " << dLNdX << std::endl;
       std::cout << "dLNdY " << dLNdY << std::endl;
       std::cout << "divN " << divN << std::endl;
       std::cout << "gradPhi " << gradPhi << std::endl;
       exit(0);
     }
       */


      // Haliday curvature calculation
      ////////////////////////////////

      // kappa = nx*ny*(dnydx + dnxdy) -  nx*nx*dnydy - ny*ny*dnxdx;
      realT dnxdx = 0;
      realT dnxdy = 0;
      realT dnydx = 0;
      realT dnydy = 0;
      R1Tensor n;
      for( size_t a=0 ; a<nodeList.size(); ++a ){
        const localIndex nd = nodeList[a];
        R1Tensor nodePos; domain.m_feNodeManager.GetPosition(nd, nodePos);
        R1Tensor l = nodePos-center;
        R1Tensor vN = node_gradPhis[nd];
        gradPhi += vN;
        vN.Normalize();

        n += vN;
        dnxdx += vN[0]/l[0];
        dnydx += vN[1]/l[0];
        dnxdy += vN[0]/l[1];
        dnydy += vN[1]/l[1];
      }

      n = 0.25*n;
      gradPhi = 0.25*gradPhi;

      dnxdx *= 0.25;
      dnydx *= 0.25;
      dnxdy *= 0.25;
      dnydy *= 0.25;

      gradPhis[fc] = gradPhi;

      realT nx = n[0];
      realT ny = n[1];
      curvatures[fc] = nx*ny*(dnydx + dnxdy) -  nx*nx*dnydy - ny*ny*dnxdx;


    }

  }

  // synchronize fields
  if(partition)
    partition->SynchronizeFields(syncedImmiscibleFields, CommRegistry::immiscibleFluidSolver);

}


void ImmiscibleFluidSolverImplicit::OverwriteOldGeometricQuantities( PhysicalDomainT& domain){

  Array1dT<R1Tensor>& faceCenter_new = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  rArray1d& fluidVolume_new  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  rArray1d& aperture_new = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr);

  Array1dT<R1Tensor>& faceCenter_old = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr+oldStr );
  rArray1d& fluidVolume_old  = domain.m_feFaceManager.GetFieldData<realT>("Volume_old");
  rArray1d& aperture_old = domain.m_feFaceManager.GetFieldData<realT>(ApertureStr+oldStr);

  // update face quantities
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf ) {
    faceCenter_old[*kf] = faceCenter_new[*kf];
    aperture_old[*kf] = aperture_new[*kf];
    fluidVolume_old[*kf] = fluidVolume_new[*kf];
  }

  // update edge properties
  Array1dT<R1Tensor>& edgeCenter_new = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr );
  Array1dT<realT>& edgeLength_new = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr );
  Array1dT<R1Tensor>& edgeCenter_old = domain.m_feEdgeManager.GetFieldData<R1Tensor>( EdgeCenterStr +oldStr );
  Array1dT<realT>& edgeLength_old = domain.m_feEdgeManager.GetFieldData<realT>( EdgeLengthStr +oldStr);

  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    localIndex eg = itr->first;
    edgeCenter_old[eg] = edgeCenter_new[eg];
    edgeLength_old[eg] = edgeLength_new[eg];
    
  }
}



void ImmiscibleFluidSolverImplicit::CalculateMassRate( PhysicalDomainT& domain,
                                                         SpatialPartition& partition, realT time, realT dt )
{

  SetupSystem (domain,partition, time);

  Assemble    (domain,partition, time, dt);

  Solve       (domain,partition);
  
}

// set initial fluid density based on known pressure;
void ImmiscibleFluidSolverImplicit::InitializeDensity( PhysicalDomainT& domain)
{
  rArray1d& mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& density  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const rArray1d& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  for( lSet::const_iterator fc=m_faceSet->begin() ; fc!=m_faceSet->end() ; ++fc ) {
    density[*fc] = rho_EOS(pressure[*fc],m_bulk_modulus,m_rho_o );
    mass[*fc] = density[*fc]*fluidVolume[*fc];
  }
}


void ImmiscibleFluidSolverImplicit::UpdateEOS( PhysicalDomainT& domain,const realT dt)
{
  rArray1d& mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();
  rArray1d& density  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  rArray1d& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  const rArray1d& massRate = domain.m_feFaceManager.GetFieldData<realT>("massRate");

  for( lSet::const_iterator fc=m_faceSet->begin() ; fc!=m_faceSet->end() ; ++fc ) {
    mass[*fc] += massRate[*fc] *dt;
    density[*fc] = mass[*fc] / fluidVolume[*fc];
    realT P = P_EOS(density[*fc],m_bulk_modulus,m_rho_o );
    pressure[*fc] = P;
    // propagate pressure to children
    lArray1d& childFaces = domain.m_feFaceManager.m_childIndices[*fc];
    for(unsigned i =0; i < childFaces.size(); ++i){
      pressure[childFaces[i]] = P;
    }
  }
}


void ImmiscibleFluidSolverImplicit::UpdateFlux( const realT time, SpatialPartition& partition,PhysicalDomainT& domain)
{
  
  const rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);
  const rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");
  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);
  volFlux = 0.0;
  
  const rArray1d& pressures    = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  
  Array1dT<R1Tensor>& flowVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>( FluidVelocityStr );
  flowVelocity = 0.0;
  
  const Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );
  
  iArray1d& is_ghost       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
  
  // debug
  //Array1dT<R1Tensor>& pressureVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>( "PressureVelocity" );
  //pressureVelocity = 0.0;
  //Array1dT<R1Tensor>& deltaPVelocity = domain.m_feFaceManager.GetFieldData<R1Tensor>( "deltaPVelocity" );
  //deltaPVelocity = 0.0;


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
      volFlux[eg] =  edgePermeabilities[eg]*(Pa+edgeDeltaP[eg]-Pb); // flux A -> B
      
      R1Tensor vecFlux;

      vecFlux = faceCenters[kfb];
      vecFlux -= faceCenters[kfa];


      vecFlux.Normalize();

      //debug
      //R1Tensor dpFlux;
      //dpFlux = vecFlux*edgePermeabilities[eg]*edgeDeltaP[eg];

      vecFlux *= volFlux[eg];

      if(is_ghost[kfa] < 0) {
        flowVelocity[kfa] += vecFlux;

        // debug
        //deltaPVelocity[kfa] += dpFlux;
        //pressureVelocity[kfa] += vecFlux-dpFlux;

      }
      if(is_ghost[kfb] < 0){
        flowVelocity[kfb] += vecFlux;

        // debug
        //deltaPVelocity[kfb] += dpFlux;
        //pressureVelocity[kfb] += vecFlux-dpFlux;
      }
      
    }  
  }
  
  ApplyBoundaryCondition<realT>(this, &ImmiscibleFluidSolverImplicit::PressureBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time );
  ApplyBoundaryCondition<realT>(this, &ImmiscibleFluidSolverImplicit::PressureBoundaryCondition_VelocityUpdate,
                                domain, domain.m_feEdgeManager, "OutflowBoundary", time );

  flowVelocity *= 0.5;
  //pressureVelocity *=0.5;
  //deltaPVelocity *= 0.5;

  // sync edge fluxes - needed for color function update? - maybe not if pressures have been updated already
  partition.SynchronizeFields(syncedFields, CommRegistry::steadyStateParallelPlateFlowSolver);
    
}

/**
 *
 *
**/
void ImmiscibleFluidSolverImplicit::UpdateColorFunction( PhysicalDomainT& domain,
                                                         SpatialPartition& partition,const realT time, const realT dt)
{
  rArray1d& phis = domain.m_feFaceManager.GetFieldData<realT>("ColorFunction");
  rArray1d& deltaPhis = domain.m_feFaceManager.GetFieldData<realT>("ColorFunctionIncrement");
  //rArray1d& node_phis = domain.m_feNodeManager.GetFieldData<realT>("ColorFunction");
  Array1dT<R1Tensor>& gradPhis = domain.m_feFaceManager.GetFieldData<R1Tensor>("ColorFunctionGradient");
//  rArray1d& node_weight = domain.m_feNodeManager.GetFieldData<realT>("NetFaceWeight");
//  rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>( ApertureStr );

 // iArray1d& edge_ghost_rank = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

  Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);

  const rArray1d& faceFluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

  Array1dT<R1Tensor>& faceCenters = domain.m_feFaceManager.GetFieldData<R1Tensor>( FaceCenterStr );




  // loop over edges - calculate color function flux
  std::map<localIndex,lArray1d>::iterator itrEnd = m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.begin(); itr!=itrEnd; ++itr )
  {
    localIndex eg = itr->first;
    int numFaces = itr->second.size();

    if( numFaces == 2 ){

      eg = itr->first;
      localIndex fca = itr->second[0];
      localIndex fcb = itr->second[1];

//      localIndex& nda = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
//      localIndex& ndb = domain.m_feEdgeManager.m_toNodesRelation(eg,1);

     // localIndex& nda = domain.m_feEdgeManager.m_edgesToNodes(eg,0);
     // localIndex& ndb = domain.m_feEdgeManager.m_edgesToNodes(eg,1);
     //   localIndex& nda = domain.m_feEdgeManager.m_toNodesRelation(eg,0);
     //   localIndex& ndb = domain.m_feEdgeManager.m_toNodesRelation(eg,1);

      realT phiA = phis[fca];
      realT phiB = phis[fcb];

      realT volA = faceFluidVolume[fca];
      realT volB = faceFluidVolume[fcb];

      R1Tensor l = faceCenters[fcb]- faceCenters[fca]; // normalized branch vector
//      realT lnorm = l.Normalize();

      R1Tensor t; // normalized edge vector
      domain.m_feEdgeManager.EdgeVector(domain.m_feNodeManager, eg, t);
//      realT tnorm = t.Normalize();

      //        realT height = 0.5*(apertures[fca] + apertures[fcb]);

      // Calculate color gradient
      realT dotlt = Dot(l,t);
      R1Tensor n = l - dotlt*t; n.Normalize();
//      realT dotln = Dot(l,n);

      //       realT dPhidt = (node_phis[ndb] - node_phis[nda])/tnorm;
      //       realT dPhidl = (phiB-phiA)/lnorm;

      //        realT dPhidn = (dPhidl - dPhidt*dotlt)/dotln;

      //R1Tensor colorGrad = dPhidn * n + dPhidt * t;
      R1Tensor colorGradA =  gradPhis[fca]; colorGradA.Normalize();
      R1Tensor colorGradB =  gradPhis[fcb]; colorGradB.Normalize();


      // Perform anti-diffusion
      realT cosThetaBetaA = Dot(colorGradA,n)*m_colorFunctionData.m_beta;
      realT cosThetaBetaB = Dot(colorGradB,n)*m_colorFunctionData.m_beta;
      realT D = m_colorFunctionData.m_D*dt;//*tnorm*height/dotln;// - fixme should we account for aperture or base on lb method (i.e. prop to total fluid)
      realT minVol = std::min(volA,volB);

      realT phiFlux = D*(phiB*(1.0-(1.0-phiB)*cosThetaBetaB)
          -phiA*(1.0+(1.0-phiA)*cosThetaBetaA))*minVol;

      // upwind advection
      if (false){
        // original (working) version
        if(volFlux[eg] > 0){
          // vol flux is from a-> b
          // phiFlux records flux into a
          phiFlux -= phiA*volFlux[eg] *dt;
        } else {
          phiFlux -= phiB*volFlux[eg] *dt;
        }
      } else {
        // w anti-diffusion on advective flux

        if(volFlux[eg] > 0){
          // vol flux is from a-> b
          // phiFlux records flux into a
          realT adTerm = 1.0+ (1.0-phiA)*cosThetaBetaA;
          phiFlux -= adTerm*phiA*volFlux[eg] *dt;
        } else {
          realT adTerm = 1.0 - (1.0-phiB)*cosThetaBetaB;
          phiFlux -= adTerm*phiB*volFlux[eg] *dt;
        }
      }


      //  phiFlux -= 0.5*(phiA+phiB)*volFlux[eg] *dt;

      // limit flux to avoid invalid color functions
      /*
        if( ( phiFlux > (1.0-phiA)*volA)  || (phiB*volB - phiFlux < 0.0) ){
          phiFlux = std::min( (1.0-phiA)*volA ,phiB*volB );
        } else if (  (-phiFlux >  (1.0-phiB)*volB) || (phiA*volA + phiFlux < 0.0) ){
          phiFlux = -std::min( (1.0-phiB)*volB,phiA*volA);
        }
       */

      deltaPhis[fca] += phiFlux/volA;
      deltaPhis[fcb] -= phiFlux/volB;

    }
  }


  ApplyBoundaryCondition<realT>(this, &ImmiscibleFluidSolverImplicit::PressureBoundaryCondition_ColorFunctionUpdate,
                                domain, domain.m_feEdgeManager, std::string(Field<FieldInfo::pressure>::Name()), time );
  // outflow boundary assumed impermeable to non-wetting phase.

  // Increment color function
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {

    phis[*kf] += deltaPhis[*kf];

    if(deltaPhis[*kf] > 0.0 && phis[*kf] < 0.999 ){
      realT maxdt = 0.01*dt*(1.0-phis[*kf])/fabs(deltaPhis[*kf]);
      m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt,maxdt);
    //  std::cout << maxdt << std::endl;
    } else if(deltaPhis[*kf] < 0.0 && phis[*kf] > 0.001){
      realT maxdt = 0.01*dt*phis[*kf]/fabs(deltaPhis[*kf]);
      m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt,maxdt);
    }
  }

  // synchronize fields
  partition.SynchronizeFields(syncedImmiscibleFields, CommRegistry::immiscibleFluidSolver);

}




/// Update the velocity fluxes on the boundary faces
void ImmiscibleFluidSolverImplicit::
       PressureBoundaryCondition_VelocityUpdate( PhysicalDomainT& domain,
                                                 ObjectDataStructureBaseT& object __attribute__((unused)),
                                                 BoundaryConditionBase* bc, const lSet& set, realT time){
  
  rArray1d& edgePermeabilities = domain.m_feEdgeManager.GetFieldData<realT>(PermeabilityStr);

  //const rArray1d& colorFunction = domain.m_feFaceManager.GetFieldData<realT>("ColorFunction");
  rArray1d& edgeDeltaP = domain.m_feEdgeManager.GetFieldData<realT>("DeltaP");
  
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

        volFlux[*eg] =  edgePermeabilities[*eg]*(Pa-Pb+edgeDeltaP[*eg]); // flux out of a into b

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

/// Implement the color function fluxes on the boundary faces
void ImmiscibleFluidSolverImplicit::
       PressureBoundaryCondition_ColorFunctionUpdate( PhysicalDomainT& domain,
                                                      ObjectDataStructureBaseT& object __attribute__((unused)),
                                                      BoundaryConditionBase* , const lSet& set, realT )
{
   const rArray1d& phis = domain.m_feFaceManager.GetFieldData<realT>("ColorFunction");
   rArray1d& deltaPhis = domain.m_feFaceManager.GetFieldData<realT>("ColorFunctionIncrement");

   const Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);

   const rArray1d& faceFluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

   const iArray1d& face_ghost_rank = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
//   const iArray1d& edge_ghost_rank = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

   for( lSet::const_iterator eg=set.begin() ; eg!=set.end() ; ++eg ) {

     std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
     if(itr != m_edgesToFaces.end() ){

       localIndex fc = itr->second[0];
       realT vol = faceFluidVolume[fc];

       if(vol > 0.0 && face_ghost_rank[fc] < 0){
         
         realT phiFlux;
         if(volFlux[*eg] > 0.0){
           // flowing out
           phiFlux = -phis[fc]*volFlux[*eg] *m_dt;
         } else {
           if(phis[fc] > 0.5){
             // round influx phi to 1
             phiFlux = -volFlux[*eg] *m_dt;
           } else {
             // round influx to 0
             phiFlux = 0;
           }
         }
         deltaPhis[fc] += phiFlux/vol;
       }
     }
   }

}

/// Implement the color function fluxes on the boundary faces
void ImmiscibleFluidSolverImplicit::
     OutflowPressureBoundaryCondition_ColorFunctionUpdate( PhysicalDomainT& domain,
                                                           ObjectDataStructureBaseT& object __attribute__((unused)),
                                                           BoundaryConditionBase* , const lSet& set, realT )
{
//   const rArray1d& phis = domain.m_feFaceManager.GetFieldData<realT>("ColorFunction");
   rArray1d& deltaPhis = domain.m_feFaceManager.GetFieldData<realT>("ColorFunctionIncrement");

   const Array1dT<realT>& volFlux =domain.m_feEdgeManager.GetFieldData<realT>(VolumetricFluxStr);

   const rArray1d& faceFluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();

   const iArray1d& face_ghost_rank = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
//   const iArray1d& edge_ghost_rank = domain.m_feEdgeManager.GetFieldData<FieldInfo::ghostRank>();

   for( lSet::const_iterator eg=set.begin() ; eg!=set.end() ; ++eg ) {

     std::map<localIndex,lArray1d>::iterator itr = m_edgesToFaces.find(*eg);
     if(itr != m_edgesToFaces.end() ){

       localIndex fc = itr->second[0];
       realT vol = faceFluidVolume[fc];

       if(vol > 0.0 && face_ghost_rank[fc] < 0){
         
         realT phiFlux;
         if(volFlux[*eg] > 0.0){
           // flowing out
           phiFlux = 0.0; // -phis[fc]*volFlux[*eg] *m_dt;
         } else {
         /*  if(phis[fc] > 0.5){
             // round influx phi to 1
             phiFlux = -volFlux[*eg] *m_dt;
           } else {
             // round influx to 0
             phiFlux = 0;
           }*/
           phiFlux = 0;
         }
         deltaPhis[fc] += phiFlux/vol;
       }
     }
   }

}



void ImmiscibleFluidSolverImplicit::SetMaxStableTimeStep( PhysicalDomainT& ,
                                                           const sArray1d& namesOfSolverRegions,
                                                           SpatialPartition& partition __attribute__((unused)) )
{
  m_stabledt.m_maxdt =  m_min_dt;

  m_stabledt.m_maxdt *= this->m_courant;

}


/// Register solver in the solver factory
REGISTER_SOLVER( ImmiscibleFluidSolverImplicit )
