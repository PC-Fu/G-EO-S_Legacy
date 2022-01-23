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
 * @file LagrangeDynamicsParallelPlateFlowExplicit_Proppant.cpp
 * @author settgast1
 * @date Feb 28, 2011
 */

#include "LagrangeDynamicsParallelPlateFlowExplicit_Proppant.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "SurfaceGeneration/FractunatorBase.h"
#include "TwoDADRSolver_Proppant.h"
#include "ObjectManagers/ProblemManagerT.h"

LagrangeDynamicsParallelPlateFlowExplicit_Proppant::LagrangeDynamicsParallelPlateFlowExplicit_Proppant( const std::string& name,
                                                                                      ProblemManagerT* const pm ):
SolverBase(name,pm),
m_ldSolve(name,pm),
m_oldPPStepTime(0),
m_newPPStepTime(0),
m_pp_dt(0)
//m_ppSolve(name,pm)
{
	  m_solverMapPtr = &(pm->m_solvers);
}

LagrangeDynamicsParallelPlateFlowExplicit_Proppant::~LagrangeDynamicsParallelPlateFlowExplicit_Proppant()
{
  // TODO Auto-generated destructor stub
}

void LagrangeDynamicsParallelPlateFlowExplicit_Proppant::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML(hdn);

  //m_ppSolve.ReadXML(hdn);
  m_FlowSolverName = hdn->GetAttributeString("flowSolver");
  if( m_FlowSolverName.size() > 0){
     std::cout << "  Parallel Plate Proppant Flow Solver " << m_FlowSolverName << std::endl;
  } else {
	 GPException("Error LagrangeDynamicsParallelPlateFlowExplicit_Proppant: No flowSolver name provided.");
  }
  m_ProppantSolverName  = hdn->GetAttributeString("proppantSolver");
  if( m_ProppantSolverName.size() > 0){
     std::cout << "  Parallel Plate Proppant Transport Solver " << m_ProppantSolverName << std::endl;
  } else {
	 std::cout << "No proppantSolver name provided." <<std::endl; // not an error for now
  }

  m_ldSolve.ReadXML(hdn);

  m_kJn = hdn->GetAttributeOrDefault<realT>("normalJointStiffness", 1.0e10);
  m_kJs = hdn->GetAttributeOrDefault<realT>("shearJointStiffness", 1.0e10);
  m_COFJ = fabs(hdn->GetAttributeOrDefault<realT>("COFJoint", 0.5));
  m_fLockedInSIF =  hdn->GetAttributeOrDefault<realT>("lockedInSIFFactor", 1.0);

  m_faceStrengthRandomFactor = hdn->GetAttributeOrDefault<realT>("faceStrengthRandomFactor", 0.0);


  m_nSubsteps = hdn->GetAttributeOrDefault<realT>("substeps", 1);

  m_doFracture =  hdn->GetAttributeOrDefault<bool>("doFracture", true);

}


void LagrangeDynamicsParallelPlateFlowExplicit_Proppant::RegisterFields( PhysicalDomainT& domain )
{
  //m_ppSolve->RegisterFields(domain);
  m_ldSolve.RegisterFields(domain);

  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("hydroForce", true, true);

  // TODO: remove the following and use interface logic
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  domain.m_feFaceManager.AddKeylessDataField<realT>("initialContactStress", true, true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("delta0N", true, false);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("gapShear0", true, false);
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("stressShear0", true, false);
  //////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  domain.m_feFaceManager.AddKeylessDataField<R1Tensor>("contactStress", false, true);
  if (m_faceStrengthRandomFactor > 1e-99)
    domain.m_feFaceManager.AddKeylessDataField<realT>("faceStrengthRandomFactor", true, true);
  /*if (m_ppSolve.m_leakoffCoef > 0.0)
  {
    domain.m_feFaceManager.AddKeylessDataField<realT>("totalLeakedVolume", true, true);
    domain.m_feFaceManager.AddKeylessDataField<realT>("initialSaturatedTime", true, true);
  }*/

  // new - but should be in base
  domain.m_feFaceManager.AddKeylessDataField<realT>("appliedPressure", true, true);

}

void LagrangeDynamicsParallelPlateFlowExplicit_Proppant::Initialize(PhysicalDomainT& domain, SpatialPartition& partition)
{
  //m_ppSolve.Initialize(domain,partition);
  m_ldSolve.Initialize(domain,partition);

  // link to ppSolver
  std::map<std::string,SolverBase*>& solverMap = *m_solverMapPtr;
  SolverBase* solverPtr = solverMap[m_FlowSolverName];
  m_ppSolve = dynamic_cast<ParallelPlateProppantFlowSolverImplicit*>( solverPtr );
  if(m_ppSolve == 0) throw GPException("Error LagrangeDynamicsParallelPlateFlowExplicit_Proppant: Flow solver" + m_FlowSolverName + " not found.");

  // link to propSolver
  if(!m_ProppantSolverName.empty()){
    SolverBase* solverPtrB = solverMap[m_ProppantSolverName];
    m_propSolve = dynamic_cast<ParallelPlateProppantSolver*>( solverPtrB );
    if(m_propSolve == 0) throw GPException("Error LagrangeDynamicsParallelPlateFlowExplicit_Proppant: Proppant solver" + m_ProppantSolverName + " not found.");
  } else {
	m_propSolve = NULL;
  }

  rArray1d& initialContactStress = domain.m_feFaceManager.GetFieldData<realT>("initialContactStress");
  rArray1d& delta0N = domain.m_feFaceManager.GetFieldData<realT>("delta0N");
  rArray1d* faceStrengthRandomFactor = domain.m_feFaceManager.GetFieldDataPointer<realT>("faceStrengthRandomFactor");
  delta0N = m_kJn;

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.m_numFaces ; ++kf )
  {
    if( domain.m_feFaceManager.m_toElementsRelation[kf].size() > 1 )
    {
      R1Tensor fc;
      R1Tensor fn;

      R2SymTensor stress0;
      R2SymTensor stress1;

      R1Tensor t0, t1;
      R1Tensor temp;

      ElementRegionT& er0 = *(domain.m_feFaceManager.m_toElementsRelation[kf][0].first);
      ElementRegionT& er1 = *(domain.m_feFaceManager.m_toElementsRelation[kf][1].first);
      const localIndex elemIndex0 = domain.m_feFaceManager.m_toElementsRelation[kf][0].second;
      const localIndex elemIndex1 = domain.m_feFaceManager.m_toElementsRelation[kf][1].second;

      er0.m_mat->StateData(elemIndex0,0)->TotalStress(stress0);
      er1.m_mat->StateData(elemIndex1,0)->TotalStress(stress1);

      // normal from away from element 0, into element 1
      domain.m_feFaceManager.FaceCenterAndNormal( domain.m_feNodeManager, kf, fc, fn );

      t0.AijBj(stress0,fn);
      t1.AijBj(stress1,fn);

      realT t0n = Dot(t0,fn);
      realT t1n = Dot(t1,fn);

      // TODO: Add normal stress initialization routine for interfaces and call here
      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////
      initialContactStress[kf] = 0.5 * fabs( t0n + t1n );
      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////
    }
  }
  if (faceStrengthRandomFactor!=NULL)
  {
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.m_numFaces ; ++kf )
    {
      localIndex iNd = domain.m_feFaceManager.m_toNodesRelation[kf][0];
      realT xFace = (*domain.m_feNodeManager.m_refposition)[iNd][0];
      srand(int(xFace/0.001));
      (*faceStrengthRandomFactor)[kf] = 1.0 + ( (rand()*1.0) / RAND_MAX - 0.5) * m_faceStrengthRandomFactor;
    }
  }
/*
  if (m_ppSolve.m_leakoffCoef > 0.0)
  {
    rArray1d& totalLeakedVolume = domain.m_feFaceManager.GetFieldData<realT>("totalLeakedVolume");
    totalLeakedVolume = 0.0;
    rArray1d& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");
    initialSaturatedTime = std::numeric_limits<realT>::max();
  }
*/

}

void LagrangeDynamicsParallelPlateFlowExplicit_Proppant::InitializeCommunications( PartitionBase& partition )
{

  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("tracer");
  partition.SetBufferSizes( syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver);
}


void LagrangeDynamicsParallelPlateFlowExplicit_Proppant::ApplyForcesFromContact(PhysicalDomainT& domain,
                                                                       StableTimeStep& timeStep,
                                                                       const realT dt )
{
  //-----------------------------
  //external faces
  //-----------------------------
  {
    Array1dT<R1Tensor>& contactForce = domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce>();
    contactForce = 0.0;
    Array1dT<R1Tensor>& decontactForce = domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::contactForce>();
    decontactForce = 0.0;

    //update nodal positions, velocities, and accelerations before updating face geometry
    //also, reset rotational and translational accelerations as well as forces and moments
    domain.m_discreteElementManager.UpdateNodalStatesZeroForcesAndAccelerations();

    //update face geometry and sort faces if necessary
    const bool resort = domain.m_externalFaces.RecalculateNeighborList(domain.m_feNodeManager,
                                                                       domain.m_discreteElementSurfaceNodes,
                                                                       domain.m_discreteElementManager,
                                                                       dt);

    //if a resort has been triggered, then you also need to update the contact manager
    if (resort)
      domain.m_contactManager.Update(domain.m_externalFaces.m_neighborList);

    {
      Array1dT<Array1dT<R1Tensor> > xs;
      xs.resize(domain.m_externalFaces.DataLengths());
      domain.m_externalFaces.UpdateGeometricContactProperties(dt, domain, xs);
#ifdef STATES_ON_CONTACTS
      domain.m_externalFaces.UpdateAndApplyContactForcesFromStatesOnCommonPlanes(dt, domain);
#else
      domain.m_externalFaces.UpdateAndApplyContactStresses(timeStep, dt, domain, xs);
#endif
    }

    //for parallel: do NOT need to synchronize nodal force fields across processes for DE nodes
    //before moving to centroid, since each process will do that calculation (redundantly) itself
    // ... there's therefore no need for explicit synchrony
    //as long as nodal states remain synchronous, everything will be fine!
#ifndef DEEFC
    domain.m_discreteElementManager.ApplyDiscreteElementContactForces(45, 0.2);
#else
    domain.m_discreteElementManager.ApplyNodalForces();
#endif
  }
}




double LagrangeDynamicsParallelPlateFlowExplicit_Proppant::TimeStep( const realT& time,
                                                          const realT& dt,
                                                          const int cycleNumber,
                                                          PhysicalDomainT& domain,
                                                          const sArray1d& namesOfSolverRegions ,
                                                          SpatialPartition& partition ,
                                                          FractunatorBase* const fractunator )
{

  if( fractunator!=NULL )
  {
    if( cycleNumber%fractunator->m_checkInterval==0 && m_doFracture )
    {
      fractunator->SeparationDriver( domain.m_feNodeManager,
                                     domain.m_feEdgeManager,
                                     domain.m_feFaceManager,
                                     domain.m_externalFaces,
                                     domain.m_feElementManager,
                                     partition, false , time);
    }
  }

  bool doParallelPlateStep = ( (time + dt) >= m_newPPStepTime);

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();
  m_ppSolve->m_stabledt.m_maxdt = 1.0; //0.9*std::numeric_limits<double>::max();

  Array1dT<R1Tensor>& hgforce = domain.m_feNodeManager.GetFieldData<FieldInfo::hgforce> ();
  hgforce = 0.0;

  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
//  Array1dT<lArray1d>& edgeToFlowFaces = domain.m_edgeManager.GetVariableOneToManyMap("edgeToFlowFaces");

  Array1dT<R1Tensor>* u0 = domain.m_feNodeManager.GetFieldDataPointer<R1Tensor>("displacement0");
  Array1dT<R1Tensor>* unet = domain.m_feNodeManager.GetFieldDataPointer<R1Tensor>("netDisplacement");


  // Parallel plate solve
  if(doParallelPlateStep){
	// get new flow surface
    m_ppSolve->DefineFlowSets( domain );
  }

  if(doParallelPlateStep){
    m_ppSolve->GenerateParallelPlateGeometricQuantities( domain,time,m_pp_dt );
    //new
    m_ppSolve->GenerateSlurryParameters( domain );
  }



  if(doParallelPlateStep){

	    if (partition.m_rank == 0 && false){
	      std::cout << "Updating parallel plate solver: t " << m_newPPStepTime << " dt " <<m_pp_dt << std::endl;
	    }
    m_ppSolve->CalculateMassRate(domain, partition, time,m_pp_dt );
  }




  // update nodes
  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart1( domain.m_feNodeManager, time, dt );

  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator elementRegionIter = domain.m_feElementManager.m_ElementRegions.begin() ;
       elementRegionIter != domain.m_feElementManager.m_ElementRegions.end() ;
       ++elementRegionIter )
  {
    ElementRegionT& elementRegion = elementRegionIter->second;

    elementRegion.CalculateVelocityGradients(domain.m_feNodeManager);

    elementRegion.MaterialUpdate(dt);

    elementRegion.CalculateNodalForces(domain.m_feNodeManager, m_stabledt, dt);
  }


  BoundaryConditionFunctions::ApplyTractionBoundaryCondition( domain ,time);

  // -- CONTACT --
  if (domain.m_externalFaces.m_contactActive)
    ApplyForcesFromContact(domain, this->m_stabledt, dt);

  m_stabledt.m_maxdt *= m_ldSolve.m_courant;



  // Update parallel plate geometric quantities

  if(doParallelPlateStep){
    m_ppSolve->GenerateParallelPlateGeometricQuantities( domain,0,0 );
    // new
    m_ppSolve->GenerateSlurryParameters( domain );
  }


  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;

    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("tracer");
    partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver );
  }

  m_ppSolve->UpdateEOS( domain, time,dt );
  if(doParallelPlateStep){
    // new
	//    m_ppSolve->UpdateEOS( domain, time,m_pp_dt, true ); // used to update mass only in pp step
    m_ppSolve->UpdateFlux( m_oldPPStepTime,m_pp_dt, domain, partition);
    m_ppSolve->OverwriteOldGeometricQuantities(domain);
  } else {
	//	m_ppSolve->UpdateEOS( domain, time,dt, false );

  }



  // update applied pressure  // new
  rArray1d& appliedPressure = domain.m_feFaceManager.GetFieldData<realT>("appliedPressure");
  const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 0 )
    {
        appliedPressure[kf] = faceFluidPressure[kf];
      //  appliedPressure[kf] += ( faceFluidPressure[kf] - appliedPressure[kf] ) / 0.5 * dt;
        if( appliedPressure[kf] < 0.0 ) appliedPressure[kf] = 0.0;
    }
  }


  // apply fluid pressure to faces
  const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
//  const rArray1d& faceFluidPressure = domain.m_faceManager.GetFieldData<FieldInfo::pressure>();
  //const rArray1d& appliedPressure = domain.m_feFaceManager.GetFieldData<realT>("appliedPressure");
  const rArray1d& initialContactStress = domain.m_feFaceManager.GetFieldData<realT>("initialContactStress");
  Array1dT<R1Tensor>& nodalForce = domain.m_feNodeManager.GetFieldData<FieldInfo::force>();
  Array1dT<R1Tensor>& hydroForce = domain.m_feNodeManager.GetFieldData<R1Tensor>("hydroForce");
  Array1dT<R1Tensor>& contactForce = domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce>();
  Array1dT<R1Tensor>& contactStress = domain.m_feFaceManager.GetFieldData<R1Tensor>("contactStress");

//  Array1dT<R1Tensor>& cohesiveForce = domain.m_feNodeManager.GetFieldData<R1Tensor>("cohesiveForce");

//  const Array1dT<R1Tensor>& cohesiveTraction = domain.m_feFaceManager.GetFieldData<R1Tensor>("cohesiveTraction");

  hydroForce = 0.0;
  contactForce = 0.0;
//  cohesiveForce = 0.0;

  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
    if( flowFaceType[kf] == 0 && childFaceIndex[kf].size() >= 1 )
    {

      realT pressure = appliedPressure[kf];

      if( pressure < initialContactStress[kf] )
      {
        pressure = initialContactStress[kf];
      }

      localIndex faceIndex[2];
      if (childFaceIndex[kf].size() == 1)
      {
        faceIndex[0] = kf;
        faceIndex[1] = childFaceIndex[kf][0];
      }
      else
      {
        faceIndex[0] = childFaceIndex[kf][0];
        faceIndex[1] = childFaceIndex[kf][1];
      }

      realT stressPen;
      R1Tensor stressShear;
      CalculateContactStress(domain, time, dt, kf, faceIndex, stressPen, stressShear);


      //////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////

      const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                              domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

      for( int side=0 ; side<2 ; ++side )
      {
        const localIndex faceID = faceIndex[side];
        const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, faceID );

        R1Tensor hForce = N[side];
        R1Tensor pForce = N[side];
        R1Tensor cForce, sForce;
        int direction = -1;
        if(side==0) direction = 1;

        hForce *= -area * pressure / domain.m_feFaceManager.m_toNodesRelation[faceID].size();
        pForce *= -area * stressPen / domain.m_feFaceManager.m_toNodesRelation[faceID].size();
        sForce = stressShear;
        contactStress[faceID] = pForce / area * domain.m_feFaceManager.m_toNodesRelation[faceID].size();
        contactStress[faceID] += sForce * direction;
        sForce *= direction * area / domain.m_feFaceManager.m_toNodesRelation[faceID].size();

//        cForce = cohesiveTraction[kf];
        cForce *= direction * area / domain.m_feFaceManager.m_toNodesRelation[faceID].size();


        for( lArray1d::const_iterator nodeID=domain.m_feFaceManager.m_toNodesRelation[faceID].begin() ;
            nodeID!=domain.m_feFaceManager.m_toNodesRelation[faceID].end() ; ++nodeID )
        {
          nodalForce[*nodeID] += hForce;
          nodalForce[*nodeID] += pForce;
          nodalForce[*nodeID] += cForce;
          nodalForce[*nodeID] += sForce;

          hydroForce[*nodeID] += hForce;
          contactForce[*nodeID] += pForce;
          contactForce[*nodeID] += sForce;

//          cohesiveForce[*nodeID] += cForce;
        }

      }
    }

  }

//  m_ldSolve.ApplyGapDamping( domain.m_nodeManager,
 //                            domain.m_faceManager, dt );







  LagrangeExplicitDynamicsFunctions::LinearPointUpdatePart2( domain.m_feNodeManager, time, dt, m_ldSolve.m_dampingM );

  if (unet != NULL && u0 != NULL)
  {
    for (localIndex i = 0; i < domain.m_feNodeManager.DataLengths(); ++i)
    {
      (*unet)[i] = (*domain.m_feNodeManager.m_displacement)[i];
      (*unet)[i] -= (*u0)[i];
    }
  }

  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;

    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("tracer");
    {
      rArray1d* initialSaturatedTime = domain.m_feFaceManager.GetFieldDataPointer<realT>("initialSaturatedTime");
      if (initialSaturatedTime != NULL )
        syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("initialSaturatedTime");
    }


    partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver );
  }


  // Advance Proppant
  if(doParallelPlateStep && m_propSolve){
	  m_propSolve->TimeStep(m_oldPPStepTime, m_pp_dt,cycleNumber, domain, namesOfSolverRegions, partition, fractunator);
  }

  // Update flow solver times
  if(doParallelPlateStep){
	    m_oldPPStepTime = m_newPPStepTime;
	    m_pp_dt = std::min( m_ppSolve->m_stabledt.m_maxdt, m_nSubsteps*dt);
	    //m_pp_dt = std::min(m_pp_dt, m_propSolve->m_stabledt.m_maxdt);
	    m_pp_dt = SetGlobalMaxTimeStep( m_pp_dt );
	    m_newPPStepTime +=  m_pp_dt;
  }


  // Update timestep
  if(  m_stabledt.m_maxdt > m_ppSolve->m_stabledt.m_maxdt )
    m_stabledt.m_maxdt = m_ppSolve->m_stabledt.m_maxdt;
 // if(  m_stabledt.m_maxdt > m_propSolve->m_stabledt.m_maxdt )
 //	    m_stabledt.m_maxdt = m_propSolve->m_stabledt.m_maxdt;

  // if (cycleNumber%10000 == 0)
  if (cycleNumber%1000 == 0)
   {
     int rank ;
     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
     if(rank == 0)
       std::cout << "t=" << time << "rank: " << rank << "Solid dt: " << m_stabledt.m_maxdt<<", Flow dt: "<<m_ppSolve->m_stabledt.m_maxdt << std::endl;
   }


  return dt;
}

realT LagrangeDynamicsParallelPlateFlowExplicit_Proppant::SetGlobalMaxTimeStep(realT local_dt){
	  realT temp = local_dt;
	  realT nextdt;
	  MPI_Allreduce (&temp,&nextdt,1,MPI_DOUBLE,MPI_MIN ,MPI_COMM_WORLD);

	  return nextdt;
}


void LagrangeDynamicsParallelPlateFlowExplicit_Proppant::CalculateContactStress(PhysicalDomainT& domain,
                                                                       const realT& time,
                                                                       const realT& dt,
                                                                       localIndex& kf,
                                                                       const localIndex faceIndex[],
                                                                       realT& stressPen,
                                                                       R1Tensor& stressShear)
{
  stressPen = 0;
  stressShear = 0;
  rArray1d& delta0N = domain.m_feFaceManager.GetFieldData<realT>("delta0N");
  Array1dT<R1Tensor>& gapShear0 = domain.m_feFaceManager.GetFieldData<R1Tensor>("gapShear0");
  Array1dT<R1Tensor>& stressShear0 = domain.m_feFaceManager.GetFieldData<R1Tensor>("stressShear0");
  rArray1d& effectiveStressN = domain.m_feFaceManager.GetFieldData<realT>("effectiveStressN");


  const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                          domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

  R1Tensor T[2];
  domain.m_feFaceManager.FaceTangential(domain.m_feNodeManager, faceIndex[0], T[0], T[1]);
  R1Tensor Nbar = N[0];
  Nbar -= N[1];
  Nbar.Normalize();

  const R1Tensor gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, kf );
  const realT gapNormal = Dot(gap, Nbar);
  R1Tensor gapShear = Dot(gap, T[0]) * T[0] + Dot(gap, T[1]) * T[1];
  R1Tensor gapShearInc(gapShear);
  gapShearInc -= gapShear0[kf];


  realT deltaN ;
  realT tBuffer = dt * 500;
  if (time == 0)
  {
    deltaN = delta0N[kf];
  }
  else if (time < tBuffer )
  {
    deltaN = delta0N[kf] * (1.0 - time / tBuffer * m_fLockedInSIF);
  }
  else deltaN = (1-m_fLockedInSIF) * delta0N[kf];

  if( -gapNormal + deltaN <= 0.0 )
  {
    stressPen = 0.0;
    stressShear = 0.0;
  }
  else
  {
    stressPen = (-gapNormal + deltaN) * m_kJn ;

    R1Tensor stressShearPrj = stressShear0[kf];
    stressShearPrj += gapShearInc * m_kJs;

    if ( pow(stressShearPrj.L2_Norm(),2) < pow(stressPen*m_COFJ, 2) )
    {
      stressShear = stressShearPrj;
    }
    else  // We have to calculate plasticity stuff
    {
      realT du[2], t[2], x[2];
      du[0] = Dot(gapShearInc, T[0]);
      du[1] = Dot(gapShearInc, T[1]);
      t[0] = Dot(stressShear0[kf], T[0]);
      t[1] = Dot(stressShear0[kf], T[1]);

      realT a, b, c;
      a = pow(t[0],2) + pow(t[1],2);
      if (a == 0) a = pow( t[0]+du[0]*m_kJs ,2) + pow(t[1]+ du[1]*m_kJs,2);

      b = -2 * ( t[0]*(t[0] + m_kJs * du[0]) + t[1]*(t[1] + m_kJs * du[1]));
      c = pow(t[0] + m_kJs * du[0], 2) + pow(t[1] + m_kJs * du[1], 2) - pow(stressPen*m_COFJ, 2);

      realT b2_4ac = pow(b,2)- 4 * a * c;

      if (b2_4ac < 0.0 || a==0)
      {
        t[0] = 0;
        t[1] = 0;
      }
      else
      {
        x[0] = (-b - pow(b2_4ac, 0.5))/2/a;
        x[1] = (-b + pow(b2_4ac, 0.5))/2/a;

        if (fabs(x[0]) > fabs(x[1])) x[0] = x[1];

        t[0] += m_kJs * du[0] - x[0] * t[0];
        t[1] += m_kJs * du[1] - x[0] * t[1];
      }
      stressShear = t[0] * T[0] + t[1] * T[1];

    }
  }
  effectiveStressN[kf] = stressPen;
  stressShear0[kf] = stressShear;
  gapShear0[kf] = gapShear;
}



void LagrangeDynamicsParallelPlateFlowExplicit_Proppant::PostProcess (PhysicalDomainT& domain,
                                                             SpatialPartition& partition,
                                                             const sArray1d& namesOfSolverRegions)
{
  m_ppSolve->PostProcess(domain, partition, namesOfSolverRegions);
  m_ldSolve.PostProcess(domain, partition, namesOfSolverRegions);
}

/// Register solver in the solver factory
REGISTER_SOLVER( LagrangeDynamicsParallelPlateFlowExplicit_Proppant )

