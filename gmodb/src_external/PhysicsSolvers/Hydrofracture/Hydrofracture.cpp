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
 * @file Hydrofracture.cpp
 * @author settgast1 
 * @date Feb 28, 2011
 */







#include "Hydrofracture.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "PhysicsSolvers/PhysicsSolverStrings.h"
#include "Utilities/StringUtilities.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "SurfaceGeneration/Fractunator2.h"
#include "SurfaceGeneration/Fractunator2D.h"
#include "SurfaceGeneration/Fractunator3.h"
#include "PhysicsSolvers/Lagrange/LagrangeHelperFunctions.h"

#include "PhysicsSolvers/FractureFlow/TwoDADRSolver_Proppant.h"


Hydrofracture::Hydrofracture(  const std::string& name,
                               ProblemManagerT* const pm ):
SolverBase(name,pm),
m_ldSolve(NULL),
m_ppSolve(NULL),
m_mfSolve(NULL),
epetra_comm(pm->m_epetraComm),
m_solvers(pm->m_solvers)
{

//  m_ppSolve = new ParallelPlateFlowSolver("hydrofracpp",pm);
//  m_ldSolve = new LagrangeLargeStrain("hydrofracld",pm);
  /* Empty */
}

Hydrofracture::~Hydrofracture()
{
}

void Hydrofracture::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML(hdn);

  /*
  m_enableTimeIntegrationOption[Implicit]= hdn->GetAttributeOrDefault<int> ("Implicit", 0);
  m_enableTimeIntegrationOption[Explicit]= hdn->GetAttributeOrDefault<int> ("Explicit", 0);

  {
    int temp = hdn->GetAttributeOrDefault<int> ("timeIntegrationOption", 0);
    m_timeIntegrationOption = IntToTimeIntegration( temp );
    m_enableTimeIntegrationOption[m_timeIntegrationOption] = 1;
  }
  */


  m_ppSolverName = hdn->GetAttributeString("ppSolverName");
  m_lgSolverName = hdn->GetAttributeString("lgSolverName");
  m_mfSolverName = hdn->GetAttributeString("mfSolverName");
  m_propSolverName = hdn->GetAttributeString("propSolverName");

  if( m_ppSolverName.empty() )
  {
    throw GPException( "Hydrofracture.cpp: No ppSolverName Specified");
  }

  if( m_lgSolverName.empty() )
  {
    throw GPException( "Hydrofracture.cpp: No lgSolverName Specified");
  }

  if( !m_mfSolverName.empty() )
  {
    std::string temp = hdn->GetAttributeString("dtMatrixFlow");
    if( temp.empty() )
    {
      throw GPException( "The time step size for matrix flow solver must be specified");
    }
    else
    {
      m_dtMatrixFlow = atof( temp.c_str() );
      m_tNextMatrixFlowUpdate = std::numeric_limits<realT>::min();
    }
  }



  if( !m_propSolverName.empty() )
  {
	  m_nSubsteps = hdn->GetAttributeOrDefault<realT>("substeps", 1);
  }


  std::string tempStr = hdn->GetAttributeString("timeIntegration");
  if( strIsInt(tempStr) ){
	int temp = fromString<int>(tempStr);
    m_timeIntegrationOption = IntToTimeIntegration( temp );
  } else {
	if( tempStr == "Implicit"){
        m_timeIntegrationOption =  Implicit;

	} else if(tempStr == "Explicit"){
  	    m_timeIntegrationOption = Explicit;

	} else {
        throw GPException("Invalid value input into Hydrofracture::IntToTimeIntegration()");
        m_timeIntegrationOption =  numTimeIntegrationEnums;

	}
  }

  // implicit
  if(m_timeIntegrationOption == Implicit)
  {
  }

  // explicit

  if(m_timeIntegrationOption == Explicit){

    m_kJn = hdn->GetAttributeOrDefault<realT>("normalJointStiffness", 1.0e10);
    m_kJs = hdn->GetAttributeOrDefault<realT>("shearJointStiffness", 1.0e10);

    m_fLockedInSIF =  hdn->GetAttributeOrDefault<realT>("lockedInSIFFactor", 1.0);

    m_faceStrengthRandomFactor = hdn->GetAttributeOrDefault<realT>("faceStrengthRandomFactor", 0.0);
    m_COFJ = fabs(hdn->GetAttributeOrDefault<realT>("COFJoint", 0.5));
    m_jointCohesion = fabs(hdn->GetAttributeOrDefault<realT>("jointCohesion", 0.0));
  }

}


void Hydrofracture::RegisterFields( PhysicalDomainT& domain )
{


  if( m_solvers.find(m_ppSolverName) != m_solvers.end() )
  {
    m_ppSolve = dynamic_cast<ParallelPlateFlowSolverBase*>(stlMapLookup( m_solvers, m_ppSolverName, "m_solvers" ));
  }
  else
  {
    throw GPException( "Hydrofracture.cpp: No ppSolverName Specified");
  }

  if( m_solvers.find(m_lgSolverName) != m_solvers.end() )
  {
    m_ldSolve = dynamic_cast<LagrangeSolverBase*>(stlMapLookup( m_solvers, m_lgSolverName, "m_solvers" ));
  }
  else
  {
    throw GPException( "Hydrofracture.cpp: No lgSolverName Specified");
  }

  if (!m_mfSolverName.empty())
  {
    if( m_solvers.find(m_mfSolverName) != m_solvers.end() )
    {
      m_mfSolve = dynamic_cast<MatrixFlowSolver*>(stlMapLookup( m_solvers, m_mfSolverName, "m_solvers" ));
    }
    else
    {
      throw GPException( "Hydrofracture.cpp: Could not find the mfSolver requested.");
    }
  } else {
	  m_mfSolve = NULL;
  }

  if (!m_propSolverName.empty())
  {
    if( m_solvers.find(m_propSolverName) != m_solvers.end() )
    {
    	m_propSolve = dynamic_cast<ParallelPlateProppantSolver*>(stlMapLookup( m_solvers, m_propSolverName, "m_solvers" ));
    }
    else
    {
      throw GPException( "Hydrofracture.cpp: Could not find the propSolver requested.");
    }
  } else {
	  m_propSolve = NULL;  // I didn't think this was necessary - but apparently the optimized version will not set the pointer to NULL by default.
  }




  // universal fields
  domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("hydroForce", true, true);


  // implicit
  if(m_timeIntegrationOption == Implicit){

    domain.m_feNodeManager.AddKeyedDataField<FieldInfo::hgforce>();
    domain.m_feNodeManager.AddKeylessDataField<R1Tensor>("effectiveNormal",true,true);

    domain.m_feFaceManager.AddKeylessDataField<realT>("mass0",false,false);
    domain.m_feFaceManager.AddKeylessDataField<realT>("p0",false,false);
    domain.m_feFaceManager.AddKeylessDataField<realT>("MassRate",true,true);

  }

  // explicit
  if(m_timeIntegrationOption == Explicit){

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


    if (m_ppSolve->m_leakoffCoef > 0.0)
    {
      domain.m_feFaceManager.AddKeylessDataField<realT>("totalLeakedVolume", true, true);
      domain.m_feFaceManager.AddKeylessDataField<realT>("initialSaturatedTime", true, true);
    }

    if (!m_mfSolverName.empty())
    {
      domain.m_feFaceManager.AddKeylessDataField<realT>("totalLeakedVolume", true, true);
    }


    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_x", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_y", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_z", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_xy", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_yz", false, true);
    domain.m_feNodeManager.AddKeylessDataField<realT>("nsigma_zx", false, true);

    rArray1d* COFJPointer = domain.m_feFaceManager.GetFieldDataPointer<realT>("COFJ");
    if (COFJPointer == NULL)
    {
      domain.m_feFaceManager.AddKeylessDataField<realT>("COFJ", true, false);
      m_COFJSetByInitialCondition = false;
    }
    else
    {
      m_COFJSetByInitialCondition = true;
    }

    rArray1d* jointCohesionPointer = domain.m_feFaceManager.GetFieldDataPointer<realT>("jointCohesion");
    if (jointCohesionPointer == NULL)
    {
      domain.m_feFaceManager.AddKeylessDataField<realT>("jointCohesion", true, false);
      m_jointCohesionSetByInitialCondition = false;
    }
    else
    {
      m_jointCohesionSetByInitialCondition = true;
    }
  }


}

void Hydrofracture::Initialize(PhysicalDomainT& domain, SpatialPartition& partition)
{





  if( m_timeIntegrationOption==Implicit){

    m_ldSolve->m_timeIntegrationOption = m_ldSolve->QuasiStatic;
    m_ldSolve->m_enableTimeIntegrationOption[m_timeIntegrationOption] = 1;

  }

  if(m_timeIntegrationOption==Explicit){

    m_ldSolve->m_timeIntegrationOption = m_ldSolve->ExplicitDynamic;
    m_ldSolve->m_enableTimeIntegrationOption[m_timeIntegrationOption] = 1;

    // explicit
    rArray1d& initialContactStress = domain.m_feFaceManager.GetFieldData<realT>("initialContactStress");
    rArray1d& delta0N = domain.m_feFaceManager.GetFieldData<realT>("delta0N");

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
    for( localIndex kf=0 ; kf<domain.m_feFaceManager.m_numFaces ; ++kf )
    {
      delta0N[kf] = m_kJn;
    }

    if (m_ppSolve->m_leakoffCoef > 0.0 && !m_mfSolverName.empty())
      throw GPException( "Both Carter's leakoff coefficient and porous media solver are defined.  You can only have one!");


    if (m_ppSolve->m_leakoffCoef > 0.0)  //If mf solver is active, totalLeakedVolume is initialized there.
    {
      rArray1d& totalLeakedVolume = domain.m_feFaceManager.GetFieldData<realT>("totalLeakedVolume");
      totalLeakedVolume = 0.0;
      rArray1d& initialSaturatedTime = domain.m_feFaceManager.GetFieldData<realT>("initialSaturatedTime");
      initialSaturatedTime = std::numeric_limits<realT>::max();
    }


    if ( !m_COFJSetByInitialCondition)
    {
      rArray1d& COFJ = domain.m_feFaceManager.GetFieldData<realT>("COFJ");
      COFJ = m_COFJ;
    }

    if (!m_jointCohesionSetByInitialCondition)
    {
      rArray1d& jointCohesion = domain.m_feFaceManager.GetFieldData<realT>("jointCohesion");
      jointCohesion = m_jointCohesion;
    }

  }
}

void Hydrofracture::InitializeCommunications( PartitionBase& partition )
{
  if( m_timeIntegrationOption==Implicit){
    m_ppSolve->InitializeCommunications( partition );
    m_ldSolve->InitializeCommunications( partition );

    m_syncedFields.clear();

    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::displacement>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_ldSolve->m_trilinosIndexStr);

    m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("Mass");
    m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(m_ppSolve->TrilinosIndexString());
    m_syncedFields[PhysicalDomainT::FiniteElementEdgeManager].push_back(PS_STR::VolumetricFluxStr);

    partition.SetBufferSizes(m_syncedFields, CommRegistry::hydrofractureSolver);
  }

  if( m_timeIntegrationOption==Explicit){

    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
    m_syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("tracer");

    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_x");
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_y");
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_z");
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_xy");
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_yz");
    m_syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_zx");


    partition.SetBufferSizes(m_syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver);
    // fixme prob should be hydro solver comm, but need to check if is calling sync within PP solver.
  }
}




double Hydrofracture::TimeStep( const realT& time,
                                      const realT& dt,
                                      const int cycleNumber,
                                      PhysicalDomainT& domain,
                                      const sArray1d& namesOfSolverRegions,
                                      SpatialPartition& partition,
                                      FractunatorBase* const fractunator )
{
 // std::cout<<"m_timeIntegrationOption = "<<m_timeIntegrationOption<<std::endl;
  if( m_timeIntegrationOption==Explicit)
  {
    TimeStepExplicit( time, dt, cycleNumber, domain, namesOfSolverRegions,
                      partition, fractunator);
  }
  else if( m_timeIntegrationOption==Implicit)
  {
    TimeStepImplicit( time, dt, cycleNumber, domain, namesOfSolverRegions,
                      partition, fractunator);

  }

  return dt;
}

void Hydrofracture::TimeStepExplicit( const realT& time,
                                      const realT& dt,
                                      const int cycleNumber,
                                      PhysicalDomainT& domain,
                                      const sArray1d& namesOfSolverRegions,
                                      SpatialPartition& partition,
                                      FractunatorBase* const fractunator )
{
  // std::cout<<"Explicit Time Step"<<std::endl;

  //////////////////////////////////////////////


  if( fractunator!=NULL )
  {
    if( cycleNumber%fractunator->m_checkInterval==0 )
    {
      fractunator->SeparationDriver( domain.m_feNodeManager,
                                     domain.m_feEdgeManager,
                                     domain.m_feFaceManager,
                                     domain.m_externalFaces,
                                     domain.m_feElementManager,
                                     partition, false , time);
    }
  }

  bool doParallelPlateStep = true;
  bool doProppantStep = false;
  if(m_propSolve){
	  doParallelPlateStep = ( (time + dt) >= m_newPPStepTime);
	  doProppantStep = doParallelPlateStep;
  }

  m_stabledt.m_maxdt = std::numeric_limits<double>::max();
  Array1dT<R1Tensor>& hgforce = domain.m_feNodeManager.GetFieldData<FieldInfo::hgforce> ();
  hgforce = 0.0;

  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
//  Array1dT<lArray1d>& edgeToFlowFaces = domain.m_edgeManager.GetVariableOneToManyMap("edgeToFlowFaces");

  Array1dT<R1Tensor>* u0 = domain.m_feNodeManager.GetFieldDataPointer<R1Tensor>("displacement0");
  Array1dT<R1Tensor>* unet = domain.m_feNodeManager.GetFieldDataPointer<R1Tensor>("netDisplacement");

  if(doParallelPlateStep && doProppantStep){
	    m_ppSolve->DefineFlowSets( domain );
  }

  if(doParallelPlateStep){
    m_ppSolve->GenerateParallelPlateGeometricQuantities( domain, time, dt );
    if(doProppantStep){
    	m_ppSolve->GenerateSlurryParameters( domain );
    	m_ppSolve->CalculateMassRate(domain, partition, time,m_pp_dt );
    } else {
      m_ppSolve->CalculateAndApplyMassFlux( dt, domain ); // fixme m_ppSolve->CalculateMassRate(domain, partition, time,m_pp_dt );
    }
  }
 
  if (m_ppSolve->m_leakoffCoef > 0.0) m_ppSolve->CalculateCarterLeakOff(time, dt, domain);

  if (!m_mfSolverName.empty()) m_ppSolve->CalculateMatrixFlowLeakOff(time, dt, domain);

  m_ppSolve->ApplyFluxBoundaryCondition(time, dt, cycleNumber, partition.m_rank, domain);

  // update nodes
  LagrangeHelperFunctions::LinearPointUpdatePart1( domain.m_feNodeManager, time, dt );
  domain.m_feNodeManager.ZeroDetachedNodeVelocity();

  m_ldSolve->ProcessElementRegions(domain.m_feNodeManager,
                                  domain.m_feElementManager,
                                  namesOfSolverRegions, dt);
  m_stabledt.m_maxdt = m_ldSolve->m_stabledt.m_maxdt;


  BoundaryConditionFunctions::ApplyTractionBoundaryCondition( domain ,time);

  // -- CONTACT --
  if (domain.m_externalFaces.m_contactActive)
    m_ldSolve->ApplyForcesFromContact(domain, m_stabledt, dt);

  m_stabledt.m_maxdt *= m_ldSolve->m_courant;

  if (cycleNumber%10000 == 0)
  {
    int rank ;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::cout << "t=" << time << "rank: " << rank << "Solid dt: " << m_stabledt.m_maxdt<<", Flow dt: "<<m_ppSolve->m_stabledt.m_maxdt << std::endl;
  }

//  std::cout<<m_stabledt.m_maxdt<<", "<<m_ppSolve->m_stabledt.m_maxdt;

  if(  m_stabledt.m_maxdt > m_ppSolve->m_stabledt.m_maxdt )
    m_stabledt.m_maxdt = m_ppSolve->m_stabledt.m_maxdt;



  // -- CONTACT --
//  if (domain.m_externalFaces.m_contactActive)
//    m_ldSolve->ApplyForcesFromContact(domain, dt);


  if(doParallelPlateStep){
    m_ppSolve->GenerateParallelPlateGeometricQuantities( domain,0,0);
    if(doProppantStep){
      m_ppSolve->GenerateSlurryParameters( domain );
    }
  }

  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;

    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::acceleration>::Name());
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::velocity>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(Field<FieldInfo::mass>::Name());
    syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("tracer");
    partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver );
  }

  m_ppSolve->UpdateEOS( time, dt, domain  );
  if(doProppantStep){
    m_ppSolve->UpdateFlux( m_oldPPStepTime, dt, domain, partition);
    m_ppSolve->OverwriteOldGeometricQuantities(domain);
  }

  if(m_propSolve){
	  // update applied pressure  // new
	  rArray1d& appliedPressure = domain.m_feFaceManager.GetFieldData<realT>("appliedPressure");
	  const rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
	  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
	  {
	    if( flowFaceType[kf] == 0 )
	    {
	        appliedPressure[kf] = faceFluidPressure[kf];
	        if( appliedPressure[kf] < 0.0 ) appliedPressure[kf] = 0.0;
	    }
	  }
  }


  if (!m_mfSolverName.empty() && (time + 1.01 * dt) >= m_tNextMatrixFlowUpdate)
  {
    std::cout << "I am updating matrix pressure at " << time + dt << std::endl;

    m_ppSolve->CalculateNodalPressure(domain, partition);

    m_mfSolve->TimeStep(time, dt, cycleNumber, domain,namesOfSolverRegions,partition, fractunator);

    m_mfSolve->CalculateLeakoff(m_dtMatrixFlow, domain);

    m_tNextMatrixFlowUpdate += m_dtMatrixFlow;

  }

  if ( !m_mfSolverName.empty() && time + dt + m_stabledt.m_maxdt >= m_tNextMatrixFlowUpdate)
  {
    //This only changes the local time step expectation.
    m_stabledt.m_maxdt = m_tNextMatrixFlowUpdate - time - dt;
    std::cout << "I am adjusting next time step to " << m_stabledt.m_maxdt << std::endl;
  }

  // apply fluid pressure to faces
  const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
//  const rArray1d& faceFluidPressure = domain.m_faceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& appliedPressure = domain.m_feFaceManager.GetFieldData<realT>("appliedPressure");
  const rArray1d& initialContactStress = domain.m_feFaceManager.GetFieldData<realT>("initialContactStress");
  Array1dT<R1Tensor>& nodalForce = domain.m_feNodeManager.GetFieldData<FieldInfo::force>();
  Array1dT<R1Tensor>& hydroForce = domain.m_feNodeManager.GetFieldData<R1Tensor>("hydroForce");
  Array1dT<R1Tensor>& contactForce = domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce>();
  Array1dT<R1Tensor>& contactStress = domain.m_feFaceManager.GetFieldData<R1Tensor>("contactStress");
  rArray1d& faceArea = domain.m_feFaceManager.GetFieldData<realT>("faceArea");

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


      const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                              domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};


      // Apply fluid pressure to solid surfaces
      for( int side=0 ; side<2 ; ++side )
      {
        const localIndex faceID = faceIndex[side];
        const realT area = faceArea[faceID];

        R1Tensor hForce = N[side];

        hForce *= -area * pressure / domain.m_feFaceManager.m_toNodesRelation[faceID].size();


        for( lArray1d::const_iterator nodeID=domain.m_feFaceManager.m_toNodesRelation[faceID].begin() ;
            nodeID!=domain.m_feFaceManager.m_toNodesRelation[faceID].end() ; ++nodeID )
        {
          nodalForce[*nodeID] += hForce;
          hydroForce[*nodeID] += hForce;
        }
      }

      // The following is how we hacked contact in the hydrofrac solver initially.
      // It is replaced by the proper implementation of the contact manager
      // But I am keeping this here as a fallback option until the contact part is more stable.
      if (!domain.m_externalFaces.m_contactActive)
      {
        realT stressPen;
        R1Tensor stressShear;
        CalculateContactStress(domain, time, dt, kf, faceIndex, stressPen, stressShear);


        //////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////

        for( int side=0 ; side<2 ; ++side )
        {
          int direction = -1;
          if(side==0) direction = 1;
          const localIndex faceID = faceIndex[side];
          const realT area = faceArea[faceID];

          R1Tensor pForce = N[side];  // normal contact
          R1Tensor  sForce;  //shear contact

//          if (childFaceIndex[kf].size() == 1)
//          {
//            direction = -1;
//            if(side==0) direction = 1;
//          }
//          else
//          {
//            R1Tensor n0, n1;
//            n0 = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, kf, true);
//            n1 = domain.m_feFaceManager.FaceNormal(domain.m_feNodeManager, faceID, true);
//            if (Dot(n0, n1) > 0.0)
//            {
//              direction = 1;
//            }
//            else
//            {
//              direction = -1;
//            }
//
//          }


          pForce *= -area * stressPen / domain.m_feFaceManager.m_toNodesRelation[faceID].size();
          sForce = stressShear;
          contactStress[faceID] = pForce / area * domain.m_feFaceManager.m_toNodesRelation[faceID].size();
          contactStress[faceID] += sForce * direction;
          sForce *= direction * area / domain.m_feFaceManager.m_toNodesRelation[faceID].size();

          for( lArray1d::const_iterator nodeID=domain.m_feFaceManager.m_toNodesRelation[faceID].begin() ;
              nodeID!=domain.m_feFaceManager.m_toNodesRelation[faceID].end() ; ++nodeID )
          {
            nodalForce[*nodeID] += pForce;
            nodalForce[*nodeID] += sForce;

            contactForce[*nodeID] += pForce;
            contactForce[*nodeID] += sForce;

            //          cohesiveForce[*nodeID] += cForce;
          }
        }
      }
    }
  }



  LagrangeHelperFunctions::LinearPointUpdatePart2( domain.m_feNodeManager, time, dt, m_ldSolve->m_gravityVector, m_ldSolve->DampingM() );

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
  if(doParallelPlateStep && m_propSolve){
	    m_oldPPStepTime = m_newPPStepTime;
	    m_pp_dt = std::min( m_ppSolve->m_stabledt.m_maxdt, m_nSubsteps*dt);
	    //m_pp_dt = std::min(m_pp_dt, m_propSolve->m_stabledt.m_maxdt);
	    m_pp_dt = SetGlobalMaxTimeStep( m_pp_dt );
	    m_newPPStepTime +=  m_pp_dt;
  }



  if (time > m_ldSolve->m_timeToSnapshotDisp)
  {
    m_ldSolve->SnapshotNodalDisplacement(domain.m_feNodeManager);
    m_ldSolve->m_timeToSnapshotDisp = std::numeric_limits<realT>::max();
  }



}

//void Hydrofracture::ApplyForcesFromContactExplicit(PhysicalDomainT& domain,
//                                                   StableTimeStep& timeStep,
//                                                   const realT dt )
//{
//  //-----------------------------
//  //external faces
//  //-----------------------------
//  {
//    Array1dT<R1Tensor>& contactForce = domain.m_feNodeManager.GetFieldData<FieldInfo::contactForce>();
//    contactForce = 0.0;
//    Array1dT<R1Tensor>& decontactForce = domain.m_discreteElementSurfaceNodes.GetFieldData<FieldInfo::contactForce>();
//    decontactForce = 0.0;
//
//    //update nodal positions, velocities, and accelerations before updating face geometry
//    //also, reset rotational and translational accelerations as well as forces and moments
//    domain.m_discreteElementManager.UpdateNodalStatesZeroForcesAndAccelerations();
//
//    //update face geometry and sort faces if necessary
//    const bool resort = domain.m_externalFaces.RecalculateNeighborList(domain.m_feNodeManager,
//                                                                       domain.m_discreteElementSurfaceNodes,
//                                                                       domain.m_discreteElementManager,
//                                                                       dt);
//
//    //if a resort has been triggered, then you also need to update the contact manager
//    if (resort)
//      domain.m_contactManager.Update(domain.m_externalFaces.m_neighborList);
//
//    {
//      Array1dT<Array1dT<R1Tensor> > xs;
//      xs.resize(domain.m_externalFaces.DataLengths());
//      domain.m_externalFaces.UpdateGeometricContactProperties(dt, domain, xs);
//#ifdef STATES_ON_CONTACTS
//      domain.m_externalFaces.UpdateAndApplyContactForcesFromStatesOnCommonPlanes(dt, domain);
//#else
//      domain.m_externalFaces.UpdateAndApplyContactStresses(timeStep, dt, domain, xs);
//#endif
//    }
//
//    //for parallel: do NOT need to synchronize nodal force fields across processes for DE nodes
//    //before moving to centroid, since each process will do that calculation (redundantly) itself
//    // ... there's therefore no need for explicit synchrony
//    //as long as nodal states remain synchronous, everything will be fine!
//#ifndef DEEFC
//    domain.m_discreteElementManager.ApplyDiscreteElementContactForces(45, 0.2);
//#else
//    domain.m_discreteElementManager.ApplyNodalForces();
//#endif
//  }
//}
//



void Hydrofracture::TimeStepImplicit( const realT& time,
                                      const realT& dt,
                                      const int cycleNumber,
                                      PhysicalDomainT& domain,
                                      const sArray1d& namesOfSolverRegions,
                                      SpatialPartition& partition,
                                      FractunatorBase* const fractunator )
{
#if 0
  m_ppSolve->TimeStep( time, dt, domain, namesOfSolverRegions, partition, NULL );

#else
#define FLOW 1


  m_ppSolve->m_boundPhysicalAperture = 0;


  Array1dT<R1Tensor>& effectiveNodalNormal = domain.m_feNodeManager.GetFieldData<R1Tensor>( "effectiveNormal");
  for( localIndex a=0 ; a<domain.m_feNodeManager.DataLengths() ; ++a )
  {
    domain.m_feNodeManager.CalculateEffectiveNormal( a, domain.m_feFaceManager, effectiveNodalNormal[a] );
  }


  SetupSystem( domain, partition, time );


  m_ppSolve->RegisterTemporaryFields( domain );
  m_ldSolve->RegisterTemporaryFields( domain );

  m_ppSolve->GenerateParallelPlateGeometricQuantities( domain,time,dt );



  rArray1d& mass     = domain.m_feFaceManager.GetFieldData<FieldInfo::mass>();

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if( rank==0 )
  {
    mass[4] = 1.0e-2 * time;
    if( mass[4] > 0.1) mass[4] = 0.1;
  }

  m_ppSolve->UpdateEOS( time, dt, domain );

  const rArray1d& density  = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();
  const rArray1d& pressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& fluidVolume  = domain.m_feFaceManager.GetFieldData<FieldInfo::volume>();
  const rArray1d& aperture  = domain.m_feFaceManager.GetFieldData<realT>( "Aperture" );



  realT stepDt = dt;
  int converged = 0;
  realT lastNorm = 0;
  realT thisNorm = 0;

  m_ppSolve->FillTemporaryFields( domain );
  m_ldSolve->FillTemporaryFields( domain );


  while( converged == 0)
  {

    for( int iter=0 ; iter<30 ; ++iter )
    {
      partition.SynchronizeFields(m_syncedFields, CommRegistry::hydrofractureSolver);
      thisNorm = Assemble( domain, partition, time, stepDt);

      if( iter > 1 && lastNorm > 0.0 )
      {
        realT scaleFactor = -1.0;

        if( thisNorm > lastNorm && thisNorm > 0.0 )
        for( int i=0 ; i<4 ; ++i )
        {
          scaleFactor *= 0.5;
          this->PropagateSolutionBlock(  scaleFactor, domain );
          partition.SynchronizeFields(m_syncedFields, CommRegistry::hydrofractureSolver);
          thisNorm = Assemble( domain, partition, time, stepDt);

          if( thisNorm < lastNorm )
            break;
        }
        if( thisNorm > lastNorm )
        {
          std::cout<<"Line Search Failed thisNorm>lastNorm ("<<thisNorm<<">"<<lastNorm<<std::endl;
          break;
        }

      }


      bool print = false;

      if( iter==0 && rank==0 )
      {
//        printf(" rank iter                mass          aperture                volume             density            pressure            residual \n");
      }
      if( size==1 )
      {
        printf("iter, residual = %4d, %16.12e\n", iter, thisNorm );
        if( print )
        {
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[4],  aperture[4],  fluidVolume[4],  density[4],  pressure[4] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[9],  aperture[9],  fluidVolume[9],  density[9],  pressure[9] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[23], aperture[23], fluidVolume[23], density[23], pressure[23]);
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[27], aperture[27], fluidVolume[27], density[27], pressure[27]);
        }
      }
      else if( rank==0 )
      {
        printf("iter, residual = %4d, %16.12e\n", iter, thisNorm );
        if( print )
        {
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[4],  aperture[4],  fluidVolume[4],  density[4],  pressure[4] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[14], aperture[14], fluidVolume[14], density[14], pressure[14] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[23], aperture[23], fluidVolume[23], density[23], pressure[23] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[31], aperture[31], fluidVolume[31], density[31], pressure[31] );
        }
      }
      else if( rank==1 )
      {
        if( print )
        {
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[15], aperture[15], fluidVolume[15], density[15], pressure[15] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[4],  aperture[4],  fluidVolume[4],  density[4],  pressure[4] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[23], aperture[23], fluidVolume[23], density[23], pressure[23] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[31], aperture[31], fluidVolume[31], density[31], pressure[31] );
        }
      }
      else if( rank==2 )
      {
        if( print )
        {
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[14], aperture[14], fluidVolume[14], density[14], pressure[14] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[23], aperture[23], fluidVolume[23], density[23], pressure[23] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[4],  aperture[4],  fluidVolume[4],  density[4],  pressure[4] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[31], aperture[31], fluidVolume[31], density[31], pressure[31] );
        }
      }
      else if( rank==3 )
      {
        if( print )
        {
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[15], aperture[15], fluidVolume[15], density[15], pressure[15] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[24], aperture[24], fluidVolume[24], density[24], pressure[24] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[32], aperture[32], fluidVolume[32], density[32], pressure[32] );
          printf("%4d %4d  % 16.12e % 16.12e % 16.12e % 16.12e % 16.12e \n",rank, iter, mass[4],  aperture[4],  fluidVolume[4],  density[4],  pressure[4] );
        }
      }


      if( iter>1 && thisNorm < 1.0e-5 )
      {
        converged = 1;
        break;
      }

#if 1
      this->SolveBlock( domain, partition, m_epetraBlockSystem, time, stepDt );
#else
      this->m_matrix = m_epetraBlockSystem.m_matrix[0][0];
      this->m_rhs = m_epetraBlockSystem.m_rhs[0];
      this->m_solution = m_epetraBlockSystem.m_solution[0];
      this->Solve( domain, partition, time, stepDt );
#endif

#if 0


      std::ofstream mat_output;
      std::ofstream rhs_output;
      std::ofstream sol_output;



      char mat00Namei[100];
      char mat01Namei[100];
      char mat10Namei[100];
      char mat11Namei[100];
      char rhs0Namei[100];
      char rhs1Namei[100];
      char sol0Namei[100];
      char sol1Namei[100];
      sprintf(mat00Namei, "mat00_%02d", iter );
      sprintf(mat01Namei, "mat01_%02d", iter );
      sprintf(mat10Namei, "mat10_%02d", iter );
      sprintf(mat11Namei, "mat11_%02d", iter );
      sprintf(rhs0Namei, "rhs0_%02d", iter );
      sprintf(rhs1Namei, "rhs1_%02d", iter );
      sprintf(sol0Namei, "sol0_%02d", iter );
      sprintf(sol1Namei, "sol1_%02d", iter );


      EpetraExt::RowMatrixToMatlabFile(mat00Namei,*m_epetraBlockSystem.m_matrix[0][0]);
      EpetraExt::RowMatrixToMatlabFile(mat01Namei,*m_epetraBlockSystem.m_matrix[0][1]);
      EpetraExt::RowMatrixToMatlabFile(mat10Namei,*m_epetraBlockSystem.m_matrix[1][0]);
      EpetraExt::RowMatrixToMatlabFile(mat11Namei,*m_epetraBlockSystem.m_matrix[1][1]);

      EpetraExt::MultiVectorToMatlabFile(rhs0Namei,*m_epetraBlockSystem.m_rhs[0]);
      EpetraExt::MultiVectorToMatlabFile(rhs1Namei,*m_epetraBlockSystem.m_rhs[1]);

      EpetraExt::MultiVectorToMatlabFile(sol0Namei,*m_epetraBlockSystem.m_solution[0]);
      EpetraExt::MultiVectorToMatlabFile(sol1Namei,*m_epetraBlockSystem.m_solution[1]);



      /*
      sprintf(mat00Namei, "mat00_%02d_rank%02d", iter, rank );
      sprintf(mat01Namei, "mat01_%02d_rank%02d", iter, rank );
      sprintf(mat10Namei, "mat10_%02d_rank%02d", iter, rank );
      sprintf(mat11Namei, "mat11_%02d_rank%02d", iter, rank );
      sprintf(rhs0Namei, "rhs0_%02d_rank%02d", iter, rank );
      sprintf(rhs1Namei, "rhs1_%02d_rank%02d", iter, rank );
      sprintf(sol0Namei, "sol0_%02d_rank%02d", iter, rank );
      sprintf(sol1Namei, "sol1_%02d_rank%02d", iter, rank );

      mat_output.precision(12);
      rhs_output.precision(12);
      sol_output.precision(12);

      mat_output.open( mat00Namei );
      m_epetraBlockSystem.m_matrix[0][0]->Print( mat_output );
      mat_output.close();

      mat_output.open( mat01Namei );
      m_epetraBlockSystem.m_matrix[0][1]->Print( mat_output );
      mat_output.close();

      mat_output.open( mat10Namei );
      m_epetraBlockSystem.m_matrix[1][0]->Print( mat_output );
      mat_output.close();

      mat_output.open( mat11Namei );
      m_epetraBlockSystem.m_matrix[1][1]->Print( mat_output );
      mat_output.close();



      rhs_output.open( rhs0Namei );
      m_epetraBlockSystem.m_rhs[0]->Print( rhs_output );
      rhs_output.close();

      rhs_output.open( rhs1Namei );
      m_epetraBlockSystem.m_rhs[1]->Print( rhs_output );
      rhs_output.close();

      sol_output.open( sol0Namei );
      m_epetraBlockSystem.m_solution[0]->Print( sol_output );
      sol_output.close();

      sol_output.open( sol1Namei );
      m_epetraBlockSystem.m_solution[1]->Print( sol_output );
      sol_output.close();
      */
#endif

      lastNorm = thisNorm;
    }
    if( converged==0 )
    {
      stepDt *= 0.5;
      this->m_stabledt.m_maxdt = stepDt * 1.1;
      m_ppSolve->OverwriteFieldsWithTemporaryFields(domain);
      m_ldSolve->OverwriteFieldsWithTemporaryFields(domain);
      if( rank==0 )
        std::cout<<"No newton convergence...substepping at dt="<<stepDt<<" for this step "<<std::endl;
    }
    else
    {
      if( rank==0 )
        std::cout<<"Newton convergence...residual="<<thisNorm<<std::endl;
      this->m_stabledt.m_maxdt = stepDt * 1.1;
    }
  }
#endif

  m_ppSolve->UpdateFlux( time, stepDt, domain, partition );

  rArray1d& massRate     = domain.m_feFaceManager.GetFieldData<realT>("MassRate");
  massRate = 0;
  const rArray1d& mass_n = domain.m_feFaceManager.GetFieldData<realT>("FluidMass_n");

  massRate = mass;
  massRate -= mass_n;
  massRate /= stepDt;



  m_ppSolve->DeregisterTemporaryFields( domain );
  m_ldSolve->DeregisterTemporaryFields( domain );

}

void Hydrofracture::SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time)
{

  m_ldSolve->SetNumRowsAndTrilinosIndices( domain, partition,
                                          m_deformationVariablesLocalRows,
                                          m_deformationVariablesGlobalRows );

  m_ppSolve->DefineFlowSets( domain );
#if FLOW==1

  m_ppSolve->SetNumRowsAndTrilinosIndices( domain, partition,
                                          m_flowVariablesLocalRows,
                                          m_flowVariablesGlobalRows );
#endif



  int dim = 3;

  m_deformationVariablesLocalRows  *= dim;
  m_deformationVariablesGlobalRows *= dim;


#if 1


#if USECPP11==1
  std::shared_ptr<Epetra_Map>& dispRowMap = m_epetraBlockSystem.GetRowMap(EpetraBlock::displacement);
  std::shared_ptr<Epetra_Map>& flowRowMap = m_epetraBlockSystem.GetRowMap(EpetraBlock::fluidMass);

  dispRowMap = std::make_shared<Epetra_Map>( m_deformationVariablesGlobalRows, m_deformationVariablesLocalRows, 0, epetra_comm );
  flowRowMap = std::make_shared<Epetra_Map>( m_flowVariablesGlobalRows, m_flowVariablesLocalRows, 0, epetra_comm );


  std::shared_ptr<Epetra_FECrsGraph>& sparsity_00 = m_epetraBlockSystem.GetSparsity( EpetraBlock::displacement, EpetraBlock::displacement );
  std::shared_ptr<Epetra_FECrsGraph>& sparsity_01 = m_epetraBlockSystem.GetSparsity( EpetraBlock::displacement, EpetraBlock::fluidMass );
  std::shared_ptr<Epetra_FECrsGraph>& sparsity_10 = m_epetraBlockSystem.GetSparsity( EpetraBlock::fluidMass, EpetraBlock::displacement );
  std::shared_ptr<Epetra_FECrsGraph>& sparsity_11 = m_epetraBlockSystem.GetSparsity( EpetraBlock::fluidMass, EpetraBlock::fluidMass );

  sparsity_00 = std::make_shared<Epetra_FECrsGraph>( Copy, *dispRowMap, 0 );
  sparsity_01 = std::make_shared<Epetra_FECrsGraph>( Copy, *dispRowMap, 0 );
  sparsity_11 = std::make_shared<Epetra_FECrsGraph>( Copy, *flowRowMap, 0 );
  sparsity_10 = std::make_shared<Epetra_FECrsGraph>( Copy, *flowRowMap, 0 );
#else
  Epetra_Map*& dispRowMap = m_epetraBlockSystem.GetRowMap(EpetraBlock::displacement);
  Epetra_Map*& flowRowMap = m_epetraBlockSystem.GetRowMap(EpetraBlock::fluidMass);

  delete dispRowMap;
  delete flowRowMap;

  dispRowMap = new Epetra_Map( m_deformationVariablesGlobalRows, m_deformationVariablesLocalRows, 0, epetra_comm );
  flowRowMap = new Epetra_Map( m_flowVariablesGlobalRows, m_flowVariablesLocalRows, 0, epetra_comm );

  Epetra_FECrsGraph*& sparsity_00 = m_epetraBlockSystem.GetSparsity( EpetraBlock::displacement, EpetraBlock::displacement );
  Epetra_FECrsGraph*& sparsity_01 = m_epetraBlockSystem.GetSparsity( EpetraBlock::displacement, EpetraBlock::fluidMass );
  Epetra_FECrsGraph*& sparsity_10 = m_epetraBlockSystem.GetSparsity( EpetraBlock::fluidMass, EpetraBlock::displacement );
  Epetra_FECrsGraph*& sparsity_11 = m_epetraBlockSystem.GetSparsity( EpetraBlock::fluidMass, EpetraBlock::fluidMass );

  delete sparsity_00;
  delete sparsity_01;
  delete sparsity_10;
  delete sparsity_11;

  sparsity_00 = new Epetra_FECrsGraph( Copy, *dispRowMap, 0 );
  sparsity_01 = new Epetra_FECrsGraph( Copy, *dispRowMap, 0 );
  sparsity_10 = new Epetra_FECrsGraph( Copy, *flowRowMap, 0 );
  sparsity_11 = new Epetra_FECrsGraph( Copy, *flowRowMap, 0 );
#endif



  const iArray1d& deformationTrilinosIndex = domain.m_feNodeManager.GetFieldData<int>( m_ldSolve->m_trilinosIndexStr);
  const iArray1d& flowTrilinosIndex = domain.m_feFaceManager.GetFieldData<int>( m_ppSolve->TrilinosIndexString() );



  LagrangeSolverBase::RegionMap::const_iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();


  // Loop over regular Elements
  for(; region != end_region; ++region)
  {
    const ElementRegionT& elemRegion = region->second;
    const unsigned numNodesPerElement = elemRegion.m_numNodesPerElem;

    for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
    {
      const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[element];

      iArray1d nonEmptyDOF;

      for(unsigned i=0; i<numNodesPerElement; ++i)
      {
        const localIndex localNodeIndex = localNodeIndices[i];

        const long long dof =  dim*deformationTrilinosIndex[localNodeIndex];

        for( int d=0 ; d<dim ; ++d )
        {
          nonEmptyDOF.push_back(dof+d);
        }
      }
      sparsity_00->InsertGlobalIndices( nonEmptyDOF.size(),
                                        nonEmptyDOF.data(),
                                        nonEmptyDOF.size(),
                                        nonEmptyDOF.data() );

    }
  }



  // Loop over Flow Faces
  const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );

  std::map<localIndex,lArray1d>::iterator itrEnd = m_ppSolve->m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_ppSolve->m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {

    iArray1d nonEmptyDispDOF;
    iArray1d nonEmptyFlowDOF;

    const lArray1d& facelist = itr->second;
    const unsigned int numFaces = facelist.size();

    for( unsigned int kr=0 ; kr<numFaces ; ++kr )
    {
      const unsigned int r = facelist[kr];
      const long long rfaceDOF = flowTrilinosIndex[r] ;
      nonEmptyFlowDOF.push_back( rfaceDOF );

      for( localIndex a=0 ; a<domain.m_feFaceManager.m_toNodesRelation[r].size() ; ++a )
      {
        const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[r][a];
        const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[childFaceIndex[r][0]][a];

        const long long nodeDOF1 = dim*deformationTrilinosIndex[localNodeIndex1];
        const long long nodeDOF2 = dim*deformationTrilinosIndex[localNodeIndex2];

        for( int d=0 ; d<dim ; ++d )
        {
          nonEmptyDispDOF.push_back( nodeDOF1+d );
          nonEmptyDispDOF.push_back( nodeDOF2+d );
        }
      }
    }

    sparsity_00->InsertGlobalIndices( nonEmptyDispDOF.size(),
                                      nonEmptyDispDOF.data(),
                                      nonEmptyDispDOF.size(),
                                      nonEmptyDispDOF.data() );

    sparsity_11->InsertGlobalIndices( nonEmptyFlowDOF.size(),
                                      nonEmptyFlowDOF.data(),
                                      nonEmptyFlowDOF.size(),
                                      nonEmptyFlowDOF.data() );

    sparsity_01->InsertGlobalIndices( nonEmptyDispDOF.size(),
                                      nonEmptyDispDOF.data(),
                                      nonEmptyFlowDOF.size(),
                                      nonEmptyFlowDOF.data() );

    sparsity_10->InsertGlobalIndices( nonEmptyFlowDOF.size(),
                                      nonEmptyFlowDOF.data(),
                                      nonEmptyDispDOF.size(),
                                      nonEmptyDispDOF.data() );



  } // edge loop


  for( int n=0 ; n<2 ; ++n )
  {
    for( int m=0 ; m<2 ; ++m )
    {

// JOSH: Note m and n are switched here.  Need to figure out how to properly determine domain map in parallel.

//      m_epetraBlockSystem.m_sparsity[n][m]->GlobalAssemble(*m_epetraBlockSystem.m_rowMap[m],*m_epetraBlockSystem.m_rowMap[n]);
      m_epetraBlockSystem.m_sparsity[n][m]->GlobalAssemble(*m_epetraBlockSystem.m_rowMap[m],*m_epetraBlockSystem.m_rowMap[n]);
      m_epetraBlockSystem.m_sparsity[n][m]->OptimizeStorage();
//      m_epetraBlockSystem.m_sparsity[n][m]->Print(std::cout);
#if USECPP11==1
      m_epetraBlockSystem.m_matrix[n][m] = std::make_shared<Epetra_FECrsMatrix>(Copy,*(m_epetraBlockSystem.m_sparsity[n][m]) );
#else
      delete m_epetraBlockSystem.m_matrix[n][m];
      m_epetraBlockSystem.m_matrix[n][m] = new Epetra_FECrsMatrix(Copy,*(m_epetraBlockSystem.m_sparsity[n][m]) );
#endif

//      Epetra_FECrsMatrix* junk = &(*(m_epetraBlockSystem.m_matrix[n][m]));
//      std::cout<<junk->NumGlobalRows()<<", "<<junk->NumGlobalCols()<<std::endl;
    }
#if USECPP11==1
    m_epetraBlockSystem.m_solution[n]  = std::make_shared<Epetra_FEVector>( *(m_epetraBlockSystem.m_rowMap[n]) );
    m_epetraBlockSystem.m_rhs[n]       = std::make_shared<Epetra_FEVector>( *(m_epetraBlockSystem.m_rowMap[n]) );
#else

    delete m_epetraBlockSystem.m_solution[n];
    delete m_epetraBlockSystem.m_rhs[n];

    m_epetraBlockSystem.m_solution[n]  = new Epetra_FEVector( *(m_epetraBlockSystem.m_rowMap[n]) );
    m_epetraBlockSystem.m_rhs[n]       = new Epetra_FEVector( *(m_epetraBlockSystem.m_rowMap[n]) );
#endif
  }


#else
  m_rowMap = std::make_shared<Epetra_Map>(m_deformationVariablesGlobalRows + m_flowVariablesGlobalRows,
                                          m_deformationVariablesLocalRows + m_flowVariablesLocalRows,
                                          0, epetra_comm );

  m_sparsity = std::make_shared<Epetra_FECrsGraph>( Copy, *m_rowMap, 0 );


  const iArray1d& deformationTrilinosIndex = domain.m_feNodeManager.GetFieldData<int>( m_ldSolve->m_trilinosIndexStr);
  const iArray1d& flowTrilinosIndex = domain.m_feFaceManager.GetFieldData<int>( m_ppSolve->TrilinosIndexString() );


  LagrangeSolverBase<3>::RegionMap::const_iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();

  iArray1d nonEmptyDOF;

  for(; region != end_region; ++region)
  {
    const ElementRegionT& elemRegion = region->second;
    const unsigned numNodesPerElement = elemRegion.m_numNodesPerElem;

//    iArray1d nonEmptyDOF (dim*numNodesPerElement);

    for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
    {
      const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[element];

      for(unsigned i=0; i<numNodesPerElement; ++i)
      {
        const localIndex localNodeIndex = localNodeIndices[i];

        const long long dof =  dim*deformationTrilinosIndex[localNodeIndex];

        for( int d=0 ; d<dim ; ++d )
        {
          nonEmptyDOF.push_back(dof+d);
        }
      }

    }

  }


  const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );

  std::map<localIndex,lArray1d>::iterator itrEnd = m_ppSolve->m_edgesToFaces.end();
  for( std::map<localIndex,lArray1d>::iterator itr = m_ppSolve->m_edgesToFaces.begin(); itr!=itrEnd  ; ++itr )
  {
    const lArray1d& facelist = itr->second;
    const unsigned int numFaces = facelist.size();

    for( unsigned int kr=0 ; kr<numFaces ; ++kr )
    {
      const unsigned int r = facelist[kr];
      const long long rfaceDOF = this->m_deformationVariablesGlobalRows + flowTrilinosIndex[r] ;
      nonEmptyDOF.push_back( rfaceDOF );

      for( localIndex a=0 ; a<domain.m_feFaceManager.m_toNodesRelation[r].size() ; ++a )
      {
        const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[r][a];
        const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[childFaceIndex[r][0]][a];

        const long long nodeDOF1 = dim*deformationTrilinosIndex[localNodeIndex1];
        const long long nodeDOF2 = dim*deformationTrilinosIndex[localNodeIndex2];

        for( int d=0 ; d<dim ; ++d )
        {
          nonEmptyDOF.push_back( nodeDOF1+d );
          nonEmptyDOF.push_back( nodeDOF2+d );
        }
      }
    }

  } // edge loop

  m_sparsity->InsertGlobalIndices(nonEmptyDOF.size(),
                                  nonEmptyDOF.data(),
                                  nonEmptyDOF.size(),
                                  nonEmptyDOF.data() );

  m_sparsity->GlobalAssemble();
  m_sparsity->OptimizeStorage();

#endif


}

realT Hydrofracture::Assemble( PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt)
{
  m_ppSolve->GenerateParallelPlateGeometricQuantities( domain, time,dt );
  m_ppSolve->UpdateEOS( time, dt, domain );


  for( int n=0 ; n<2 ; ++n )
  {
    for( int m=0 ; m<2 ; ++m )
    {
      m_epetraBlockSystem.m_matrix[n][m]->Scale(0.0);
    }
    m_epetraBlockSystem.m_rhs[n]->Scale(0.0);
  }

  m_ldSolve->SetMatrixPtr( m_epetraBlockSystem.m_matrix[0][0] );
  m_ldSolve->SetRhsPtr( m_epetraBlockSystem.m_rhs[0] );
  m_ldSolve->SetRowMapPtr( m_epetraBlockSystem.m_rowMap[0] );

  realT nodalForceScale0 = m_ldSolve->Assemble( domain, partition, time+dt );

  m_ppSolve->SetMatrixPtr( m_epetraBlockSystem.m_matrix[1][1] );
  m_ppSolve->SetRhsPtr( m_epetraBlockSystem.m_rhs[1] );
  m_ppSolve->SetRowMapPtr( m_epetraBlockSystem.m_rowMap[1] );

#if USECPP11==1
//  Epetra_FECrsMatrix* junk = (m_epetraBlockSystem.m_matrix[1][0]).get();
#else
//  Epetra_FECrsMatrix* junk = m_epetraBlockSystem.m_matrix[1][0];
#endif

  realT fluidMassScale = m_ppSolve->Assemble( domain, m_epetraBlockSystem, time, dt );


  realT nodalForceScale1 = AssembleCouplingTerms( domain,
                                                  m_ppSolve->m_dwdu,
                                                  m_ppSolve->m_dwdw );
  realT residualNormForce = 0;
  realT residualNormMass = 0;
  int dummy;
  double* residual = NULL;

  m_epetraBlockSystem.m_rhs[0]->ExtractView(&residual,&dummy);
  realT nodalForceScale = std::max( nodalForceScale0 , nodalForceScale1 );

  for( int i=0 ; i<m_deformationVariablesLocalRows ; ++i )
  {
    residualNormForce += residual[i]*residual[i]/this->m_deformationVariablesGlobalRows;
  }



  m_epetraBlockSystem.m_rhs[1]->ExtractView(&residual,&dummy);
  for( int i=0 ; i<m_flowVariablesLocalRows ; ++i )
  {
    residualNormMass += residual[i]*residual[i]/this->m_flowVariablesGlobalRows;
  }


  realT localSumData[2] = { residualNormMass, residualNormForce };
  realT globalSumData[2] = {0.0, 0.0};
  MPI_Allreduce (localSumData,globalSumData,2,MPI_DOUBLE,MPI_SUM ,MPI_COMM_WORLD);

  realT localMaxData[2] = { fluidMassScale, nodalForceScale };
  realT globalMaxData[2] = {0.0, 0.0};
  MPI_Allreduce (localMaxData,globalMaxData,2,MPI_DOUBLE,MPI_MAX ,MPI_COMM_WORLD);

  residualNormMass = globalSumData[0] ;
  residualNormForce = globalSumData[1] ;


  /*
  if( nodalForceScale > 1.0 )
    residualNormForce = globalSumData[1] / globalMaxData[1]*globalMaxData[1];

  if( fluidMassScale > 0.0 )
    residualNormMass = globalSumData[0] / globalMaxData[0]*globalMaxData[0];
*/


  return ( sqrt( residualNormForce+residualNormMass ) );
}



realT Hydrofracture::AssembleCouplingTerms( const PhysicalDomainT& domain,
                                            const Array1dT< rArray1d >& dwdu,
                                            const rArray1d& dwdw )
{
  const int dim = 3;
  const iArray1d& deformationTrilinosIndex = domain.m_feNodeManager.GetFieldData<int>( m_ldSolve->m_trilinosIndexStr);
  const iArray1d& flowTrilinosIndex = domain.m_feFaceManager.GetFieldData<int>( m_ppSolve->TrilinosIndexString() );

  realT maxForce = 0.0 ;

  const rArray1d* const ppFlowPressure = domain.m_feFaceManager.GetFieldDataPointer<FieldInfo::pressure>();

  if( ppFlowPressure!=NULL )
  {
    const rArray1d& density_np1 = domain.m_feFaceManager.GetFieldData<FieldInfo::density>();

    Epetra_IntSerialDenseVector  dispDofIndex;
    Epetra_IntSerialDenseVector  flowDofIndex;

    Epetra_SerialDenseVector     face_rhs;
    Epetra_SerialDenseMatrix     matrix_00;
    Epetra_SerialDenseMatrix     matrix_01;
    Epetra_SerialDenseMatrix     matrix_10;


    // determine ghost elements
    const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();
    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
    const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");
    const rArray1d& dPdM = domain.m_feFaceManager.GetFieldData<realT>("dPdM");


    for(localIndex kf = 0 ; kf < domain.m_feFaceManager.m_numFaces ; ++kf)
    {
      if( isGhost[kf] < 0 && flowFaceType[kf]==0 )
      {

        const localIndex faceIndex[2] = { kf, childFaceIndex[kf][0] };

        const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                                domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

        R1Tensor Nbar = N[0];
        Nbar -= N[1];
        Nbar.Normalize();


        const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf].size();

        {
          const int nDispDOF = 2*dim*numNodes;
          const int nFlowDOF = 1;
          dispDofIndex.Resize(nDispDOF);
          flowDofIndex.Resize(nFlowDOF);

          face_rhs.Resize(nDispDOF);
          matrix_00.Reshape(nDispDOF,nDispDOF);
          matrix_01.Reshape(nDispDOF,nFlowDOF);
          matrix_10.Reshape(nFlowDOF,nDispDOF);

          face_rhs.Scale(0.0);
          matrix_00.Scale(0.0);
          matrix_01.Scale(0.0);
          matrix_10.Scale(0.0);


          for( localIndex a=0 ; a<numNodes ; ++a )
          {
            const localIndex aa = a == 0 ? a : numNodes - a;
            const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[0]][a];
            const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[1]][aa];

            const int aDof = 2*a*dim;
            for( int i=0 ; i<dim ; ++i )
            {
              dispDofIndex[aDof+i]       = dim*deformationTrilinosIndex[localNodeIndex1]+i;
              dispDofIndex[aDof+dim+i]   = dim*deformationTrilinosIndex[localNodeIndex2]+i;
            }
          }
          flowDofIndex[0] = flowTrilinosIndex[kf];

          const realT dPdV = -dPdM[kf]*density_np1[kf];

          // todo use correct area
          const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, kf );

          const realT hydroForceMag = (*ppFlowPressure)[kf] * area;
          const realT hydroStiffnessMag = dPdM[kf] * area;


          for( localIndex a=0 ; a<numNodes ; ++a )
          {
            // rhs and mass derivative

            const int aDof = 2*a*dim;
            for( int i=0 ; i<dim ; ++i )
            {
              face_rhs[aDof+i]     +=  hydroForceMag * Nbar[i] / numNodes ;
              face_rhs[aDof+dim+i] += -hydroForceMag * Nbar[i] / numNodes ;

              maxForce = std::max( maxForce, fabs(hydroForceMag * Nbar[i] / numNodes) );

              matrix_01(aDof+i,0)     += -hydroStiffnessMag * Nbar[i] / numNodes;
              matrix_01(aDof+dim+i,0) +=  hydroStiffnessMag * Nbar[i] / numNodes;
            }



            // displacement derivative
            for( localIndex b=0 ; b<numNodes ; ++b)
            {
              const int bDof = 2*b*dim;
              for( int i=0 ; i<dim ; ++i )
              {
                for( int j=0 ; j<dim ; ++j )
                {
                  matrix_00(aDof+i,bDof+j)          -= dPdV*area*dwdw(kf)*dwdu(kf)(bDof+j)*Nbar[i]     * ( area / numNodes );
                  matrix_00(aDof+i,bDof+dim+j)      -= dPdV*area*dwdw(kf)*dwdu(kf)(bDof+dim+j)*Nbar[i] * ( area / numNodes );
                  matrix_00(aDof+dim+i,bDof+j)      += dPdV*area*dwdw(kf)*dwdu(kf)(bDof+j)*Nbar[i]     * ( area / numNodes );
                  matrix_00(aDof+dim+i,bDof+dim+j)  += dPdV*area*dwdw(kf)*dwdu(kf)(bDof+dim+j)*Nbar[i] * ( area / numNodes );

                }
              }
            }


          }


          m_epetraBlockSystem.m_matrix[0][0]->SumIntoGlobalValues( dispDofIndex,
                                                                   dispDofIndex,
                                                                   matrix_00 );


          m_epetraBlockSystem.m_matrix[0][1]->SumIntoGlobalValues( dispDofIndex,
                                                                   flowDofIndex,
                                                                   matrix_01 );


          m_epetraBlockSystem.m_rhs[0]->SumIntoGlobalValues(dispDofIndex,
                                                            face_rhs);


        }
      }
    }
  }  // end element

  m_epetraBlockSystem.m_matrix[0][0]->GlobalAssemble();
  m_epetraBlockSystem.m_matrix[0][1]->GlobalAssemble( *(m_epetraBlockSystem.m_rowMap[1]),
                                                      *(m_epetraBlockSystem.m_rowMap[0]) );
  m_epetraBlockSystem.m_rhs[0]->GlobalAssemble();

  return maxForce;
}


realT Hydrofracture::CheckSolutionBlock( const PhysicalDomainT& domain )
{
//  return   m_ppSolve->CheckSolution( local_solution, domain, this->m_deformationVariablesGlobalRows );
  return 1.0;
}


void Hydrofracture::PropagateSolutionBlock( const realT scalingFactor,
                                            PhysicalDomainT& domain )
{
  int dummy;
  double* solution = NULL;

  m_epetraBlockSystem.m_solution[0]->ExtractView(&solution,&dummy);
  m_ldSolve->PropagateSolution( solution, scalingFactor, domain, 0 );

  m_epetraBlockSystem.m_solution[1]->ExtractView(&solution,&dummy);
  m_ppSolve->PropagateSolution( solution, scalingFactor, domain, 0 );

}


void Hydrofracture::PostSyncConsistency( PhysicalDomainT& domain, SpatialPartition& partition )
{
  m_ldSolve->PostSyncConsistency( domain, partition );
}


void Hydrofracture::PostProcess (PhysicalDomainT& domain,
                                                             SpatialPartition& partition,
                                                             const sArray1d& namesOfSolverRegions)
{
  if( m_timeIntegrationOption==Explicit)
  {
  m_ppSolve->PostProcess(domain, partition, namesOfSolverRegions);
  m_ldSolve->PostProcess(domain, partition, namesOfSolverRegions);
  if (!m_mfSolverName.empty()) m_mfSolve->PostProcess(domain, partition, namesOfSolverRegions);
  }


  rArray1d& sigma_x = domain.m_feNodeManager.GetFieldData<realT>("nsigma_x");
  rArray1d& sigma_y = domain.m_feNodeManager.GetFieldData<realT>("nsigma_y");
  rArray1d& sigma_z = domain.m_feNodeManager.GetFieldData<realT>("nsigma_z");
  rArray1d& sigma_xy = domain.m_feNodeManager.GetFieldData<realT>("nsigma_xy");
  rArray1d& sigma_yz = domain.m_feNodeManager.GetFieldData<realT>("nsigma_yz");
  rArray1d& sigma_zx = domain.m_feNodeManager.GetFieldData<realT>("nsigma_zx");

  sigma_x = 0.0;
  sigma_y = 0.0;
  sigma_z = 0.0;
  sigma_xy = 0.0;
  sigma_yz = 0.0;
  sigma_zx = 0.0;

  for( sArray1d::const_iterator regionName = namesOfSolverRegions.begin() ;
      regionName != namesOfSolverRegions.end() ; ++regionName )
  {
    std::map<std::string, ElementRegionT>::iterator iter = domain.m_feElementManager.m_ElementRegions.find(*regionName);
    if(iter == domain.m_feElementManager.m_ElementRegions.end())
      continue;

    ElementRegionT& elementRegion = iter->second;


    for( localIndex k=0 ; k<elementRegion.m_numElems ; ++k )
    {
      R2SymTensor s;
      realT pressure = 0.0;
      realT ex, ey, ez, exy, eyz, ezx;
      for( localIndex a=0 ; a<elementRegion.m_numIntegrationPointsPerElem ; ++a )
      {
        const MaterialBaseStateData& state = *(elementRegion.m_mat->StateData(k,a));
        s += state.devStress;
        pressure += state.pressure;
      }
      s /= elementRegion.m_numIntegrationPointsPerElem;
      pressure /= elementRegion.m_numIntegrationPointsPerElem;
      ex =  s(0,0) + pressure;
      ey =  s(1,1) + pressure;
      ez =  s(2,2) + pressure;
      exy =  s(0,1);
      eyz =  s(1,2);
      ezx =  s(0,2);

      for (localIndex j = 0; j < elementRegion.m_toNodesRelation.Dimension(1); ++j)
      {
        localIndex ind = elementRegion.m_toNodesRelation[k][j];
        sigma_x[ind] +=ex;
        sigma_y[ind] +=ey;
        sigma_z[ind] +=ez;
        sigma_xy[ind] +=exy;
        sigma_yz[ind] +=eyz;
        sigma_zx[ind] +=ezx;
      }
    }
  }

  for (localIndex i=0; i<domain.m_feNodeManager.DataLengths(); ++i)
  {
    if (domain.m_feNodeManager.m_toElementsRelation[i].size() > 0)
    {
      sigma_x[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_y[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_z[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_xy[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_yz[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
      sigma_zx[i] /= domain.m_feNodeManager.m_toElementsRelation[i].size();
    }
  }

  {
    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_x");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_y");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_z");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_xy");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_yz");
    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("nsigma_zx");

    partition.SynchronizeFields( syncedFields, CommRegistry::lagrangeParallelPlateFlowSolver);
  }

}


void Hydrofracture::CalculateContactStress(PhysicalDomainT& domain,
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
  rArray1d& COFJ = domain.m_feFaceManager.GetFieldData<realT>("COFJ");
  rArray1d& jointCohesion = domain.m_feFaceManager.GetFieldData<realT>("jointCohesion");

  rArray1d* proppedWidth = domain.m_feFaceManager.GetFieldDataPointer<realT>("proppedWidth");
  bool propped(true);
  if (proppedWidth == NULL) propped = false;



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

  if( (!propped && -gapNormal + deltaN <= 0.0) || (propped && -gapNormal + deltaN + (*proppedWidth)[kf]<= 0.0) )
  {
    stressPen = 0.0;
    stressShear = 0.0;
  }
  else
  {
    if (!propped)
    {
      stressPen = (-gapNormal + deltaN) * m_kJn ;
    }
    else
    {
      stressPen = (-gapNormal + deltaN + (*proppedWidth)[kf]) * m_kJn ;
    }

    R1Tensor stressShearPrj = stressShear0[kf];
    stressShearPrj += gapShearInc * m_kJs;

    if ( pow(stressShearPrj.L2_Norm(),2) < pow(stressPen*COFJ[kf] + jointCohesion[kf], 2) )
    {
      stressShear = stressShearPrj;
    }
    else  // We have to calculate plasticity stuff
    {
      jointCohesion[kf] = 0.0;
      realT du[2], t[2], x[2];
      du[0] = Dot(gapShearInc, T[0]);
      du[1] = Dot(gapShearInc, T[1]);
      t[0] = Dot(stressShear0[kf], T[0]);
      t[1] = Dot(stressShear0[kf], T[1]);

      realT a, b, c;
      a = pow(t[0],2) + pow(t[1],2);
      if (a == 0) a = pow( t[0]+du[0]*m_kJs ,2) + pow(t[1]+ du[1]*m_kJs,2);

      b = -2 * ( t[0]*(t[0] + m_kJs * du[0]) + t[1]*(t[1] + m_kJs * du[1]));
      c = pow(t[0] + m_kJs * du[0], 2) + pow(t[1] + m_kJs * du[1], 2) - pow(stressPen*COFJ[kf], 2);

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

realT Hydrofracture::SetGlobalMaxTimeStep(realT local_dt){
	  realT temp = local_dt;
	  realT nextdt;
	  MPI_Allreduce (&temp,&nextdt,1,MPI_DOUBLE,MPI_MIN ,MPI_COMM_WORLD);

	  return nextdt;
}


/// Register solver in the solver factory
REGISTER_SOLVER( Hydrofracture )

