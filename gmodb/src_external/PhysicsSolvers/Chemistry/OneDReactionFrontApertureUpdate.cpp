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
 * @file TwoDSteadyStateParallelPlateFlowSolver.cpp
 * @author walsh24
 * @date June 1, 2011
 */

#include "OneDReactionFrontApertureUpdate.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"

//#include "Constitutive/Material/Materials.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"


using namespace PS_STR;

namespace{
	const realT TINY = 1e-64;	
}


ReactionFrontApertureUpdate::ReactionFrontApertureUpdate( const std::string& name, ProblemManagerT* const pm):
SolverBase(name,pm)
{
  //ReadXML(hdn);
}

ReactionFrontApertureUpdate::~ReactionFrontApertureUpdate()
{
  // TODO Auto-generated destructor stub
}

void ReactionFrontApertureUpdate::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{  

  m_faceSetName = hdn->GetAttributeString("faceset");

  Ku = hdn->GetAttributeOrDefault("Ku","60 MPa");
  Kam = hdn->GetAttributeOrDefault<realT>("Kam",Ku);
  Kd = hdn->GetAttributeOrDefault<realT>("Kd",Ku);

  Kamp = hdn->GetAttributeOrDefault<realT>("Kamp",Kam*0.22); // plastic stiffness K''
  Kdp = hdn->GetAttributeOrDefault<realT>("Kdp",Kd*0.30);

  // slider yield strengths
  sigmaYam = hdn->GetAttributeOrDefault("sigmaYam","2 MPa");
  sigmaYd = hdn->GetAttributeOrDefault("sigmaYd","3 MPa");


  pillarRadius = hdn->GetAttributeOrDefault("pillarRadius","1 mm");
  rr = pillarRadius*pillarRadius;
  
  bmin = hdn->GetAttributeOrDefault<realT>("MinimumAperture",0);
  	
  
}


void ReactionFrontApertureUpdate::RegisterFields( PhysicalDomainT& domain )
{

  domain.m_feFaceManager.AddKeyedDataField<FieldInfo::pressure>();

  domain.m_feFaceManager.AddKeylessDataField<realT>("InitialAperture" ,true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>(ApertureStr ,true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("NormalStress" ,true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("x_rf_0" ,true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>("x_rf_2" ,true,true);
}

void ReactionFrontApertureUpdate::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
    
  FaceManagerT& faceManager = domain.m_feFaceManager;
  m_faceSet = &(faceManager.GetSet(m_faceSetName));


}


void ReactionFrontApertureUpdate::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back(ApertureStr);
  syncedFields[PhysicalDomainT::FiniteElementFaceManager].push_back("InitialAperture");
  // partition.SetBufferSizes(syncedFields, CommRegistry::reactionFrontApertureUpdate); // fixme reintroduce once branch has been merged with trunk and tests updated
  partition.SetBufferSizes(syncedFields, CommRegistry::immiscibleFluidSolver);


}





double ReactionFrontApertureUpdate::TimeStep( const realT& time ,
                                            const realT& dt ,
                                            const int cycleNumber,
                                            PhysicalDomainT& domain,
                                            const sArray1d& namesOfSolverRegions,
                                            SpatialPartition& partition ,
                                            FractunatorBase* const fractunator )
{
	

  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();

  //FaceManagerT& faceManager = domain.m_feFaceManager;
  rArray1d& apertures = domain.m_feFaceManager.GetFieldData<realT>("Aperture");
  const rArray1d& initialApertures = domain.m_feFaceManager.GetFieldData<realT>("InitialAperture");
  const rArray1d& porePressures   = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  const rArray1d& normalStresses = domain.m_feFaceManager.GetFieldData<realT>("NormalStress");

  const rArray1d& amorphousExtents = domain.m_feFaceManager.GetFieldData<realT>("x_rf_0");
  const rArray1d& depletedExtents  = domain.m_feFaceManager.GetFieldData<realT>("x_rf_2");

     
  //loop over faces
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
    const localIndex& fc =  *kf;
    realT bu = initialApertures[fc];
    realT b = apertures[fc];
    realT sigmaEff = normalStresses[fc] - porePressures[fc];
    realT amSiExtent = amorphousExtents[fc];
    realT depletedExtent = depletedExtents[fc];
    apertures[fc] = CalculateNewAperture(bu, b, sigmaEff,amSiExtent,depletedExtent);
  }

  // sync new apertures

 // partition.SynchronizeFields(syncedFields, CommRegistry::reactionFrontApertureUpdate); // fixme - reintroduce once branch has been merged with trunk and test updated
  partition.SynchronizeFields(syncedFields, CommRegistry::immiscibleFluidSolver);


  return dt;
}


realT ReactionFrontApertureUpdate::CalculateNewAperture(realT bu, realT b,realT sigmaEff,realT amSiExtent, realT depletedExtent){

	realT bo = bu/(1-sigmaEff/Ku);
	realT eta = (bo-b)/bo;
  
	realT sigmaA = FindSpringSliderYieldStress(Kam,Kamp,sigmaYam,eta); // loading only - will need to be modified for unloading also
	realT sigmaD = FindSpringSliderYieldStress(Kd,Kdp,sigmaYd,eta);
	realT sigmaU = eta*Ku;

	realT aRatio,dRatio,uRatio;
	CalculateRelativeAreas(amSiExtent, depletedExtent, aRatio,dRatio, uRatio);

	realT sigmaStar = aRatio*sigmaA + dRatio*sigmaD + uRatio*sigmaU; // current estimate of net effective stress
	realT dsigma = sigmaEff- sigmaStar;

	realT KKam = FindSpringSliderK(Kam,Kamp,sigmaA,sigmaYam);
	realT KKd  = FindSpringSliderK(Kd,Kdp,sigmaD,sigmaYd);
	realT KKu = Ku; // might make this a spring slider later

	realT   Keff = aRatio*KKam + dRatio*KKd + uRatio*KKu;
	realT   dEta = dsigma/(Keff);

	eta += dEta;

	b = bo*(1-eta);

	if (b < bmin) b = bmin;
  	
	return b;
}



/// Register solver in the solver factory
REGISTER_SOLVER( ReactionFrontApertureUpdate )
