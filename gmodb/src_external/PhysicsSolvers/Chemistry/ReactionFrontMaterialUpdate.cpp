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

#include "ReactionFrontMaterialUpdate.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"

using namespace PS_STR;

namespace{
	const realT TINY = 1e-64;	
}


ReactionFrontMaterialUpdate::ReactionFrontMaterialUpdate( const std::string& name,
                                                          ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

ReactionFrontMaterialUpdate::~ReactionFrontMaterialUpdate()
{
  // TODO Auto-generated destructor stub
}

void ReactionFrontMaterialUpdate::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{  
  SolverBase::ReadXML( hdn );

  m_frontNames = hdn->GetStringVector("fronts");
  m_bulkModuli = hdn->GetAttributeVector<realT>("bulk_moduli");
  m_shearModuli = hdn->GetAttributeVector<realT>("shear_moduli");
  m_faceSetName = hdn->GetAttributeString("faceset");
  m_numFronts = m_frontNames.size();
  
  if( m_bulkModuli.size()!= m_shearModuli.size() )
        throw GPException("Number of bulk moduli does not match number of shear moduli. \n");
  	
  if(m_numFronts != m_bulkModuli.size())
        throw GPException("Number of bulk moduli does not match number of fronts. \n");
  	
  
}


void ReactionFrontMaterialUpdate::RegisterFields( PhysicalDomainT& domain )
{
	
  for(unsigned i =0; i < m_frontNames.size(); ++i){
    domain.m_feFaceManager.AddKeylessDataField<realT>(m_frontNames[i],true,true);
  }
}

void ReactionFrontMaterialUpdate::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
    
  FaceManagerT& faceManager = domain.m_feFaceManager;
  m_faceSet = &(faceManager.GetSet(m_faceSetName));
}





double ReactionFrontMaterialUpdate::TimeStep( const realT& time ,
                                            const realT& dt ,
                                            const int cycleNumber,
                                            PhysicalDomainT& domain,
                                            const sArray1d& namesOfSolverRegions,
                                            SpatialPartition& partition ,
                                            FractunatorBase* const fractunator )
{
	
  ElementRegionT& elementRegion = domain.m_feElementManager.m_ElementRegions[namesOfSolverRegions[0]];
  const NodeManagerT& nodeManager = domain.m_feNodeManager;
  const FaceManagerT& faceManager = domain.m_feFaceManager;
	
  // get front pointers
  std::vector< rArray1d* > frontPtrs(m_numFronts);
  m_frontDistances.resize(m_numFronts);
  for(localIndex i=0; i < m_numFronts; ++i){
    frontPtrs[i] = &(domain.m_feFaceManager.GetFieldData<realT>(m_frontNames[i]));
  }    
     
  //loop over faces
  for( lSet::const_iterator kf=m_faceSet->begin() ; kf!=m_faceSet->end() ; ++kf )
  {
  	localIndex fc = *kf;
    for(localIndex i =0; i < m_numFronts; ++i){
    	m_frontDistances(i) = (*frontPtrs[i])[fc];
    }
    realT maxDistance = m_frontDistances(m_numFronts-1);
    
    R1Tensor faceCenter;
    faceManager.FaceCenter( nodeManager, fc, faceCenter );
    
    // collect list of elements 
    lSet checkedElements;
    lSet elementsToCheck;
    
    ElementIdPair eid = faceManager.m_toElementsRelation[fc][0];
    // check first element
    if( elementRegion.ContainsElement(eid) ){
      localIndex el = eid.second;
      R1Tensor elementCenter = elementRegion.GetElementCenter(el, nodeManager); // FIXME should be getting quadrature point locations not element centers
    	
      realT dist = (elementCenter-faceCenter).L2_Norm(); // FIXME - use face normal ???
      if(dist < maxDistance){
      	
      	UpdateElementModuli(elementRegion,el,dist);
      
        // find neighboring elements
        elementRegion.GetElementNeighbors(el,faceManager,elementsToCheck);
        checkedElements.insert(el); 	
      }
    }
    
    while(elementsToCheck.size() > 0){
      localIndex el = *(elementsToCheck.begin());
      checkedElements.insert(el);
      
      R1Tensor elementCenter = elementRegion.GetElementCenter(el, nodeManager);
      realT dist = (elementCenter-faceCenter).L2_Norm(); // FIXME - use face normal??
      
      if(dist < maxDistance){
      	
      	UpdateElementModuli(elementRegion,el,dist);
      	
        // find nbrs
        std::set<localIndex> nbrs;
        elementRegion.GetElementNeighbors(el,faceManager,nbrs);
        
        // add unchecked neighbors to list to check
        set_difference (nbrs.begin(), nbrs.end(), 
                        checkedElements.begin(), checkedElements.end(), 
                        std::inserter(elementsToCheck, elementsToCheck.end()) );
        
      }
      elementsToCheck.erase(el);
    }
    
  }
 return dt;
}


void ReactionFrontMaterialUpdate::UpdateElementModuli(ElementRegionT& elementRegion, localIndex element, realT distance){
  
  bool doUpdate = true;
    
  realT K = elementRegion.m_mat->StateData(element,0)->BulkModulus;
  unsigned i = 0;
  for(; i < m_numFronts; ++i){
    if( isEqual(K,m_bulkModuli[i]) )
    {
      doUpdate = false;
      break;
    }
    if(distance < m_frontDistances[i]) break;
  }
  if(doUpdate){
    elementRegion.m_mat->StateData(element,0)->BulkModulus = m_bulkModuli[i];
    elementRegion.m_mat->StateData(element,0)->ShearModulus = m_shearModuli[i];
  }
  	
}



/// Register solver in the solver factory
REGISTER_SOLVER( ReactionFrontMaterialUpdate )
