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
 * @file GroupLabeler.cpp
 * @author walsh24
 * @date Nov 21, 2013
 */

#include "GroupLabeler.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"

using namespace PS_STR;

namespace{
	const realT TINY = 1e-64;	
}


GroupLabeler::GroupLabeler( const std::string& name,
                                                          ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

GroupLabeler::~GroupLabeler()
{
  // TODO Auto-generated destructor stub
}

void GroupLabeler::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{  
  SolverBase::ReadXML( hdn );

  m_groupIdFieldName = "groupID";
  m_isMemberFieldName = "memberField";
  
}


void GroupLabeler::RegisterFields( PhysicalDomainT& domain )
{
	
  domain.m_feFaceManager.AddKeylessDataField<int>(m_groupIdFieldName,true,true);
  domain.m_feFaceManager.AddKeylessDataField<realT>(m_isMemberFieldName,true,true);

}

void GroupLabeler::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{
    
}





double GroupLabeler::TimeStep( const realT& time ,
                                            const realT& dt ,
                                            const int cycleNumber,
                                            PhysicalDomainT& domain,
                                            const sArray1d& namesOfSolverRegions,
                                            SpatialPartition& partition ,
                                            FractunatorBase* const fractunator )
{
	  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();

	// currently connect faces by edges - ultimately want to connect object a by object b

  rArray1d& isMemberField = domain.m_feFaceManager.GetFieldData<realT>(m_isMemberFieldName);
  iArray1d& groupIdField = domain.m_feFaceManager.GetFieldData<int>(m_groupIdFieldName);
  groupIdField = -1;  // -1 indicates not in any group (failed member test)
	
  //ElementRegionT& elementRegion = domain.m_feElementManager.m_ElementRegions[namesOfSolverRegions[0]];

  //const FaceManagerT& faceManager = domain.m_feFaceManager;
  UnorderedVariableOneToManyRelation& edgesToFaces = domain.m_feEdgeManager.m_toFacesRelation;


  std::map<localIndex, DSN > faceNodeMap;

  // loop over faces (member object), build nodes
  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf )
  {
  	if( isMemberField[kf] > 0  ){  // passes test
		faceNodeMap[kf] = DSN(-1,0); // -1 = node info - not used here
  	}
  }
     
  //loop over edges ( connecting object)
  for( localIndex ke=0 ; ke<domain.m_feEdgeManager.DataLengths() ; ++ke )
  {
	const lSet& faces = edgesToFaces[ke];
	//localIndex numFaces = faces.size();
	lSet connectedFaces;
	// find member faces and build nodes
    for(lSet::iterator fc = faces.begin(); fc != faces.end(); ++fc){
    	localIndex kf = *fc;
    	if( isMemberField[kf] > 0  ){  // passes test
    		connectedFaces.insert(kf);
    	}
    }
	// build and connect all valid local faces
    for(lSet::iterator fca = connectedFaces.begin(); fca != connectedFaces.end(); ++fca){
    	localIndex kfa = *fca;
        DSN& nodeA = faceNodeMap[kfa];

    	lSet::iterator fcb = fca; ++fcb;
        for(; fcb != connectedFaces.end(); ++fcb){
        	localIndex kfb = *fcb;
            DSN& nodeB = faceNodeMap[kfb];
            nodeA.Union(nodeB);
        }
    }
  }

  // build set of unique roots
  std::map<localIndex, DSN >::iterator faceNode_itr = faceNodeMap.begin();
  std::map<localIndex, DSN >::iterator faceNode_itr_end = faceNodeMap.end();
  std::map<DSN*,int> rootMap;
  faceNode_itr = faceNodeMap.begin();
  for(  ; faceNode_itr != faceNode_itr_end; ++faceNode_itr){
	  DSN& node = faceNode_itr->second;
	  DSN* root = node.FindRoot();
      rootMap[root]=0;
  }

  // number roots
  int idOffset = 0; // placeholder for parallel version
  std::map<DSN*,int>::iterator rootMap_itr_end = rootMap.end();
  int regionCount = 0;
  for(  std::map<DSN*,int>::iterator rootMap_itr = rootMap.begin(); rootMap_itr != rootMap_itr_end; ++rootMap_itr){
	rootMap_itr->second = regionCount + idOffset;  // group id
    ++regionCount;
  }

  // record groups
  faceNode_itr = faceNodeMap.begin();
  for(  ;faceNode_itr != faceNode_itr_end; ++faceNode_itr){
	  localIndex fc = faceNode_itr->first;
	  DSN& node = faceNode_itr->second;
	  DSN* root = node.FindRoot();
	  groupIdField[fc] = rootMap[root];
  }

  return dt;
}



/// Register solver in the solver factory
REGISTER_SOLVER( GroupLabeler )
