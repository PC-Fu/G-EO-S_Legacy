/**
 * @file Fractunator3.cpp
 * @author settgast1
 * @date Jul 14, 2011
 */

#include "Fractunator3.h"
#include <limits.h>
#include "IO/BinStream.h"

#include "MPI_Communications/SpatialPartition.h"

#include "SurfaceGeneration/FractunatorFactory.h"

static localIndex GetOtherFaceEdge( const std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                                    const localIndex thisFace, const localIndex thisEdge )
{
  localIndex nextEdge = LOCALINDEX_MAX;

  const std::pair<localIndex,localIndex>& faceToEdges = stlMapLookup( localFacesToEdges, thisFace );
  if( faceToEdges.first == thisEdge )
  {
    nextEdge = faceToEdges.second;
  }
  else if( faceToEdges.second == thisEdge )
  {
    nextEdge = faceToEdges.first;
  }
  else
  {
    std::cout<<"";
    throw GPException("Fractunator3::Couldn't find thisEdge in localFacesToEdges[thisFace]");
  }
  return nextEdge;
}

static void CheckForAndRemoveDeadEndPath( const localIndex edgeIndex,
                                          const iArray1d& isEdgeExternal,
                                          std::map< localIndex, std::set<localIndex> >& edgesToRuptureReadyFaces,
                                          std::map< localIndex, std::pair<localIndex,localIndex> >& localVFacesToVEdges,
                                          lSet& nodeToRuptureReadyFaces)
{


  localIndex thisEdge = edgeIndex;

  // if the edge is internal and the edge is only attached to one ruptured face...
  while( isEdgeExternal[thisEdge]!=1 )
  {

    //    std::set<localIndex>& edgeToRuptureReadyFaces = stlMapLookup(edgesToRuptureReadyFaces,thisEdge);
    std::set<localIndex>& edgeToRuptureReadyFaces = edgesToRuptureReadyFaces[thisEdge];

    if( edgeToRuptureReadyFaces.size()!=1 )
      break;

    // then the index for the face that is a "dead end"
    localIndex deadEndFace = *(edgeToRuptureReadyFaces.begin());


    std::pair<localIndex,localIndex>& localVFaceToVEdges = stlMapLookup(localVFacesToVEdges,deadEndFace);

    // get the edge on the other side of the "dead end" face
    localIndex nextEdge;
    if( localVFaceToVEdges.first == thisEdge )
      nextEdge = localVFaceToVEdges.second;
    else if( localVFaceToVEdges.second == thisEdge )
      nextEdge = localVFaceToVEdges.first;
    else
    {
      throw GPException("Fractunator3::FindFracturePlanes: Could not find the next edge when removing dead end faces.");
    }

    // delete the face from the working arrays
    edgeToRuptureReadyFaces.erase( deadEndFace );
    edgesToRuptureReadyFaces[nextEdge].erase(deadEndFace);
    nodeToRuptureReadyFaces.erase(deadEndFace);

    // if all the faces have been deleted, then go ahead and delete the top level entry
    if( edgeToRuptureReadyFaces.empty() )
      edgesToRuptureReadyFaces.erase(thisEdge);
    if( edgesToRuptureReadyFaces[nextEdge].empty() )
      edgesToRuptureReadyFaces.erase(nextEdge);

    // now increment the "thisEdge" to point to the other edge on the face that was just deleted
    thisEdge = nextEdge;
  }

}


Fractunator3::Fractunator3():
        FractunatorBase()
{
  // TODO Auto-generated constructor stub


  //  m_virtualNodes.AddMap< Array1dT<lSet> >("usedFaces");
}

Fractunator3::~Fractunator3()
{
  // TODO Auto-generated destructor stub
}


void Fractunator3::RegisterFieldsAndMaps( NodeManagerT& nodeManager,
                                          EdgeManagerT& edgeManager,
                                          FaceManagerT& faceManager )
{


  // the virtual FaceManager's rutpureState will be used with the following definitions:
  //   ruptureState = 0 means not reached rupture criteria
  //                = 1 means has reached rupture criteria
  //                = 2 means that face has already been split...but it can still
  //                  be part of a new separation path!!!
  m_virtualFaces.AddKeylessDataField<int>( "ruptureState",true, true );
  m_virtualFaces.AddKeylessDataField<realT>( "separationCoeff",true, true );

  faceManager.AddKeylessDataField<realT>( "stressNOnFace",true, true );
  faceManager.AddKeylessDataField<R1Tensor>( "stressTOnFace",true, true );
  rArray1d& stressNOnFace = faceManager.GetFieldData<realT>("stressNOnFace");
  stressNOnFace = -std::numeric_limits<realT>::max();


  FractunatorBase::RegisterFieldsAndMaps( nodeManager, edgeManager, faceManager );

  //FIXME hack to make 3D hydrofrac work (Fu)
  faceManager.AddKeylessDataField<realT>( "effectiveStressN",true, true );
  //
  if (m_failCriterion > 0)
  {
    edgeManager.AddKeylessDataField<realT>("SIF_I", true, true );
    edgeManager.AddKeylessDataField<realT>("SIF_II", true, true );
    edgeManager.AddKeylessDataField<realT>("SIF_III", true, true );
    faceManager.AddKeylessDataField<realT>("SIFonFace", true, true );
    faceManager.AddKeylessDataField<localIndex>("primaryCandidateFace", true, false );
    edgeManager.AddKeylessDataField<int>("LayersFromDomainBoundary", true, true );
  }


  nodeManager.AddKeylessDataField<R1Tensor>("displacement0", true, false );
  nodeManager.AddKeylessDataField<R1Tensor>("netDisplacement", false, true );

  edgeManager.AddKeylessDataField<realT>("kinkAngle", false, true );


}


void Fractunator3::ReadXML( TICPP::HierarchicalDataNode& hdn )
{
  FractunatorBase::ReadXML(hdn);

  m_insituStress3D = 0.0;
  m_insituStress3D = hdn.GetAttributeVector<realT>("insituStress3D", ",");
  m_failCriterion = hdn.GetAttributeOrDefault<int>("failCriterion",0);
  m_markExtendedLayer =  hdn.GetAttributeOrDefault<int>("markExtendedLayer",1);

  m_saturationPressureCuttoff = hdn.GetAttributeOrDefault<realT>("saturationPressureCutoff",0.01);
  m_faceToEdgeProjectionTol = hdn.GetAttributeOrDefault<realT>("faceToEdgeProjectionTol",0.3);

}



void Fractunator3::Initialize( NodeManagerT& nodeManager,
                               EdgeManagerT& edgeManager,
                               FaceManagerT& faceManager,
                               ElementManagerT& elementManager)
{

  //All nodes are separable unless the input file specifies an explicit separable set
  FractunatorBase::Initialize( nodeManager, edgeManager, faceManager, elementManager );


  if (m_failCriterion >0)  edgeManager.SetLayersFromDomainBoundary( nodeManager );



  // set up the virtual objects used for the fracture path searching routines
  m_virtualNodes.resize( nodeManager.DataLengths() );
  m_virtualNodes.m_nodeToEdgeMap = nodeManager.m_nodeToEdgeMap;
  m_virtualNodes.m_nodeToFaceMap = nodeManager.m_nodeToFaceMap;

  m_virtualEdges.resize( edgeManager.DataLengths() );
  m_virtualEdges.m_toNodesRelation = edgeManager.m_toNodesRelation;
  m_virtualEdges.m_toFacesRelation = edgeManager.m_toFacesRelation;
  m_virtualEdges.m_toFacesRelation = edgeManager.m_toFacesRelation;
  m_virtualEdges.m_isExternal = edgeManager.m_isExternal;



  m_virtualFaces.resize(faceManager.DataLengths());
  m_virtualFaces.m_toEdgesRelation = faceManager.m_toEdgesRelation;
  m_virtualFaces.m_toNodesRelation = faceManager.m_toNodesRelation;
  m_virtualFaces.m_toElementsRelation = faceManager.m_toElementsRelation;
  m_virtualFaces.m_Sets = faceManager.m_Sets;

  //Move reference location to create preexisting stress field.
  //TODO Maybe later we should move this and its counterpart in Fractunator2D to hydrofrac solver
  //Move reference location to create preexisting stress field.
  Array1dT<R1Tensor>& u = nodeManager.GetFieldData<FieldInfo::displacement>();
  Array1dT<R1Tensor>& x0 = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  Array1dT<R1Tensor>& du = nodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();

  Array1dT<R1Tensor>& u0 = nodeManager.GetFieldData<R1Tensor>("displacement0");


  //Move the ref position of nodes to generate the desired stress field.
  // This only works for uniform material properties

  if (!m_insituStress3D.empty())
  {
    const ElementRegionT& elemRegion = (elementManager.m_ElementRegions.begin())->second;
    const MaterialBaseParameterData& parameter = *(elemRegion.m_mat->ParameterData(0));
    const realT E = parameter.E;
    const realT v = parameter.Nu;
    rArray1d strain(6);

    strain[0] = m_insituStress3D[0] - m_insituStress3D[1] * v - m_insituStress3D[2] * v;
    strain[1] = -m_insituStress3D[0] * v + m_insituStress3D[1] - m_insituStress3D[2] * v;
    strain[2] = -m_insituStress3D[0] *v - m_insituStress3D[1] * v + m_insituStress3D[2];
    strain[3] = (1+v) * m_insituStress3D[3];
    strain[4] = (1+v) * m_insituStress3D[4];
    strain[5] = (1+v) * m_insituStress3D[5];

    strain /= E;

    for (localIndex i = 0; i != nodeManager.DataLengths(); ++i)
    {
      R1Tensor x_current = x0[i];
      x0[i][0] = x_current[0] / (1 + strain[0]);
      x0[i][1] = x_current[1] / (1 + strain[1]);
      x0[i][2] = x_current[2] / (1 + strain[2]);
      x0[i][0] = (x0[i][0] - x0[i][1] * strain[5]) / (1 - 0.25 * strain[5] * strain[5]);
      x0[i][0] = (x0[i][0] - x0[i][2] * strain[4]) / (1 - 0.25 * strain[4] * strain[4]);
      x0[i][1] = (x0[i][1] - x0[i][0] * strain[5]) / (1 - 0.25 * strain[5] * strain[5]);
      x0[i][1] = (x0[i][1] - x0[i][2] * strain[3]) / (1 - 0.25 * strain[3] * strain[3]);
      x0[i][2] = (x0[i][2] - x0[i][0] * strain[4]) / (1 - 0.25 * strain[4] * strain[4]);
      x0[i][2] = (x0[i][2] - x0[i][1] * strain[3]) / (1 - 0.25 * strain[3] * strain[3]);

      u[i] = x_current;
      u[i] -= x0[i];
      u0[i] = u[i];
      du[i] = u[i];
    }


    for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator elementRegionIter = elementManager.m_ElementRegions.begin() ;
        elementRegionIter != elementManager.m_ElementRegions.end() ;
        ++elementRegionIter )
    {
      ElementRegionT& elementRegion = elementRegionIter->second;

      elementRegion.CalculateVelocityGradients(nodeManager);

      elementRegion.MaterialUpdate(0.0);

    }


  }
}






void Fractunator3::UpdateRuptureStates( NodeManagerT& nodeManager,
                                        EdgeManagerT& edgeManager,
                                        FaceManagerT& faceManager,
                                        ElementManagerT& elementManager,
                                        Array1dT<lSet>& nodesToRupturedFaces,
                                        Array1dT<lSet>& edgesToRupturedFaces,
                                        const bool prefrac )
{

  if( !prefrac )
  {
    //    m_virtualFaces.UpdateRuptureStates( elementManager, nodeManager, m_separableSet, this->m_failstress );
    faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, this->m_failstress );
  }

  PostUpdateRuptureStates( nodeManager,
                           edgeManager,
                           faceManager,
                           elementManager,
                           nodesToRupturedFaces,
                           edgesToRupturedFaces);

}

void Fractunator3::PostUpdateRuptureStates( NodeManagerT& nodeManager,
                                            EdgeManagerT& edgeManager,
                                            FaceManagerT& faceManager,
                                            ElementManagerT& elementManager,
                                            Array1dT<lSet>& nodesToRupturedFaces,
                                            Array1dT<lSet>& edgesToRupturedFaces)

{
  nodesToRupturedFaces.resize(m_virtualNodes.DataLengths());
  edgesToRupturedFaces.resize(m_virtualEdges.DataLengths());

  iArray1d& faceRuptureState = faceManager.GetFieldData<int>( "ruptureState" );
  iArray1d& vfaceRuptureState = m_virtualFaces.GetFieldData<int>( "ruptureState" );

  // assign the values of the nodeToRupturedFaces and edgeToRupturedFaces arrays.
  for( localIndex kf=0 ; kf<m_virtualFaces.DataLengths() ; ++kf )
  {
    vfaceRuptureState[kf] = faceRuptureState[kf];
    if( vfaceRuptureState[kf] >0 )
    {
      for( localIndex a=0 ; a<m_virtualFaces.m_toNodesRelation[kf].size() ; ++a )
      {
        const localIndex nodeIndex = m_virtualFaces.m_toNodesRelation[kf][a];
        nodesToRupturedFaces[nodeIndex].insert( kf );
      }

      for( localIndex a=0 ; a<m_virtualFaces.m_toEdgesRelation[kf].size() ; ++a )
      {
        const localIndex edgeIndex = m_virtualFaces.m_toEdgesRelation[kf][a];
        edgesToRupturedFaces[edgeIndex].insert( kf );
      }
    }
  }

  for( iArray1d::iterator i=edgeManager.m_isExternal.begin() ; i!=edgeManager.m_isExternal.end() ; ++i )
  {
    if( *i == -1 )
    {
      *i = 1;
      throw GPException("edgeManager.m_isExternal=-1. Call Pengcheng if you see this error");
      // I couldn't figure out why we need this loop. I am putting this exception here so that when it is actually invoked, we will know why we use it.

    }
  }


  iArray1d numberOfRupturedFaces = nodeManager.GetFieldData<int>("numberOfRupturedFaces");

  for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
  {
    const localIndex parentIndex = nodeManager.GetParentIndex( a );
    numberOfRupturedFaces[a] = nodesToRupturedFaces[parentIndex].size();
  }
}

int Fractunator3::CheckEdgeSplitability( const localIndex edgeID,
                                         NodeManagerT& nodeManager,
                                         FaceManagerT& faceManager,
                                         EdgeManagerT& edgeManager,
                                         const bool prefrac)
{
  //     Return value = -1, this edge won't split for sure, don't do any more work;
  //                  = 0, edge is along a tip, but the fracture connected to it is not saturated yet.  We will only calculate SIF but will not perform splitting.
  //                  = 1, edge is along a tip and the adjacent fracture is saturated, more work to be done; or this is a dry simulation
  //                  = 2, this is a singular edge, we need split it.
  //                  = 3, this is an eligible kink, we need to process it as a kink


  int isSplitable = -1;
  const rArray1d* faceFluidPressure = faceManager.GetFieldDataPointer<FieldInfo::pressure>();
  const iArray1d* flowFaceType = faceManager.GetFieldDataPointer<int>("flowFaceType");

  if (edgeManager.m_isExternal[edgeID] == 0)
  {
    isSplitable = -1;
    return (isSplitable);
  }

  // We first count the external faces connected to this edge;
  int nExternalFaces = 0;
  lArray1d faceInvolved;
  for( lSet::const_iterator iface=edgeManager.m_toFacesRelation[edgeID].begin() ;
      iface!=edgeManager.m_toFacesRelation[edgeID].end() ; ++iface )
  {
    if (faceManager.m_isExternal[*iface] == 1)
    {
      nExternalFaces++;
      faceInvolved.push_back(*iface);
    }
  }

  if (nExternalFaces%2 == 1)
  {
    //    char msg[200];
    //    sprintf(msg, "Error! Edge %d has an odd number of external faces.", int(edgeID));
    //    throw GPException(msg);
    //    std::cout << "Error! Edge " << int(edgeID) << " has an odd number of external faces. "
    //        << (*nodeManager.m_refposition)[edgeManager.m_toNodesRelation[edgeID][0]][0] << " ,"
    //        << (*nodeManager.m_refposition)[edgeManager.m_toNodesRelation[edgeID][0]][1] << " ,"
    //        << (*nodeManager.m_refposition)[edgeManager.m_toNodesRelation[edgeID][0]][2] << " ,";
    //    isSplitable = -1;
    return(isSplitable);
  }

  if (nExternalFaces >= 4)
  {
    isSplitable = 2;
    return (isSplitable);
  }

  if (nExternalFaces == 2)
  {
    localIndex parentFace = LOCALINDEX_MAX;
    if (faceManager.m_parentIndex[faceInvolved[0]] == LOCALINDEX_MAX && faceManager.m_parentIndex[faceInvolved[1]] == faceInvolved[0])
    {
      parentFace = faceInvolved[0];
    }
    else if (faceManager.m_parentIndex[faceInvolved[1]] == LOCALINDEX_MAX && faceManager.m_parentIndex[faceInvolved[0]] == faceInvolved[1])
    {
      parentFace = faceInvolved[1];
    }

    if (parentFace == LOCALINDEX_MAX)
    {
      isSplitable = -1;
    }
    else
    {
      if ( faceFluidPressure == NULL || flowFaceType == NULL || prefrac || m_allowVacuumFrac == 1)  //This is a dry simulation.  No fluid involved.
      {
        isSplitable = 1;
      }
      else
      {
        if ((*flowFaceType)[parentFace] > -1 && (*faceFluidPressure)[parentFace] > m_saturationPressureCuttoff)
        {
          isSplitable = 1;
        }
        else
        {
          isSplitable = 0;
        }
      }
    }


  }

  return (isSplitable);
}

int Fractunator3::CheckNodeSplitability( const localIndex nodeID,
                                         NodeManagerT& nodeManager,
                                         FaceManagerT& faceManager,
                                         EdgeManagerT& edgeManager,
                                         const bool prefrac)
{
  return (1);
  //
  //  if (prefrac)
  //  {
  //    return (1);
  //  }
  //  else if (m_failCriterion == 0)
  //  {
  //    return (1);
  //  }
  //  else if (nodeManager.m_isExternal[nodeID] == 1)
  //  {
  //    return (1);
  //  }
  //  else
  //  {
  //    return (0);
  //  }
}

bool Fractunator3::FindFracturePlanes( const localIndex nodeID,
                                       const NodeManagerT& nodeManager,
                                       const EdgeManagerT& edgeManager,
                                       const FaceManagerT& faceManager,
                                       const Array1dT<lSet>& nodesToRupturedFaces,
                                       const Array1dT<lSet>& edgesToRupturedFaces,
                                       lSet& separationPathFaces,
                                       std::map<localIndex,int>& edgeLocations,
                                       std::map<localIndex,int>& faceLocations,
                                       std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{
  const localIndex parentNodeIndex = nodeManager.GetParentIndex(nodeID);

  const lSet& vNodeToRupturedFaces = nodesToRupturedFaces[parentNodeIndex];
  const lSet& vNodeToVFaces = m_virtualNodes.m_nodeToFaceMap[parentNodeIndex];
  const lSet& vNodeToVEdges = m_virtualNodes.m_nodeToEdgeMap[parentNodeIndex];

  const lSet& nodeToEdges = nodeManager.m_nodeToEdgeMap[nodeID];
  const lSet& nodeToFaces = nodeManager.m_nodeToFaceMap[nodeID];
  //  const iArray1d& faceRuptureState = m_virtualFaces.GetFieldData<int>( "ruptureState" );

  const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements = nodeManager.m_toElementsRelation[nodeID] ;


  const iArray1d& isEdgeExternal = m_virtualEdges.m_isExternal;

  const Array1dT<lArray1d>& vFaceToVEdges = m_virtualFaces.m_toEdgesRelation;


  //const lSet& usedFaces = m_virtualNodes.GetUnorderedVariableOneToManyMap("usedFaces")[parentNodeIndex];
  const lSet& usedFaces = nodeManager.GetUnorderedVariableOneToManyMap("usedFaces")[nodeID];

  // **** local working arrays *****************************************************************************************

  // array to hold the faces ready for rupture. It is filled with the intersection of the virtual parent faces associated
  // with all faces attached to the node, and all ruptured virtual faces attached to the virtual parent node.
  lSet vNodeToRuptureReadyVFaces;
  for( lSet::const_iterator i=nodeToFaces.begin() ; i!=nodeToFaces.end() ; ++i )
  {
    const localIndex parentFaceIndex = ( faceManager.m_parentIndex[*i] == LOCALINDEX_MAX ) ? *i : faceManager.m_parentIndex[*i];

    if( vNodeToRupturedFaces.count(parentFaceIndex) > 0 )
    {
      vNodeToRuptureReadyVFaces.insert(parentFaceIndex);
    }
  }


  // local map to hold the edgesToRuptureReadyFaces
  std::map< localIndex, std::set<localIndex> > edgesToRuptureReadyFaces;
  for( lSet::const_iterator edgeIndex=vNodeToVEdges.begin() ; edgeIndex!=vNodeToVEdges.end() ; ++edgeIndex )
  {
    if( !(edgesToRupturedFaces[*edgeIndex].empty()) )
      edgesToRuptureReadyFaces[*edgeIndex].insert( edgesToRupturedFaces[*edgeIndex].begin(), edgesToRupturedFaces[*edgeIndex].end() );
  }


  // need a map from faces to edges that are attached to the node
  std::map< localIndex, std::pair<localIndex,localIndex> > localVFacesToVEdges;
  for( lSet::const_iterator kf=vNodeToVFaces.begin() ; kf!=vNodeToVFaces.end() ; ++kf )
  {
    localIndex edge[2] = { INT_MAX,INT_MAX };
    int count = 0;
    for( lArray1d::const_iterator ke=vFaceToVEdges[*kf].begin() ; ke!=vFaceToVEdges[*kf].end() ; ++ke )
    {
      if( m_virtualEdges.hasNode( *ke, parentNodeIndex ) )
      {
        edge[count++] = *ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      throw GPException("Fractunator3::FindFracturePlanes: invalid edge.");
    }


    localVFacesToVEdges[*kf] = std::make_pair(edge[0],edge[1]);

    if( m_verbose ==2 )
      std::cout<<"localFacesToEdges["<<*kf<<"] = ( "<<localVFacesToVEdges[*kf].first<<", "<<localVFacesToVEdges[*kf].second<<" )"<<std::endl;
  }


  // ***** remove dead end paths ***************************************************************************************
  // if the edge is not external, and the size of edgesToRupturedFaces is less than 2, then the edge is a dead-end
  // as far as a rupture plane is concerned. The face associated with the edge should be removed from the working
  // list of ruptured faces.

  // loop over all the edges
  for( lSet::const_iterator edgeIndex=vNodeToVEdges.begin() ; edgeIndex!=vNodeToVEdges.end() ; ++edgeIndex )
    //  for( std::map< localIndex, std::set<localIndex> >::const_iterator edgeIndex=edgesToRuptureReadyFaces.begin() ; edgeIndex!=edgesToRuptureReadyFaces.end() ; ++edgeIndex )
  {

    CheckForAndRemoveDeadEndPath( *edgeIndex,
                                  isEdgeExternal,
                                  edgesToRuptureReadyFaces,
                                  localVFacesToVEdges,
                                  vNodeToRuptureReadyVFaces);

  }

  // if there are no ruptured faces attached to the node, then we are done.
  // or if there are no faces that have not been used in a rupture path for this node...we are done.
  if( vNodeToRuptureReadyVFaces.empty() )//|| nodeToRuptureReadyFaces.size() == usedFaces.size() )
  {
    return false;
  }

  // ***** find separation path ****************************************************************************************

  // ***** find starting face *****
  // We need to find a starting point for the path. The path must have a face that does has not been used in a previous
  // path for this node...otherwise it is the same path as used previously.
  localIndex startingEdge = INT_MAX;
  localIndex startingFace = INT_MAX;
  bool startingEdgeExternal = false;

  for( lSet::const_iterator i=vNodeToRuptureReadyVFaces.begin() ; i!=vNodeToRuptureReadyVFaces.end() ; ++i )
  {
    // check to see if this face has been used to split this node as part of a previously used path
    if( usedFaces.count(*i)==0 )
    {
      // great! It hasn't. It's on like Donkey Kong.
      startingFace = *i;

      if( isEdgeExternal[localVFacesToVEdges[startingFace].first] )
      {
        startingEdge = localVFacesToVEdges[startingFace].first;
        startingEdgeExternal = true;
        break;
      }
      else if( isEdgeExternal[localVFacesToVEdges[startingFace].second] )
      {
        startingEdge = localVFacesToVEdges[startingFace].second;
        startingEdgeExternal = true;
        break;
      }
      else
      {
        startingEdge = localVFacesToVEdges[startingFace].first;
      }
    }
  }

  // if the starting face was not set, then we don't have a rupture surface....so just quit.
  if( startingFace==INT_MAX || startingEdge==INT_MAX )
  {
    return false;
    //    throw GPException("Fracturantor3::FindFracturePlanes: couldn't set starting face/edge");
  }



  // so now the working arrays have been purged of any faces that are on a dead-end path. All remaining faces
  // are part of a separation plane...of course, there can be more than one...which is bad. We will just take the first
  // path we find, and call this function again after the selected path is processed. Since the ruptureState of a face
  // is set to 2 after it is ruptured, if we enforce that candidate paths must have a face with a ruptureState of 1, then
  // everything will work out. Also since the new nodes that are created will have higher node indices than the
  // current node, they will be checked for separation prior to completion of the separation driver.



  // We now have to define the separation plane over which a node/face/edge will be split, and all elements on one side
  // of the plane get one set of objects, and all elements on the other side get the other set.



  {
    // now we start the process of setting the separation path. Begin by
    localIndex thisEdge = startingEdge;
    localIndex thisFace = startingFace;

    localIndex nextEdge = INT_MAX;
    localIndex nextFace = INT_MAX;

    //localIndex lastEdge = INT_MAX;
    //localIndex lastFace = INT_MAX;

    // the seprationPath is used to hold combinations of edge and face
    std::map<localIndex,int> facesInPath;
    std::map<localIndex,int> edgesInPath;

    int numFacesInPath = 0;
    edgesInPath[thisEdge] = numFacesInPath;
    facesInPath[thisFace] = numFacesInPath++;

    lArray1d facePath;
    lArray1d edgePath;

    facePath.push_back(thisFace);
    edgePath.push_back(thisEdge);

    // now walk from face->edge->face->edge etc. until we get to an external edge, or back to the startingEdge.
    // the breakFlag indicates that we have found a complete separation path
    bool breakFlag = false;
    while ( !breakFlag )
    {

      // get the next edge in the path...it is on the other side of "thisFace", so assign the other edge on the face as
      // the next edge

      nextEdge = GetOtherFaceEdge( localVFacesToVEdges, thisFace,  thisEdge );


      // if the nextEdge has already been used in the path, and the nextEdge is not the starting edge, then we have
      // to take a step back and try a different path
      if( edgesInPath.count(nextEdge)==1 && nextEdge!=startingEdge )
      {
        // first check to see if we can use the path without the preceding
        return false;
      }

      // if we have reached an external face, or the edge is already in the path, then we are done
      if( (isEdgeExternal[nextEdge]==1 && startingEdgeExternal ) || edgesInPath.count(nextEdge)==1 )
      {
        // check to see if nextEdge is the startingEdge. If not, then all faces must that are before the nextEdge must
        // NOT be included in the path!!!
        if( nextEdge!=startingEdge && !(isEdgeExternal[nextEdge]==1 && startingEdgeExternal ) )
        {
          std::cout<<std::endl;


          std::cout<<"  NodeID, ParentID = "<<nodeID<<", "<<parentNodeIndex<<std::endl;
          std::cout<<"  Starting Edge/Face = "<<startingEdge<<", "<<startingFace<<std::endl;
          std::cout<<"  Face Separation Path = ";
          for( lArray1d::const_iterator kf=facePath.begin() ; kf!=facePath.end() ; ++kf )
          {
            std::cout<<*kf<<", ";
          }
          std::cout<<std::endl;

          std::cout<<"  Edge Separation Path = ";
          for( lArray1d::const_iterator kf=edgePath.begin() ; kf!=edgePath.end() ; ++kf )
          {
            std::cout<<*kf<<", ";
          }
          std::cout<<std::endl;


          throw GPException("crap");
        }

        // add faces in the path to separationPathFaces
        for( std::map<localIndex,int>::const_iterator kf=facesInPath.begin() ; kf!=facesInPath.end() ; ++kf )
        {
          separationPathFaces.insert( kf->first );
        }

        // break out of the while loop
        breakFlag = true;
      }
      else
      {
        // if the previous if statement is false, then what if we have reached an external edge, but the starting edge
        // was not external?? This means that we must continue the process from the edge opposite the startingEdge on the
        // startingFace....which is hard-coded as the second entry in localFacesToEdges.
        if( isEdgeExternal[nextEdge]==1 )
        {
          nextEdge = localVFacesToVEdges[startingFace].second;
        }

        // I sure hope that this is true!!
        if( edgesToRuptureReadyFaces[nextEdge].size() > 1 )
        {
          // we need to pick another face attached to the "next edge"
          // increment the face and edge, and add to the separationPathFaces


          {
            // OK...so we have an iterator that points to a candidate face. We prefer to move towards a face that is
            // ruptureState 1, so that we can get as much splitting done in this event. So we will loop over all the
            // faces attached to the edge, and pick one with ruptureState==1, otherwise just pick any one.
            bool pathFound = false;
            std::pair<ElementRegionT*,localIndex> thisElem0 = this->m_virtualFaces.m_toElementsRelation[thisFace][0];
            std::pair<ElementRegionT*,localIndex> thisElem1 = this->m_virtualFaces.m_toElementsRelation[thisFace][1];


            // nextFaceQuality is intended to keep how desirable a face is for the rupture path.
            // A value of:
            //    0 -> the face is kind of creppy
            //    1 -> the face is does not turn a corner around the elements surrounding thisFace
            //    2 -> the face has not been used in a separation path
            //    3 -> a combination of 1 and 2.
            //    4 -> other edge on the face is the startingEdge.
            //
            int nextFaceQuality = -1;

            for( std::set<localIndex>::const_iterator iter_edgeToFace = edgesToRuptureReadyFaces[nextEdge].begin() ;
                iter_edgeToFace!=edgesToRuptureReadyFaces[nextEdge].end() ; ++iter_edgeToFace )
            {
              if( *iter_edgeToFace != thisFace )
              {
                pathFound = true;




                const localIndex candidateFaceIndex = *iter_edgeToFace;
                int candidateFaceQuality = 0;


                localIndex candidateEdgeIndex = GetOtherFaceEdge( localVFacesToVEdges, candidateFaceIndex,  nextEdge );
                if( candidateEdgeIndex == startingEdge )
                {
                  nextFace = candidateFaceIndex;
                  break;
                }


                std::pair<ElementRegionT*,localIndex> nextElem0 = this->m_virtualFaces.m_toElementsRelation[candidateFaceIndex][0];
                std::pair<ElementRegionT*,localIndex> nextElem1 = this->m_virtualFaces.m_toElementsRelation[candidateFaceIndex][1];

                if( thisElem0 != nextElem0 && thisElem0 != nextElem1 &&
                    thisElem1 != nextElem0 && thisElem1 != nextElem1 )
                {
                  candidateFaceQuality += 1;
                }

                if( usedFaces.count(candidateFaceIndex) == 0 )
                {
                  candidateFaceQuality += 2;
                }


                if( candidateFaceQuality > nextFaceQuality )
                {
                  nextFace = candidateFaceIndex;
                  nextFaceQuality = candidateFaceQuality;
                }

                if( candidateFaceQuality == 3 )
                {
                  break;
                }
              }
            }
            if( pathFound == false )
              throw GPException("Fractunator3::FindFracturePlanes: couldn't find the next face in the rupture path");
          }

          //        lastEdge = thisEdge;
          //        lastFace = thisFace;

          thisEdge = nextEdge;
          thisFace = nextFace;
          //      separationPathFaces.insert( thisFace );
          edgesInPath[thisEdge] = numFacesInPath;
          facesInPath[thisFace] = numFacesInPath++;

          facePath.push_back(thisFace);
          edgePath.push_back(thisEdge);

        }
        else
        {
          throw GPException("Fractunator3::next edge in separation path is apparently  connected to less than 2 ruptured face");
        }

      }
    }
  }


  //***** SET LOCATIONS ************************************************************************************************



  // need a map from faces to edges that are attached to the node
  std::map< localIndex, std::pair<localIndex,localIndex> > localFacesToEdges;
  for( lSet::const_iterator kf=nodeToFaces.begin() ; kf!=nodeToFaces.end() ; ++kf )
  {
    localIndex edge[2] = { INT_MAX,INT_MAX };
    int count = 0;
    for( lArray1d::const_iterator ke=faceManager.m_toEdgesRelation[*kf].begin() ; ke!=faceManager.m_toEdgesRelation[*kf].end() ; ++ke )
    {
      if( edgeManager.hasNode( *ke, nodeID ) )
      {
        edge[count++] = *ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      throw GPException("Fractunator3::FindFracturePlanes: invalid edge.");
    }


    localFacesToEdges[*kf] = std::make_pair(edge[0],edge[1]);

  }


  // now we want to identify the objects on either side of the separation plane. First we assign an array to indicate
  // whether a face/edge is on the fracture plane.

  for( lSet::const_iterator kf=nodeToFaces.begin() ; kf!=nodeToFaces.end() ; ++kf )
  {
    // iff the face is being split NOW, the set the faceLocation = -1.
    const localIndex virtualFaceIndex = ( faceManager.m_parentIndex[*kf] == LOCALINDEX_MAX ) ? *kf : faceManager.m_parentIndex[*kf];
    if( *kf == virtualFaceIndex && faceManager.m_childIndices[*kf].empty() && separationPathFaces.count(*kf) )
    {
      faceLocations[*kf] = -1;
    }
    else
    {
      faceLocations[*kf] = INT_MIN;
    }

  }

  for( lSet::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
  {
    edgeLocations[*ke] = INT_MIN;
  }

  for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodesToElements.begin() ; k!=nodesToElements.end() ; ++k )
  {
    elemLocations[*k] = INT_MIN;
  }



  /*
  SetLocations( 0, separationPathFaces, faceManager, nodesToElements, localFacesToEdges, //nodeToEdges,
                edgeLocations, faceLocations, elemLocations );

  if( !(SetLocations( 1, separationPathFaces, faceManager, nodesToElements, localFacesToEdges, //nodeToEdges,
                      edgeLocations, faceLocations, elemLocations )) )
  {
    return false;
  }*/

  SetLocations( separationPathFaces, faceManager, nodesToElements, localFacesToEdges,
                edgeLocations, faceLocations, elemLocations );



  bool fail = false;

  for( lSet::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
  {
    if( edgeLocations[*ke] == INT_MIN )
    {
      fail = true;
    }
  }
  for( lSet::const_iterator ke=nodeToFaces.begin() ; ke!=nodeToFaces.end() ; ++ke )
  {
    if( faceLocations[*ke] == INT_MIN )
    {
      fail = true;
    }
  }
  /*
    std::cout<<"  NodeID, ParentID = "<<nodeID<<", "<<parentNodeIndex<<std::endl;
    std::cout<<"  separation path = ";
    for( lSet::const_iterator kf=separationPathFaces.begin() ; kf!=separationPathFaces.end() ; ++kf )
    {
      std::cout<<*kf<<", ";
    }
    std::cout<<std::endl;

    std::cout<<"  Starting Edge/Face = "<<startingEdge<<", "<<startingFace<<std::endl;
    for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodesToElements.begin() ; k!=nodesToElements.end() ; ++k )
    {
      std::cout<<"  elemLocations["<<k->second<<"] = "<<elemLocations[*k]<<std::endl;
    }

    for( lSet::const_iterator ke=nodeToFaces.begin() ; ke!=nodeToFaces.end() ; ++ke )
    {
      std::cout<<"  faceLocations["<<*ke<<"] = "<<faceLocations[*ke]<<std::endl;
    }

    for( lSet::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
    {
      std::cout<<"  edgeLocations["<<*ke<<"] = "<<edgeLocations[*ke]<<std::endl;
    }
   */
  if( fail )
  {

    //    throw GPException("Fractunator3::FindFracturePlanes: unset element,face, or edge");
    return false;
  }
  return true;
}






bool Fractunator3::SetLocations( const lSet& separationPathFaces,
                                 const FaceManagerT& faceManager,
                                 const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements,
                                 const std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                                 std::map<localIndex,int>& edgeLocations,
                                 std::map<localIndex,int>& faceLocations,
                                 std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{
  bool rval = true;
  //  const localIndex separationFace = *(separationPathFaces.begin());

  // insert an element attached to the separation face
  //  std::pair<ElementRegionT*,localIndex> elem0 = m_virtualFaces.m_FaceToElementMap[separationFace][0] ;

  std::pair<ElementRegionT*,localIndex> elem0 = *(nodesToElements.begin()) ;


  SetElemLocations( 0,
                    elem0,
                    separationPathFaces,
                    faceManager,
                    nodesToElements,
                    localFacesToEdges,
                    edgeLocations,
                    faceLocations,
                    elemLocations );

  return rval;
}

bool Fractunator3::SetElemLocations( const int location,
                                     const std::pair< ElementRegionT*, localIndex >& k,
                                     const lSet& separationPathFaces,
                                     const FaceManagerT& faceManager,
                                     const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements,
                                     const std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                                     std::map<localIndex,int>& edgeLocations,
                                     std::map<localIndex,int>& faceLocations,
                                     std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{

  const int otherlocation = (location==0) ? 1 : 0;

  elemLocations[k] = location;


  // loop over all faces on the element
  for( localIndex kf=0 ; kf<k.first->m_toFacesRelation.Dimension(1) ; ++kf )
  {

    // define the actual face index, and the virtual face index
    const localIndex faceIndex = k.first->m_toFacesRelation(k.second,kf);
    const localIndex virtualFaceIndex = ( faceManager.m_parentIndex[faceIndex] == LOCALINDEX_MAX ) ? faceIndex : faceManager.m_parentIndex[faceIndex];


    // see if we can find the face in the faceLocations array.
    std::map<localIndex,int>::iterator iterFace = faceLocations.find(faceIndex);
    // if we can find the face in the faceLocations array, then we must process the face, otherwise it is not
    // connected to the node, so we do nothing.
    if( iterFace != faceLocations.end() )
    {

      if( faceLocations[faceIndex]==otherlocation )
        faceLocations[faceIndex] = -1;
      else if( faceLocations[faceIndex] == INT_MIN )
        faceLocations[faceIndex] = location;

      std::map< localIndex, std::pair<localIndex,localIndex> >::const_iterator iterF2E = localFacesToEdges.find(faceIndex);

      if( iterF2E != localFacesToEdges.end() )
      {
        const localIndex edge0 = (iterF2E->second).first;
        const localIndex edge1 = (iterF2E->second).second;

        if( edgeLocations[edge0]==otherlocation )
          edgeLocations[edge0] = -1;
        else if( edgeLocations[edge0] == INT_MIN )
          edgeLocations[edge0] = location;

        if( edgeLocations[edge1]==otherlocation )
          edgeLocations[edge1] = -1;
        else if( edgeLocations[edge1] == INT_MIN )
          edgeLocations[edge1] = location;

      }




      // now we add the element that is a neighbor to the face
      // of course, this only happens if there are more than one element
      // attached to the face.
      if( m_virtualFaces.m_toElementsRelation[virtualFaceIndex].size() > 1 )
      {


        const std::pair<ElementRegionT*,localIndex>& elemIndex0 = m_virtualFaces.m_toElementsRelation[virtualFaceIndex][0];
        const std::pair<ElementRegionT*,localIndex>& elemIndex1 = m_virtualFaces.m_toElementsRelation[virtualFaceIndex][1];

        const std::pair<ElementRegionT*,localIndex>& nextElem = ( elemIndex0 == k ) ? elemIndex1 : elemIndex0;
        const int nextLocation = (separationPathFaces.count(virtualFaceIndex)==0) ? location : otherlocation;

        // if the first element is the one we are on, and the element is attached
        // to the splitting node, then add the second element to the list.
        if( nodesToElements.find(nextElem)!=nodesToElements.end() )
        {
          if( elemLocations[nextElem]==INT_MIN )
          {
            SetElemLocations( nextLocation,
                              nextElem,
                              separationPathFaces,
                              faceManager,
                              nodesToElements,
                              localFacesToEdges,
                              edgeLocations,
                              faceLocations,
                              elemLocations );
          }
        }
      }
    }
  }

  return true;
}

void Fractunator3::PerformFracture( const localIndex nodeID,
                                    NodeManagerT& nodeManager,
                                    EdgeManagerT& edgeManager,
                                    FaceManagerT& faceManager,
                                    ExternalFaceManagerT& externalFaceManager,
                                    ElementManagerT& elementManager,
                                    ModifiedObjectLists& modifiedObjects,
                                    Array1dT<lSet>& nodesToRupturedFaces ,
                                    Array1dT<lSet>& edgesToRupturedFaces ,
                                    const lSet& separationPathFaces,
                                    const std::map<localIndex,int>& edgeLocations,
                                    const std::map<localIndex,int>& faceLocations,
                                    const std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{

  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );



  iArray1d* flowEdgeType = edgeManager.GetFieldDataPointer<int>("flowEdgeType");
  iArray1d* flowFaceType = faceManager.GetFieldDataPointer<int>("flowFaceType");
  rArray1d& stressNOnFace = faceManager.GetFieldData<realT>("stressNOnFace");
  Array1dT<R1Tensor>& stressTOnFace = faceManager.GetFieldData<R1Tensor>("stressTOnFace");
  Array1dT<lSet>* edgeToFlowFaces = edgeManager.GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");

  //  const Array1dT<lArray1d>& childEdgeIndex = edgeManager.GetVariableOneToManyMap( "childIndices" );
  const OrderedVariableOneToManyRelation& childFaceIndex = faceManager.GetVariableOneToManyMap( "childIndices" );

  //  const localIndex parentNodeIndex = nodeManager.GetParentIndex( nodeID );


  //  lSet& usedFaces = m_virtualNodes.GetUnorderedVariableOneToManyMap("usedFaces")[parentNodeIndex];
  //  lSet& usedFaces = nodeManager.GetUnorderedVariableOneToManyMap("usedFaces")[nodeID];
  Array1dT<lSet>& usedFaces = nodeManager.GetUnorderedVariableOneToManyMap("usedFaces");
  usedFaces[nodeID].insert( separationPathFaces.begin(), separationPathFaces.end() );


  rArray1d* delta0N = faceManager.GetFieldDataPointer<realT>("delta0N");
  Array1dT<R1Tensor>* stressShear0 = faceManager.GetFieldDataPointer<R1Tensor>("stressShear0");


  // ***** split all the objects first *****

  // Split the node into two, using the original index, and a new one.
  localIndex newNodeIndex;
  if( m_verbose )
  {
    std::cout<<"\nSplitting node "<<nodeID<<" along separation plane faces ";
    for( lSet::const_iterator i=separationPathFaces.begin() ; i!=separationPathFaces.end() ; ++i )
    {
      std::cout<<*i<<", ";
    }
    std::cout<<std::endl;
  }


  nodeManager.SplitObject(nodeID, rank, newNodeIndex);
  modifiedObjects.newNodes.insert( newNodeIndex );
  modifiedObjects.modifiedNodes.insert( nodeID );

  //TODO HACK...should recalculate mass
  const realT newMass = 0.5 * (*nodeManager.m_mass)[nodeID];
  (*nodeManager.m_mass)[nodeID] = newMass;
  (*nodeManager.m_mass)[newNodeIndex] = newMass;

  lSet& usedFacesNew = nodeManager.GetUnorderedVariableOneToManyMap("usedFaces")[newNodeIndex];
  usedFacesNew = usedFaces[nodeID];


  if( m_verbose ) std::cout<<"\nDone splitting node "<<nodeID<<" into nodes "<<nodeID<<" and "<<newNodeIndex<<std::endl;

  // HACK...the node to element map is a bastard-child and is not managed by the
  // database
  nodeManager.m_toElementsRelation.resize( nodeManager.m_numNodes );

  // split edges
  std::map<localIndex,localIndex> splitEdges;
  // loop over all edges connected to the node
  for( std::map<localIndex,int>::const_iterator iter_edge=edgeLocations.begin() ; iter_edge!=edgeLocations.end() ; ++iter_edge )
  {
    const localIndex& parentEdgeIndex = iter_edge->first;
    const int& location = iter_edge->second;

    // if the edge is on the separation plane, then split it
    if( location == -1  )
    {
      localIndex newEdgeIndex;

      edgeManager.SplitObject(parentEdgeIndex, rank, newEdgeIndex);

      if( m_verbose ) std::cout<<"  Split edge "<<parentEdgeIndex<<" into edges "<<parentEdgeIndex<<" and "<<newEdgeIndex<<std::endl;

      splitEdges[parentEdgeIndex] = newEdgeIndex;
      modifiedObjects.newEdges.insert( newEdgeIndex );
      modifiedObjects.modifiedEdges.insert( parentEdgeIndex );


      for( int a=0 ; a<2 ; ++a )
      {
        edgeManager.m_toNodesRelation(newEdgeIndex,a) = edgeManager.m_toNodesRelation(parentEdgeIndex,a);
      }

      if( flowEdgeType != NULL )
      {
        (*flowEdgeType)[newEdgeIndex] = -1;
        (*flowEdgeType)[parentEdgeIndex] = 0;
      }


    } //    if( location == -1  )
  } // for( std::map<localIndex,int>::const_iterator iter_edge...





  // split the faces
  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  iArray1d& vruptureState = m_virtualFaces.GetFieldData<int>("ruptureState");
  std::map<localIndex,localIndex> splitFaces;

  // loop over all faces attached to the nodeID
  for( std::map<localIndex,int>::const_iterator iter_face=faceLocations.begin() ; iter_face!=faceLocations.end() ; ++iter_face )
  {
    const localIndex faceIndex = iter_face->first;
    const localIndex parentFaceIndex = faceManager.GetParentIndex( faceIndex );
    const int location = iter_face->second;
    // if the face is on the separation plane, then split it
    if( location == -1 && vruptureState[parentFaceIndex] == 1 )
    {
      localIndex newFaceIndex;

      if( faceManager.SplitObject( parentFaceIndex, rank, newFaceIndex ) )
      {
        if (delta0N != NULL)
        {
          if (stressNOnFace[parentFaceIndex] < 0)
          {
            // We pass the contact stiffness in from the solver via field dleta0 itself.
            (*delta0N)[parentFaceIndex] = -stressNOnFace[parentFaceIndex] /(*delta0N)[parentFaceIndex];
            (*stressShear0)[parentFaceIndex] = stressTOnFace[parentFaceIndex];
          }
          else
          {
            (*delta0N)[parentFaceIndex] = 0.0;
            (*stressShear0)[parentFaceIndex] = 0.0;
          }
        }


        if( m_verbose ) std::cout<<"  Split face "<<parentFaceIndex<<" into faces "<<parentFaceIndex<<" and "<<newFaceIndex<<std::endl;

        splitFaces[parentFaceIndex] = newFaceIndex;
        modifiedObjects.newFaces.insert( newFaceIndex );
        modifiedObjects.modifiedFaces.insert( parentFaceIndex );

        vruptureState[parentFaceIndex] = 2;
        ruptureState[parentFaceIndex] = 2;


        faceManager.m_toEdgesRelation[newFaceIndex] = faceManager.m_toEdgesRelation[parentFaceIndex];

        faceManager.m_toNodesRelation[newFaceIndex] = faceManager.m_toNodesRelation[parentFaceIndex];

        faceManager.m_externalFaces.insert(newFaceIndex );
        faceManager.m_externalFaces.insert(parentFaceIndex );

        if( flowFaceType != NULL )
        {
          (*flowFaceType)[newFaceIndex] = -1;
          (*flowFaceType)[parentFaceIndex] = 0;
        }

        if( edgeToFlowFaces!= NULL )
        {
          for( lArray1d::const_iterator ke=m_virtualFaces.m_toEdgesRelation[parentFaceIndex].begin() ; ke!=m_virtualFaces.m_toEdgesRelation[parentFaceIndex].end() ; ++ke )
          {
            //std::cout<<*ke<<", "<<parentFaceIndex<<std::endl;
            (*edgeToFlowFaces)[*ke].insert(parentFaceIndex);
          }
        }

        // Fu: All edges of the parent face should be external now.
        // We have to do the following because isExternal attribute of the tip edge is not handled by the splitter.
        for( lArray1d::iterator j = faceManager.m_toEdgesRelation[parentFaceIndex].begin() ;
            j!=faceManager.m_toEdgesRelation[parentFaceIndex].end() ; ++j )
        {
          edgeManager.m_isExternal[*j] = 1;
        }
        for( lArray1d::iterator j = faceManager.m_toNodesRelation[parentFaceIndex].begin() ;
            j!=faceManager.m_toNodesRelation[parentFaceIndex].end() ; ++j )
        {
          nodeManager.m_isExternal[*j] = 1;
        }

        externalFaceManager.SplitFace(parentFaceIndex, newFaceIndex, nodeManager);

      } // if( faceManager.SplitObject( faceIndex, newFaceIndex ) )
    } // if( location == -1 )
  } // for( std::map<localIndex,int>::const_iterator iter_face

  // HACK...the node to element map is a bastard-child and is not managed by the
  // database
  faceManager.m_toElementsRelation.resize( faceManager.m_numFaces );



  // ***** now correct all the relations between the objects *****

  /* To accomplish this annoying yet exceedingly important task, we will take a "top down"
   * approach. Note that this is a two way correction, i.e. if we are correcting
   * elementToNodes, we also correct nodesToElements. This is summarized as:
   * 1) Loop over elements attached to the split node.
   *     2a) correct all relations between the single  element and the nodes.
   *     2b) Loop over all faces on the element
   *         3a) For each face, correct the face relations with the element
   *         3b) For each face, correct the face relations with the nodes
   *         3c) Loop over all edges on the face
   *             4a) For each edge, correct the face relations
   *             4b) for each edge, correct the node relations
   *
   *  The element location will define which side of the rupture everything
   *  is on.
   *  - location 0 gets the original node,edge,face.
   *  - location 1 gets the new node,edge,face.
   */

  // 1) loop over all elements attached to the nodeID
  for( std::map<std::pair<ElementRegionT*, localIndex>, int>::const_iterator iter_elem =
      elemLocations.begin() ; iter_elem != elemLocations.end() ; ++iter_elem )
  {
    const int& location = iter_elem->second;

    if( location==1 )
    {
      const std::pair< ElementRegionT*, localIndex >& elem = iter_elem->first;

      ElementRegionT& elemRegion = *(elem.first);
      const localIndex elemIndex = elem.second;

      modifiedObjects.modifiedElements[elemRegion.m_regionName].insert(elemIndex);

      if( m_verbose ) std::cout<<"Element "<<elemIndex<<std::endl;

      // 2a) correct elementToNode and nodeToElement
      if( m_verbose )  std::cout<<"  Looping over all nodes on element, and correcting node<->element maps:"<<std::endl;

      {
        // loop over all nodes on element
        if( m_verbose ) std::cout<<"    m_ElementToNodeMap = ( ";
        localIndex* const elementToNodeMap = elemRegion.m_toNodesRelation[elemIndex];
        for( localIndex a=0 ; a<elemRegion.m_toNodesRelation.Dimension(1) ; ++a )
        {
          // if the node was just split
          if( elementToNodeMap[a] == nodeID )
          {

            if( m_verbose )
              std::cout<<elementToNodeMap[a]<<"->"<<newNodeIndex<<", ";

            elementToNodeMap[a] = newNodeIndex;

            nodeManager.m_toElementsRelation[newNodeIndex].insert(elem);
            nodeManager.m_toElementsRelation[nodeID].erase(elem);

          }
          else
            if( m_verbose ) std::cout<<elementToNodeMap[a]<<", ";
        }
        if( m_verbose ) std::cout<<")"<<std::endl;

        if( m_verbose )
        {
          for( localIndex a=0 ; a<elemRegion.m_toNodesRelation.Dimension(1) ; ++a )
          {
            if( m_verbose ) std::cout<<"    nodeToElemMaps["<<elementToNodeMap[a]<<"] = ( ";
            for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[elementToNodeMap[a]].begin() ;
                k!=nodeManager.m_toElementsRelation[elementToNodeMap[a]].end() ; ++k )
            {
              std::cout<<k->second<<", ";
            }
            std::cout<<" )"<<std::endl;

          }
        }
      }



      // 2b) loop over all faces on element.
      if( m_verbose )
      {
        std::cout<<"  Looping over all faces on element (parent and child):"<<std::endl;
      }


      localIndex* const elemToFaces = elemRegion.m_toFacesRelation[elemIndex];
      // we need to build a list of faces that is elemToFaces FOLLOWED by any
      // parent face of those indicated in elemToFaces

      // Now we do a loop over the facelist and process all the faces
      for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
      {

        // set both faceID and newFaceID to the parent face.
        const localIndex faceIndex = elemToFaces[kf];
        const bool isNewFace = (splitFaces.count(faceIndex)>0) ? true : false;
        const localIndex newFaceIndex = isNewFace ? childFaceIndex[faceIndex][0] : faceIndex;

        //        std::map<localIndex,localIndex>::iterator iterSplitFace = splitFaces.find(faceIndex);

        //        const bool isNewFace = (splitFaces.count(faceIndex)>0) ? true : false;


        // 3a) check to see if the face was split. If so, then we will need
        // to alter the face relation with the elements in both directions.
        if( isNewFace )
        {

          // replace the parent face with the child face in elementToFace. Now
          // faceID is the parent face, and newFaceID is the child face.
          elemToFaces[kf] = childFaceIndex[ faceIndex ][0];

          // add the element to the child faceToElem
          faceManager.m_toElementsRelation[newFaceIndex].push_back( elem );

          // remove the element from the parent face
          if( faceManager.m_toElementsRelation[faceIndex][0] == elem )
          {
            faceManager.m_toElementsRelation[faceIndex].erase( faceManager.m_toElementsRelation[faceIndex].begin() );
          }
          else if( faceManager.m_toElementsRelation[faceIndex][1] == elem )
          {
            faceManager.m_toElementsRelation[faceIndex].erase( faceManager.m_toElementsRelation[faceIndex].begin()+1 );
          }


          faceManager.SortFaceNodes( nodeManager, faceIndex );
          faceManager.SortFaceNodes( nodeManager, newFaceIndex );



        } // if( splitFaces.count( faceID ) > 0 )

        modifiedObjects.modifiedFaces.insert( faceIndex );




        // 3b) correct faceToNodes and nodeToFaces

        if( m_verbose )
        {
          const localIndex parentFace = faceManager.m_parentIndex[newFaceIndex];
          if( parentFace!=LOCALINDEX_MAX )
          {
            std::cout<<"    m_FaceToNodeMap["<<parentFace<<"->"<<newFaceIndex<<"] = ( ";
          }
          else
          {
            std::cout<<"    m_FaceToNodeMap["<<newFaceIndex<<"] = ( ";
          }
        }

        // loop over all nodes on the face.
        for( lArray1d::iterator nodeIndex=faceManager.m_toNodesRelation[newFaceIndex].begin() ;
            nodeIndex!=faceManager.m_toNodesRelation[newFaceIndex].end() ; ++nodeIndex )
        {
          if( m_verbose ) std::cout<<*nodeIndex;

          // if the facenode is the one that is being split
          if( *nodeIndex == nodeID )
          {
            *nodeIndex = newNodeIndex;

            // if it is not a new face.
            if( !isNewFace )
            {
              // remove the face from the nodeToFaceMap of the parent node.
              nodeManager.m_nodeToFaceMap[nodeID].erase(faceIndex);

              // add the face to the nodeToFaceMap of the new node.
              nodeManager.m_nodeToFaceMap[*nodeIndex].insert(faceIndex);

            }
            else
            {
              // it is a new face

              // insert the newFace into the nodeToFaceMap of the newNode
              nodeManager.m_nodeToFaceMap[*nodeIndex].insert(newFaceIndex);


            }

            if( m_verbose ) std::cout<<"->"<<*nodeIndex<<", ";

          }
          else // the node is not being split
          {
            nodeManager.m_nodeToFaceMap[*nodeIndex].insert(newFaceIndex);

            //            if( faceRuptureState[newFaceIndex ] )
            //              nodeToRupturedFaces[*nodeIndex].push_back(newFaceIndex);

            if( m_verbose ) std::cout<<", ";
          }

        }
        if( m_verbose ) std::cout<<")"<<std::endl;




        // faceToEdges
        if( m_verbose )
        {
          const localIndex parentFace = faceManager.m_parentIndex[newFaceIndex];
          if( parentFace!=LOCALINDEX_MAX )
          {
            std::cout<<"    m_FaceToEdgeMap["<<parentFace<<"->"<<newFaceIndex<<"] = ( ";
          }
          else
          {
            std::cout<<"    m_FaceToEdgeMap["<<newFaceIndex<<"] = ( ";
          }
        }
        // loop over all edges on face
        for( lArray1d::iterator edgeIndex=faceManager.m_toEdgesRelation[newFaceIndex].begin() ; edgeIndex!=faceManager.m_toEdgesRelation[newFaceIndex].end() ; ++edgeIndex )
        {

          // if the edge was just split
          if( splitEdges.count( *edgeIndex ) > 0 )
          {
            if( faceIndex == newFaceIndex )
              edgeManager.m_toFacesRelation[*edgeIndex].erase(faceIndex);

            *edgeIndex = splitEdges[*edgeIndex];
          }
          edgeManager.m_toFacesRelation[*edgeIndex].insert(newFaceIndex);
          modifiedObjects.modifiedEdges.insert( *edgeIndex );

          if( m_verbose ) std::cout<<*edgeIndex;



          //edgesToNodes
          if( m_verbose )
          {
            std::cout<<"(";
          }

          {
            localIndex* const nodeIndex = edgeManager.m_toNodesRelation[*edgeIndex];
            for( unsigned int a=0 ; a<edgeManager.m_toNodesRelation.Dimension(1) ; ++a )
            {
              if( nodeIndex[a] == nodeID )
              {

                if( m_verbose ) std::cout<<nodeIndex[a];

                nodeIndex[a] = newNodeIndex;
                nodeManager.m_nodeToEdgeMap[nodeID].erase(*edgeIndex);

                if( m_verbose ) std::cout<<"->"<<nodeIndex[a]<<", ";

              }
              else
                if( m_verbose ) std::cout<<nodeIndex[a]<<", ";

              nodeManager.m_nodeToEdgeMap[nodeIndex[a]].insert(*edgeIndex);
            }
            if( m_verbose ) std::cout<<")";
          }
          if( m_verbose ) std::cout<<", ";
        }
        if( m_verbose ) std::cout<<")"<<std::endl;
      } // for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
    } // if( location==1 )
  } // for( std::map<std::pair<ElementRegionT*, localIndex>, int>::const_iterator iter_elem = elemLocations.begin()



  //**************************************************************************
  // THIS IS ALL JUST CONSISTENCY CHECKING
  //**************************************************************************


  if( m_verbose == 1 )
  {
    std::cout<<"CONSISTENCY CHECKING OF THE MAPS"<<std::endl;

    for( std::map< std::pair< ElementRegionT*, localIndex >, int>::const_iterator iter_elem=elemLocations.begin() ; iter_elem!=elemLocations.end() ; ++iter_elem )
    {
      const std::pair< ElementRegionT*, localIndex >& elem = iter_elem->first;

      ElementRegionT& elemRegion = *(elem.first);
      const localIndex elemIndex = elem.second;


      lSet elemNodes;


      std::cout<<"Element "<<elemIndex<<"\n";
      std::cout<<" elementToNodes = ";
      for( int a=0; a<8 ; ++a )
      {
        elemNodes.insert(elemRegion.m_toNodesRelation(elemIndex,a));
        std::cout<<elemRegion.m_toNodesRelation(elemIndex,a)<<", ";
      }
      std::cout<<std::endl;

      std::cout<<" elementToFaces->edges->nodes = ";


      // Now we do a loop over the facelist and process all the faces
      for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
      {
        lSet faceNodes;

        localIndex faceIndex  = elemRegion.m_toFacesRelation(elemIndex,kf);

        if( kf>0 )
          std::cout<<"                              = ";


        std::cout<<faceIndex<<"( ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex faceNodeID = faceManager.m_toNodesRelation[faceIndex][b];
          faceNodes.insert(faceNodeID);
          if( elemNodes.count(faceNodeID) == 0 && kf<elemRegion.m_numFacesPerElement )
            std::cout<<"*";
          std::cout<<faceNodeID<<",";
        }
        std::cout<<" )      ";



        std::cout<<faceIndex<<"[ ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex edgeIndex = faceManager.m_toEdgesRelation[faceIndex][b];
          std::cout<<edgeIndex<<"( ";
          for( int c=0 ; c<2 ; ++c )
          {
            localIndex edgeNodeID = edgeManager.m_toNodesRelation(edgeIndex,c);
            if( elemNodes.count(edgeNodeID) == 0  && kf<elemRegion.m_numFacesPerElement )
              std::cout<<"*";
            if( faceNodes.count(edgeNodeID) == 0 )
              std::cout<<"#";
            std::cout<<edgeNodeID<<",";
          }
          std::cout<<" ), ";
        }
        std::cout<<" ] \n";

      }
      std::cout<<std::endl;

    }

  }

  if( m_verbose == 2 )
  {
    std::cout<<" elementToFaces->edges->nodes = ";
    //      for( int a=0; a<6 ; ++a )
    {
      lSet faceNodes;
      localIndex faceIndex = 73847;

      //        if( a>0 )
      //          std::cout<<"                              = ";

      std::cout<<faceIndex<<"( ";
      for( int b=0 ; b<4 ; ++b )
      {
        localIndex faceNodeID = faceManager.m_toNodesRelation[faceIndex][b];
        faceNodes.insert(faceNodeID);
        //          if( elemNodes.count(nodeID) == 0 )
        //            std::cout<<"*";
        std::cout<<faceNodeID<<",";
      }
      std::cout<<" )      ";



      std::cout<<faceIndex<<"[ ";
      for( int b=0 ; b<4 ; ++b )
      {
        localIndex edgeIndex = faceManager.m_toEdgesRelation[faceIndex][b];
        std::cout<<edgeIndex<<"( ";
        for( int c=0 ; c<2 ; ++c )
        {
          localIndex edgeNodeID = edgeManager.m_toNodesRelation(edgeIndex,c);
          //            if( elemNodes.count(nodeID) == 0 )
          //              std::cout<<"*";
          if( faceNodes.count(edgeNodeID) == 0 )
            std::cout<<"#";
          std::cout<<edgeNodeID<<",";
        }
        std::cout<<" ), ";
      }
      std::cout<<" ] \n";

    }
    std::cout<<std::endl;

  }


  if( m_verbose == 2 )
  {
    // nodeToEdge
    Array1dT<lSet> tempNodesToEdges( nodeManager.m_numNodes );

    for( localIndex ke=0 ; ke<edgeManager.DataLengths() ; ++ke )
    {
      for( localIndex b= 0 ; b<edgeManager.m_toNodesRelation.Dimension(1) ; ++b )
      {
        localIndex nodeIndex = edgeManager.m_toNodesRelation(ke,b);
        tempNodesToEdges[nodeIndex].insert(ke);
      }
    }
    std::cout<<"Check NodeToEdge "<<std::endl;
    for( localIndex a=0 ; a<nodeManager.m_numNodes ; ++a )
    {
      std::cout<<"m_nodesToEdges["<<a<<"] = ( ";
      for( lSet::const_iterator iedge=nodeManager.m_nodeToEdgeMap[a].begin() ;
          iedge!=nodeManager.m_nodeToEdgeMap[a].end() ; ++iedge )
      {
        if( tempNodesToEdges[a].count(*iedge) == 0 )
          std::cout<<"*";

        std::cout<<*iedge<<", ";
      }
      std::cout<<")    (";

      for( lSet::const_iterator iedge=tempNodesToEdges[a].begin() ;
          iedge!=tempNodesToEdges[a].end() ; ++iedge )
      {
        if( nodeManager.m_nodeToEdgeMap[a].count(*iedge) == 0 )
          std::cout<<"*";

        std::cout<<*iedge<<", ";
      }
      std::cout<<")"<<std::endl;
    }


  }

  if( m_verbose == 2 )
  {
    // nodeToFace
    Array1dT<lSet> tempNodesToFaces( nodeManager.m_numNodes );
    for( localIndex kf=0 ; kf<faceManager.m_numFaces ; ++kf )
    {
      for( lArray1d::const_iterator b=faceManager.m_toNodesRelation[kf].begin() ;
          b!=faceManager.m_toNodesRelation[kf].end() ; ++b )
      {
        tempNodesToFaces[*b].insert(kf);
      }
    }
    std::cout<<"Check NodeToFace "<<std::endl;
    for( localIndex a=0 ; a<nodeManager.m_numNodes ; ++a )
    {
      std::cout<<"m_nodeToFaceMap["<<a<<"] = ( ";
      for( lSet::const_iterator iface=nodeManager.m_nodeToFaceMap[a].begin() ;
          iface!=nodeManager.m_nodeToFaceMap[a].end() ; ++iface )
      {
        if( tempNodesToFaces[a].count(*iface) == 0 )
          std::cout<<"*";

        std::cout<<*iface<<", ";
      }
      std::cout<<")    (";

      for( lSet::const_iterator iface=tempNodesToFaces[a].begin() ;
          iface!=tempNodesToFaces[a].end() ; ++iface )
      {
        if( nodeManager.m_nodeToFaceMap[a].count(*iface) == 0 )
          std::cout<<"*";

        std::cout<<*iface<<", ";
      }
      std::cout<<")"<<std::endl;
    }

  }



  if( m_verbose == 2 )
  {


    // nodeToElement
    Array1dT<std::set<std::pair< ElementRegionT*, localIndex > > > tempNodesToElems( nodeManager.m_numNodes );
    for( std::map< std::string, ElementRegionT >::iterator ielem=elementManager.m_ElementRegions.begin() ;
        ielem!=elementManager.m_ElementRegions.end() ; ++ielem )
    {
      ElementRegionT& elemRegion = ielem->second;
      for( localIndex k=0 ; k<elemRegion.m_numElems ; ++k )
      {
        std::pair< ElementRegionT*, localIndex > elem = std::make_pair(&elemRegion,k);

        for( localIndex a=0 ; a<elemRegion.m_toNodesRelation.Dimension(1) ; ++a )
        {
          tempNodesToElems[elemRegion.m_toNodesRelation(k,a)].insert(elem);
        }
      }
    }
    std::cout<<"Check NodeToElem "<<std::endl;
    for( localIndex a=0 ; a<nodeManager.m_numNodes ; ++a )
    {
      std::cout<<"m_NodeToElementMap["<<a<<"] = ( ";
      for( std::set<std::pair< ElementRegionT*, localIndex > >::const_iterator ielem=nodeManager.m_toElementsRelation[a].begin() ;
          ielem!=nodeManager.m_toElementsRelation[a].end() ; ++ielem )
      {
        if( tempNodesToElems[a].count(*ielem) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")    (";

      for( std::set<std::pair< ElementRegionT*, localIndex > >::const_iterator ielem=tempNodesToElems[a].begin() ;
          ielem!=tempNodesToElems[a].end() ; ++ielem )
      {
        if( nodeManager.m_toElementsRelation[a].count(*ielem) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")"<<std::endl;
    }


    // edgeToFace
    Array1dT<lSet> tempEdgeToFaces( edgeManager.DataLengths() );
    for( localIndex kf=0 ; kf<faceManager.m_numFaces ; ++kf )
    {
      for( lArray1d::const_iterator b=faceManager.m_toEdgesRelation[kf].begin() ;
          b!=faceManager.m_toEdgesRelation[kf].end() ; ++b )
      {
        tempEdgeToFaces[*b].insert(kf);
      }
    }
    std::cout<<"Check EdgeToFace "<<std::endl;
    for( localIndex ke=0 ; ke<edgeManager.DataLengths() ; ++ke )
    {
      std::cout<<"m_edgesToFaces["<<ke<<"] = ( ";
      for( lSet::const_iterator iface=edgeManager.m_toFacesRelation[ke].begin() ;
          iface!=edgeManager.m_toFacesRelation[ke].end() ; ++iface )
      {
        if( tempEdgeToFaces[ke].count(*iface) == 0 )
          std::cout<<"*";

        std::cout<<*iface<<", ";
      }
      std::cout<<")    (";

      for( lSet::const_iterator iface=tempEdgeToFaces[ke].begin() ;
          iface!=tempEdgeToFaces[ke].end() ; ++iface )
      {
        if( edgeManager.m_toFacesRelation[ke].count(*iface) == 0 )
          std::cout<<"*";

        std::cout<<*iface<<", ";
      }
      std::cout<<")"<<std::endl;
    }

    // faceToElement
    OneToOneRelation& parentFaceIndex = faceManager.GetOneToOneMap("parentIndex");
    Array1dT<std::set<std::pair< ElementRegionT*, localIndex > > > tempFacesToElems( faceManager.m_numFaces );
    for( std::map< std::string, ElementRegionT >::iterator ielem=elementManager.m_ElementRegions.begin() ;
        ielem!=elementManager.m_ElementRegions.end() ; ++ielem )
    {
      ElementRegionT& elemRegion = ielem->second;
      for( localIndex k=0 ; k<elemRegion.m_numElems ; ++k )
      {
        std::pair< ElementRegionT*, localIndex > elem = std::make_pair(&elemRegion,k);

        for( localIndex a=0 ; a<elemRegion.m_toFacesRelation.Dimension(1) ; ++a )
        {
          const localIndex faceID = elemRegion.m_toFacesRelation(k,a);
          tempFacesToElems[faceID].insert(elem);

          if( parentFaceIndex[faceID] != LOCALINDEX_MAX )
          {
            tempFacesToElems[parentFaceIndex[faceID]].insert(elem);
          }
        }
      }
    }
    std::cout<<"Check FacesToElem "<<std::endl;
    for( localIndex a=0 ; a<faceManager.m_numFaces ; ++a )
    {
      std::cout<<"m_FaceToElementMap["<<a<<"] = ( ";

      for( Array1dT<std::pair< ElementRegionT*, localIndex > >::const_iterator ielem=faceManager.m_toElementsRelation[a].begin() ;
          ielem!=faceManager.m_toElementsRelation[a].end() ; ++ielem )
      {
        if( tempFacesToElems[a].count(*ielem) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")    (";

      for( std::set<std::pair< ElementRegionT*, localIndex > >::const_iterator ielem=tempFacesToElems[a].begin() ;
          ielem!=tempFacesToElems[a].end() ; ++ielem )
      {

        if( faceManager.m_toElementsRelation[a].size() == 2 )
        {
          if( (faceManager.m_toElementsRelation[a][0] != *ielem) && (faceManager.m_toElementsRelation[a][1] != *ielem) )
            std::cout<<"*";
        }
        else if ( faceManager.m_toElementsRelation[a].size() )
        {
          if( (faceManager.m_toElementsRelation[a][0] != *ielem)  )
            std::cout<<"*";
        }
        else
        {
          std::cout<<"****";
        }


        std::cout<<ielem->second<<", ";
      }
      std::cout<<")"<<std::endl;
    }
  }

  CorrectSplitNodalMass(nodeManager, nodeID, nodeManager.m_childIndices[nodeID][0]);
}


void Fractunator3::WriteSiloDerived( SiloFile& siloFile,
                                     const int cycleNum,
                                     const realT problemTime,
                                     const bool isRestart )
{
  m_virtualFaces.WriteSiloMesh( siloFile, "virtual", m_virtualNodes, cycleNum, problemTime, isRestart);

  m_virtualNodes.WriteSilo( siloFile, "Fractunator/virtualNodes", "virtual", DB_NODECENT, cycleNum, problemTime, isRestart);
  m_virtualEdges.WriteSilo( siloFile, "Fractunator/virtualEdges", "virtual", DB_EDGECENT, cycleNum, problemTime, isRestart);
  m_virtualFaces.WriteSilo( siloFile, "Fractunator/virtualFaces", "virtual", DB_ZONECENT, cycleNum, problemTime, isRestart);
}

void Fractunator3::ReadSiloDerived( const SiloFile& siloFile,
                                    const int cycleNum,
                                    const realT problemTime,
                                    const bool isRestart )
{
  m_virtualNodes.ReadSilo( siloFile, "virtualNodes", "virtual", DB_NODECENT, cycleNum, problemTime, isRestart);
  m_virtualEdges.ReadSilo( siloFile, "virtualEdges", "virtual", DB_EDGECENT, cycleNum, problemTime, isRestart);
  m_virtualFaces.ReadSilo( siloFile, "virtualFaces", "virtual", DB_ZONECENT, cycleNum, problemTime, isRestart);

}

void Fractunator3::PreexistingFracture2D( NodeManagerT& nodeManager,
                                          EdgeManagerT& edgeManager,
                                          FaceManagerT& faceManager,
                                          ExternalFaceManagerT& externalFaceManager,
                                          ElementManagerT& elementManager,
                                          SpatialPartition& partition,
                                          const bool prefrac)
{
  // Do nothing.
}

int Fractunator3::SeparationDriver( NodeManagerT& nodeManager,
                                    EdgeManagerT& edgeManager,
                                    FaceManagerT& faceManager,
                                    ExternalFaceManagerT& externalFaceManager,
                                    ElementManagerT& elementManager,
                                    SpatialPartition& partition,
                                    const bool prefrac,
                                    const realT time)
{


  Array1dT<lSet> nodesToRupturedFaces;
  Array1dT<lSet> edgesToRupturedFaces;

  //  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator elementRegionIter = elementManager.m_ElementRegions.begin() ;
  //      elementRegionIter != elementManager.m_ElementRegions.end() ;
  //      ++elementRegionIter )
  //  {
  //    ElementRegionT& elementRegion = elementRegionIter->second;
  //
  //    elementRegion.CalculateVelocityGradients(nodeManager);
  //  }
  //  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");

  if (!prefrac)
  {

    if (m_failCriterion >0 )  // Stress intensity factor based criterion and mixed criterion.
    {
      if (m_failCriterion == 1)
      {
        faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, std::numeric_limits<realT>::max()); //this->m_failstress );
      }
      else
      {
        faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, m_failstress);
      }

      rArray1d& SIFonFace = faceManager.GetFieldData<realT>("SIFonFace");
      SIFonFace = std::numeric_limits<double>::min();

      if (m_failCriterion == 1)
      {
        iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
        for (localIndex a=0; a<faceManager.DataLengths(); ++a)
        {
          if (faceManager.m_toElementsRelation[a].size() == 2) ruptureState[a]=0;
        }
      }


      IdentifyRupturedFaces( nodeManager,
                             edgeManager,
                             faceManager,
                             elementManager,
                             partition,
                             prefrac );

    }

    else
    {
      UpdateRuptureStates( nodeManager,
                           edgeManager,
                           faceManager,
                           elementManager,
                           nodesToRupturedFaces,
                           edgesToRupturedFaces,
                           prefrac );

    }
  }
  else  // In the prefrac call, we need this to get the stressNOnFace, which will be used in the initialization of contacts for preexisting fractures.
  {
    for (localIndex kf = 0; kf < faceManager.DataLengths(); ++kf) faceManager.CalculateStressOnFace(elementManager, nodeManager, kf);
  }


  if (prefrac)
  {
    ModifiedObjectLists modifiedObjects;
    CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, prefrac);
  }

  // We do this here to get the nodesToRupturedFaces etc.
  // The fail stress check inside has been disabled
  PostUpdateRuptureStates( nodeManager,
                           edgeManager,
                           faceManager,
                           elementManager,
                           nodesToRupturedFaces,
                           edgesToRupturedFaces);

  int rval = 0;
  int rank ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //  Array1dT<MaterialBaseStateDataT*>&  temp = elementManager.m_ElementRegions["PM1"].m_materialStates;

  const iArray1d& isNodeGhost = nodeManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d& isSeparable = nodeManager.GetFieldData<int>("isSeparable");
  const iArray1d& layersFromDomainBoundary = nodeManager.GetFieldData<int>("LayersFromDomainBoundary");
  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  const iArray1d& isFaceGhost = faceManager.GetFieldData<FieldInfo::ghostRank>();



  // process nodes on the interior
  {
    ModifiedObjectLists modifiedObjects;
    for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
    {
      //      const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);
      if( layersFromDomainBoundary[a]>1 &&
          (isSeparable[a] || prefrac)&&
          isNodeGhost[a]<0 &&
          nodeManager.m_toElementsRelation[a].size()>1 &&
          CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac) > 0 ) //&&
        //          nodesToRupturedFaces[a].size()>0 )
      {
        rval += ProcessNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager, modifiedObjects, prefrac ) ;
      }
    }
    if (m_failCriterion == 1)
    {
      const lArray1d& primaryCandidateFace = faceManager.GetFieldData<localIndex>("primaryCandidateFace");


      for (localIndex a=0; a<faceManager.DataLengths(); ++a)
      {
        if (isFaceGhost[a]<0 && ruptureState[a] == -1 && ruptureState[primaryCandidateFace[a]] !=2 )
        {
          ruptureState[a] = 1;
        }
      }

      for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
      {

        //      const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);

        if( layersFromDomainBoundary[a]>1 &&
            (isSeparable[a] || prefrac) &&
            isNodeGhost[a]<0 &&
            nodeManager.m_toElementsRelation[a].size()>1 &&
            CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac) > 0 ) //&&
          //            nodesToRupturedFaces[a].size()>0 )
        {
          rval += ProcessNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager, modifiedObjects, prefrac ) ;
        }
      }

    }

    CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, false);

    MarkBirthTime(faceManager, modifiedObjects, time);
  }


  for( int color=0 ; color<partition.NumColor() ; ++color )
  {
    ModifiedObjectLists modifiedObjects;
    if( partition.Color() == color )
    {

      // process "near-boundary" nodes
      for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
      {

        //        const localIndex parentNodeIndex = nodeManager.GetParentIndex(a);

        if( layersFromDomainBoundary[a]<=1 &&
            (isSeparable[a] || prefrac) &&
            isNodeGhost[a]<0 &&
            nodeManager.m_toElementsRelation[a].size()>1 &&
            CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac) > 0 )// &&
          //           nodesToRupturedFaces[a].size()>0 )
        {
          rval += ProcessNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager, modifiedObjects, prefrac ) ;
        }
      }

      CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, false);

      MarkBirthTime(faceManager, modifiedObjects, time);

    }




    // TODO need to add to rval as a result of this communication
    partition.ModifyGhostsAndNeighborLists( modifiedObjects );

    // If a face is split by a domains that does not own this face, the rupture state for the virtual face will not be communicated to the owner.
    // The following is to fix this problem.
    iArray1d* vfaceRuptureState = m_virtualFaces.GetFieldDataPointer<int>( "ruptureState" );
    if (vfaceRuptureState != NULL)
    {
      const iArray1d& faceRuptureState = faceManager.GetFieldData<int>( "ruptureState" );
      for( localIndex kf=0 ; kf<(*vfaceRuptureState).size() ; ++kf )
      {
        if( faceRuptureState[kf]==2 )
        {
          (*vfaceRuptureState)[kf]=2;
        }
      }

    }


  }

  if (m_failCriterion == 1)
  {
    const lArray1d& primaryCandidateFace = faceManager.GetFieldData<localIndex>("primaryCandidateFace");

    {
      ModifiedObjectLists modifiedObjects;

      // Turn on rupture state for secondary fracture faces
      for (localIndex a=0; a<faceManager.DataLengths(); ++a)
      {
        if (isFaceGhost[a] < 0 && ruptureState[a] == -1 && ruptureState[primaryCandidateFace[a]] !=2 )
        {
          ruptureState[a] = 1;
          modifiedObjects.modifiedFaces.insert(a);
        }
      }
      partition.ModifyGhostsAndNeighborLists( modifiedObjects );
    }

    for( int color=0 ; color<partition.NumColor() ; ++color )
    {
      ModifiedObjectLists modifiedObjects;
      if( partition.Color() == color )
      {

        // process "near-boundary" nodes
        for( localIndex a=0 ; a<nodeManager.DataLengths() ; ++a )
        {

          if( layersFromDomainBoundary[a]<=1 &&
              (isSeparable[a] || prefrac) &&
              isNodeGhost[a]<0 &&
              nodeManager.m_toElementsRelation[a].size()>1 &&
              CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac) > 0 ) //&&
            //              nodesToRupturedFaces[a].size()>0 )
          {
            rval += ProcessNode( a, nodeManager, edgeManager, faceManager, externalFaceManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager, modifiedObjects, prefrac ) ;
          }
        }
        CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, false);

        MarkBirthTime(faceManager, modifiedObjects, time);
      }



      // TODO need to add to rval as a result of this communication
      partition.ModifyGhostsAndNeighborLists( modifiedObjects );

      // If a face is split by a domains that does not own this face, the rupture state for the virtual face will not be communicated to the owner.
      // The following is to fix this problem.
      iArray1d* vfaceRuptureState = m_virtualFaces.GetFieldDataPointer<int>( "ruptureState" );
      if (vfaceRuptureState != NULL)
      {
        const iArray1d& faceRuptureState = faceManager.GetFieldData<int>( "ruptureState" );
        for( localIndex kf=0 ; kf<(*vfaceRuptureState).size() ; ++kf )
        {
          if( faceRuptureState[kf]==2 )
          {
            (*vfaceRuptureState)[kf]=2;
          }
        }

      }
    }
  }




  /*
  for( std::map< std::string, ElementRegionT >::iterator i=elementManager.m_ElementRegions.begin() ;
       i != elementManager.m_ElementRegions.end() ; ++i )
  {
    i->second.CalculateNodalMasses( nodeManager ) ;
  }
   */

  return rval;
}



void Fractunator3::IdentifyRupturedFaces( NodeManagerT& nodeManager,
                                          EdgeManagerT& edgeManager,
                                          FaceManagerT& faceManager,
                                          ElementManagerT& elementManager,
                                          SpatialPartition& partition,
                                          const bool prefrac  )
{
  const iArray1d& isEdgeGhost = edgeManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d& layersEdgeFromBoundary = edgeManager.GetFieldData<int>("LayersFromDomainBoundary");
  lArray1d& primaryCandidateFace = faceManager.GetFieldData<localIndex>("primaryCandidateFace");
  primaryCandidateFace = std::numeric_limits<localIndex>::max();

  int rank ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //"Heal" faces that were marked but not split.
  if (m_failCriterion >0)
  {
    iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
    for (localIndex iFace = 0; iFace < faceManager.DataLengths(); ++iFace)
    {
      if (faceManager.m_isExternal[iFace] == 0)
      {
        ruptureState[iFace] = 0;
      }
    }
  }

  // We use the color map scheme because we can mark a face to be rupture ready from a partition where the face is a ghost.

  // Process interior edges
  {
    ModifiedObjectLists modifiedObjects;

    for (localIndex iEdge = 0; iEdge != edgeManager.DataLengths(); ++iEdge)
    {
      if(isEdgeGhost[iEdge] < 0 && layersEdgeFromBoundary[iEdge]>1 )
      {
        int edgeMode = CheckEdgeSplitability(iEdge,
                                             nodeManager,
                                             faceManager,
                                             edgeManager,
                                             prefrac);
        if (edgeMode == 0 || edgeMode == 1) // We need to calculate SIF
        {
          R1Tensor vecTipNorm, vecTip;
          realT SIF = CalculateEdgeSIF(iEdge,
                                       nodeManager,
                                       edgeManager,
                                       faceManager,
                                       elementManager,
                                       vecTipNorm, vecTip);

          if ( SIF > MinimumToughnessOnEdge(iEdge, edgeManager, faceManager) * m_thresholdForEvaluateFace )
          {
            MarkRuptureFaceFromEdge(iEdge,
                                    nodeManager,
                                    edgeManager,
                                    faceManager,
                                    vecTipNorm,
                                    vecTip,
                                    modifiedObjects,
                                    edgeMode);
          }
        }
      }
    }
  }

  // Process near boundary edges
  {
    for( int color=0 ; color<partition.NumColor() ; ++color )
    {
      ModifiedObjectLists modifiedObjects;
      if( partition.Color() == color )
      {
        for (localIndex iEdge = 0; iEdge != edgeManager.DataLengths(); ++iEdge)
        {

          if(isEdgeGhost[iEdge] < 0 && layersEdgeFromBoundary[iEdge]<=1 )
          {
            int edgeMode = CheckEdgeSplitability(iEdge,
                                                 nodeManager,
                                                 faceManager,
                                                 edgeManager,
                                                 prefrac);
            if (edgeMode == 0 || edgeMode == 1) // We need to calculate SIF
            {
              R1Tensor vecTipNorm, vecTip, vecEdge;
              realT SIF = CalculateEdgeSIF(iEdge,
                                           nodeManager,
                                           edgeManager,
                                           faceManager,
                                           elementManager,
                                           vecTipNorm, vecTip);

              if (SIF >  MinimumToughnessOnEdge(iEdge, edgeManager, faceManager) * m_thresholdForEvaluateFace) // && edgeMode == 1)
              {
                MarkRuptureFaceFromEdge(iEdge,
                                        nodeManager,
                                        edgeManager,
                                        faceManager,
                                        vecTipNorm,
                                        vecTip,
                                        modifiedObjects,
                                        edgeMode);
              }
            }
          }
        }
        partition.ModifyGhostsAndNeighborLists( modifiedObjects );
      }
    }
  }

}

realT Fractunator3::CalculateEdgeSIF( const localIndex edgeID,
                                      NodeManagerT& nodeManager,
                                      EdgeManagerT& edgeManager,
                                      FaceManagerT& faceManager,
                                      ElementManagerT& elementManager,
                                      R1Tensor& vecTipNorm, R1Tensor& vecTip)
{
  realT rval;
  localIndex nExternalFaces = 0;
  lArray1d faceInvolved;
  rArray1d& SIF_I = edgeManager.GetFieldData<realT>("SIF_I");
  rArray1d& SIF_II = edgeManager.GetFieldData<realT>("SIF_II");
  rArray1d& SIF_III = edgeManager.GetFieldData<realT>("SIF_III");

  SIF_I[edgeID] = 0.0;
  SIF_II[edgeID] = 0.0;
  SIF_III[edgeID] = 0.0;

  for( lSet::const_iterator iface=edgeManager.m_toFacesRelation[edgeID].begin() ;
      iface!=edgeManager.m_toFacesRelation[edgeID].end() ; ++iface )
  {
    if (faceManager.m_isExternal[*iface] == 1)
    {
      nExternalFaces++;
      faceInvolved.push_back(*iface);
    }
  }
  if (nExternalFaces > 2)
  {
    throw GPException("Error! This is a singular edge, not a tip.  This should not happen!");
  }

  localIndex faceA, faceAp;
  if ( (faceManager.m_parentIndex[faceInvolved[0]] == LOCALINDEX_MAX && faceManager.m_parentIndex[faceInvolved[1]] == faceInvolved[0]) ||
      (faceManager.m_parentIndex[faceInvolved[1]] == LOCALINDEX_MAX && faceManager.m_parentIndex[faceInvolved[0]] == faceInvolved[1]) )
  {
    faceA = faceInvolved[0];
    faceAp = faceInvolved[1];
  }
  else
  {
    char msg[200];
    sprintf(msg, "Error! Edge %d has two external faces, but the parent-child relationship is wrong.", int(edgeID));
    throw GPException(msg);
  }


  // We define three unit vectors
  // vecEdge: pointing from node 0 to node 1 along the tip edge
  // vecTip: pointing from the opening into the solid
  // vecTipNorm: normal of the one of the fracture faces;  vecTip X vecTipNorm should point to the direction of vecEdge
  // These new definitions are consistent with those in 2D.  See the illustration in Fractunator2D.

  vecTipNorm = faceManager.FaceNormal(nodeManager, faceA);
  vecTipNorm -= faceManager.FaceNormal(nodeManager, faceAp);
  vecTipNorm.Normalize();

  R1Tensor vecEdge;
  edgeManager.EdgeVector(nodeManager, edgeID, vecEdge);
  vecEdge.Normalize();

  vecTip.Cross(vecTipNorm, vecEdge);
  vecTip.Normalize();
  R1Tensor v0, v1;
  edgeManager.EdgeCenter(nodeManager, edgeID, v0);
  faceManager.FaceCenter(nodeManager, faceA, v1);
  v0 -= v1;

  if (Dot(v0,vecTip) < 0) vecTip *= -1.0;
  if (Dot(Cross(vecTip, vecTipNorm), vecEdge) < 0)
  {
    vecTipNorm *= -1;
    faceA = faceInvolved[1];
    faceAp = faceInvolved[0];
  }


  //Now we need to figure out if a special situation applies to this edge
  // where the fracture face is a quad and three of the nodes are still pinched
  // We use a different algorithm

  bool threeNodesPinched(false);
  lArray1d openNodeID;

  if (faceManager.m_toNodesRelation[faceA].size() == 4)  // Only quads have this problem
  {
    int numSharedNodes = 2;
    lArray1d lNodeFaceA, lNodeFaceAp;

    lNodeFaceA = faceManager.m_toNodesRelation[faceA];
    lNodeFaceAp = faceManager.m_toNodesRelation[faceAp];


    //We remove all the shared nodes and the one remains should be the open one.
    lNodeFaceAp.erase(std::find(lNodeFaceAp.begin(),lNodeFaceAp.end(), edgeManager.m_toNodesRelation[edgeID][0]));
    lNodeFaceAp.erase(std::find(lNodeFaceAp.begin(),lNodeFaceAp.end(), edgeManager.m_toNodesRelation[edgeID][1]));
    lNodeFaceA.erase(std::find(lNodeFaceA.begin(),lNodeFaceA.end(), edgeManager.m_toNodesRelation[edgeID][0]));
    lNodeFaceA.erase(std::find(lNodeFaceA.begin(),lNodeFaceA.end(), edgeManager.m_toNodesRelation[edgeID][1]));

    for( lArray1d::iterator j = faceManager.m_toNodesRelation[faceA].begin() ;
        j!=faceManager.m_toNodesRelation[faceA].end() ; ++j )
    {
      localIndex iNd = *j;
       if (iNd != edgeManager.m_toNodesRelation[edgeID][0] && iNd != edgeManager.m_toNodesRelation[edgeID][1])
       {
         if (std::find(faceManager.m_toNodesRelation[faceAp].begin(), faceManager.m_toNodesRelation[faceAp].end(), iNd) != faceManager.m_toNodesRelation[faceAp].end())
         {
           numSharedNodes++;
           lNodeFaceA.erase(std::find(lNodeFaceA.begin(),lNodeFaceA.end(), iNd));
           lNodeFaceAp.erase(std::find(lNodeFaceAp.begin(),lNodeFaceAp.end(), iNd));
         }
       }
    }

    if (numSharedNodes == 4)
    {
      throw GPException("Error.  The fracture face has four shared nodes with its child.  This should not happen.");
    }
    else if (numSharedNodes == 3)
    {
      threeNodesPinched = true;
      if (lNodeFaceA.size() != 1 || lNodeFaceAp.size() != 1)
      {
        throw GPException("Error. These two faces share three nodes but the number of remaining nodes is not one.  Something is wrong");
      }
      else
      {
        openNodeID.push_back(lNodeFaceA[0]);
        openNodeID.push_back(lNodeFaceAp[0]);
      }
    }
  }

  localIndex convexCorner(std::numeric_limits<localIndex>::max());
  if (threeNodesPinched) // Now we need to identify which node on the edge is the convex point and which one is the concave corner.  The convex node must share an edge with the open node.
  {
    localIndex iNd, jNd;
    iNd = edgeManager.m_toNodesRelation[edgeID][0];
    jNd = edgeManager.m_toNodesRelation[edgeID][1];
    for( lArray1d::iterator j = faceManager.m_toEdgesRelation[faceA].begin() ;
        j!=faceManager.m_toEdgesRelation[faceA].end() ; ++j )
    {
      localIndex edge = *j;
      if ((openNodeID[0] == edgeManager.m_toNodesRelation[edge][0] && iNd == edgeManager.m_toNodesRelation[edge][1]) ||
          (openNodeID[0] == edgeManager.m_toNodesRelation[edge][1] && iNd == edgeManager.m_toNodesRelation[edge][0])
         )
      {
        convexCorner = iNd;
        break;
      }
      if ((openNodeID[0] == edgeManager.m_toNodesRelation[edge][0] && jNd == edgeManager.m_toNodesRelation[edge][1]) ||
          (openNodeID[0] == edgeManager.m_toNodesRelation[edge][1] && jNd == edgeManager.m_toNodesRelation[edge][0])
         )
      {
        convexCorner = jNd;
        break;
      }
    }

    if (convexCorner == std::numeric_limits<localIndex>::max())
      throw GPException("Error.  This is a three-node-pinched edge but I cannot find the convex corner");

  }



  // Calculate element forces acting on this edge.  Need to add nodal forces from two nodes up.
  //An element has to be within the range of this edge to be included.
  //For the threeNodesPinched case, we only use the force on the node at the convex point, not the concave point.  The force at the former is ususally greater, so we just pick the great one instead of doing a geometrical check.

  localIndex nElemEachSide[2], nGhostElem;
  R1Tensor fNodeO = static_cast<R1Tensor>(0.0);
  nElemEachSide[0] = 0;
  nElemEachSide[1] = 0;
  R1Tensor xEdge;
  edgeManager.EdgeCenter(nodeManager, edgeID, xEdge);


  if (!threeNodesPinched)
  {
    for( unsigned int a=0 ; a<edgeManager.m_toNodesRelation.Dimension(1) ; ++a ) // Loop through the two nodes
    {
      localIndex nodeID = edgeManager.m_toNodesRelation(edgeID,a);

      for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodeID].begin() ;
          k!=nodeManager.m_toElementsRelation[nodeID].end() ; ++k )
      {
        ElementRegionT* elementRegion = k->first;
        localIndex iEle = k->second;
        R1Tensor fN, xEle;


        xEle = elementRegion->GetElementCenter(iEle, nodeManager);

        realT ndist, udist, segmentLength;
        R1Tensor ptPrj;
        GeometryUtilities::ProjectPointToLineSegment( (*nodeManager.m_refposition)[nodeID],
                                                      (*nodeManager.m_refposition)[edgeManager.m_toNodesRelation(edgeID,1-a)],
                                                      xEle,
                                                      ndist, udist, segmentLength,
                                                      ptPrj);
        if (udist <= edgeManager.EdgeLength(nodeManager, edgeID) && udist > 0.0)
        {
          elementRegion->CalculateNodalForcesFromOneElement( nodeID, iEle, nodeManager, fN);

          xEle -= xEdge;

          if (Dot(xEle, vecTipNorm) > 0)
          {
            nElemEachSide[0] += 1;
            fNodeO += fN;
          }
          else
          {
            nElemEachSide[1] +=1;
            fNodeO -= fN;
          }
        }
      } // Loop through all elements connected to this node

    }  // loop over two nodes on the edge
  }
  else
  {

    localIndex nodeID = convexCorner;
    fNodeO = 0.0;

    for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodeID].begin() ;
        k!=nodeManager.m_toElementsRelation[nodeID].end() ; ++k )
    {
      ElementRegionT* elementRegion = k->first;
      localIndex iEle = k->second;
      R1Tensor fN, xEle;

      xEle = elementRegion->GetElementCenter(iEle, nodeManager);

      {
        elementRegion->CalculateNodalForcesFromOneElement( nodeID, iEle, nodeManager, fN);

        xEle -= xEdge;

        if (Dot(xEle, vecTipNorm) > 0)
        {
          nElemEachSide[0] += 1;
          fNodeO += fN;
        }
        else
        {
          nElemEachSide[1] +=1;
          fNodeO -= fN;
        }
      }
    } // Loop through all elements connected to this node
  }

  if (nElemEachSide[0]>=1 && nElemEachSide[1]>=1) fNodeO /= 2.0;
    //We have contributions from both sides.  The two sizes are the two sides of the fracture plane.  If the fracture face is on domain boundary, it's possible to have just one side.

  localIndex tipFaces[2];
  tipFaces[0] = faceA;
  tipFaces[1] = faceAp;

  // We have to subtract the nodal force at other nodes on these two open faces to take into account the effects of surface traction along the fracture.
  R1Tensor fFaceA[2];
  // fFaceA is actually the average nodal force on each face.  Assuming homogeneous meshing.

  for (localIndex i=0; i<2; ++i)
  {
    localIndex faceID = tipFaces[i];
    nGhostElem = 0;
    fFaceA[i] = 0.0;

    for( lArray1d::iterator j = faceManager.m_toNodesRelation[faceID].begin() ;
        j!=faceManager.m_toNodesRelation[faceID].end() ; ++j )
    {
      localIndex iNd = *j;
      if (iNd != edgeManager.m_toNodesRelation(edgeID,0) && iNd != edgeManager.m_toNodesRelation(edgeID,1) )
      {
        nGhostElem = 0;
        R1Tensor fNode = static_cast<R1Tensor>(0.0);

        for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[iNd].begin() ;
            k!=nodeManager.m_toElementsRelation[iNd].end() ; ++k )
        {
          ElementRegionT* elementRegion = k->first;
          localIndex iEle = k->second;
          iArray1d& elem_is_ghost = elementRegion->GetFieldData<FieldInfo::ghostRank>();
          R1Tensor fN;

          elementRegion->CalculateNodalForcesFromOneElement( iNd, iEle, nodeManager, fN);

          fNode += fN;

          if (elem_is_ghost[iEle] >= 0) nGhostElem +=1;
        }

        if (nGhostElem == nodeManager.m_toElementsRelation.size())
        {
          // All elements connected to this node are ghost elements.  This implies that half of the elements are in the next partition.
          fNode *= 2.0;
        }

        fNode /= ( faceManager.m_toNodesRelation[faceID].size() - 2 );
        fFaceA[i] += fNode;
      }

    }

  }

  R1Tensor tipForce;
  tipForce[0] = Dot(fNodeO, vecTipNorm) + Dot(fFaceA[0], vecTipNorm) / 2.0 - Dot(fFaceA[1], vecTipNorm) /2.0;
  tipForce[1] = Dot(fNodeO, vecTip) + Dot(fFaceA[0], vecTip) / 2.0 - Dot(fFaceA[1], vecTip) /2.0;
  tipForce[2] = Dot(fNodeO, vecEdge) + Dot(fFaceA[0], vecEdge) / 2.0 - Dot(fFaceA[1], vecEdge) /2.0;

  R1Tensor tipDisplacement, tipOpening, tipFaceDisplacement[2];

  if( !threeNodesPinched )
  {
    for (localIndex i=0; i<2; ++i)
    {
      localIndex faceID = tipFaces[i];
      tipFaceDisplacement[i] = 0.0;

      for( lArray1d::iterator j = faceManager.m_toNodesRelation[faceID].begin() ;
          j!=faceManager.m_toNodesRelation[faceID].end() ; ++j )
      {
        localIndex iNd = *j;
        if (iNd != edgeManager.m_toNodesRelation(edgeID,0) && iNd != edgeManager.m_toNodesRelation(edgeID,1) )
        {
          tipFaceDisplacement[i] += (*nodeManager.m_displacement)[iNd];
        }
      }

      tipFaceDisplacement[i] /= (faceManager.m_toNodesRelation[faceID].size() - 2);
    }
    tipDisplacement = tipFaceDisplacement[1];
    tipDisplacement -= tipFaceDisplacement[0];
  }
  else
  {
    tipDisplacement = (*nodeManager.m_displacement)[openNodeID[1]];
    tipDisplacement -= (*nodeManager.m_displacement)[openNodeID[0]];


// This is Randy's old algorithm to handle the three-pinched situation.  I am keeping it because it's quite clever in how it matched the node pairs across the two fracture faces.
//    const localIndex numNodes = faceManager.m_toNodesRelation[tipFaces[0]].size();
//
//    for( localIndex a=0 ; a<numNodes ; ++a )
//    {
//      const localIndex aa = a == 0 ? a : numNodes - a;
//
//      const localIndex nodeID0 = faceManager.m_toNodesRelation[tipFaces[0]][a];
//      const localIndex nodeID1 = faceManager.m_toNodesRelation[tipFaces[1]][aa];
//
//      if ( nodeID0 != edgeManager.m_toNodesRelation(edgeID,0) && nodeID0 != edgeManager.m_toNodesRelation(edgeID,1) &&
//           nodeID1 != edgeManager.m_toNodesRelation(edgeID,0) && nodeID1 != edgeManager.m_toNodesRelation(edgeID,1) )
//      {
//        tipDisplacement = (*nodeManager.m_displacement)[nodeID1];
//        tipDisplacement -= (*nodeManager.m_displacement)[nodeID0];
//
//        if( tipDisplacement.L2_Norm() > maxTipDisplacement.L2_Norm() )
//        {
//          maxTipDisplacement = tipDisplacement;
//        }
//      }
//    }
  }


  tipOpening[0] = Dot(tipDisplacement, vecTipNorm);
  tipOpening[1] = Dot(tipDisplacement, vecTip);
  tipOpening[2] = Dot(tipDisplacement, vecEdge);

  //if (Dot(Cross(vecTip, vecTipNorm),vecEdge) < 0.0) tipOpening[1] *= -1.0;



  realT tipArea;
  tipArea = faceManager.SurfaceArea(nodeManager, faceA, true);
  if (faceManager.m_toNodesRelation[faceA].size() == 3)
  {
    tipArea *= 2.0;
  }

  SIF_I[edgeID] = pow(fabs(tipForce[0] * tipOpening[0] / 2.0 / tipArea), 0.5);
  SIF_II[edgeID] = pow(fabs(tipForce[1] * tipOpening[1] / 2.0 / tipArea), 0.5);
  SIF_III[edgeID] = pow(fabs(tipForce[2] * tipOpening[2] / 2.0 / tipArea), 0.5);

  if (tipForce[1] < 0.0) SIF_II[edgeID] *= -1.0;

  if (tipOpening[0] < 0) SIF_I(edgeID) *= -1.0;

  if (SIF_I[edgeID] > 0.0)
  {
    rval = pow(SIF_I[edgeID]*SIF_I[edgeID]+SIF_II[edgeID]*SIF_II[edgeID]+SIF_III[edgeID]*SIF_III[edgeID], 0.5);
  }
  else
  {
    rval = -1.0;
  }

  return rval;
}

void Fractunator3::MarkRuptureFaceFromEdge ( const localIndex edgeID,
                                             NodeManagerT& nodeManager,
                                             EdgeManagerT& edgeManager,
                                             FaceManagerT& faceManager,
                                             R1Tensor& vecTipNorm,
                                             R1Tensor& vecTip,
                                             ModifiedObjectLists& modifiedObjects,
                                             const int edgeMode)
{

  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  rArray1d& stressNOnFace = faceManager.GetFieldData<realT>("stressNOnFace");
  rArray1d& SIFonFace = faceManager.GetFieldData<realT>("SIFonFace");
  rArray1d& SIF_I = edgeManager.GetFieldData<realT>("SIF_I");
  rArray1d& SIF_II = edgeManager.GetFieldData<realT>("SIF_II");
  lArray1d& primaryCandidateFace = faceManager.GetFieldData<localIndex>("primaryCandidateFace");
  iArray1d& isFaceSeparable = faceManager.GetFieldData<int>("isSeparable");
  rArray1d& faceToughness = faceManager.GetFieldData<realT>("faceToughness");


  iArray1d eligibleFaces;
  realT lowestSIF = std::numeric_limits<realT>::max();
  realT highestSIF = std::numeric_limits<realT>::min();
  realT lowestStress = std::numeric_limits<realT>::max();
  realT lowestScore = std::numeric_limits<realT>::max();
  realT highestScore = std::numeric_limits<realT>::min();
  realT secondScore = std::numeric_limits<realT>::min();
  localIndex faceWithHighestScore = std::numeric_limits<localIndex>::max();
  localIndex faceWithSecondScore = std::numeric_limits<localIndex>::max();

  R1Tensor vecEdge;
  edgeManager.EdgeVector(nodeManager, edgeID, vecEdge);
  vecEdge.Normalize();

  for( lSet::const_iterator a=edgeManager.m_toFacesRelation[edgeID].begin() ;
      a!=edgeManager.m_toFacesRelation[edgeID].end() ; ++a )
  {
    localIndex iface = *a;

    if (faceManager.m_toElementsRelation[iface].size() == 2  &&
        CheckOrphanElement(faceManager, iface) == 0)
    {
      R1Tensor fc, fn, vecFace;
      faceManager.FaceCenterAndNormal( nodeManager, iface, fc, fn );
      faceManager.InFaceVectorNormalToEdge(nodeManager,
                                           edgeManager,
                                           iface, edgeID,
                                           vecFace);
      if ( Dot(vecTip,vecFace) > cos(m_maxTurnAngle))
      {
        eligibleFaces.push_back(iface);
        realT thetaFace = acos(Dot(vecTip,vecFace)*0.999999);  // We multiply this by 0.9999999 to avoid an exception caused by acos a number slightly larger than 1.

        if (Dot(Cross(vecTip,vecFace),vecEdge) < 0.0)
        {
          thetaFace *= -1.0;
        }

        SIFonFace[iface] = cos(thetaFace / 2.0) * ( SIF_I[edgeID] * cos(thetaFace / 2.0) * cos(thetaFace / 2.0) - 1.5 * SIF_II[edgeID] * sin(thetaFace) );

        highestSIF = std::max(highestSIF, SIFonFace[iface]/faceToughness[iface]);
        lowestSIF = std::min(lowestSIF, SIFonFace[iface]/faceToughness[iface]);
        lowestStress = std::min(lowestStress, stressNOnFace[iface]);

      }
    }
  }


  iArray1d pickedFaces;
  if (eligibleFaces.size() >=1)
  {
    realT lengthscale = edgeManager.EdgeLength(nodeManager, edgeID);

    for (localIndex i = 0; i < eligibleFaces.size(); ++i )
    {
      localIndex iface = eligibleFaces[i];
      realT splitabilityScore = SIFonFace[iface] - lowestSIF * faceToughness[iface] + (stressNOnFace[iface] - lowestStress) * sqrt(lengthscale);
      lowestScore = std::min(lowestScore, splitabilityScore);

      if (faceWithHighestScore == std::numeric_limits<localIndex>::max())
      {
        faceWithHighestScore = iface;
        highestScore = splitabilityScore;
      }
      else if (splitabilityScore > highestScore)
      {
        faceWithSecondScore = faceWithHighestScore;
        secondScore = highestScore;
        faceWithHighestScore = iface;
        highestScore = splitabilityScore;
      }
      else if (splitabilityScore > secondScore)
      {
        faceWithSecondScore = iface;
        secondScore = splitabilityScore;
      }
    }

    pickedFaces.push_back(faceWithHighestScore);

    if ( eligibleFaces.size() >= 3 && (highestScore - secondScore) < 0.1 * (highestScore - lowestScore))
    {
      pickedFaces.push_back(faceWithSecondScore);
    }

  }

  for (localIndex i = 0; i < pickedFaces.size(); ++i)
  {
    localIndex pickedFace = pickedFaces[i];

    if (highestSIF > 1.0 && edgeMode == 1 && i == 0 && isFaceSeparable[pickedFace] == 1)
    {
      ruptureState[pickedFace] = 1;
      modifiedObjects.modifiedFaces.insert(pickedFace);
    }
    else if (highestSIF > 1.0 && edgeMode == 1 && i == 1 && isFaceSeparable[pickedFace] == 1)
    {
      ruptureState[pickedFace] = -1;
      modifiedObjects.modifiedFaces.insert(pickedFace);
      primaryCandidateFace[pickedFace] = faceWithHighestScore;
    }


    // We didn't really need to do this unless the criterion above has been satisfied.
    // We are calculating this regardless the criterion for debugging purpose.
    if (m_markExtendedLayer == 1 && highestSIF > 1.0 && edgeMode == 1 )
    {
      // Next we mark the faces that are 1) connected to this face, and 2) attached to one node of the edge (implicitly satisfied), and 3) almost co-plane with this face
      for( lArray1d::iterator j = faceManager.m_toEdgesRelation[pickedFace].begin() ;
          j!=faceManager.m_toEdgesRelation[pickedFace].end() ; ++j )
      {
        localIndex iedge = *j;
        if (iedge != edgeID)
        {
          for( lSet::const_iterator a=edgeManager.m_toFacesRelation[iedge].begin() ;
              a!=edgeManager.m_toFacesRelation[iedge].end() ; ++a )
          {
            localIndex iface = *a;
            if (iface != pickedFace && isFaceSeparable[iface] == 1 &&
                ( faceManager.IsNodeOnFace(iface, edgeManager.m_toNodesRelation[edgeID][0]) ||
                    faceManager.IsNodeOnFace(iface, edgeManager.m_toNodesRelation[edgeID][1])))
            {
              R1Tensor fc, fn, vecFace, fn0, fc0, ptPrj;
              realT nDist, uDist, segmentLength;
              faceManager.FaceCenterAndNormal( nodeManager, iface, fc, fn );
              faceManager.FaceCenterAndNormal( nodeManager, pickedFace, fc0, fn0 );
              faceManager.InFaceVectorNormalToEdge(nodeManager,
                                                   edgeManager,
                                                   iface, edgeID,
                                                   vecFace);
              GeometryUtilities::ProjectPointToLineSegment( (*nodeManager.m_refposition)[edgeManager.m_toNodesRelation(edgeID,0)],
                                                            (*nodeManager.m_refposition)[edgeManager.m_toNodesRelation(edgeID,1)],
                                                            fc,
                                                            nDist, uDist, segmentLength, ptPrj);

              // thetaFace does not strictly speaking apply to this face since the tip edge is not a edge of this face.
              // We calculate it as if this face is coplane with the master face
              realT thetaFace = acos(Dot(vecTip,vecFace)*0.999999);
              if (Dot(Cross(vecTip,vecFace),vecEdge) < 0.0)
              {
                thetaFace *= -1.0;
              }

              if ( Dot(vecTip,vecFace) > cos(m_maxTurnAngle) &&
                  uDist / segmentLength > -m_faceToEdgeProjectionTol &&
                  uDist / segmentLength < 1 + m_faceToEdgeProjectionTol &&
                  fabs(Dot(vecEdge, fn)) < 0.3 &&  // this face is kind of parallel to the tip edge.
                  fabs(Dot(fn0, fn)) > 0.95)  // co-plane
              {
                // Calculate SIFonFace is not really necessary but we keep it here for now for debugging.
                // Marking of the extended layer is purely based on geometrical and topological criteria.
                SIFonFace[iface] = std::max(SIFonFace[iface],
                                            SIF_I[edgeID] * cos(thetaFace / 2.0) * cos(thetaFace / 2.0) - 1.5 * SIF_II[edgeID] * sin(thetaFace) );

                if (highestSIF > 1.0 && edgeMode == 1)
                {
                  ruptureState[iface] = ruptureState[pickedFace];
                  modifiedObjects.modifiedFaces.insert(iface);
                  primaryCandidateFace[iface] = primaryCandidateFace[pickedFace];
                }
              }
            }
          }
        }
      }
    }
  }
}

realT Fractunator3::CalculateKinkAngle (const localIndex edgeID,
                                        const NodeManagerT& nodeManager,
                                        EdgeManagerT& edgeManager,
                                        FaceManagerT& faceManager)
{
  lArray1d faces;
  realT kinkAngle;

  for( lSet::const_iterator iface=edgeManager.m_toFacesRelation[edgeID].begin() ;
      iface!=edgeManager.m_toFacesRelation[edgeID].end() ; ++iface )
  {
    if (faceManager.m_isExternal[*iface] == 1) faces.push_back(*iface);
  }

  if (faces.size() != 2)
  {
    return(-1.0);
  }
  else
  {
    // First check if the two faces are parent-child pairs
    if (faceManager.m_parentIndex[faces[0]]==faces[1] || faceManager.m_parentIndex[faces[1]]==faces[0] )
    {
      return(0.0);
    }

    R1Tensor vecFace[3];
    faceManager.InFaceVectorNormalToEdge(nodeManager, edgeManager, faces[0], edgeID, vecFace[0]);
    faceManager.InFaceVectorNormalToEdge(nodeManager, edgeManager, faces[1], edgeID, vecFace[1]);
    vecFace[2] = vecFace[0];
    vecFace[2] += vecFace[1];
    vecFace[2] /= 2.0;

    kinkAngle = acos(Dot(vecFace[0],vecFace[1])*0.999999) / 3.141592653589793238462 * 180.0;

    R1Tensor vecFaceNorm;
    vecFaceNorm = faceManager.FaceNormal(nodeManager, faces[0]);
    vecFaceNorm  += faceManager.FaceNormal(nodeManager, faces[1]);
    vecFaceNorm /= 2.0;

    if (Dot(vecFace[2], vecFaceNorm) < 0.0)
      kinkAngle = 360.0 - kinkAngle;

    return(kinkAngle);

  }
}

void Fractunator3::CalculateKinkAngles (FaceManagerT& faceManager,
                                        EdgeManagerT& edgeManager,
                                        NodeManagerT& nodeManager,
                                        ModifiedObjectLists& modifiedObjects,
                                        const bool prefrac)
{
  rArray1d& kinkAngle = edgeManager.GetFieldData<realT>("kinkAngle");

  if (prefrac)
  {
    for (localIndex edgeID = 0; edgeID < edgeManager.DataLengths(); ++edgeID)
    {
      kinkAngle[edgeID] = CalculateKinkAngle(edgeID, nodeManager, edgeManager, faceManager);
    }
  }
  else
  {
    for( lSet::const_iterator i=modifiedObjects.newEdges.begin() ; i!=modifiedObjects.newEdges.end() ; ++i )
    {
      kinkAngle[*i] = CalculateKinkAngle(*i, nodeManager, edgeManager, faceManager);
    }
    for( lSet::const_iterator i=modifiedObjects.modifiedEdges.begin() ; i!=modifiedObjects.modifiedEdges.end() ; ++i )
    {
      kinkAngle[*i] = CalculateKinkAngle(*i, nodeManager, edgeManager, faceManager);
    }
  }
}

/// Register solver in the solver factory
REGISTER_FRACTUNATOR( Fractunator3)
