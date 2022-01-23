/**
 * @file Fractunator2.cpp
 * @author settgast1
 * @date Jul 14, 2011
 */

#include "Fractunator2.h"
#include <limits.h>

Fractunator2::Fractunator2():
FractunatorBase()
{
  // TODO Auto-generated constructor stub

}

Fractunator2::~Fractunator2()
{
  // TODO Auto-generated destructor stub
}







bool Fractunator2::FindFracturePlanes( const localIndex nodeID,
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

/*
  std::cout<<"checking node "<<nodeID<<std::endl;

  if( nodeID==33 )
    std::cout<<""<<std::endl;
*/
  const lSet& nodeToRupturedFaces = nodesToRupturedFaces[nodeID];
  const lSet& nodeToFaces = nodeManager.GetUnorderedVariableOneToManyMap( "nodeToFaceMap" )[nodeID];
  const lSet& nodeToEdges = nodeManager.GetUnorderedVariableOneToManyMap( "nodeToEdgeMap" )[nodeID];


  const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements = nodeManager.m_toElementsRelation[nodeID] ;


  const iArray1d& isEdgeExternal = edgeManager.m_isExternal;

  const Array1dT<lArray1d>& faceToEdges = faceManager.m_toEdgesRelation;


  // **** local working arrays ****

  // array to hold the faces ready for rupture
  std::set<localIndex> nodeToRuptureReadyFaces;
  for( lSet::const_iterator i=nodeToRupturedFaces.begin() ;
       i!=nodeToRupturedFaces.end() ; ++i )
  {
      nodeToRuptureReadyFaces.insert(*i);
  }



  if( nodeToRuptureReadyFaces.size() == 0 )
    return false;






  //
  std::map< localIndex, std::set<localIndex> > edgesToRuptureReadyFaces;
  for( lSet::const_iterator edgeIndex=nodeToEdges.begin() ; edgeIndex!=nodeToEdges.end() ; ++edgeIndex )
  {
    if( !(edgesToRupturedFaces[*edgeIndex].empty()) )
      edgesToRuptureReadyFaces[*edgeIndex].insert( edgesToRupturedFaces[*edgeIndex].begin(), edgesToRupturedFaces[*edgeIndex].end() );
  }

  // need a map from faces to edges attached to the node
  std::map< localIndex, std::pair<localIndex,localIndex> > localFacesToEdges;


  for( lSet::const_iterator kf=nodeToFaces.begin() ; kf!=nodeToFaces.end() ; ++kf )
  {
    localIndex edge[2] = { INT_MAX,INT_MAX };
    int count = 0;
    for( lArray1d::const_iterator ke=faceToEdges[*kf].begin() ; ke!=faceToEdges[*kf].end() ; ++ke )
    {
      if( edgeManager.hasNode( *ke, nodeID ) )
      {
        edge[count++] = *ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
//      throw GPException("Fractunator2::FindFracturePlanes: invalid edge.");
    }


    localFacesToEdges[*kf] = std::make_pair(edge[0],edge[1]);

    if( m_verbose ==2 )
      std::cout<<"localFacesToEdges["<<*kf<<"] = ( "<<localFacesToEdges[*kf].first<<", "<<localFacesToEdges[*kf].second<<" )"<<std::endl;
  }


  // ****
  // if the edge is not external, and the size of edgesToRupturedFaces is less than 2, then the edge is a dead-end
  // as far as a rupture plane is concerned. The face associated with the edge should be removed from the working
  // list of ruptured faces.

  // loop over all the edges
  for( lSet::const_iterator edgeIndex=nodeToEdges.begin() ; edgeIndex!=nodeToEdges.end() ; ++edgeIndex )
  {

    localIndex thisEdge = *edgeIndex;

//    std::cout<<thisEdge<<std::endl;
    // if the edge is internal and the edge is only attached to one ruptured face
    while( isEdgeExternal[thisEdge]!=1 && edgesToRuptureReadyFaces[thisEdge].size()==1 )
    {
      // the index for the face that is a "dead end"
      localIndex deadEndFace = *(edgesToRuptureReadyFaces[thisEdge].begin());

      // get the edge on the other side of the face
      localIndex nextEdge;
      if( localFacesToEdges[deadEndFace].first == thisEdge )
        nextEdge = localFacesToEdges[deadEndFace].second;
      else if( localFacesToEdges[deadEndFace].second == thisEdge )
        nextEdge = localFacesToEdges[deadEndFace].first;
      else
      {
//        return false;
        std::cout<<"nodeID, thisEdge = "<<nodeID<<", "<<thisEdge<<std::endl;
        std::cout<<"deadEndFace, localFacesToEdges[deadEndFace] = "<<", "<<deadEndFace<<", ( "<<localFacesToEdges[deadEndFace].first<<", "<<localFacesToEdges[deadEndFace].second<<" )"<<std::endl;

        std::cout<<"  nodeToEdges["<<nodeID<<"] = ( ";
        for( lSet::const_iterator i=nodeToEdges.begin() ; i!=nodeToEdges.end() ; ++i )
        {
          std::cout<<*i<<", ";
        }
        std::cout<<" )"<<std::endl;

        std::cout<<"  nodeToFaces["<<nodeID<<"] = ( ";
        for( lSet::const_iterator i=nodeToFaces.begin() ; i!=nodeToFaces.end() ; ++i )
        {
          std::cout<<*i<<", ";
        }
        std::cout<<" )"<<std::endl;

        std::cout<<"  faceToEdges["<<deadEndFace<<"] = ( ";
        for( lArray1d::const_iterator i=faceToEdges[deadEndFace].begin() ; i!=faceToEdges[deadEndFace].end() ; ++i )
        {
          std::cout<<*i<<", ";
        }
        std::cout<<" )"<<std::endl;

//        std::cout<<"  edgeToNodes[76922] = ( "<<edgeManager.m_edgesToNodes(76922,0)<<", "<<edgeManager.m_edgesToNodes(76922,1)<<" )"<<std::endl;


//        throw GPException("Fractunator2::FindFracturePlanes: Could not find the next edge when removing dead end faces.");
      }

      // delete the face from the working arrays
      edgesToRuptureReadyFaces[*edgeIndex].erase( deadEndFace );
      edgesToRuptureReadyFaces[nextEdge].erase(deadEndFace);
      nodeToRuptureReadyFaces.erase(deadEndFace);

      // if all the faces have been deleted, then go ahead and delete the top level entry
      if( edgesToRuptureReadyFaces[*edgeIndex].empty() )
        edgesToRuptureReadyFaces.erase(*edgeIndex);
      if( edgesToRuptureReadyFaces[nextEdge].empty() )
        edgesToRuptureReadyFaces.erase(nextEdge);

      // now increment the "thisEdge" to point to the other edge on the face that was just deleted
      thisEdge = nextEdge;
    }
  }

  if( nodeToRuptureReadyFaces.empty() )
  {
    return false;
  }






  // so now the working arrays have been purged of any faces that are on a dead-end path. All remaining faces
  // are part of a separation plane...of course, there can be more than one...which is bad. We will just take the first
  // path we find, and call this function again after the selected path is processed. That will happen automatically
  // since the new nodes that are created will have higher node indices than the current node, and will be checked for
  // separation prior to completion of the separation driver.



  // We now have to define the separation plane over which a node/face/edge will be split, and all elements on one side
  // of the plane get one set of objects, and all elements on the other side get the other set.


  // these are the edge and face where the separation surface starts.
  localIndex startingEdge = INT_MAX;
  localIndex startingFace = INT_MAX;


  // the startingEdge and startingFace needs to be set. It is best to start with an external face if we have one.
  // loop over all edges that have an entry in edgesToRuptureReadyFaces (i.e. all edges that are attached to a ruptured
  // node.

  for( std::map< localIndex, std::set<localIndex> >::const_iterator ke=edgesToRuptureReadyFaces.begin() ;
       ke!=edgesToRuptureReadyFaces.end() ; ++ke )
  {

    // make sure there is a face still attached to the edge, as it could have been removed when we got rid of dead ends
    // ...actually, this shouldn't ever happen as we have already removed such edges from the map.
    if( ke->second.size() > 0 )
    {

      startingEdge = ke->first;
      startingFace = *(ke->second.begin());
    }
    // if the size is 1, then the edge is only attached to one ruptured face, which means that it is external. This is
    // the case that we want. The starting edge and face already were set above, so just break at this point
    if( ke->second.size() == 1 && isEdgeExternal[ke->first]==1 )
    {
      break;
    }
  }

  // if the starting face was not set, then we don't have a rupture surface....so just quit.
  if( startingFace==INT_MAX )
    return false;


  // now we start the process of setting the separation path. Begin by
  localIndex thisEdge = startingEdge;
  localIndex thisFace = startingFace;

  localIndex nextEdge;

//  separationPathFaces.insert( startingFace );

  // the seprationPath is used to hold combinations of edge and face
  std::map<localIndex,int> facesInPath;
  std::map<localIndex,int> edgesInPath;

  int numFacesInPath = 0;
  edgesInPath[thisEdge] = numFacesInPath;
  facesInPath[thisFace] = numFacesInPath++;

  // now walk from face->edge->face->edge etc. until we get to another external edge, or back to the startingEdge.
  bool breakFlag = false;
  while ( !breakFlag )
  {
    // assign the other edge on the face as the next edge
    if( localFacesToEdges[thisFace].first == thisEdge )
    {
      nextEdge = localFacesToEdges[thisFace].second;
    }
    else if( localFacesToEdges[thisFace].second == thisEdge )
    {
      nextEdge = localFacesToEdges[thisFace].first;
    }
    else
    {
      throw GPException("Fractunator2::FindFracturePlanes breakpoint 2");
    }

    // if we have reached an external face, or the edge we started with, then we are done
    if( isEdgeExternal[nextEdge]==1 || edgesInPath.count(nextEdge)==1 )
    {
      const int startingIndex = edgesInPath[nextEdge];
      for( std::map<localIndex,int>::const_iterator kf=facesInPath.begin() ; kf!=facesInPath.end() ; ++kf )
      {
//        std::cout<<kf->first<<", "<<kf->second<<std::endl;
        if( kf->second >= startingIndex )
        {
          separationPathFaces.insert( kf->first );
        }
      }
      breakFlag = true;
    }
    else if( edgesToRuptureReadyFaces[nextEdge].size() )
    {
      // we need to pick another face attached to the "next edge"
      // increment the face and edge, and add to the separationPathFaces
      std::set<localIndex>::const_iterator iter_edgeToFace = edgesToRuptureReadyFaces[nextEdge].begin();
      if( *iter_edgeToFace == thisFace )
      {
        thisFace=*(++iter_edgeToFace);
      }
      else
      {
        bool pathFound = false;
        for( ; iter_edgeToFace!=edgesToRuptureReadyFaces[nextEdge].end() ; ++iter_edgeToFace )
        {
          if( *iter_edgeToFace == thisFace )
          {
            thisFace=*(--iter_edgeToFace);
            pathFound = true;
            break;
          }
        }
        if( pathFound == false )
          throw GPException("Fractunator2::FindFracturePlanes breakpoint 3");
      }

      thisEdge = nextEdge;
//      separationPathFaces.insert( thisFace );
      edgesInPath[thisEdge] = numFacesInPath;
      facesInPath[thisFace] = numFacesInPath++;
    }
    else
    {
      std::cout<<"you are screwed"<<std::endl;
      throw GPException("Fractunator2::FindFracturePlanes breakpoint 3b");
    }

  }





  // now we want to identify the objects on either side of the separation plane. First we assign an array to indicate
  // whether a face/edge is on the fracture plane.
  for( lSet::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
  {
    edgeLocations[*ke] = INT_MIN;
  }

  for( lSet::const_iterator ke=nodeToFaces.begin() ; ke!=nodeToFaces.end() ; ++ke )
  {
    faceLocations[*ke] = INT_MIN;
  }

  // now loop over the separation plane, and set the location arrays.
  for( lSet::const_iterator kf=separationPathFaces.begin() ; kf!=separationPathFaces.end() ; ++kf )
  {
    faceLocations[*kf] = -1;
    edgeLocations[localFacesToEdges[*kf].first] = -1;
    edgeLocations[localFacesToEdges[*kf].second] = -1;
  }


  for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=nodesToElements.begin() ; k!=nodesToElements.end() ; ++k )
  {
    elemLocations[*k] = INT_MIN;
  }

  SetLocations( 0, separationPathFaces, faceManager, nodesToElements, localFacesToEdges,
                edgeLocations, faceLocations, elemLocations );

  if( !(SetLocations( 1, separationPathFaces, faceManager, nodesToElements, localFacesToEdges,
                    edgeLocations, faceLocations, elemLocations )) )
  {
    return false;
  }



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

  if( fail )
  {
    return false;
    std::cout<<"Unset face or edge\n";
    std::cout<<"Separating node "<<nodeID<<" on separation plane consisting of faces: ";
    for( lSet::const_iterator i=separationPathFaces.begin() ; i!=separationPathFaces.end() ; ++i )
    {
      std::cout<<*i<<", ";
    }
    std::cout<<std::endl;


    for( std::map< std::pair< ElementRegionT*, localIndex >, int>::const_iterator iter_elem=elemLocations.begin() ; iter_elem!=elemLocations.end() ; ++iter_elem )
    {
      const std::pair< ElementRegionT*, localIndex >& elem = iter_elem->first;

      ElementRegionT& elemRegion = *(elem.first);
      const localIndex elemIndex = elem.second;

      const int location = iter_elem->second;


      std::cout<<"Element "<<elemIndex<<" at location "<<location<<"\n";


       std::cout<<" faces->edges->nodes = ";
       for( int a=0; a<6 ; ++a )
       {
         localIndex faceIndex = elemRegion.m_toFacesRelation(elemIndex,a);

         if( a>0 )
           std::cout<<"                      = ";



         if( faceManager.m_parentIndex[faceIndex] != LOCALINDEX_MAX )
         {
           std::cout<<faceManager.m_parentIndex[faceIndex]<<"->";
         }
         std::cout<<faceIndex<<"[ ";
         for( int b=0 ; b<4 ; ++b )
         {
           localIndex edgeIndex = faceManager.m_toEdgesRelation[faceIndex][b];
           std::cout<<edgeIndex<<", ";
         }
         std::cout<<" ] \n";

       }
       std::cout<<std::endl;

    }

    for( std::map< std::pair< ElementRegionT*, localIndex >, int>::const_iterator i=elemLocations.begin() ; i!=elemLocations.end() ; ++i )
    {
      std::cout<<"( "<<i->first.second<<","<<i->second<<") , ";
    }
    std::cout<<std::endl;


    std::cout<<"faceLocations = ";
    for( std::map<localIndex,int>::const_iterator i=faceLocations.begin() ; i!=faceLocations.end() ; ++i )
    {
      std::cout<<"( "<<i->first<<","<<i->second<<") , ";
    }
    std::cout<<std::endl;

    std::cout<<"edgeLocations = ";
    for( std::map<localIndex,int>::const_iterator i=edgeLocations.begin() ; i!=edgeLocations.end() ; ++i )
    {
      std::cout<<"( "<<i->first<<","<<i->second<<") , ";
    }
    std::cout<<std::endl;
    throw GPException("Fractunator2::FindFracturePlanes breakpoint 6...smart guy");
  }

return true;
}



bool Fractunator2::SetLocations( const int location,
                                const lSet& separationPathFaces,
                                const FaceManagerT& faceManager,
                                const std::set< std::pair<ElementRegionT*,localIndex> >& nodesToElements,
                                std::map< localIndex, std::pair<localIndex,localIndex> >& localFacesToEdges,
                                std::map<localIndex,int>& edgeLocations,
                                std::map<localIndex,int>& faceLocations,
                                std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{
  bool rval = true;
  const localIndex separationFace = *(separationPathFaces.begin());
  std::set< std::pair<ElementRegionT*,localIndex> > elem0 ;
  std::set< std::pair<ElementRegionT*,localIndex> > processedElements ;

  // insert an element attached to the separation face
  elem0.insert( faceManager.m_toElementsRelation[separationFace][location] );
//  elemLocations[faceManager.m_FaceToElementMap[separationFace][location]] = location;

  if( m_verbose ) std::cout<<"Setting Location "<<location<<std::endl;


  bool addedElem = false;
  do
  {
    addedElem = false;
    for( std::set< std::pair<ElementRegionT*,localIndex> >::iterator k=elem0.begin() ; k!=elem0.end() ; ++k )
    {
      // make sure that we have not already processed the element
      if( processedElements.count(*k)==0 )
      {
        processedElements.insert(*k);
  //      if( m_verbose ) std::cout<<"  processing Element "<<k->second<<std::endl;
        for( localIndex kf=0 ; kf<k->first->m_toFacesRelation.Dimension(1) ; ++kf )
        {

          const localIndex faceIndex = k->first->m_toFacesRelation(k->second,kf);

  //        if( m_verbose ) std::cout<<"    processing Face "<<faceIndex<<std::endl;

          std::map<localIndex,int>::iterator iterFace = faceLocations.find(faceIndex);

          if( iterFace != faceLocations.end() )
          {
            // make sure the face is not on the separation plane
            if( iterFace->second != -1 )
            {
              iterFace->second = location;

              if( edgeLocations[localFacesToEdges[faceIndex].first] == INT_MIN )
                edgeLocations[localFacesToEdges[faceIndex].first] =location;
              if( edgeLocations[localFacesToEdges[faceIndex].second] == INT_MIN )
                edgeLocations[localFacesToEdges[faceIndex].second] = location;


              // now we add the element that is a neighbor to the "linkingFace"
              // of course, this only happens if there are more than one element
              // attached to the face.
              if( faceManager.m_toElementsRelation[faceIndex].size() > 1 )
              {
                const std::pair<ElementRegionT*,localIndex>& elemIndex0 = faceManager.m_toElementsRelation[faceIndex][0];
                const std::pair<ElementRegionT*,localIndex>& elemIndex1 = faceManager.m_toElementsRelation[faceIndex][1];

                // if the first element is the one we are on, and the element is attached
                // to the splitting node, then add the second element to the list.
                if( ( elemIndex0 == *k ) && ( nodesToElements.find(elemIndex1)!=nodesToElements.end() ) )
                {
                  elem0.insert(elemIndex1);
                  addedElem = true;
                }
                // if the second element is the one we are on, and the element is attached
                // to the splitting node, then add the first element to the list.
                else if( ( elemIndex1 == *k ) && ( nodesToElements.find(elemIndex0)!=nodesToElements.end() ) )
                {
                  elem0.insert(elemIndex0);
                  addedElem = true;
                }
              }
            }
          }


        }
        if( addedElem == true)
        {
          break;
        }
      }
    }
  }while(addedElem==true);


  for( std::set< std::pair<ElementRegionT*,localIndex> >::const_iterator k=elem0.begin() ; k!=elem0.end() ; ++k  )
  {
    if( m_verbose ) std::cout<<"  Setting Element "<<k->second<<" to location "<<location<<std::endl;

    if( elemLocations.find(*k) != elemLocations.end() )
    {
//      std::cout<<k->second<<", "<<elemLocations[*k]<<std::endl;
      if( elemLocations[*k]==INT_MIN )
      {
        elemLocations[*k] = location;
      }
      else
      {
        rval = false;
      }
    }

  }

return rval;
}

void Fractunator2::PerformFracture( const localIndex parentNodeIndex,
                                   NodeManagerT& nodeManager,
                                   EdgeManagerT& edgeManager,
                                   FaceManagerT& faceManager,
                                   ExternalFaceManagerT& externalFaceManager,
                                   ElementManagerT& elementManager,
                                   ModifiedObjectLists& modifiedObjects ,
                                   Array1dT<lSet>& nodesToRupturedFaces,
                                   Array1dT<lSet>& edgesToRupturedFaces,
                                   const lSet& separationPathFaces,
                                   const std::map<localIndex,int>& edgeLocations,
                                   const std::map<localIndex,int>& faceLocations,
                                   const std::map< std::pair< ElementRegionT*, localIndex >, int>& elemLocations )
{

  const int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, const_cast<int*>(&rank) );

  const OrderedVariableOneToManyRelation& childEdgeIndex = edgeManager.GetVariableOneToManyMap( "childIndices" );
  const OrderedVariableOneToManyRelation& childFaceIndex = faceManager.GetVariableOneToManyMap( "childIndices" );




  // ***** split all the objects first *****

  // Split the node into two, using the original index, and a new one.
  localIndex newNodeIndex;
  if( m_verbose )
  {
    std::cout<<"\nSplitting node "<<parentNodeIndex<<" along separation plane faces ";
    for( lSet::const_iterator i=separationPathFaces.begin() ; i!=separationPathFaces.end() ; ++i )
    {
      std::cout<<*i<<", ";
    }
    std::cout<<std::endl;
  }


  nodeManager.SplitObject(parentNodeIndex, rank, newNodeIndex);

  if( m_verbose ) std::cout<<"\nDone splitting node "<<parentNodeIndex<<" into nodes "<<parentNodeIndex<<" and "<<newNodeIndex<<std::endl;

  // HACK...the node to element map is a bastard-child and is not managed by the
  // database
  nodeManager.m_toElementsRelation.resize( nodeManager.m_numNodes );
  nodesToRupturedFaces.push_back( nodesToRupturedFaces.back() );

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

      for( int a=0 ; a<2 ; ++a )
      {
        edgeManager.m_toNodesRelation(newEdgeIndex,a) = edgeManager.m_toNodesRelation(parentEdgeIndex,a);
      }

      edgesToRupturedFaces.push_back( edgesToRupturedFaces.back() );

    } //    if( location == -1  )
  } // for( std::map<localIndex,int>::const_iterator iter_edge...





  // split the faces
  iArray1d& ruptureState = faceManager.GetFieldData<int>("ruptureState");
  std::map<localIndex,localIndex> splitFaces;

  // loop over all faces attached to the nodeID
  for( std::map<localIndex,int>::const_iterator iter_face=faceLocations.begin() ; iter_face!=faceLocations.end() ; ++iter_face )
  {
    const localIndex& parentFaceIndex = iter_face->first;
    const int& location = iter_face->second;
    // if the face is on the separation plane, then split it
    if( location == -1 )
    {
      localIndex newFaceIndex;

      if( faceManager.SplitObject( parentFaceIndex, rank, newFaceIndex ) )
      {
        if( m_verbose ) std::cout<<"  Split face "<<parentFaceIndex<<" into faces "<<parentFaceIndex<<" and "<<newFaceIndex<<std::endl;

        splitFaces[parentFaceIndex] = newFaceIndex;

        ruptureState[newFaceIndex] = 0;


        faceManager.m_toEdgesRelation[newFaceIndex] = faceManager.m_toEdgesRelation[parentFaceIndex];

        faceManager.m_toNodesRelation[newFaceIndex] = faceManager.m_toNodesRelation[parentFaceIndex];

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
          if( elementToNodeMap[a] == parentNodeIndex )
          {

            if( m_verbose )
              std::cout<<elementToNodeMap[a]<<"->"<<newNodeIndex<<", ";

            elementToNodeMap[a] = newNodeIndex;

            nodeManager.m_toElementsRelation[newNodeIndex].insert(elem);
            nodeManager.m_toElementsRelation[parentNodeIndex].erase(elem);

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
            faceManager.m_toElementsRelation[faceIndex].erase( faceManager.m_toElementsRelation[faceIndex].begin() );
          else if( faceManager.m_toElementsRelation[faceIndex][1] == elem )
            faceManager.m_toElementsRelation[faceIndex].erase( faceManager.m_toElementsRelation[faceIndex].begin()+1 );

        } // if( splitFaces.count( faceID ) > 0 )




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
          if( *nodeIndex == parentNodeIndex )
          {
            *nodeIndex = newNodeIndex;

            // if it is not a new face.
            if( !isNewFace )
            {
              // remove the face from the nodeToFaceMap of the parent node.
              nodeManager.m_nodeToFaceMap[parentNodeIndex].erase(faceIndex);

              // add the face to the nodeToFaceMap of the new node.
              nodeManager.m_nodeToFaceMap[*nodeIndex].insert(faceIndex);

              // if the face is ruptured, then it should be added to the
              // nodeToRupturedFaces map.
//              if( faceRuptureState[faceIndex ] )
//              nodesToRupturedFaces[*nodeIndex].erase(faceIndex);
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

          if( isNewFace )
          {
            nodesToRupturedFaces[parentNodeIndex].erase(faceIndex);
            nodesToRupturedFaces[parentNodeIndex].erase(newFaceIndex);
            nodesToRupturedFaces[newNodeIndex].erase(faceIndex);
            nodesToRupturedFaces[newNodeIndex].erase(newFaceIndex);
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

          if( isNewFace )
          {
            edgesToRupturedFaces[*edgeIndex].erase(faceIndex);
            edgesToRupturedFaces[*edgeIndex].erase(newFaceIndex);
            if( splitEdges.count( *edgeIndex ) > 0 )
            {
              edgesToRupturedFaces[childEdgeIndex[*edgeIndex][0]].erase(faceIndex);
              edgesToRupturedFaces[childEdgeIndex[*edgeIndex][0]].erase(newFaceIndex);
            }
          }



          // if the edge was just split
          if( splitEdges.count( *edgeIndex ) > 0 )
          {
            if( faceIndex == newFaceIndex )
              edgeManager.m_toFacesRelation[*edgeIndex].erase(faceIndex);
            *edgeIndex = childEdgeIndex[*edgeIndex][0];
          }
          edgeManager.m_toFacesRelation[*edgeIndex].insert(newFaceIndex);
          if( m_verbose ) std::cout<<*edgeIndex;

//          if( faceRuptureState[newFaceIndex ] )
            edgesToRupturedFaces[*edgeIndex].erase(faceIndex);




          //edgesToNodes
          if( m_verbose )
          {
            std::cout<<"(";
          }

          {
            localIndex* const nodeIndex = edgeManager.m_toNodesRelation[*edgeIndex];
            for( unsigned int a=0 ; a<edgeManager.m_toNodesRelation.Dimension(1) ; ++a )
            {
              if( nodeIndex[a] == parentNodeIndex )
              {

                if( m_verbose ) std::cout<<nodeIndex[a];

                nodeIndex[a] = newNodeIndex;
                nodeManager.m_nodeToEdgeMap[parentNodeIndex].erase(*edgeIndex);

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

      localIndex* const elemToFaces = elemRegion.m_toFacesRelation[elemIndex];

      lArray1d facelist;
      // first fill the facelist with the faces in elemToFaces
      for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
      {
        facelist.push_back(elemToFaces[kf]);
      }
      for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
      {
        const localIndex parentFaceIndex = faceManager.m_parentIndex[elemToFaces[kf]];
        if( parentFaceIndex!=LOCALINDEX_MAX )
        {
          facelist.push_back(parentFaceIndex);
        }
      }

      // Now we do a loop over the facelist and process all the faces
      for( int kf=0 ; kf<int(facelist.size()) ; ++kf )
      {
        lSet faceNodes;

        localIndex faceIndex ;
        if( kf< elemRegion.m_numFacesPerElement )
          faceIndex = elemRegion.m_toFacesRelation(elemIndex,kf);
        else
          faceIndex = facelist[kf];

        if( kf>0 )
          std::cout<<"                              = ";

        if( kf>=elemRegion.m_numFacesPerElement )
          std::cout<<"P-";

        std::cout<<faceIndex<<"( ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex nodeID = faceManager.m_toNodesRelation[faceIndex][b];
          faceNodes.insert(nodeID);
          if( elemNodes.count(nodeID) == 0 && kf<elemRegion.m_numFacesPerElement )
            std::cout<<"*";
          std::cout<<nodeID<<",";
        }
        std::cout<<" )      ";



        std::cout<<faceIndex<<"[ ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex edgeIndex = faceManager.m_toEdgesRelation[faceIndex][b];
          std::cout<<edgeIndex<<"( ";
          for( int c=0 ; c<2 ; ++c )
          {
            localIndex nodeID = edgeManager.m_toNodesRelation(edgeIndex,c);
            if( elemNodes.count(nodeID) == 0  && kf<elemRegion.m_numFacesPerElement )
              std::cout<<"*";
            if( faceNodes.count(nodeID) == 0 )
              std::cout<<"#";
            std::cout<<nodeID<<",";
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
        localIndex nodeID = faceManager.m_toNodesRelation[faceIndex][b];
        faceNodes.insert(nodeID);
//          if( elemNodes.count(nodeID) == 0 )
//            std::cout<<"*";
        std::cout<<nodeID<<",";
      }
      std::cout<<" )      ";



      std::cout<<faceIndex<<"[ ";
      for( int b=0 ; b<4 ; ++b )
      {
        localIndex edgeIndex = faceManager.m_toEdgesRelation[faceIndex][b];
        std::cout<<edgeIndex<<"( ";
        for( int c=0 ; c<2 ; ++c )
        {
          localIndex nodeID = edgeManager.m_toNodesRelation(edgeIndex,c);
//            if( elemNodes.count(nodeID) == 0 )
//              std::cout<<"*";
          if( faceNodes.count(nodeID) == 0 )
            std::cout<<"#";
          std::cout<<nodeID<<",";
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


}

void Fractunator2::PreexistingFracture2D( NodeManagerT& nodeManager,
                                EdgeManagerT& edgeManager,
                                FaceManagerT& faceManager,
                                ExternalFaceManagerT& externalFaceManager,
                                ElementManagerT& elementManager,
                                SpatialPartition& partition,
                                const bool prefrac)
{
// Do nothing.
}
