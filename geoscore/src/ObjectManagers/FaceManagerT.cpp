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
 * @file FaceManagerT.cpp
 * @author settgast1
 * @date Feb 4, 2011
 */


#include "FaceManagerT.h"
#include "ExternalFaceManagerT.h"
#include "ElementManagerT.h"
#include "PhysicalDomainT.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "IO/BinStream.h"
#include "Utilities/Utilities.h"
#include "ArrayT/Array3dT.h"

#include "BoundaryConditions/BoundaryConditions.h"

#include "ExternalFaceStructs.h"

#include "DataStructures/VectorFields/TempODS.h"

#include "Constitutive/CohesiveZone/CohesiveZoneFactory.h"
#include "Constitutive/Material/MaterialFactory.h"

#include <limits.h>

/**
 *
 * @return
 */
FaceManagerT::FaceManagerT()
:ObjectDataStructureBaseT( ObjectDataStructureBaseT::FaceManager ),
m_numFaces(this->m_DataLengths),
m_toNodesRelation(m_VariableOneToManyMaps["FaceToNodeMap"]),
m_toEdgesRelation(m_VariableOneToManyMaps["FaceToEdgeMap"]),
m_externalFaces(m_Sets["ExternalFaces"]),
#if USECPP11==1
m_cohesiveZone(),
#else
m_cohesiveZone(NULL),
#endif
m_elemManagerHACK(NULL)
{
  //0-based; note that the following field is ALSO 0
  //for faces that are not external faces, so check isExternal before using
  this->AddKeylessDataField<localIndex>("externalFaceIndex", true, true);

  this->AddKeylessDataField<R1Tensor>("FaceCenter",true,true);

//  //@annavarapusr1: Adding stabilization and weighting fields for Nitsche's method
//  this->AddKeylessDataField<realT> ( "nitscheStab_n" , false, false );
//  this->AddKeylessDataField<realT> ( "nitscheStab_t1" , false, false );
//  this->AddKeylessDataField<realT> ( "nitscheStab_t2" , false, false );
//  this->AddKeylessDataField<realT> ( "nitscheGamma" , false, false );

}

/**
 *
 * @return
 */
FaceManagerT::~FaceManagerT()
{
#if USECPP11!=1
  if( m_cohesiveZone )
    delete m_cohesiveZone;
#endif

}

void FaceManagerT::BuildFaces( const NodeManagerT& nodeManager, const ElementManagerT& elementManager )
{

  lArray1d tempNodeList;
  Array1dT<lArray1d> tempFaceToNodeMap;

  localIndex numFaces = 0;
  Array1dT<lArray1d> facesByLowestNode;



  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::const_iterator elementRegionIter = elementManager.m_ElementRegions.begin() ;
       elementRegionIter != elementManager.m_ElementRegions.end() ;
       ++elementRegionIter )
  {
    const ElementRegionT& elementRegion = elementRegionIter->second;

//    if( elementRegion.m_ElementDimension == 3 )
    {
      for( localIndex k=0 ; k<elementRegion.m_numElems ; ++k )
      {
        // kelf = k'th element local face index
        for( int kelf=0 ; kelf<elementRegion.m_numFacesPerElement ; ++kelf )
        {
          // get the nodes associated with the local face
          elementRegion.GetFaceNodes( k, kelf, tempNodeList );

          // sort the nodes
          std::sort(tempNodeList.begin(), tempNodeList.end() );

          // get the lowest node index from the list for simplicity
          const localIndex& lowNode = tempNodeList[0];

          // now check to see if the lowest node index has an entry in the facesByLowestNode vector
          if( facesByLowestNode.size() < (lowNode+1) )
          {
            // the node has not been entered, so add it.
            facesByLowestNode.resize(lowNode+1);

            // this a new face, so add it,
            AddNewFace( k, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, elementRegion );
          }
          else
          {
            // does the node have an entry? If not, then this has to be a new face.
            if( facesByLowestNode[lowNode].empty() )
            {
              // this a new face, so add it,
              AddNewFace( k, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, elementRegion );
            }
            else
            {
              // the node does have an entry, so it is possible that the facet has already be assigned a number

              // make a flag to indicate whether the face is a duplicate...assume that it isn't unless this is disproved.
              bool duplicate = false;


              // there are faces in facesByLowestNode, so lets loop over them and check for duplicates
              for( lvector::iterator existingFaceIndex = facesByLowestNode[lowNode].begin() ;
                  existingFaceIndex != facesByLowestNode[lowNode].end() ; ++existingFaceIndex )
              {
                // this is the nodelist of the face that we are testing agains
                const lArray1d& existingFaceNodelist = tempFaceToNodeMap[*existingFaceIndex];

                // test to see if the size of the nodelists are the same....
                if( existingFaceNodelist.size() == tempNodeList.size() )
                {
                  // since the size is the same, then we should test the nodes...they are sorted, so
                  // the std::equal() algorithm will work for this.
                  if( std::equal( existingFaceNodelist.begin(), existingFaceNodelist.end(), tempNodeList.begin() ) )
                  {
                    // they are equal!
                    duplicate = true;

                    // add the element to the faceToElement map
                    m_toElementsRelation[*existingFaceIndex].push_back( std::pair<ElementRegionT*, localIndex>( const_cast<ElementRegionT*>(&elementRegion), k) );

                    // add the face to the elementToFaceMap for the element region.
                    elementRegion.m_toFacesRelation(k,kelf) = *existingFaceIndex;

                    // now remove the entry from the face that we were checking against from the facesByLowestNode list...
                    // because it is no longer possible that it will have another element that has this face.
		    facesByLowestNode[lowNode].erase( existingFaceIndex );

                    // break the loop
                    break;
                  }
                }
              }
              if( !duplicate )
              {
                // the face is not a duplicate of any in the facesByLowestNode list, so we need to add a new face.
                AddNewFace( k, kelf, numFaces, facesByLowestNode, tempNodeList, tempFaceToNodeMap, elementRegion );
              }
            }
          }
        }
      }
    }
  }

  // resize the data vectors according to the number of faces indicated by the size of tempFaceToNodeMap
  this->resize(tempFaceToNodeMap.size());

  // set m_FaceToNodeMap
  this->m_toNodesRelation = tempFaceToNodeMap;


  // make sets from nodesets
  for( std::map< std::string, lSet >::const_iterator i = nodeManager.m_Sets.begin() ;
       i != nodeManager.m_Sets.end() ; ++i )
  {
    const std::string& setname = i->first;
    const lSet& set = i->second;
    this->ConstructSetFromSetAndMap( set, this->m_toNodesRelation, setname );
  }

  
  // sort the face node lists
  SortAllFaceNodes(nodeManager);


  Array1dT<R1Tensor>& faceCenter = this->GetFieldData<R1Tensor>( "FaceCenter" );
  for( localIndex k=0 ; k<m_numFaces ; ++k )
  {
    FaceCenter( nodeManager, k , faceCenter[k] );

  }
}




void FaceManagerT::AddNewFace( const localIndex& k,
                               const localIndex& kelf,
                               localIndex& numFaces,
                               Array1dT<lArray1d>& facesByLowestNode,
                               lArray1d& tempNodeList,
                               Array1dT<lArray1d>& tempFaceToNodeMap,
                               const ElementRegionT& elementRegion )
{
  Array1dT< std::pair< ElementRegionT*, localIndex > > tempFaceToElemEntry;

  // and add the face to facesByLowestNode[]
  facesByLowestNode[tempNodeList[0]].push_back(numFaces);


  // add the face to the elementToFaceMap
  elementRegion.m_toFacesRelation(k,kelf) = numFaces;

  // add the nodes to the faceToNodeMap
  tempFaceToNodeMap.push_back(tempNodeList);

  // add to the element information to the faceToElementMap
  tempFaceToElemEntry.push_back( std::pair<ElementRegionT*, localIndex>( const_cast<ElementRegionT*>(&elementRegion), k) );
  m_toElementsRelation.push_back(tempFaceToElemEntry);

  // now increment numFaces to reflect the number of faces rather than the index of the new face
  ++numFaces;
}

void FaceManagerT::SetDomainBoundaryObjects(const ObjectDataStructureBaseT* const referenceObject  )
{
  // assume that all faces will be on a domain boundary. Set it to zero if it is found to have two elements that it
  // is connected to.
  iArray1d& isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();
  isDomainBoundary = 1;

  for( localIndex kf=0 ; kf<m_numFaces ; ++kf )
  {
    if( m_toElementsRelation(kf).size() == 2 )
    {
      isDomainBoundary(kf) = 0;
    }
  }

}


void FaceManagerT::SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  )
{
  const iArray1d& isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();
  m_isExternal = 0;
  for( localIndex k=0 ; k<m_numFaces ; ++k )
  {
    if( isDomainBoundary[k]==1 && ( m_matchedBoundaryFaces.find(k)==m_matchedBoundaryFaces.end() ) )
    {
      m_isExternal[k] = 1;
      m_externalFaces.insert(k);
    }
  }
}

















/*
 * @author Scott Johnson
 * @brief Returns the centroid of the face 
 * @param[in] nodeManager Node manager that owns the face manager's nodes
 * @param[in] faceIndex Local index of the face to query
 * @param[out] center Center of the face
 */
realT FaceManagerT::FaceCenter( const NodeManagerT& nodeManager, const localIndex faceIndex, R1Tensor& center ) const
{
  center = 0.0;
  const lArray1d& nodeList = m_toNodesRelation[faceIndex];
  return GeometryUtilities::Centroid_3DPolygon(nodeList,
                                        *(nodeManager.m_refposition),
                                        *(nodeManager.m_displacement),
                                        center);
}

/*
 * @author Scott Johnson
 * @brief Returns the centroid of the face
 * @param[in] nodeManager Node manager that owns the face manager's nodes
 * @param[in] faceIndex Local index of the face to query
 * @param[out] center Center of the face
 * @param[out] normal Normal to the face
 */
realT FaceManagerT::FaceCenterAndNormal( const NodeManagerT& nodeManager,
                                         const localIndex faceIndex,
                                         R1Tensor& center,
                                         R1Tensor& normal,
                                         const bool referenceState ) const
{
  center = 0.0;
  normal = 0.0;
  if(referenceState)
  {
    return GeometryUtilities::Centroid_3DPolygon(m_toNodesRelation[faceIndex],
                                                 *(nodeManager.m_refposition),
                                                 center,
                                                 normal );
  }
  else
  {
    return GeometryUtilities::Centroid_3DPolygon(m_toNodesRelation[faceIndex],
                                                 *(nodeManager.m_refposition),
                                                 *(nodeManager.m_displacement),
                                                 center, normal);
  }
}

/*
 * @author Scott Johnson
 * @brief Returns the normal of the face
 * @param[in] nodeManager Node manager that owns the face manager's nodes
 * @param[in] faceIndex Local index of the face to query
 * @return Normal to the face
 */
R1Tensor FaceManagerT::FaceNormal( const NodeManagerT& nodeManager, const localIndex faceIndex, const bool referenceFlag) const
{
  R1Tensor center(0.0);
  R1Tensor normal(0.0);

  if( referenceFlag )
  {
    GeometryUtilities::FaceNormal3DPolygonReferenceConfig( m_toNodesRelation[faceIndex] ,
                                                          *(nodeManager.m_refposition),
                                                          normal );
  }
  else
  {
    GeometryUtilities::Centroid_3DPolygon( m_toNodesRelation[faceIndex],
                                          *(nodeManager.m_refposition),
                                          *(nodeManager.m_displacement),
                                          center,
                                          normal );
  }

  return normal;
}

/*
 * @author Scott Johnson
 * @brief Returns the normal of the face
 * @param[in] nodeManager Node manager that owns the face manager's nodes
 * @param[in] faceIndex Local index of the face to query
 * @param[out] normal Normal to the face
 */
realT FaceManagerT::FaceNormal( const NodeManagerT& nodeManager, const localIndex faceIndex, R1Tensor& normal) const
{
  R1Tensor center(0.0);
  normal = 0.0;
  return GeometryUtilities::Centroid_3DPolygon(m_toNodesRelation[faceIndex],
                                               *(nodeManager.m_refposition),
                                               *(nodeManager.m_displacement),
                                               center, normal);
}

/*
 * @author Pengcheng Fu
 * @brief Returns a pair of normalized tangential vectors
 * @param[in] nodeManager Node manager that owns the face manager's nodes
 * @param[in] faceIndex Local index of the face to query
 * @param[out] vectors (perpendicular to each other) along the face; non unique.
 */
void FaceManagerT::FaceTangential( const NodeManagerT& nodeManager, const localIndex faceIndex, R1Tensor& tanA, R1Tensor& tanB) const
{
  tanA = (*nodeManager.m_refposition)[m_toNodesRelation[faceIndex][1]];
  tanA += (*nodeManager.m_displacement)[m_toNodesRelation[faceIndex][1]];
  tanA -= (*nodeManager.m_refposition)[m_toNodesRelation[faceIndex][0]];
  tanA -= (*nodeManager.m_displacement)[m_toNodesRelation[faceIndex][0]];
  tanA.Normalize();

  if (m_toNodesRelation[faceIndex].size() <= 2)  //2D
  {
    tanB[0] = 0.0;
    tanB[1] = 0.0;
    tanB[2] = 1.0;
  }
  else
  {
    R1Tensor norm = FaceNormal(nodeManager, faceIndex);
    tanB.Cross(norm, tanA);
    tanB.Normalize();
  }
}

/*
 * @author Scott Johnson
 *
 * @brief Surface area by summing areas of triangular facets on the face.
 *
 *
 *
 */
realT FaceManagerT::SurfaceArea( const NodeManagerT& nodeManager,
                                 const localIndex faceIndex,
                                 const bool referenceFlag ) const
{
  const lArray1d& nodeList = m_toNodesRelation[faceIndex];

  R1Tensor junk1, junk2;

  if( referenceFlag )
  {
    return GeometryUtilities::Centroid_3DPolygon( nodeList ,
                                                  *(nodeManager.m_refposition),
                                                  junk1,
                                                  junk2 );
  }
  else
  {
    return GeometryUtilities::Centroid_3DPolygon(nodeList,
                                                 *(nodeManager.m_refposition),
                                                 *(nodeManager.m_displacement),
                                                 junk1,
                                                 junk2 );
  }
}

/*
 * @author Scott Johnson
 *
 * Returns the area of a face projected onto a surface with the given normal
 *
 */
realT FaceManagerT::ProjectedArea( const NodeManagerT& nodeManager, const localIndex faceIndex, const R1Tensor& norm) const
{
  const lArray1d& nodeList = m_toNodesRelation[faceIndex];
  Array1dT<R1TensorT<2> > positionsProjected;
  R1TensorT<2> min, max;
  R1Tensor e1, e2;
  if(!GeometryUtilities::VectorsInPlane(norm, e1, e2))
    throw GPException("Normal is pathologically 0");
  GeometryUtilities::CartesianPointsProjectedToPlanarPoints(nodeList,
                                                            *(nodeManager.m_refposition),
                                                            *(nodeManager.m_displacement),
                                                            e1, e2,
                                                            positionsProjected,
                                                            min, max);
  return GeometryUtilities::Area_2DPolygon(positionsProjected);
}

/**
 * @author Scott Johnson
 * @brief Get the bounding sphere center and radius for the given face
 * @param[in] nodeManager Node manager associated with the current FaceManager
 * @param[in] faceIndex Local face index to query
 * @param[out] center Center of the face in global frame
 * @param[out] radius Radius of the bounding sphere about the face center
 */
void FaceManagerT::FaceBoundingSphere( const NodeManagerT& nodeManager, const localIndex faceIndex, R1Tensor& center, realT& radius) const
{
  FaceCenter(nodeManager, faceIndex, center);
  radius = 0.;
  R1Tensor rr(0.);
  for( localIndex a=0 ; a<m_toNodesRelation[faceIndex].size() ; ++a )
  {
    const localIndex nodeIndex = m_toNodesRelation[faceIndex][a];
    rr = center;
    rr -= (*nodeManager.m_refposition)[nodeIndex];
    rr -= (*nodeManager.m_displacement)[nodeIndex];
    const realT tmp = rr.L2_Norm();
    radius = tmp > radius ? tmp : radius;
  }
}

/**
 * @brief Get center, velocity, radius
 * @author Scott Johnson
 * @param[in] nodeManager Node object manager
 * @param[in] faceIndex Index of the face
 * @param[in] dt Timestep
 * @param[out] xfc Center of the face
 * @param[out] dxfc Velocity of the face at the mid-step
 * @param[out] radius Radius of the bounding sphere including displacement over the step
*/
void FaceManagerT::FaceBoundingSphere(
    const NodeManagerT& nodeManager,
    const localIndex faceIndex,
    const realT dt,
    R1Tensor& center,
    R1Tensor& velocity,
    realT& radius,
    R1Tensor& normal,
    const bool referenceState) const
{
  FaceCenterAndNormal(nodeManager, faceIndex, center, normal);

  radius = 0.0;
  R1Tensor rr(0.0);
  velocity = 0.0;

  const lArray1d& nodeList = m_toNodesRelation[faceIndex];
  for( localIndex a=0 ; a<nodeList.size() ; ++a )
  {
    const localIndex nodeIndex = nodeList[a];
    velocity += (*nodeManager.m_velocity)[nodeIndex];
    rr = center;
    rr -= (*nodeManager.m_refposition)[nodeIndex];
    if(!referenceState)
      rr -= (*nodeManager.m_displacement)[nodeIndex];
    const realT rtmp = Dot(rr,rr);
    radius = rtmp > radius ? rtmp : radius;
  }
  radius = sqrt(radius);
  velocity *= 1.0/nodeList.size();
  if(!referenceState)
    radius += 0.5 * dt * velocity.L2_Norm();
}

/**
 * @brief Nodal positions on the face
 * @author Scott Johnson
 * @param[in] faceIndex Index of the face
 * @param[out] xs List of face node global coordinates
*/
void FaceManagerT::NodalPositions(  const NodeManagerT& nodeManager,
                                    const localIndex faceIndex,
                                    Array1dT<R1Tensor>& xs) const
{
  xs.clear();

  const lArray1d& nodeList = m_toNodesRelation[faceIndex];
  for( localIndex a=0 ; a<nodeList.size() ; ++a )
  {
    const localIndex nodeIndex = nodeList[a];
    R1Tensor tmp = (*nodeManager.m_refposition)[nodeIndex];
    tmp += (*nodeManager.m_displacement)[nodeIndex];
    //add the position to the nodal position list
    xs.push_back(tmp);
  }
}

/**
 * @brief Get area, length, center, velocity, normal, etc of face
 * @author Scott Johnson
 * @param[in] nodeManager Node object manager
 * @param[in] faceIndex Index of the face
 * @param[in] dt Timestep
 * @param[out] xfc Center of the face
 * @param[out] dxfc Velocity of the face at the mid-step
 * @param[out] xmin AABB lower coordinate
 * @param[out] xmax AABB upper coordinate
 * @param[out] area Area of the face
 * @param[out] length Length of the face
 * @param[out] dxs List of face node velocities at the mid-step
*/
void FaceManagerT::FaceProperties( const NodeManagerT& nodeManager,
                                   const localIndex faceIndex,
                                   const realT dt,
                                   R1Tensor & xfc,
                                   R1Tensor & dxfc,
                                   ExternalFaceStruct& efs) const
{
  efs.xs.clear();
  efs.dxs.clear();

  FaceCenter(nodeManager, faceIndex, xfc);
  efs.area = SurfaceArea(nodeManager, faceIndex);

  efs.xmin = std::numeric_limits<realT>::max();
  efs.xmax = efs.xmin;
  efs.xmax *= -1.0;
  dxfc = 0.;

  const lArray1d& nodeList = m_toNodesRelation[faceIndex];
  efs.xs.reserve(nodeList.size());
  efs.dxs.reserve(nodeList.size());
  for( localIndex a=0 ; a<nodeList.size() ; ++a )
  {
    const localIndex nodeIndex = nodeList[a];
    R1Tensor tmp = (*nodeManager.m_refposition)[nodeIndex];
    tmp += (*nodeManager.m_displacement)[nodeIndex];
    efs.xmin.SetMin(tmp);
    efs.xmax.SetMax(tmp);
    //add the position to the nodal position list
    efs.xs.push_back(tmp);
    //add the velocity pushed forward by 1/2 step to the nodal velocity list
    tmp = (*nodeManager.m_acceleration)[nodeIndex];
    tmp *= 0.5 * dt;
    tmp += (*nodeManager.m_velocity)[nodeIndex];
    efs.dxs.push_back(tmp);
    //add the velocity contribution to the face
    dxfc += tmp;
  }
  dxfc *= 1.0/nodeList.size();
}

/**
 * @brief Get area, length, center, velocity, normal, etc of face
 * @author Scott Johnson
 * @param[in] nodeManager Node object manager
 * @param[in] faceIndex Index of the face
 * @param[in] dt Timestep
 * @param[out] xmin AABB lower coordinate
 * @param[out] xmax AABB upper coordinate
 * @param[out] area Area of the face
 * @param[out] length Length of the face
 * @param[out] xs List of face node global coordinates
 * @param[out] dxs List of face node velocities at the mid-step
*/
void FaceManagerT::FaceProperties( const NodeManagerT& nodeManager,
                                   const localIndex faceIndex,
                                   const realT dt,
                                   ExternalFaceStruct& efs) const
{
  efs.dxs.clear();
  efs.xs.clear();

  efs.area = SurfaceArea(nodeManager, faceIndex);

  efs.xmin = std::numeric_limits<realT>::max();
  efs.xmax = efs.xmin;
  efs.xmax *= -1.0;

  const lArray1d& nodeList = m_toNodesRelation[faceIndex];
  efs.dxs.reserve(nodeList.size());
  efs.xs.reserve(nodeList.size());
  for( localIndex a=0 ; a<nodeList.size() ; ++a )
  {
    const localIndex nodeIndex = nodeList[a];
    R1Tensor tmp = (*nodeManager.m_refposition)[nodeIndex];
    tmp += (*nodeManager.m_displacement)[nodeIndex];
    efs.xmin.SetMin(tmp);
    efs.xmax.SetMax(tmp);
    efs.xs.push_back(tmp);
    //add the velocity pushed forward by 1/2 step to the nodal velocity list
    tmp = (*nodeManager.m_acceleration)[nodeIndex];
    tmp *= 0.5 * dt;
    tmp += (*nodeManager.m_velocity)[nodeIndex];
    efs.dxs.push_back(tmp);
  }
}















/*
 * @author walsh24
 *
 * Arrange the nodes in counter-clockwise order for all of the faces in the face manager
 */
void FaceManagerT::SortAllFaceNodes(const NodeManagerT& nodeManager){
  for(localIndex f =0; f < m_numFaces; ++f)
  {
    SortFaceNodes(nodeManager, f );
  }
}

/*
 * @author walsh24
 *
 * Arranges face nodes in counter-clockwise order.
 * Sets approximate face normal = element[0].center - face center
 * Then sorts the node order counter-clockwise around the element center
*/
void FaceManagerT::SortFaceNodes(const NodeManagerT& nodeManager, const localIndex faceIndex )
{

  Array1dT<localIndex>& faceNodes = m_toNodesRelation[faceIndex];
  const localIndex firstNodeIndex = faceNodes[0];
  const unsigned int numFaceNodes = faceNodes.size();

  // get face center (average vertex location) and store node coordinates
  Array1dT<R1Tensor> faceCoords(numFaceNodes);
  R1Tensor fc;
  for( unsigned int n =0; n < numFaceNodes; ++n){
  	localIndex nd = faceNodes[n];
  	faceCoords[n] = (*nodeManager.m_refposition)[nd] + (*nodeManager.m_displacement)[nd];
    fc += faceCoords[n];
  }
  fc /= realT(numFaceNodes);

  // find center of element 0
  R1Tensor ec;
  ElementIdPair eid = m_toElementsRelation[faceIndex][0];
  ec = eid.first->GetElementCenter(eid.second,nodeManager);


  R1Tensor ex, ey, ez;
  // Approximate face normal direction (unscaled)
  if (numFaceNodes == 2)  //2D only.
  {
    ex = (*nodeManager.m_refposition)[faceNodes[1]] + (*nodeManager.m_displacement)[faceNodes[1]];
    ex -= (*nodeManager.m_refposition)[faceNodes[0]] + (*nodeManager.m_displacement)[faceNodes[0]];
    ey = ec;
    ey -= fc;

  }
  else if (eid.first->m_numFacesPerElement == 1)
  {
    //  The original/default algorithm does not work for shell elements where the face is the element itself
    //  In the new algorithm, we construct ez based on two vectors ex and ey in the plane.
    //Because ex and ey are generally not perpendicular to each other, we have to replace ey after we get ez.

      ex = faceCoords[0];
      ex -= fc;
      ex /= ex.L2_Norm();

      ey = faceCoords[1];
      ey -= fc;
      ey /= ey.L2_Norm();
      ez.Cross(ex, ey);
      if (ez.L2_Norm() < 0.01)  // Node 0, 1, and face center are roughly along one straight line.  We use another node to construct the vectors.
      {
        ey = faceCoords[2]-fc;
        ey /= ey.L2_Norm();
      }

      ez.Cross(ex, ey); ez /= ez.L2_Norm();
      ey.Cross(ez,ex); ey /= ey.L2_Norm();
  }
  else
  {
     ez = fc;
     ez -=ec;

    /// Approximate in-plane axis
     ex = faceCoords[0]-fc; ex /= ex.L2_Norm();
     ey.Cross(ez,ex); ey /= ey.L2_Norm();
  }


  if (numFaceNodes > 2)
  {
    /// Sort nodes counterclockwise around face center
    Array1dT< std::pair<realT,int> > thetaOrder(numFaceNodes);
    for( unsigned int n =0; n < numFaceNodes; ++n){
    R1Tensor v = faceCoords[n]-fc;
    thetaOrder[n] = std::pair<realT,int>(atan2(v*ey,v*ex),faceNodes[n]);
    }

    sort(thetaOrder.begin(), thetaOrder.end());

    // Reorder nodes on face
    for( unsigned int n =0; n < numFaceNodes; ++n)
    {
      faceNodes[n] = thetaOrder[n].second;
    }

    lArray1d tempFaceNodes(numFaceNodes);
    localIndex firstIndexIndex = 0;
    for( unsigned int n =0; n < numFaceNodes; ++n)
    {
      tempFaceNodes[n] = thetaOrder[n].second;
      if( tempFaceNodes[n] == firstNodeIndex )
      {
        firstIndexIndex = n;
      }
    }
    for( unsigned int n=0; n < numFaceNodes; ++n)
    {
      const localIndex index = firstIndexIndex+n < numFaceNodes ? firstIndexIndex+n : firstIndexIndex+n-numFaceNodes;
      faceNodes[n] = tempFaceNodes[index];
    }

  }
  else
  {
    ez.Cross(ex, ey);
    // The element should be on the right hand side of the vector from node 0 to node 1.
    // This ensure that the normal vector of an external face points to outside the element.
    if (ez[2] > 0)
    {
      localIndex itemp = faceNodes[0];
      faceNodes[0] = faceNodes[1];
      faceNodes[1] = itemp;

    }
  }

}

void FaceManagerT::EdgeVectors( const NodeManagerT& nodeManager, const localIndex faceIndex, Array1dT<R1Tensor>& edgeVectors )
{
  edgeVectors.resize(4);
  for( localIndex a=1 ; a<m_toNodesRelation[faceIndex].size() ; ++a )
  {
    const localIndex nodeIndex0 = m_toNodesRelation[faceIndex][a-1];
    const localIndex nodeIndex1 = m_toNodesRelation[faceIndex][a];

    edgeVectors[a-1]  = (*nodeManager.m_refposition)[nodeIndex1];
    edgeVectors[a-1] += (*nodeManager.m_displacement)[nodeIndex1];

    edgeVectors[a-1] -= (*nodeManager.m_refposition)[nodeIndex0];
    edgeVectors[a-1] -= (*nodeManager.m_displacement)[nodeIndex0];
  }
}

// Calculate a unit vector in the face and normal to the edge given
void FaceManagerT::InFaceVectorNormalToEdge(const NodeManagerT& nodeManager,
                                            const EdgeManagerT& edgeManager,
                                            localIndex iFace,
                                            localIndex iEdge,
                                            R1Tensor& v)
{
  R1Tensor faceCenter, xEdge[2], prj;
  FaceCenter(nodeManager, iFace, faceCenter);
  xEdge[0] = (*nodeManager.m_refposition)[edgeManager.m_toNodesRelation(iEdge,0)];
  xEdge[1] = (*nodeManager.m_refposition)[edgeManager.m_toNodesRelation(iEdge,1)];

  realT ndist, udist, segmentLength;
  GeometryUtilities::ProjectPointToLineSegment( xEdge[0], xEdge[1], faceCenter, ndist, udist, segmentLength, prj);
  v = faceCenter;
  v -=prj;
  v.Normalize();
}


/// Calculates the vector from node 0 to node 1
void FaceManagerT::FaceVector2D(const NodeManagerT& nodeManager, localIndex iFace, localIndex iNd, R1Tensor& v)
{
  const localIndex& node1 = m_toNodesRelation[iFace][1];
  const localIndex& node0 = m_toNodesRelation[iFace][0];
  if (iNd == node0)
  {
  v =  (*nodeManager.m_refposition)[node1];
  //v += (*nodeManager.m_displacement)[node1];
  v -= (*nodeManager.m_refposition)[node0];
  //v -= (*nodeManager.m_displacement)[node0];
  }
  else
  {
    v =  (*nodeManager.m_refposition)[node0];
    v -= (*nodeManager.m_refposition)[node1];
  }

}

int FaceManagerT::ParentNormalDirection( const int lfn )
{
  int dir;

  if( nsdof == 2 )
  {
    if( lfn == 0 )
      dir = 2;
    else if( lfn == 1 )
      dir = -1;
    else if( lfn == 2 )
      dir = -2;
    else if( lfn == 3 )
      dir = 1;
    else
      dir = INT_MIN;
  }
  else
  {
    if( lfn == 0 )
      dir = 1;
    else if( lfn == 1 )
      dir = -1;
    else if( lfn == 2 )
      dir = 2;
    else if( lfn == 3 )
      dir = -2;
    else if( lfn == 4 )
      dir = 3;
    else if( lfn == 5 )
      dir = -3;
    else
      dir = INT_MIN;
  }
  return dir;
}

#if 0
void FaceManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                                     Array1dT<gArray1d>& objectToCompositionObject )
{

  compositionObjectManager.CheckObjectType( ObjectDataStructureBaseT::NodeManager );

  iArray1d& isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();

  int kf_ext = 0;
  for( localIndex kf=0 ; kf<m_numFaces ; ++kf )
  {

    if( isDomainBoundary(kf) != 0 )
    {

      for( lArray1d::const_iterator lnode=m_toNodesRelation(kf).begin() ; lnode!=m_toNodesRelation(kf).end() ; ++lnode )
      {
        globalIndex gnode = compositionObjectManager.m_localToGlobalMap(*lnode);
        objectToCompositionObject(kf_ext).push_back( gnode );
      }
      std::sort( objectToCompositionObject(kf_ext).begin(), objectToCompositionObject(kf_ext).end() );
      ++kf_ext;
    }
  }
}
#else
void FaceManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                                     Array1dT<gArray1d>& objectToCompositionObject )
{

  compositionObjectManager.CheckObjectType( ObjectDataStructureBaseT::NodeManager );

  iArray1d& isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();


  objectToCompositionObject.clear();
  for( localIndex kf=0 ; kf<m_numFaces ; ++kf )
  {

    if( isDomainBoundary(kf) != 0 )
    {
      gArray1d temp;

      for( lArray1d::const_iterator lnode=m_toNodesRelation(kf).begin() ; lnode!=m_toNodesRelation(kf).end() ; ++lnode )
      {
        const globalIndex gnode = compositionObjectManager.m_localToGlobalMap(*lnode);
        temp.push_back( gnode );
      }
      std::sort( temp.begin(), temp.end() );
      temp.insert( temp.begin(), this->m_localToGlobalMap[kf] );
      objectToCompositionObject.push_back(temp);
    }
  }
}
#endif


template< typename T_indices >
unsigned int FaceManagerT::PackFaces( const T_indices& sendfaces,
                                      const NodeManagerT& nodeManager,
                                      const EdgeManagerT* const edgeManager,
                                      bufvector& buffer,
                                      const bool packConnectivityToGlobal,
                                      const bool packFields,
                                      const bool packMaps,
                                      const bool packSets ) const
{

  unsigned int sizeOfPacked = 0;

  const std::string label = "FaceData";
  sizeOfPacked += buffer.Pack(label);

  // pack data in the base
  sizeOfPacked += ObjectDataStructureBaseT::PackBaseObjectData( buffer, sendfaces, packFields, packMaps, packSets, packConnectivityToGlobal );





  // pack the face specific data
  for( typename T_indices::const_iterator faceIndex=sendfaces.begin() ; faceIndex!=sendfaces.end() ; ++faceIndex )
  {
    const localIndex* const nodelist = m_toNodesRelation[*faceIndex].data();

    lArray1d::size_type numnodes = m_toNodesRelation[*faceIndex].size();
    sizeOfPacked += buffer.Pack(numnodes);

    for( lArray1d::size_type a=0 ; a<numnodes ; ++a )
    {
      globalIndex gnode = GLOBALINDEX_MAX;
      if( packConnectivityToGlobal )
      {
        gnode = nodeManager.m_localToGlobalMap(nodelist[a]) ;
      }
      else
      {
        gnode = nodelist[a] ;
      }
      sizeOfPacked += buffer.Pack(gnode);
    }


    if( edgeManager != NULL )
    {
      lArray1d::size_type numedges = m_toEdgesRelation[*faceIndex].size();
      sizeOfPacked += buffer.Pack(numedges);
      for( lArray1d::const_iterator ke=m_toEdgesRelation[*faceIndex].begin() ;
           ke!=m_toEdgesRelation[*faceIndex].end() ; ++ke )
      {
        globalIndex gedge = GLOBALINDEX_MAX;
        if( packConnectivityToGlobal )
        {
          gedge = edgeManager->m_localToGlobalMap( *ke );
        }
        else
        {
          gedge = *ke;
        }

        sizeOfPacked += buffer.Pack(gedge);
      }
    }
  }


  return sizeOfPacked;
}
template unsigned int FaceManagerT::PackFaces( const lSet& sendfaces,const NodeManagerT& nodeManager,const EdgeManagerT* const edgeManager,bufvector& buffer, const bool, const bool, const bool, const bool ) const;
template unsigned int FaceManagerT::PackFaces( const lArray1d& sendfaces,const NodeManagerT& nodeManager,const EdgeManagerT* const edgeManager,bufvector& buffer, const bool, const bool, const bool, const bool ) const;




unsigned int FaceManagerT::UnpackFaces( const char*& buffer,
                                        const NodeManagerT& nodeManager,
                                        const EdgeManagerT* const edgeManager,
                                        lArray1d& faceReceiveLocalIndices,
                                        const bool unpackConnectivityToLocal,
                                        const bool unpackFields,
                                        const bool unpackMaps,
                                        const bool unpackSets  )
{
  unsigned int sizeOfUnpacked = 0;



  const std::string label = "FaceData";
  std::string temp;

  sizeOfUnpacked += bufvector::Unpack( buffer, temp );
  if( label.compare(temp)!=0 )
  {
    throw GPException("FaceManagerT::UnpackFaces: buffer location incorrect\n");
  }


  // unpack data from base object
  sizeOfUnpacked += ObjectDataStructureBaseT::UnpackBaseObjectData( buffer, faceReceiveLocalIndices, unpackFields, unpackMaps, unpackSets, unpackConnectivityToLocal );

  const lArray1d::size_type numReceivedFaces = faceReceiveLocalIndices.size();


  // unpack face specific data
  for( lArray1d::size_type kf=0 ; kf<numReceivedFaces ; ++kf )
  {
    // unpack nodes that make up faces
    lArray1d::size_type numnodes;
    sizeOfUnpacked += bufvector::Unpack( buffer, numnodes );

    lArray1d& faceToNodeMap = m_toNodesRelation(faceReceiveLocalIndices(kf));
    faceToNodeMap.resize(numnodes);
    for( lArray1d::size_type a=0 ; a<numnodes ; ++a )
    {
      globalIndex gnode;
      sizeOfUnpacked += bufvector::Unpack( buffer, gnode );

      if( unpackConnectivityToLocal )
      {
        const localIndex lnode = stlMapLookup( nodeManager.m_globalToLocalMap, gnode );
        faceToNodeMap[a] = lnode;
      }
      else
      {
        faceToNodeMap[a] = gnode;
      }

    }

    // unpack edges that make up faces
    if( edgeManager != NULL )
    {
      lArray1d::size_type numedges;
      sizeOfUnpacked += bufvector::Unpack( buffer, numedges );

      lArray1d& faceToEdgeMap = m_toEdgesRelation(faceReceiveLocalIndices(kf));
      faceToEdgeMap.resize( numedges );
      for( lArray1d::size_type a=0 ; a<numedges ; ++a )
      {
        globalIndex gedge;
        sizeOfUnpacked += bufvector::Unpack( buffer, gedge );
        if( unpackConnectivityToLocal )
        {
          const localIndex ledge = stlMapLookup( edgeManager->m_globalToLocalMap, gedge );
          faceToEdgeMap[a] = ledge;
        }
        else
        {
          faceToEdgeMap[a] = gedge;
        }
      }
    }
  }


  return sizeOfUnpacked;

}

void FaceManagerT::ConnectivityFromGlobalToLocal( const lSet& indices,
                                                  const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                                  const std::map<globalIndex,localIndex>& edgeGlobalToLocal )
{
  for( lSet::const_iterator kf=indices.begin() ; kf!=indices.end() ; ++kf )
  {
    for( lArray1d::iterator iter_nodeIndex=m_toNodesRelation[*kf].begin() ; iter_nodeIndex!=m_toNodesRelation[*kf].end() ; ++iter_nodeIndex )
    {
      const globalIndex gnode = *iter_nodeIndex;
      const localIndex lnode = stlMapLookup( nodeGlobalToLocal, gnode );
      *iter_nodeIndex = lnode;

    }

    for( lArray1d::iterator iter_edgeIndex=m_toEdgesRelation[*kf].begin() ; iter_edgeIndex!=m_toEdgesRelation[*kf].end() ; ++iter_edgeIndex )
    {
      const globalIndex gnode = *iter_edgeIndex;
      const localIndex lnode = stlMapLookup( edgeGlobalToLocal, gnode );
      *iter_edgeIndex = lnode;

    }
  }
}




void FaceManagerT::AddToFaceToElementMap( const ElementManagerT& elementManager,
                                          const std::map<std::string,lArray1d>& newElementIndices )
{
  // because the faceToElementMap is an odd creature, it is not managed by ObjectDataStructureBaseT...so we must
  // resize.
  m_toElementsRelation.resize(m_numFaces);

  // iterate over all element regions
  for( std::map<ElementManagerT::RegKeyType, ElementRegionT>::const_iterator ielemRegion=elementManager.m_ElementRegions.begin() ;
       ielemRegion!=elementManager.m_ElementRegions.end() ; ++ielemRegion )
  {
    // the element region is the mapped value of the iterator
    const std::string& regionName = ielemRegion->first;
    const ElementRegionT& elemRegion = ielemRegion->second;

    std::map<std::string,lArray1d>::const_iterator i=newElementIndices.find( regionName );
    if( i!= newElementIndices.end() )
    {
      const lArray1d& elementIndices = i->second;

      // loop over all elements in list
      for( size_t k=0 ; k<elementIndices.size(); ++k)
      {
        const localIndex elemIndex = elementIndices(k);

        // get the elementToNodeMap for element k
        const localIndex* const elementToFaceMap = elemRegion.m_toFacesRelation[elemIndex];

        // loop over all nodes in elementToNodeMap
        for( int a=0 ; a<elemRegion.m_numFacesPerElement ; ++a )
        {
          // get local index of the node from elementToNodeMap
          const localIndex localFaceIndex = elementToFaceMap[a];

          // now set the NodeToElementMap as a combination of element region and element index.
          m_toElementsRelation(localFaceIndex).push_back( std::make_pair(const_cast<ElementRegionT*>(&elemRegion),elemIndex));
        }
      }
    }
  }
}







void FaceManagerT::PreSeparateFaces( const std::string& setname, const int setState )
{

  // preset separation states for specific facesets
  const lSet& sepset=this->GetSet(setname);
  iArray1d& ruptureState = this->GetFieldData<int>("ruptureState");
  rArray1d& separationCoeff = this->GetFieldData<realT>("separationCoeff");


  for( lSet::const_iterator faceIndex=sepset.begin() ; faceIndex!=sepset.end() ; ++faceIndex )
  {
    if( m_toElementsRelation[*faceIndex].size() > 1 )
    {
      ruptureState[*faceIndex] = setState;
      separationCoeff[*faceIndex] = 0.5;
    }
  }

}

// PFu note: Looks like elementManager is passed in but not used.
void FaceManagerT::UpdateRuptureStates( const ElementManagerT& elementManager ,
                                        const NodeManagerT& nodeManager,
                                        const std::string& separableSet,
                                        const realT failval  )
{

//  Array1dT<R1Tensor>& maxTraction = this->GetFieldData<R1Tensor>("maxTraction");
  rArray1d* failStressPointer = this->GetFieldDataPointer<realT>("faceFailStress");

  std::map< std::string, lSet >::const_iterator setMap = this->m_Sets.find( separableSet );

  if( setMap==this->m_Sets.end() )
  {
    for( localIndex kf=0 ; kf<m_numFaces ; ++kf )
    {
      if (failStressPointer == NULL)
      {
        UpdateRuptureState( elementManager,nodeManager, kf, failval );
      }
      else
      {
        UpdateRuptureState( elementManager,nodeManager, kf, (*failStressPointer)[kf] );
      }
    }
  }
  else
  {
    for( lSet::const_iterator kf=setMap->second.begin() ; kf!=setMap->second.end() ; ++kf )
    {
      if (failStressPointer == NULL)
      {
        UpdateRuptureState( elementManager,nodeManager, *kf, failval );
      }
      else
      {
        UpdateRuptureState( elementManager,nodeManager, *kf, (*failStressPointer)[*kf] );
      }
    }

  }

}

void FaceManagerT::UpdateRuptureState( const ElementManagerT& elementManager,
                                       const NodeManagerT& nodeManager,
                                       const localIndex kf,
                                       const realT failval )
{
  R1Tensor fc;
  R1Tensor fn, ft0, ft1;



  R2SymTensor stress0;
  R2SymTensor stress1;

  R1Tensor t0, t1;
  R1Tensor temp;

  iArray1d& ruptureState = this->GetFieldData<int>("ruptureState");
  rArray1d* stressNOnFace = this->GetFieldDataPointer<realT>("stressNOnFace");
  Array1dT<R1Tensor>* stressTOnFace = this->GetFieldDataPointer<R1Tensor>("stressTOnFace");
  rArray1d* faceStrengthRandomFactor =  this->GetFieldDataPointer<realT>("faceStrengthRandomFactor");

  if( m_toElementsRelation[kf].size() > 1 )
    //Warning (fu): With the introduction of flow face regions, this is not safe anymore.  I need to modify then when I get a chance.

  {
//      unsigned long long Size = m_FaceToElementMap[kf].size();

    ElementRegionT& er0 = *(m_toElementsRelation[kf][0].first);
    ElementRegionT& er1 = *(m_toElementsRelation[kf][1].first);
    const localIndex elemIndex0 = m_toElementsRelation[kf][0].second;
    const localIndex elemIndex1 = m_toElementsRelation[kf][1].second;

    stress0 = 0;
    realT pressure0;
    {
      er0.m_mat->MeanPressureDevStress(elemIndex0,pressure0, stress0);
      stress0.PlusIdentity(pressure0);
    }

    stress1 = 0;
    realT pressure1;
    {
      er1.m_mat->MeanPressureDevStress(elemIndex1,pressure1, stress1);
      stress1.PlusIdentity(pressure1);
    }

    // normal from away from element 0, into element 1
    FaceCenterAndNormal( nodeManager, kf, fc, fn );
    FaceTangential(nodeManager, kf, ft0, ft1);


    t0.AijBj(stress0,fn);
    t1.AijBj(stress1,fn);


    realT t0n = Dot(t0,fn);
    realT t1n = Dot(t1,fn);
    R1Tensor t0t = Dot(t0,ft0) * ft0 + Dot(t0, ft1) * ft1;
    R1Tensor t1t = Dot(t1,ft0) * ft0 + Dot(t1, ft1) * ft1;
    if (stressNOnFace != NULL)
    {
      (*stressNOnFace)[kf] = 0.5*(t0n+t1n);
    }
    if (stressTOnFace != NULL)
    {
      (*stressTOnFace)[kf] = 0.5*(t0t+t1t);
    }


    if (faceStrengthRandomFactor != NULL)
    {
      if( 0.5*(t0n+t1n) > failval * (*faceStrengthRandomFactor)[kf] && ruptureState[kf] == 0) // && fabs(fn[2]) > 0.9)  //HACK
      {
        ruptureState[kf] = 1;
      }
      else if (0.5*(t0n+t1n) <= failval * (*faceStrengthRandomFactor)[kf] && ruptureState[kf]==1 )
      {
        ruptureState[kf] = 0;
      }

    }
    else
    {
      if( 0.5*(t0n+t1n) > failval && ruptureState[kf] == 0) // && fabs(fn[2]) > 0.9)  //HACK
      {
        ruptureState[kf] = 1;
      }
      else if (0.5*(t0n+t1n) <= failval && ruptureState[kf]==1 )
      {
        ruptureState[kf] = 0;
      }
    }
  }
  else if( m_toElementsRelation[kf].size() == 1 && stressNOnFace != NULL)
  {
    //We calculate stress for external faces only for visualization purpose.
    ElementRegionT& er0 = *(m_toElementsRelation[kf][0].first);
    const localIndex elemIndex0 = m_toElementsRelation[kf][0].second;

    stress0 = 0;
    realT pressure0;
    {
      er0.m_mat->MeanPressureDevStress(elemIndex0,pressure0, stress0);
      stress0.PlusIdentity(pressure0);
    }

    // normal from away from element 0, into element 1
    FaceCenterAndNormal( nodeManager, kf, fc, fn );
    FaceTangential(nodeManager, kf, ft0, ft1);


    t0.AijBj(stress0,fn);

    realT t0n = Dot(t0,fn);
    R1Tensor t0t = Dot(t0,ft0) * ft0 + Dot(t0, ft1) * ft1;


    if (stressNOnFace != NULL) (*stressNOnFace)[kf] = t0n;
    if (stressTOnFace != NULL) (*stressTOnFace)[kf] = t0t;

  }
}

void FaceManagerT::CalculateStressOnFace( const ElementManagerT& elementManager,
                                          const NodeManagerT& nodeManager,
                                          const localIndex kf )
{

  R1Tensor fc;
  R1Tensor fn, ft0, ft1;

  R2SymTensor stress0;
  R2SymTensor stress1;

  R1Tensor t0, t1;
  R1Tensor temp;

  rArray1d* stressNOnFace = this->GetFieldDataPointer<realT>("stressNOnFace");
  Array1dT<R1Tensor>* stressTOnFace = this->GetFieldDataPointer<R1Tensor>("stressTOnFace");

  if( m_toElementsRelation[kf].size() > 1 )
  {
//      unsigned long long Size = m_FaceToElementMap[kf].size();

    ElementRegionT& er0 = *(m_toElementsRelation[kf][0].first);
    ElementRegionT& er1 = *(m_toElementsRelation[kf][1].first);
    const localIndex elemIndex0 = m_toElementsRelation[kf][0].second;
    const localIndex elemIndex1 = m_toElementsRelation[kf][1].second;

    stress0 = 0;
    realT pressure0;
    {
      er0.m_mat->MeanPressureDevStress(elemIndex0,pressure0, stress0);
      stress0.PlusIdentity(pressure0);
    }

    stress1 = 0;
    realT pressure1;
    {
      er1.m_mat->MeanPressureDevStress(elemIndex1,pressure1, stress1);
      stress1.PlusIdentity(pressure1);
    }

    // normal from away from element 0, into element 1
    FaceCenterAndNormal( nodeManager, kf, fc, fn );
    FaceTangential(nodeManager, kf, ft0, ft1);


    t0.AijBj(stress0,fn);
    t1.AijBj(stress1,fn);


    realT t0n = Dot(t0,fn);
    realT t1n = Dot(t1,fn);
    R1Tensor t0t = Dot(t0,ft0) * ft0 + Dot(t0, ft1) * ft1;
    R1Tensor t1t = Dot(t1,ft0) * ft0 + Dot(t1, ft1) * ft1;
    if (stressNOnFace != NULL)
    {
      (*stressNOnFace)[kf] = 0.5*(t0n+t1n);
    }
    if (stressTOnFace != NULL)
    {
      (*stressTOnFace)[kf] = 0.5*(t0t+t1t);
    }
  }

  else if( m_toElementsRelation[kf].size() == 1 && stressNOnFace != NULL)
  {
    //We calculate stress for external faces only for visualization purpose.
    ElementRegionT& er0 = *(m_toElementsRelation[kf][0].first);
    const localIndex elemIndex0 = m_toElementsRelation[kf][0].second;

    stress0 = 0;
    realT pressure0;
    {
      er0.m_mat->MeanPressureDevStress(elemIndex0,pressure0, stress0);
      stress0.PlusIdentity(pressure0);
    }

    // normal from away from element 0, into element 1
    FaceCenterAndNormal( nodeManager, kf, fc, fn );
    FaceTangential(nodeManager, kf, ft0, ft1);


    t0.AijBj(stress0,fn);

    realT t0n = Dot(t0,fn);
    R1Tensor t0t = Dot(t0,ft0) * ft0 + Dot(t0, ft1) * ft1;


    if (stressNOnFace != NULL) (*stressNOnFace)[kf] = t0n;
    if (stressTOnFace != NULL) (*stressTOnFace)[kf] = t0t;

  }

}



#if 0
void FaceManagerT::SplitFace( const localIndex indexToSplit,
                              const localIndex parentNodeIndex,
                              const localIndex childNodeIndex[2],
                              const lSet& splitEdges,
                              const Array1dT<lArray1d>& childEdges,
                              Array1dT<lSet>& nodesToFaces,
                              Array1dT<lSet>& edgesToFaces )
{
  localIndex newFaceIndex[2] ;
  const bool didSplit = SplitObject( indexToSplit, newFaceIndex );


  // if the face was not already split, then we need to allocate and copy the face-> maps.
  if( didSplit )
  {
    // have to do this manually!!!
    m_toElementsRelation.resize( this->m_numFaces );

    // resize the faceToEdge map for each new face to be the same size as the parent face
    m_toEdgesRelation[newFaceIndex[0]] = m_toEdgesRelation[indexToSplit];
    m_toEdgesRelation[newFaceIndex[1]] = m_toEdgesRelation[indexToSplit];

    m_toNodesRelation[newFaceIndex[0]] = m_toNodesRelation[indexToSplit];
    m_toNodesRelation[newFaceIndex[1]] = m_toNodesRelation[indexToSplit];

  }




  // modify the face/node maps

  // loop over nodes on the parent face
  for( Array1dT<lArray1d>::size_type a=0 ; a<m_toNodesRelation[indexToSplit].size() ; ++a )
  {
    const localIndex& nodeIndex = m_toNodesRelation[indexToSplit][a];

    // if the node is the one that was just split, then we have to modify the relation
    if( nodeIndex == parentNodeIndex )
    {
      // set the child nodes to the child edges faceToNodeMap
      m_toNodesRelation[newFaceIndex[0]][a] = childNodeIndex[0];
      m_toNodesRelation[newFaceIndex[1]][a] = childNodeIndex[1];

      // add the child edges to the child nodes nodeToFaceMap
      nodesToFaces[childNodeIndex[0]].insert( newFaceIndex[0] );
      nodesToFaces[childNodeIndex[1]].insert( newFaceIndex[1] );

      // erase the new faces from the parent nodes. This will be necessary if the face is already split.
      nodesToFaces[parentNodeIndex].erase( newFaceIndex[0] );
      nodesToFaces[parentNodeIndex].erase( newFaceIndex[1] );
    }
    else
    {
      if( didSplit )
      {
        // add the child faces to the child nodes nodeToFaceMap
        nodesToFaces[nodeIndex].insert( newFaceIndex[0] );
        nodesToFaces[nodeIndex].insert( newFaceIndex[1] );
      }
    }
  }







  // Fill the face to edge map for the new faces, and modify the edgeToFaces map to include
  // the new faces.


  // now loop over each edge in the parent face
  for( lArray1d::size_type ke=0 ; ke<m_toEdgesRelation[indexToSplit].size() ; ++ke )
  {
    const localIndex& edgeIndex = m_toEdgesRelation[indexToSplit][ke];

    // if the edge was just split, then we will have to modify both relations
    if( splitEdges.count(edgeIndex) > 0 )
    {
      // add the child edges to the faceToEdge map
      m_toEdgesRelation[newFaceIndex[0]][ke] = childEdges[edgeIndex][0];
      m_toEdgesRelation[newFaceIndex[1]][ke] = childEdges[edgeIndex][1];

      std::cout<<"    m_FaceToEdgeMap["<<newFaceIndex[0]<<"]["<<ke<<"] = "<<childEdges[edgeIndex][0]<<std::endl;
      std::cout<<"    m_FaceToEdgeMap["<<newFaceIndex[1]<<"]["<<ke<<"] = "<<childEdges[edgeIndex][1]<<std::endl;

      // add the new faces to the edgetoFace map entries for the child edges
      edgesToFaces[childEdges[edgeIndex][0]].insert(newFaceIndex[0]);
      edgesToFaces[childEdges[edgeIndex][1]].insert(newFaceIndex[1]);

//      std::cout<<"    edgesToFaces["<<childEdges[edgeIndex][0]<<"].insert("<<newFaceIndex[0]<<")"<<std::endl;
//      std::cout<<"    edgesToFaces["<<childEdges[edgeIndex][1]<<"].insert("<<newFaceIndex[0]<<")"<<std::endl;

      edgesToFaces[edgeIndex].erase(newFaceIndex[0]);
      edgesToFaces[edgeIndex].erase(newFaceIndex[1]);

    }
    else
    {
      if( didSplit )
      {
        edgesToFaces[edgeIndex].insert( newFaceIndex[0] );
        edgesToFaces[edgeIndex].insert( newFaceIndex[1] );
      }
    }

  }
}
#endif


void FaceManagerT::ModifyToFaceMapsFromSplit( const lSet& newFaces,
                                              const lSet& modifiedFaces,
                                              NodeManagerT& nodeManager,
                                              EdgeManagerT& edgeManager )
{

  lSet allFaces;
  allFaces.insert( newFaces.begin(), newFaces.end() );
  allFaces.insert( modifiedFaces.begin(), modifiedFaces.end() );

  // update external face set
  m_externalFaces.insert( newFaces.begin(), newFaces.end() );
  m_externalFaces.insert( modifiedFaces.begin(), modifiedFaces.end() );


  for( lSet::const_iterator faceIndex=allFaces.begin() ; faceIndex!=allFaces.end() ; ++faceIndex )
  {
    // wipe the face from all toFace
    for( lArray1d::const_iterator iter_node=this->m_toNodesRelation[*faceIndex].begin() ; iter_node!=this->m_toNodesRelation[*faceIndex].end() ; ++iter_node )
    {
      localIndex nodeIndex = *iter_node;
      localIndex parentNodeIndex = nodeManager.m_parentIndex[nodeIndex];

      while( parentNodeIndex != LOCALINDEX_MAX )
      {
        nodeManager.m_nodeToFaceMap[parentNodeIndex].erase(*faceIndex);

        nodeIndex = parentNodeIndex;
        parentNodeIndex = nodeManager.m_parentIndex[nodeIndex];
      }
    }

    for( lArray1d::const_iterator iter_edge=this->m_toEdgesRelation[*faceIndex].begin() ; iter_edge!=this->m_toEdgesRelation[*faceIndex].end() ; ++iter_edge )
    {
      localIndex edgeIndex = *iter_edge;
      localIndex parentEdgeIndex = edgeManager.m_parentIndex[edgeIndex];

      while( parentEdgeIndex != LOCALINDEX_MAX )
      {
        edgeManager.m_toFacesRelation[parentEdgeIndex].erase(*faceIndex);

        edgeIndex = parentEdgeIndex;
        parentEdgeIndex = edgeManager.m_parentIndex[edgeIndex];
      }
    }

    // add the face from toFaces
    for( lArray1d::const_iterator iter_node=this->m_toNodesRelation[*faceIndex].begin() ; iter_node!=this->m_toNodesRelation[*faceIndex].end() ; ++iter_node )
    {
      localIndex nodeIndex = *iter_node;
      nodeManager.m_nodeToFaceMap[nodeIndex].insert(*faceIndex);
    }

    for( lArray1d::const_iterator iter_edge=this->m_toEdgesRelation[*faceIndex].begin() ; iter_edge!=this->m_toEdgesRelation[*faceIndex].end() ; ++iter_edge )
    {
      localIndex edgeIndex = *iter_edge;
      edgeManager.m_toFacesRelation[edgeIndex].insert(*faceIndex);
    }
  }
}


R1Tensor FaceManagerT::CalculateGapVector( const NodeManagerT& nodeManager, const localIndex faceIndex ) const
{
  R1Tensor fc0, fc1;
  R1Tensor v;

  if( m_childIndices[faceIndex].size() == 2 )
  {
    FaceCenter( nodeManager, this->m_childIndices[faceIndex][0] , fc0 );
    FaceCenter( nodeManager, this->m_childIndices[faceIndex][1] , fc1 );
    v = fc1;
    v -= fc0 ;
  }
  else if( m_childIndices[faceIndex].size() == 1 )
  {
    FaceCenter( nodeManager, faceIndex , fc0 );
    FaceCenter( nodeManager, this->m_childIndices[faceIndex][0] , fc1 );
    v = fc1;
    v -= fc0 ;
  }

  else
  {
    v = 0.0;
  }


  return v;
}

R1Tensor FaceManagerT::CalculateGapRateVector( const NodeManagerT& nodeManager, const localIndex faceIndex ) const
{

  localIndex faceID[2];
  if (m_childIndices[faceIndex].size() == 2)
  {
    faceID[0] = m_childIndices[faceIndex][0];
    faceID[1] = m_childIndices[faceIndex][1];
  }
  else
  {
    faceID[0] = faceIndex;
    faceID[1] = m_childIndices[faceIndex][0];
  }

  R1Tensor fv[2];

  for (localIndex i=0; i<2; ++i)
  {
    fv[i] = 0;
    for (localIndex j=0; j<m_toNodesRelation[faceID[i]].size(); ++j)
    {
      fv[i] += (*nodeManager.m_velocity)[m_toNodesRelation[faceID[i]][j]];
    }
    fv[i] /= m_toNodesRelation[faceID[i]].size();
  }

  R1Tensor v;
  v = fv[1];
  v -= fv[0];

  return v;
}


//FIXME: This function does not seem to be complete.
void FaceManagerT::CalculateGapVectorDerivative( const NodeManagerT& nodeManager,
                                                 const localIndex faceIndex,
                                                 Array1dT< Array1dT<R1Tensor> >& gapDerivative,
                                                 Array1dT<iArray1d>& gapDerivativeIndices ) const
{
  localIndex faceIndices[2];
  if( m_childIndices[faceIndex].size() == 2 )
  {
    faceIndices[0] = m_childIndices[faceIndex][0];
    faceIndices[1] = m_childIndices[faceIndex][1];
  }
  else if( m_childIndices[faceIndex].size() == 1 )
  {
    faceIndices[0] = faceIndex;
    faceIndices[1] = m_childIndices[faceIndex][0];
  }
  else
  {
    faceIndices[0] = faceIndex;
    faceIndices[1] = faceIndex;
  }



  const localIndex numNodes = this->m_toNodesRelation[m_childIndices[faceIndex][0]].size();
  gapDerivative.resize( numNodes );
  gapDerivativeIndices.resize( numNodes );

  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    gapDerivative(a).resize(2);
    gapDerivative(a)[0] = -1.0/numNodes;
    gapDerivative(a)[1] =  1.0/numNodes;

    gapDerivativeIndices[a].resize(2);
    gapDerivativeIndices[a][0] = m_toNodesRelation[faceIndices[0]][a];
    gapDerivativeIndices[a][1] = m_toNodesRelation[faceIndices[1]][a];

  }

}


void FaceManagerT::SetApertureFromRigidWall( const NodeManagerT& nodeManager )
{

  // iterate over all boundary conditions.
  for( Array1dT<BoundaryConditionBase*>::const_iterator bcItr=nodeManager.m_bcData.begin() ; bcItr!=nodeManager.m_bcData.end() ; ++ bcItr )
  {

    // check to if the requested field has a wall boundary condition applied to it.
    WallBoundaryCondition* bc = dynamic_cast<WallBoundaryCondition*>(*bcItr);
    if( bc )
    {
      const R1Tensor& b = bc->m_position;
      const R1Tensor& n = bc->GetDirection(0.0); // fixme need time
      rArray1d& apertures = GetFieldData<realT>( "Aperture" );

      for(localIndex i =0; i < bc->m_setNames.size(); ++i){
        std::map< std::string, lSet >::iterator setMap = m_Sets.find( bc->m_setNames[i] );
        if( setMap != m_Sets.end() )
        {

          const lSet& set = setMap->second;

          for ( lSet::const_iterator k=set.begin() ; k!=set.end() ; ++k )
//        for( localIndex k=0 ; k<apertures.size() ; ++k )
          {
            R1Tensor fc;
            FaceCenter(nodeManager,*k,fc);

            fc -=b;
            apertures[*k] = Dot(fc,n);
//          if( apertures[*k] < 0.0 ) apertures[*k] = 0.0;
          }
        }
      }
    }
  }
}

void FaceManagerT::WriteSiloMesh( SiloFile& siloFile,
                                  const std::string& meshname,
                                  const NodeManagerT& nodeManager,
                                  const int cycleNum,
                                  const realT problemTime,
                                  const bool isRestart )
{

  realT* coords[3];
  rArray1d xcoords(nodeManager.m_numNodes);
  rArray1d ycoords(nodeManager.m_numNodes);
  rArray1d zcoords(nodeManager.m_numNodes);
  for (localIndex a = 0; a < nodeManager.m_numNodes; ++a)
  {
    R1Tensor nodePosition;
    nodePosition = (*nodeManager.m_refposition)[a];
    nodePosition += (*nodeManager.m_displacement)[a];

    xcoords[a] = nodePosition(0);
    ycoords[a] = nodePosition(1);
    zcoords[a] = nodePosition(2);
  }

  coords[0] = xcoords.data();
  coords[1] = ycoords.data();
  coords[2] = zcoords.data();

  const int numFaceTypes = 1;
  int numNodesPerFace = m_toNodesRelation[0].size(); // TODO assumes all faces have same number of nodes
  int dbZoneType = DB_ZONETYPE_POLYGON;
  if(numNodesPerFace == 3) {
    dbZoneType = DB_ZONETYPE_TRIANGLE;
  }else if(numNodesPerFace == 4){
    dbZoneType = DB_ZONETYPE_QUAD;
  }else if(numNodesPerFace == 2){
    dbZoneType = DB_ZONETYPE_BEAM;
  }

  Array1dT<localIndex*> faceConnectivity(numFaceTypes);
  Array1dT<globalIndex*> globalFaceNumbers(numFaceTypes);
  ivector fshapecnt(numFaceTypes);
  ivector fshapetype(numFaceTypes);
  ivector fshapesize(numFaceTypes);

  Array1dT<lArray2d> faceToNodeMap(numFaceTypes);


  for (int faceType = 0; faceType < numFaceTypes; ++faceType)
  {
    faceToNodeMap[faceType].resize2(m_numFaces, numNodesPerFace);

    for (localIndex k = 0; k < m_numFaces; ++k)
    {
      for (int a = 0; a < numNodesPerFace; ++a)
      {
        faceToNodeMap[faceType][k][a] = m_toNodesRelation[k][a];
      }
    }

    faceConnectivity[faceType] = faceToNodeMap[faceType].data();

    globalFaceNumbers[faceType] = m_localToGlobalMap.data();
    fshapecnt[faceType] = m_numFaces;
    fshapetype[faceType] = dbZoneType;
    fshapesize[faceType] = numNodesPerFace;
  }

  siloFile.WriteMeshObject(meshname, nodeManager.m_numNodes, coords,
                           m_localToGlobalMap.data(), numFaceTypes,
                           fshapecnt.data(), faceConnectivity.data(), globalFaceNumbers.data(),
                           NULL, fshapetype.data(), fshapesize.data(), cycleNum, problemTime);


//  m_feFaceManager.WriteSilo( siloFile, "FaceFields",meshname, DB_ZONECENT, cycleNum, problemTime, isRestart );


}


bool FaceManagerT::IsNodeOnFace(const localIndex faceIndex,
                                const localIndex nodeIndex )
{
  for( lArray1d::iterator i = m_toNodesRelation[faceIndex].begin() ;
       i!=m_toNodesRelation[faceIndex].end() ; ++i )
  {
    if (*i == nodeIndex) return (true);
  }

  return (false);
}

void FaceManagerT::WriteNonManagedDataMembersToSilo( SiloFile& siloFile,
                                                     const std::string& siloDirName,
                                                     const std::string& meshname,
                                                     const int centering,
                                                     const int cycleNum,
                                                     const realT problemTime,
                                                     const bool isRestart,
                                                     const std::string& multiRoot,
                                                     const std::string& regionName,
                                                     const lArray1d& mask)
{
  if( isRestart )
  {
    siloFile.DBWriteWrapper("m_matchedBoundaryFaces",m_matchedBoundaryFaces);

    std::map<std::string,unsigned int> regionNamesReverse;

    {
      char dir[10000];
      sArray1d regionNames;

      DBGetDir(siloFile.m_dbFilePtr,dir);
      DBSetDir(siloFile.m_dbFilePtr, "..");
      siloFile.DBReadWrapper( "regionNames", regionNames );
      DBSetDir(siloFile.m_dbFilePtr, dir);

      for( unsigned int i=0 ; i<regionNames.size() ; ++i )
      {
        regionNamesReverse[regionNames[i]] = i;
      }
    }


    Array1dT< iArray1d > toElementsRegion(m_numFaces);
    Array1dT< lArray1d > toElementsIndex(m_numFaces);


    for( localIndex k=0 ; k<m_numFaces ; ++k )
    {
      toElementsRegion[k].resize( m_toElementsRelation[k].size() );
      toElementsIndex[k].resize( m_toElementsRelation[k].size() );
      for( localIndex a=0 ; a<m_toElementsRelation[k].size() ; ++a )
      {
        toElementsRegion[k][a] = regionNamesReverse[m_toElementsRelation[k][a].first->m_regionName ];
        toElementsIndex[k][a] = m_toElementsRelation[k][a].second;
      }
    }

    siloFile.DBWriteWrapper("toElementsRegion",toElementsRegion);
    siloFile.DBWriteWrapper("toElementsIndex",toElementsIndex);
  }



  if( m_cohesiveZone )
  {
  TempODS temp;
  temp.resize( this->DataLengths() );

  sArray1d intVarNames;
  sArray1d realVarNames;
  sArray1d R1TensorVarNames;
  sArray1d R2TensorVarNames;
  sArray1d R2SymTensorVarNames;

  Array1dT<iArray1d*> intVars;
  Array1dT<rArray1d*> realVars;
  Array1dT<Array1dT<R1Tensor>*> R1Vars;
  Array1dT<Array1dT<R2Tensor>*> R2Vars;
  Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

  m_cohesiveZone->GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  temp.AllocateDummyFields( intVarNames, intVars );
  temp.AllocateDummyFields( realVarNames, realVars );
  temp.AllocateDummyFields( R1TensorVarNames, R1Vars );
  temp.AllocateDummyFields( R2TensorVarNames, R2Vars );
  temp.AllocateDummyFields( R2SymTensorVarNames, R2SymVars );

  m_cohesiveZone->Serialize( intVars, realVars, R1Vars, R2Vars, R2SymVars);

  std::string czname = siloDirName + "/CohesiveZones";
  temp.WriteSilo( siloFile, czname, meshname, centering, cycleNum, problemTime, isRestart, regionName, mask );


  temp.DeallocateDummyFields<int>(intVarNames);
  temp.DeallocateDummyFields<realT>(realVarNames);
  temp.DeallocateDummyFields<R1Tensor>(R1TensorVarNames);
  temp.DeallocateDummyFields<R2Tensor>(R2TensorVarNames);
  temp.DeallocateDummyFields<R2SymTensor>(R2SymTensorVarNames);

  }

}

void FaceManagerT::ReadNonManagedDataMembersFromSilo( const SiloFile& siloFile,
                                                      const std::string& siloDirName,
                                                      const std::string& meshname,
                                                      const int centering,
                                                      const int cycleNum,
                                                      const realT problemTime,
                                                      const bool isRestart,
                                                      const std::string& regionName,
                                                      const lArray1d& mask)
{

  if( isRestart )
  {
    siloFile.DBReadWrapper("m_matchedBoundaryFaces",m_matchedBoundaryFaces);

    sArray1d regionNames;

    {
      char dir[10000];
      DBGetDir(siloFile.m_dbFilePtr,dir);

      for( int i=0 ; i<10 ; ++i )
      {
        DBSetDir(siloFile.m_dbFilePtr, "..");
        siloFile.DBReadWrapper( "regionNames", regionNames );
        if( !(regionNames.empty()))
          break;
      }
      DBSetDir(siloFile.m_dbFilePtr, dir);
    }


    Array1dT< iArray1d > toElementsRegion(m_numFaces);
    Array1dT< lArray1d > toElementsIndex(m_numFaces);


    siloFile.DBReadWrapper("toElementsRegion",toElementsRegion);
    siloFile.DBReadWrapper("toElementsIndex",toElementsIndex);

    m_toElementsRelation.resize( m_numFaces );
    for( localIndex k=0 ; k<m_numFaces ; ++k )
    {
      m_toElementsRelation[k].resize( toElementsRegion[k].size() );
      for( localIndex a=0 ; a<m_toElementsRelation[k].size() ; ++a )
      {
        if( m_elemManagerHACK != NULL )
        {
          m_toElementsRelation[k][a].first = &(m_elemManagerHACK->m_ElementRegions[ regionNames[toElementsRegion[k][a] ] ]);
        }
        else
        {
          m_toElementsRelation[k][a].first = NULL;
        }

        m_toElementsRelation[k][a].second =  toElementsIndex[k][a];
      }
    }
  }

  if (m_cohesiveZone)
  {
    m_cohesiveZone->resize(DataLengths());

    TempODS temp;
    temp.resize(this->DataLengths());

    sArray1d intVarNames;
    sArray1d realVarNames;
    sArray1d R1TensorVarNames;
    sArray1d R2TensorVarNames;
    sArray1d R2SymTensorVarNames;

    Array1dT<iArray1d*> intVars;
    Array1dT<rArray1d*> realVars;
    Array1dT<Array1dT<R1Tensor>*> R1Vars;
    Array1dT<Array1dT<R2Tensor>*> R2Vars;
    Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

    m_cohesiveZone->GetVariableNames(intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames,
                                     R2SymTensorVarNames);

    temp.AllocateDummyFields(intVarNames, intVars);
    temp.AllocateDummyFields(realVarNames, realVars);
    temp.AllocateDummyFields(R1TensorVarNames, R1Vars);
    temp.AllocateDummyFields(R2TensorVarNames, R2Vars);
    temp.AllocateDummyFields(R2SymTensorVarNames, R2SymVars);

//  std::string czname = siloDirName + "/CohesiveZones";
    std::string czname = "CohesiveZones";
    const int err = temp.ReadSilo(siloFile, czname, meshname, centering, cycleNum, problemTime,
                                  isRestart, regionName, mask);
    if(err)
      return;

    m_cohesiveZone->Deserialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

    temp.DeallocateDummyFields<int>(intVarNames);
    temp.DeallocateDummyFields<realT>(realVarNames);
    temp.DeallocateDummyFields<R1Tensor>(R1TensorVarNames);
    temp.DeallocateDummyFields<R2Tensor>(R2TensorVarNames);
    temp.DeallocateDummyFields<R2SymTensor>(R2SymTensorVarNames);
  }
}
