/*
 * MeshGenerator.h
 *
 *  Created on: Nov 19, 2012
 *      Author: settgast1
 */

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include "IO/ticpp/TinyXMLParser.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "MPI_Communications/SpatialPartition.h"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */




class MeshGenerator
{
public:
  MeshGenerator();
  virtual ~MeshGenerator();

  void ReadXML( TICPP::HierarchicalDataNode& hdn );

  void GenerateElementRegions( PhysicalDomainT& domain );

  void GenerateMesh( SpatialPartition& partition,
                     PhysicalDomainT& domain );

  void GenerateNodesets( TICPP::HierarchicalDataNode& hdn,
                         NodeManagerT& nodeManager );

  void GetElemToNodesRelationInBox ( const std::string& elementType,
                                     const int index[],
                                     const int& iEle,
                                     int nodeIDInBox[],
                                     const int size);
private:


  int m_dim;
  rArray1d m_vertices[3];
  iArray1d m_nElems[3];
  rArray1d m_nElemScaling[3];

  sArray1d m_regionNames;

  realT m_min[3]; // Minimum extent of mesh dimensions
  realT m_max[3]; // Maximum extent of mesh dimensions

  //int m_numElems[3];
  iArray1d m_firstElemIndexForBlock[3];
  iArray1d m_lastElemIndexForBlock[3];



//  realT m_wExtensionMin[3];
//  realT m_wExtensionMax[3];
//  int m_nExtensionLayersMin[3];
//  int m_nExtensionLayersMax[3];
//  realT m_extendedMin[3];
//  realT m_extendedMax[3]; // This is the domain size after we apply n layers of elements which are of the same size as the core elements.  We will move these nodes to where they should be later when we finish the meshing.
  int m_numElemsTotal[3];
//  realT m_commonRatioMin[3];
//  realT m_commonRatioMax[3];



  sArray1d m_elementType;

  iArray1d m_numElePerBox;

  int m_trianglePattern  ; // In pattern 0, half nodes have 4 edges and the other half have 8; for Pattern 1, every node has 6.

  realT m_fPerturb;
  int m_randSeed;

  int m_mapToRadial;
  int meshAxis;
  float meshTheta;
  float meshPhi;
  float meshRout;
  float meshRact;

  inline globalIndex NodeGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0]*(m_numElemsTotal[1]+1)*(m_numElemsTotal[2]+1) + index[1]*(m_numElemsTotal[2]+1) + index[2];
    return rval;
  }

  inline globalIndex ElemGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0]*m_numElemsTotal[1]*m_numElemsTotal[2] + index[1]*m_numElemsTotal[2] + index[2];
    return rval;
  }

  inline R1Tensor NodePosition( const int a[3] )
  {
    R1Tensor X;

    for( int i=0 ; i<3 ; ++i )
    {

      int startingIndex = 0;
      int endingIndex = 0;
      unsigned int block = 0;
      for( block=0 ; block<m_nElems[i].size() ; ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[i][block];
        if( a[i]>=startingIndex && a[i]<=endingIndex )
        {
          break;
        }
      }
      realT min = m_vertices[i][block];
      realT max = m_vertices[i][block+1];

      X[i] = min + (max-min) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );
//      X[i] = m_min[i] + (m_max[i]-m_min[i]) * ( double(a[i]) / m_numElemsTotal[i] );
    }

    return X;
  }

  inline R1Tensor ElemCenterPosition( const int k[3] )
  {
    R1Tensor X;

    for( int i=0 ; i<3 ; ++i )
    {
      X[i] = m_min[i] + (m_max[i]-m_min[i]) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
    }

    return X;
  }

public:
  inline bool isRadial()
  {
    bool rval = (m_mapToRadial > 0);
    return rval;
  }

};

#endif /* MESHGENERATOR_H_ */
