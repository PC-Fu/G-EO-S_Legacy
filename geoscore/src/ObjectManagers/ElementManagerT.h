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
 * @file ElementManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 14, 2010
 */

#ifndef ELEMENTMANAGERT_H_
#define ELEMENTMANAGERT_H_

#include "Common.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "DataStructures/VectorFields/ElementRegionT.h"


/**
 * Class to manage the data stored at the element level.
 */
class ElementManagerT : public ObjectDataStructureBaseT
{
public:
  ElementManagerT();
  virtual ~ElementManagerT();

  void Initialize(  ){}

  void SetIsAttachedToSendingGhostNode( const NodeManagerT& nodeManager )
  {
    for( std::map< RegKeyType, ElementRegionT >::iterator i=m_ElementRegions.begin() ; i!=m_ElementRegions.end(); ++i )
    {
      ElementRegionT& elementRegion = i->second;
      elementRegion.SetIsAttachedToSendingGhostNode(nodeManager);
    }
  }


  using ObjectDataStructureBaseT::resize;
  globalIndex resize( const lvector& numElements,
                      const sArray1d& elementRegionNames,
                      const sArray1d& elementTypes );


  void HACKInitialConditions();

  void SetDomainBoundaryObjects(const ObjectDataStructureBaseT* const referenceObject );
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL) {}
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager ,
                                                         Array1dT<gArray1d>& objectToCompositionObject  )
  { throw GPException("ElementManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n"); }

  //void ConstructListOfBoundaryObjects( gArray1d& objectList ) const ;


  int AddElementToRegion( const int regionNumber, const int* const elemsToNodes );

  void ConstructListOfIndexesFromMap( const Array1dT< std::set< std::pair<ElementRegionT*,localIndex> > >& toElementMap,
                                      const lArray1d& nodeList,
                                      std::map< std::string, lArray1d>& localIndexes,
                                      const int depth );

  void ModifyToElementMapsFromSplit( const std::map< std::string, lSet>& modifiedElements ,
                                     NodeManagerT& nodeManager,
                                     FaceManagerT& faceManager );

  void UpdateExternalityFromSplit( const std::map< std::string, lSet>& modifiedElements ,
                                     NodeManagerT& nodeManager,
                                     EdgeManagerT& edgeManager,
                                     FaceManagerT& faceManager );

  void InitializeFlowFaceRegion();
  void GenerateFlowFaceRegion(FaceManagerT& faceManager);

  using ObjectDataStructureBaseT::WriteSilo;
  using ObjectDataStructureBaseT::ReadSilo;

  void WriteSilo( SiloFile& siloFile,
                  const std::string& meshname,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart,
                  const std::string& regionName = "none",
                  const lArray1d& mask = lArray1d() );

  void ReadSilo( const SiloFile& siloFile,
                 const std::string& meshname,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart );


  void ResetGlobalToLocalMap();

  template< typename T_indices >
  unsigned int PackElements( bufvector& buffer,
                             lSet& sendnodes,
                             lSet& sendfaces,
                             const std::map<std::string,T_indices>& elementList,
                             const NodeManagerT& nodeManager,
                             const FaceManagerT& faceManager,
                             const bool packConnectivityToGlobal,
                             const bool packFields,
                             const bool packMaps,
                             const bool packSets ) const;

  unsigned int UnpackElements( const bufvector& buffer,
                               const NodeManagerT& nodeManager,
                               const FaceManagerT& faceManager,
                               std::map< std::string, lArray1d>& elementRegionReceiveLocalIndices,
                               const bool unpackConnectivityToLocal,
                               const bool unpackFields,
                               const bool unpackMaps,
                               const bool unpackSets );

  unsigned int UnpackElements( const char*& pbuffer,
                               const NodeManagerT& nodeManager,
                               const FaceManagerT& faceManager,
                               std::map< std::string, lArray1d>& elementRegionReceiveLocalIndices,
                               const bool unpackConnectivityToLocal,
                               const bool unpackFields,
                               const bool unpackMaps,
                               const bool unpackSets );

  void ConnectivityFromGlobalToLocal( const std::map< std::string, lSet>& allReceivedElements,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );

  localIndex* ElementToNodeMap( const localIndex elemIndex )
  {
    const RegKeyType regionKey = m_regionIndexToKey[m_ElementIndexToRegionIndex[elemIndex]];
    return m_ElementRegions[regionKey].m_toNodesRelation[m_ElementIndexToRegionLocalIndex[elemIndex]];
  }

  void ElementToNodeMap( const localIndex elemIndex, lArray1d& indices )
  {
    indices.clear();
    const RegKeyType regionKey = m_regionIndexToKey[m_ElementIndexToRegionIndex[elemIndex]];
    const localIndex* iptr = m_ElementRegions[regionKey].m_toNodesRelation[m_ElementIndexToRegionLocalIndex[elemIndex]];
    indices.resize(m_ElementRegions[regionKey].m_numNodesPerElem);
    localIndex a = 0;
    for(lArray1d::iterator it = indices.begin(); it != indices.end(); ++it, ++a)
      *it = iptr[a];
  }

  int m_numElems;

  Array1dT<lArray1d>& m_ElementToNodeMap;
  Array1dT<lArray1d>& m_ElementToFaceMap;
  Array1dT<lArray1d>& m_ElementToElementMap;

  iArray1d& m_ElementIndexToRegionIndex;
  lArray1d& m_ElementIndexToRegionLocalIndex;

  sArray1d m_regionIndexToKey;


  typedef std::string RegKeyType;
  std::map< RegKeyType, ElementRegionT > m_ElementRegions;

private:
  ElementManagerT( const ElementManagerT& );
  ElementManagerT& operator=( const ElementManagerT&);


};

#endif /* ELEMENTMANAGERT_H_ */
