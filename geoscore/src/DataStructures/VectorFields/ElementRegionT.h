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

#ifndef ELEMENTOBJECTT_H_
#define ELEMENTOBJECTT_H_

#include "Common/Common.h"
#include "ObjectDataStructureBaseT.h"
#include "StableTimeStep.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "ObjectManagers/EnergyT.h"

class IntegrationRuleT;
class MaterialBaseParameterDataT;
class MaterialBaseStateDataT;
class MaterialBaseT;
class NodeManagerT;
class FaceManagerT;
class EdgeManagerT;
class ElementRegionT;
class PhysicalDomainT;

class FiniteElementBase;
class Quadrature ;
class Basis;

class MaterialBase;


/**
 * Class to manage the data stored at the element level.
 */
class ElementRegionT : public ObjectDataStructureBaseT
{
public:



  static const char ElementObjectToElementManager[];
  static const char ElementToNode[];
  static const char ElementToFace[];
  static const char ElementToEdge[];


  ElementRegionT();


  ElementRegionT(const ElementRegionT& init);
  

  virtual ~ElementRegionT();

  void Initialize( const PhysicalDomainT& domain );

  virtual void DeserializeObjectField(const std::string& ifield, const rArray1d& field);
  virtual void DeserializeObjectFields(const sArray1d& names, const Array1dT<rArray1d>& fields);
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn,
                const bool isRestart );

  globalIndex resize( const localIndex size, const bool assignGlobals=false );

  void SetGeometryBasedVariables();

  void AllocateElementLibrary( const int basis,
                               const int quadrature );

  void Initialize();

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL);
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL) {}
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager ,
                                                         Array1dT<gArray1d>& objectToCompositionObject  )
  { throw GPException("ElementRegionT::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n"); }


  void ModifyToElementMapsFromSplit( const lSet& modifiedElements ,
                                     NodeManagerT& nodeManager,
                                     FaceManagerT& faceManager );

  void UpdateExternalityFromSplit( const lSet& modifiedElements ,
                                   NodeManagerT& nodeManager,
                                   EdgeManagerT& edgeManager,
                                   FaceManagerT& faceManager );



  int CalculateShapeFunctionDerivatives( const NodeManagerT& nodeManager );


  int CalculateVelocityGradients( const NodeManagerT& nodeManager, const int calcGroup = 0 );

  int MaterialUpdate(const realT dt);

  int CalculateSmallDeformationNodalForces( NodeManagerT& nodeManager ,
                                                            StableTimeStep& timeStep,
                                                            const realT dt );
  int CalculateNodalForces( NodeManagerT& nodeManager ,
                            StableTimeStep& timeStep,
                            const realT dt );

  int CalculateNodalForcesFromOneElement(const localIndex nodeID,
                                         const localIndex elemID,
                                         NodeManagerT& nodeManager,
                                         R1Tensor& fNode );

  void CalculateNodalForceFromStress(const localIndex elemID,
                                     const NodeManagerT& nodeManager,
                                     R2SymTensor& stress,
                                     Array1dT<R1Tensor>& fNode);

  int ProcessElements(  NodeManagerT& nodeManager,
                        StableTimeStep& timeStep,
                        const realT dt  );


  int CalculateNodalMasses( NodeManagerT& nodeManager );


  void SetIsAttachedToSendingGhostNode( const NodeManagerT& nodeManager );

  template< typename T_indices >
  unsigned int PackElements( bufvector& buffer,
                             lSet& sendnodes,
                             lSet& sendfaces,
                             const T_indices& elementList,
                             const NodeManagerT& nodeManager,
                             const FaceManagerT& faceManager,
                             const bool packConnectivityToGlobal,
                             const bool packFields,
                             const bool packMaps,
                             const bool packSets  ) const;

  unsigned int UnpackElements( const char*& buffer,
                               const NodeManagerT& nodeManager,
                               const FaceManagerT& faceManager,
                               lArray1d& elementRegionReceiveLocalIndices,
                               const bool unpackConnectivityToLocal,
                               const bool unpackFields,
                               const bool unpackMaps,
                               const bool unpackSets  );

  void ConnectivityFromGlobalToLocal( const lSet& list,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );


  template< typename T_indices >
  unsigned int PackFieldsIntoBuffer( bufvector& buffer,
                                     const sArray1d& fieldNames,
                                     const T_indices& localIndices,
                                     const bool doBufferPacking=true ) const;

  template< typename T_indices >
  unsigned int PackFieldsIntoBuffer( char*& buffer,
                                     const sArray1d& fieldNames,
                                     const T_indices& localIndices,
                                     const bool doBufferPacking=true ) const;

  unsigned int UnpackFieldsFromBuffer( const char*& buffer,
                                       const sArray1d& fieldNames,
                                       const lArray1d& localIndices );

  template< typename T_indices >
  unsigned int PackAllFieldsIntoBuffer( bufvector& buffer,
                                        const T_indices& localIndices ) const;

  unsigned int UnpackAllFieldsFromBuffer( const char*& buffer,
                                          const lArray1d& localIndices );

  void UpdateElementFieldsWithGaussPointData();

  lArray2d& ElementToNodeMap();
  localIndex* ElementToNodeMap( const int elemNum ) ;
  const localIndex* ElementToNodeMap( const int elemNum ) const;
  localIndex ElementToNodeMap( const int elemNum, const int nodeNum ) const ;

  void GetFaceNodes( const localIndex elementIndex,
                     const localIndex localFaceIndex,
                     lArray1d& nodeIndicies ) const;

  void GetElementNeighbors(localIndex el, 
                           const FaceManagerT& faceManager, 
                           std::set<localIndex>& neighbors) const;
                     
  R1Tensor GetElementCenter(localIndex k, const NodeManagerT& nodeManager) const;


  void WriteSiloRegionMesh( SiloFile& siloFile,
                            const std::string& meshname,
                            const int cycleNum,
                            const realT problemTime,
                            const bool isRestart,
                            const std::string& regionName );
  
  void ReadSiloRegionMesh( const SiloFile& siloFile,
                           const std::string& meshname,
                           const int cycleNum,
                           const realT problemTime,
                           const bool isRestart,
                           const std::string& regionName );

  /// returns true if element id pair is inside region
  bool ContainsElement(const ElementIdPair& ep);

  iArray1d SiloNodeOrdering();

  std::string m_regionName;
  int m_regionNumber;
  const localIndex& m_numElems;
  unsigned int m_numNodesPerElem;
  unsigned int m_numIntegrationPointsPerElem;

  std::string m_elementType;
  std::string m_elementGeometryID;
  std::vector<std::string> m_parentFaceSetNames;

  unsigned int m_ElementDimension;


//  lArray1d& m_ElementObjectToElementManagerMap;
  FixedOneToManyRelation& m_toNodesRelation;
  FixedOneToManyRelation& m_toFacesRelation;
  //FixedOneToManyRelation& m_toEdgesRelation;  //This is neither filled nor used.


  Array1dT< Array2dT<R1Tensor> > m_dNdX;
  rArray2d m_detJ;
  rArray2d m_detJ_n;
  rArray2d m_detJ_np1;
  Array2dT< R2Tensor >  m_dUdX;
  Array2dT< R2Tensor > m_Finv;





  Array1dT< Array1dT<R2SymTensor> > m_Dadt;
  Array1dT< Array1dT<R2Tensor> >  m_Rot;


  Array1dT< rArray2d > m_Ke;
  Array1dT< rArray2d > m_matrixB;
  Array1dT< rArray2d > m_matrixE;



  unsigned int m_basis;
  unsigned int m_quadrature;

  FiniteElementBase* m_finiteElement;
  Quadrature* m_elementQuadrature;
  Basis* m_elementBasis;

//  MaterialManagerT m_material;

//  MaterialBaseT* m_materialComputations;


  //const int m_numFacesPerElement;
  //const int m_numNodesPerFace;
  int m_numFacesPerElement;
  int m_numNodesPerFace;

  EnergyT m_energy;

  realT m_hgDamp;
  realT m_hgStiff;


#if USECPP11==1
  std::unique_ptr<MaterialBase> m_mat;
#else
  MaterialBase* m_mat;
#endif


private:
  ElementRegionT& operator=(const ElementRegionT& rhs);



};


inline lArray2d& ElementRegionT::ElementToNodeMap()
{
  return m_toNodesRelation;
}


inline localIndex* ElementRegionT::ElementToNodeMap( const int elemNum )
{
  return &(m_toNodesRelation(elemNum,0));
}

inline const localIndex* ElementRegionT::ElementToNodeMap( const int elemNum ) const
{
  return &(m_toNodesRelation(elemNum,0));
}

inline localIndex ElementRegionT::ElementToNodeMap( const int elemNum, const int nodeNum ) const
{
  return m_toNodesRelation(elemNum,nodeNum);
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////










#endif /* ELEMENTOBJECTT_H_ */
