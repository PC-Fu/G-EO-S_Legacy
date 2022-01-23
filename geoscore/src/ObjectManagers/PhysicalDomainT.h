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
 * @file PhysicalDomainT.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef PHYSICAL_DOMAIN_T_H_
#define PHYSICAL_DOMAIN_T_H_

#include "DataStructures/EncapsulatedFields/NodeManager.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "ElementManagerT.h"
#include "FaceManagerT.h"
#include "ExternalFaceManagerT.h"
#include "ContactManagerT.h"
#include "CartesianGridManagerT.h"
#include "EdgeManagerT.h"
#include "EllipsoidalContactManagerT.h"
#include "DiscreteElementManagerT.h"
#include "EllipsoidalDiscreteElementManagerT.h"

#include "SurfaceGeneration/JointPopulator.h"
#include "EnergyT.h"

#ifdef SRC_EXTERNAL
#include "src_external/ObjectManagers/FaultElementManagerT.h"
#endif

class PhysicalDomainT
{
public:
  PhysicalDomainT();

  ~PhysicalDomainT();


  int Initialize();

  void SetInterObjectRelations();

  int m_globalDomainNumber;

  //note: because discrete elements have nodes with much reduced degrees-of-freedom,
  //as well as different solution requirements, we keep a separate collection specifically
  //for discrete elements, which, in turn, necessitates a separate collection of faces
  //the ExternalFaceManagerT ties the external faces of the finite elements
  //with those in the discrete elements to provide coupling between the methods

  EllipsoidalDiscreteElementManagerT m_ellipsoidalDiscreteElementManager;
  EllipsoidalContactManagerT m_ellipsoidalContactManager;

  FaceManagerT m_discreteElementSurfaceFaces;
  NodeManagerT m_discreteElementSurfaceNodes;
  DiscreteElementManagerT m_discreteElementManager;

  ContactManagerT m_contactManager;

  NodeManagerT m_feNodeManager;
  ElementManagerT m_feElementManager;
  FaceManagerT m_feFaceManager;
  EdgeManagerT m_feEdgeManager;

//  NodeManagerT m_flowNodes;
//  EdgeManagerT m_flowEdges;
//  FaceManagerT m_flowFaces;

  ExternalFaceManagerT m_externalFaces;

  CartesianGridManagerT m_cartesianGridManager;

#ifdef SRC_EXTERNAL
  FaceManagerT m_faultPatchFaces;
  NodeManagerT m_faultPatchNodes;
  CartesianGridManagerT m_faultPorePressure;
  FaultElementManagerT m_faultElementManager;
  NodeManagerT m_microseismicSourceNodes;
#endif
  JointPopulator m_jointSets;

  EnergyT m_energy;

  void WriteSilo( SiloFile& siloFile,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart,
                  const bool writeFEMMesh = true,
                  const bool writeFEMFaces = false,
                  const bool writeFEMEdges = false,
                  const bool writeDE = true,
                  const bool writeCP = true,
                  const bool writeCG = true,
                  const bool writeFaultElements = true,
                  const bool writeMicroseismicSources = true);



  void ReadSilo( const SiloFile& siloFile,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart,
                 const bool writeFEMMesh = true,
                 const bool writeFEMFaces = false,
                 const bool writeFEMEdges = false,
                 const bool writeDE = true,
                 const bool writeCP = true,
                 const bool writeCG = true,
                 const bool writeFaultElements = true,
                 const bool writeMicroseismicSources = true);

  void WriteFiniteElementMesh( SiloFile& siloFile,
                               const int cycleNum,
                               const realT problemTime,
                               const bool isRestart,
                               const bool writeFEMMesh,
                               const bool writeFEMFaces,
                               const bool writeFEMEdges);

  void ReadFiniteElementMesh( const SiloFile& siloFile,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart );

#ifdef SRC_EXTERNAL
  void WriteFaultElementMesh( SiloFile& siloFile,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart,
                              const bool writeFaultElements);
  void ReadFaultElementMesh( const SiloFile& siloFile,
                             const int cycleNum,
                             const realT problemTime,
                             const bool isRestart );
  void WriteMicroseismicSources( SiloFile& siloFile,
                                 const int cycleNum,
                                 const realT problemTime,
                                 const bool isRestart,
                                 const bool writeMSS );

  void ReadMicroseismicSources( const SiloFile& siloFile,
                                 const int cycleNum,
                                 const realT problemTime,
                                 const bool isRestart );

#endif

  void WriteDiscreteElements( SiloFile& siloFile,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart,
                              const bool writeDE );

  void ReadDiscreteElements( const SiloFile& siloFile,
                             const int cycleNum,
                             const realT problemTime,
                             const bool isRestart );

  void WriteEllipsoidalDiscreteElements( SiloFile& siloFile,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart,
                                         const bool writeEDE );

  void ReadEllipsoidalDiscreteElements( const SiloFile& siloFile,
                                        const int cycleNum,
                                        const realT problemTime,
                                        const bool isRestart );


  void WriteCommonPlanes( SiloFile& siloFile,
                          const int cycleNum,
                          const realT problemTime,
                          const bool isRestart,
                          const bool writeCP );

  void ReadCommonPlanes( const SiloFile& siloFile,
                         const int cycleNum,
                         const realT problemTime,
                         const bool isRestart );

  void WriteCartesianGrid( SiloFile& siloFile,
                          const int cycleNum,
                          const realT problemTime,
                          const bool isRestart,
                          const bool writeCG );

  void ReadCartesianGrid( const SiloFile& siloFile,
                          const int cycleNum,
                          const realT problemTime,
                          const bool isRestart );

  void RegisterBCFields();
  void RegisterBCFields(ObjectDataStructureBaseT& objectManager);

  void UpdateEnergy();

  //NOTE: the order of these is important; they must be in order of ascending dependence
  //      e.g., faces are dependent on nodes, so faces must be of greater ordinal value than nodes
  enum ObjectDataStructureKeys
  {
    EllipsoidalDiscreteElementManager  = 0,
    EllipsoidalContactManager = 1,
    DiscreteElementNodeManager = 2,
    DiscreteElementFaceManager = 3,
    DiscreteElementManager = 4,
    ContactManager = 5,
    FiniteElementNodeManager = 6,
    FiniteElementEdgeManager = 7,
    FiniteElementFaceManager = 8,
    FiniteElementElementManager = 9,
    FiniteElementElementRegion = 10,
    ExternalFaceManager = 11,
    CartesianGridManager = 12,
    FaultPatchNodeManager = 13,
    FaultPatchFaceManager = 14,
    FaultPatchElementManager = 15,
    MicroseismicSourceNodeManager = 16,
    VirtualEdgeManager = 17,
    VirtualFaceManager = 18,
    Last_ObjectDataStructureNames_Index = 18,
    numObjectDataStructureNames = 19
  };




  static const char* EllipsoidalDiscreteElementManagerStr() { return "EllipsoidalDiscreteElement"; }
  static const char* EllipsoidalContactManagerStr()         { return "EllipsoidalContactManager"; }
  static const char* EllipsoidalConditionStr()              { return "EllipsoidalDiscreteElement"; }

  static const char* DiscreteElementNodeManagerStr()        { return "DiscreteElement_NodeManager"; }
  static const char* DiscreteElementFaceManagerStr()        { return "DiscreteElement_FaceManager"; }
  static const char* DiscreteElementManagerStr()            { return "DiscreteElementManager"; }
  static const char* DiscreteElementConditionStr()          { return "DiscreteElement"; }

  static const char* ContactManagerStr()                    { return "ContactManager"; }
  static const char* CartesianGridManagerStr()              { return "CartesianGrid"; }

  static const char* FiniteElementNodeManagerStr()          { return "FiniteElement_NodeManager"; }
  static const char* FiniteElementNodeConditionStr()        { return "Node"; }
  static const char* FiniteElementEdgeManagerStr()          { return "FiniteElement_EdgeManager"; }
  static const char* FiniteElementEdgeConditionStr()        { return "Edge"; }
  static const char* FiniteElementFaceManagerStr()          { return "FiniteElement_FaceManager"; }
  static const char* FiniteElementFaceConditionStr()        { return "Face"; }
  static const char* FiniteElementElementManagerStr()       { return "FiniteElement_ElementManager"; }
  static const char* FiniteElementElementConditionStr()     { return "Element"; }
  static const char* FiniteElementElementRegionStr()        { return "FiniteElement_ElementRegion"; }

  static const char* ExternalFaceManagerStr()               { return "ExternalFaceManager"; }
  static const char* ExternalFaceConditionStr()             { return "ExternalFace"; }

  static const char* FaultNodeManagerStr()                  { return "Fault_NodeManager"; }
  static const char* FaultFaceManagerStr()                  { return "Fault_FaceManager"; }
  static const char* FaultElementManagerStr()               { return "Fault_ElementManager"; }
  static const char* FaultElementConditionStr()             { return "FaultElement"; }

  static const char* MicroseismicNodeManagerStr()           { return "Microseismic_NodeManager"; }
  static const char* MicroseismicNodeConditionStr()         { return "MicroseismicNode"; }

  static const char* VirtualEdgeManagerStr()                { return "Virtual_EdgeManager"; }
  static const char* VirtualFaceManagerStr()                { return "Virtual_FaceManager"; }


  ObjectDataStructureBaseT& GetObjectDataStructure( ObjectDataStructureKeys key,
                                                    const std::string& regionName="" );

  static ObjectDataStructureKeys GetObjectDataStructureConditionKey(const std::string name);
  static ObjectDataStructureKeys GetObjectDataStructureKey(const std::string name);
  static std::string GetObjectDataStructureName(const ObjectDataStructureKeys key);

  template< typename T_indices >
  void Pack( const ObjectDataStructureKeys name,
             const T_indices& sendIndices,
             bufvector& buffer,
             const bool packConnectivityToGlobal,
             const bool packFields,
             const bool packMaps,
             const bool packSets );

  template< typename T_indices >
  void Unpack(const ObjectDataStructureKeys name,
              const char*& pbuffer,
              T_indices& indices,
              const bool unpackConnectivityToLocal,
              const bool unpackFields,
              const bool unpackMaps,
              const bool unpackSets  );




private:
  PhysicalDomainT( const PhysicalDomainT& );
  PhysicalDomainT& operator=( const PhysicalDomainT& rhs );

};

// the following facilitates looping over object data structure keys
// nb keys must be contiguous and 0 start.
inline
PhysicalDomainT::ObjectDataStructureKeys& operator++(PhysicalDomainT::ObjectDataStructureKeys& objDSK)
{
  // if ENUMS ARE CONTIGUOUS
  int i = static_cast<int>(objDSK);
  if( ++i > PhysicalDomainT::numObjectDataStructureNames )
    i = 0;
  objDSK = static_cast<PhysicalDomainT::ObjectDataStructureKeys>(i);
  return objDSK;
}


#endif /* PHYSICAL_DOMAIN_T_H_ */
