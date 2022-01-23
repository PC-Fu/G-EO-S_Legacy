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
//  LLNL-CODE-656690
//  GMOD-B, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GMOD-B. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
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
 * @file FaultElementManagerT.h
 * @author Scott Johnson
 * @date created on Dec 2, 2012
 */

#ifndef FAULTLEMENTMANAGERT_H_
#define FAULTLEMENTMANAGERT_H_

#include "Common.h"
#include "dc3d.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "ObjectManagers/FaceManagerT.h"
#include "DataStructures/VectorFields/StableTimeStep.h"
#include "DataStructures/VectorFields/NodeManagerT.h"
#include "ObjectManagers/TableManager.h"
#include "Utilities/GeometryUtilities.h"
#include "Constitutive/Material/LinearElastic.h"
#include "Constitutive/ConstitutivePropertiesTable.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "PhysicsSolvers/Seismicity/FaultRuptureDieterich.h"
#include "DataStructures/Tables/Table.h"

#include "Constitutive/Interface/Dieterich.h"

#include "Contact/SpatialSorterFactory.h"

#if GPAC_MPI
#include <mpi.h>
#endif

/**
 * @author Scott Johnson
 * @brief Class to manager the collection of fault elements
 */
class FaultElementManagerT : public ObjectDataStructureBaseT
{
public:
  FaultElementManagerT();
  FaultElementManagerT(NodeManagerT* nm, FaceManagerT* fm);
  virtual ~FaultElementManagerT();

  void ReadXML(TableManager& tableManager, TICPP::HierarchicalDataNode* hdn);

  void Initialize( );

  virtual void DeserializeObjectField(const std::string& name, const rArray1d& field);

  void ResetStatesAndParameters( const bool resetStress = false );
  void ResetStates();

  void erase( const localIndex i );
  globalIndex resize( const localIndex size, const bool assignGlobals = false );

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject  = NULL) {}
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject  = NULL) {}
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager ,
                                                         Array1dT<gArray1d>& objectToCompositionObject  )
  { throw GPException("ExternalFaceManagerT::ExtractMapFromObjectForAssignGlobalObjectNumbers() shouldn't be called\n"); }

  realT CalculateTimestep(const bool steadyState = false);

  bool TimeStep(const realT dt, const bool steadyState = false);

  realT LastMagnitude();

  const EarthquakeSimulation::FaultRuptureData& LastEvent() const;

  void ApplyBoundaryConditions();

  void WriteSiloFaultElements( SiloFile& siloFile,
                               const std::string& siloDirName,
                               const std::string& meshname,
                               const int centering,
                               const int cycleNum,
                               const realT problemTime,
                               const bool isRestart,
                               const std::string& regionName = "none",
                               const lArray1d& mask = lArray1d() );

  void WriteSiloEvent(const int cycleNum);

  void ReadSiloFaultElements( const SiloFile& siloFile,
                              const std::string& siloDirName,
                              const std::string& meshname,
                              const int centering,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart,
                              const std::string& regionName = "none",
                              const lArray1d& mask = lArray1d() );

  inline Dieterich& DieterichObject() { return m_contact; }

  inline void SetPorePressureSteadyState()
  {
    if(m_porePressure)
      m_porePressure->SetZeroGradient(true);
  }

  inline bool HasPorePressure() const
  {
    return m_porePressure != NULL;
  }

  inline const Table<4, realT>* PorePressure() const { return m_porePressure; }

protected:
  void WriteSiloEvent( SiloFile& siloFile,
                       const std::string& siloDirName,
                       const std::string& meshname,
                       const int cycleNum,
                       const realT problemTime,
                       const std::string& regionName = "none",
                       const lArray1d& mask = lArray1d() );

  virtual globalIndex insert( const localIndex i, const bool assignGlobals = false );

  void WriteNonManagedDataMembersToSilo( SiloFile& siloFile,
                                         const std::string& siloDirName,
                                         const std::string& meshname,
                                         const int centering,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart,
                                         const std::string& multiRoot,
                                         const std::string& regionName = "none",
                                         const lArray1d& mask = lArray1d());

  void ReadNonManagedDataMembersFromSilo( const SiloFile& siloFile,
                                          const std::string& siloDirName,
                                          const std::string& meshname,
                                          const int centering,
                                          const int cycleNum,
                                          const realT problemTime,
                                          const bool isRestart,
                                          const std::string& regionName = "none",
                                          const lArray1d& mask = lArray1d());

  inline realT PorePressure(const realT time,
                            const R1Tensor& center) const
  {
    if(!m_porePressure)
      return 0.0;
    R1TensorT<4> txyz;
    txyz(0) = time;
    txyz(1) = center(0);
    txyz(2) = center(1);
    txyz(3) = center(2);
    return m_porePressure->Lookup(txyz);
  }

  inline R1Tensor MaximumHorizontalStressDirection() const
  {
    R1Tensor ret;
    ret(0) = m_transformStressFrame(0,0);
    ret(1) = m_transformStressFrame(1,0);
    ret(2) = m_transformStressFrame(2,0);
    return ret;
  }

  inline R1Tensor MinimumHorizontalStressDirection() const
  {
    R1Tensor ret;
    ret(0) = m_transformStressFrame(0,1);
    ret(1) = m_transformStressFrame(1,1);
    ret(2) = m_transformStressFrame(2,1);
    return ret;
  }

  inline R1Tensor Up() const
  {
    R1Tensor ret;
    ret(0) = m_transformStressFrame(0,2);
    ret(1) = m_transformStressFrame(1,2);
    ret(2) = m_transformStressFrame(2,2);
    return ret;
  }

  void CalculateFaceNormalsAndCenters(rArray1d& centroidsnormals);

  void SigmaTau(const R1Tensor& x, const R1Tensor& n, const R1Tensor& d,
                realT& stressNormal, realT& stressShear) const;
  static void SigmaTau(const R2Tensor& sigmaGlobal,
                       const R1Tensor& dislocation,
                       const R1Tensor& normal,
                       realT& sigma,
                       realT& tau);


  //void ReadPorePressureNUFT(const std::string& filename);

  void CalculateStiffnessMatrix();
  void Initialize(const rArray1d& centroidsnormals);

  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  //     MEADE ELEMENTS
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------


  static void CalculateMeadeTriStrains(const R1Tensor& xobs,
                                       const realT nu,
                                       const R1Tensor& slipPlaneFrame,
                                       Array1dT<R1Tensor>& xx,
                                       R2SymTensor& S);

  static void TriStrainMeade(const R1Tensor& y,
                             const realT a,
                             const realT b,
                             const realT nu,
                             const R1Tensor& B,
                             R2SymTensor& e);

  static void PlaneFrameTransform(Array1dT<R1Tensor>& xx,
                                  R2Tensor& T);

  static void DeformationToStress(const R2Tensor& strain, /* deformation tensor, see above */
                                  const realT lambda, const realT G, /* Lame's constants */
                                  R2SymTensor& sigma  /* stress tensor as a vector, see above */
                                  );

  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------
  //     OKADA ELEMENTS
  //----------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------

  static void CalculateOkadaQuadStrains(const R1Tensor& xStrike,
                                        const realT alpha,
                                        const realT depth,
                                        const realT length,
                                        const realT width,
                                        const realT dip,
                                        const R1Tensor& dislocationStrike,
                                        R1Tensor& slipStrike,
                                        R2Tensor& SStrike);

public:
  void ApplyInitialConstitutive();

  void SetInitialConstitutive(ConstitutivePropertiesTable& initialConstitutive);

private:
  //geometry
  NodeManagerT* m_nodeManager;
  FaceManagerT* m_faceManager;

  globalIndex m_maxGlobalNumber;

  localIndex m_eventIndex;
  realT m_magnitudeWriteMin;

  Array1dT< lArray1d >& m_neighborList;
  Array1dT< lSet >& m_neighborListInverse;

  realT m_contactBufferOffset;
  realT m_contactBufferFactor;

  bool m_sorted;
  SpatialSorting::SpatialSorterBase* m_sorter;

  //numerics
  bool m_useInfiniteSupport, m_isSerial;
  EarthquakeSimulation::FaultRuptureDieterich m_rupture;
  Dieterich m_contact;
  EarthquakeSimulation::BoundaryElementDataManagerT m_boundaryElementMaterialParameters;

  //material properties and boundary conditions
  LinearElasticParameterData m_properties;
  R2Tensor m_transformStressFrame;
  R1Tensor m_dislocation;

  //stress gradient along the directions of max, min horizontal stress
  //and vertical stress, respectively
  R1Tensor m_dstressDz;

  //pore pressure
  Table<4, realT>* m_porePressure;

  //ConstitutivePropertiesTable* m_initialConstitutive;
  sArray1d m_initialConstitutiveFields;
  sArray1d m_initialConstitutiveTables;
};
#endif /* FAULTLEMENTMANAGERT_H_ */
