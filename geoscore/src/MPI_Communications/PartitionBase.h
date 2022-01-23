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
 * @file PartitionBase.h
 * @author settgast1
 * @date Mar 16, 2011
 */

#ifndef PARTITIONBASE_H_
#define PARTITIONBASE_H_

#include "Common/Common.h"
#include "NeighborCommunication.h"
#include "IO/ticpp/TinyXMLParser.h"

#include <mpi.h>


class oBinStream;
class iBinStream;




class PartitionBase
{

public:

  virtual ~PartitionBase();

  virtual void Initialize() = 0;

  void SetDomain( PhysicalDomainT& domain );

  virtual bool IsCoordInPartition( const R1Tensor& elemCenter ) = 0;

  virtual bool IsCoordInContactGhostRange( const R1Tensor& elemCenter ) = 0;


  virtual void ReadXMLInput( TICPP::HierarchicalDataNode& hdn) = 0;


  virtual void AssignGlobalIndices( PhysicalDomainT& domain );

  virtual void FindMatchedBoundaryIndices( const PhysicalDomainT::ObjectDataStructureKeys key,
                                           const ObjectDataStructureBaseT& object );


  virtual void SetUpNeighborLists( PhysicalDomainT& domain,
                                   const bool contactActive );

  void SetRankOfNeighborNeighbors();

//  virtual void ResetNeighborLists( PhysicalDomainT& domain,
//                                   const int elementGhostingDepth );

  virtual void ModifyGhostsAndNeighborLists( const ModifiedObjectLists& modifiedObjects );

  template< typename T >
  void SendReceive( const Array1dT<Array1dT<T> >& sendArray, Array1dT<Array1dT<T> >& recvArray );

  void SetBufferSizes( const std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d>& fieldNames,
                       const CommRegistry::commID commID = CommRegistry::genericComm01 );

  void SynchronizeFields( const std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d>& fieldNames,
                          const CommRegistry::commID commID = CommRegistry::genericComm01 );

  void SetOwnedByRank( const std::map<PhysicalDomainT::ObjectDataStructureKeys, gArray1d>& localBoundaryGlobalIndices,
                       std::map<PhysicalDomainT::ObjectDataStructureKeys, std::map< globalIndex, int > >& boundaryOwnership);

  void SetGhostArrays( PhysicalDomainT& domain );

  lArray1d GetFaceSendIndices();

  virtual void SetContactGhostRange( const double bufferSize ) = 0;

  int NumberOfNeighbors( ) {return m_neighbors.size();}

  int m_size;
  int m_rank;

  virtual int GetColor() = 0;

  int Color() const {return m_color;}
  int NumColor() const {return m_numColors;}

  void WriteSilo( SiloFile& siloFile );

  void ReadSilo( const SiloFile& siloFile );

  int m_sizeMetis;
  void DeleteExcessNeighbors();


protected:
  PartitionBase();
  PartitionBase( const unsigned int numPartitions, const unsigned int thisPartiton );


  VectorT<NeighborCommunication> m_neighbors;

  VectorT<MPI_Request> m_mpiRequest;
  VectorT<MPI_Status> m_mpiStatus;

  R1Tensor m_contactGhostMin;
  R1Tensor m_contactGhostMax;

  int m_color;
  int m_numColors;

  PhysicalDomainT* const m_domain;

public:
  realT m_t1;
  realT m_t2;
  realT m_t3;
  realT m_t4;

  bool m_hasLocalGhosts;
  std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d> m_localGhosts;
  std::map< std::string, lArray1d> m_elementRegionsLocalGhosts;

  std::map<PhysicalDomainT::ObjectDataStructureKeys, lArray1d> m_localGhostSources;
  std::map< std::string, lArray1d> m_elementRegionsLocalGhostSources;

  int m_ghostDepth;

private:
  virtual void AssignGlobalIndices( ObjectDataStructureBaseT& object, const ObjectDataStructureBaseT& compositionObject );

  void CommunicateRequiredObjectIndices();

  virtual void WriteSiloDerived( SiloFile& siloFile ) = 0;

  virtual void ReadSiloDerived( const SiloFile& siloFile ) = 0;


};

#endif /* PARTITIONBASE_H_ */
