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
/*
 * NodeManager.cpp
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#include "NodeManager.h"
#include <fstream>

/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::NodeManager():
ObjectManager(),
m_numNodes(ObjectManager::m_numObjects)
{

  std::vector< std::pair< FieldKey, FieldType > > fields;
  fields.resize(7);

  fields[0] = std::pair< FieldKey, FieldType >(FieldInfo::referencePosition,FieldInfo::R1TensorField);
  fields[1] = std::pair< FieldKey, FieldType >(FieldInfo::displacement,FieldInfo::R1TensorField);
  fields[2] = std::pair< FieldKey, FieldType >(FieldInfo::incrementalDisplacement,FieldInfo::R1TensorField);
  fields[3] = std::pair< FieldKey, FieldType >(FieldInfo::velocity,FieldInfo::R1TensorField);
  fields[4] = std::pair< FieldKey, FieldType >(FieldInfo::acceleration,FieldInfo::R1TensorField);
  fields[5] = std::pair< FieldKey, FieldType >(FieldInfo::force,FieldInfo::R1TensorField);
  fields[6] = std::pair< FieldKey, FieldType >(FieldInfo::mass,FieldInfo::realField);

  this->RegisterFields(fields);
  this->OrganizeFields();
  this->AllocateObjectFields();
}


/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::NodeManager( const NodeManager& init ):
ObjectManager(init),
m_numNodes(ObjectManager::m_numObjects)
{
}


/**
 * @author R.R. Settgast
 * @return
 */
NodeManager::~NodeManager()
{}

/**
 * @author Settgast
 * @param geometryStream open file stream to read nodal data from
 *
 * Needs to be replaced once we settle on a file format
 */
void NodeManager::ReadAsciiNodeInput( std::ifstream& geometryStream )
{

  int globalNodeIndex;
  int junk;

  const int xOffset = m_fieldRegistry.m_fieldEnumToOffset[FieldInfo::referencePosition];

  for( std::vector<Object*>::iterator pnode = m_objects.begin() ; pnode != m_objects.end() ; ++pnode )
  {
    Object* const node = *pnode;
    R1Tensor& m_refposition = node->GetFieldFromOffset<R1Tensor>(xOffset);
    geometryStream>>globalNodeIndex>>m_refposition(0)>>m_refposition(1)>>m_refposition(2)>>junk;

  }

}



/**
 * @author R. Settgast
 * @param destination local node number of destination node
 * @param source local node number of source node
 */
void NodeManager::CopyNode( const int destination, const int source )
{
  *(m_objects[destination]) = *(m_objects[source]);

}

/**
 * @author R. Settgast
 * @param nodenum local node number of node to tranlate
 * @param offset amount to translate node
 */
void NodeManager::TranslateNode( const int nodenum, const R1Tensor& offset  )
{
  m_objects[nodenum]->GetFieldFromOffset<R1Tensor>(m_fieldRegistry.m_fieldEnumToOffset[FieldInfo::referencePosition]);

}
