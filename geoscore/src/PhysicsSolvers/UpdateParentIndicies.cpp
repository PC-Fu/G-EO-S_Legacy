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
 * @file UpdateParentIndicies.cpp
 * @author walsh24
 * @date December 3, 2013
 */

#include "SolverFactory.h"
#include "UpdateParentIndicies.h"
#include "Common/Common.h"
#include "Common/typedefs.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"



//////////////////////////////////////////////////////////////////////////////////////////

// Upate field with function

UpdateParentIndicies::UpdateParentIndicies(  const std::string& name,
                                                   ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

UpdateParentIndicies::~UpdateParentIndicies()
{
  // TODO Auto-generated destructor stub
}

/*
 * <UpdateParentIndicies name="upi"            * name of solver
 *            object="Face"  />                * object to calculate the global parent indicies on (Default=Node)
 */
void UpdateParentIndicies::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML( hdn );

  m_setNames = hdn->GetStringVector("setnames");
  if(m_setNames.empty()) m_setNames = hdn->GetStringVector("setname");
  

  std::string objectTypeStr = hdn->GetAttributeStringOrDefault("object", PhysicalDomainT::FiniteElementNodeManagerStr() );
  {
    std::string oldStr = hdn->GetAttributeStringOrDefault("objecttype", "" ); // throw error if using old syntax
    if(!oldStr.empty()) {
      throw GPException("UpdateParentIndicies: Attempting to set objecttype - use 'object' instead.");
    }
  }
  m_objectKey = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);
//=======
//  std::string objectTypeStr = hdn->GetAttributeString("objecttype");
//  if(objectTypeStr.empty())
//    throw GPException("Must define an object type for updating a field with a function");
//  m_objectKey = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);
//>>>>>>> .r1185
  
  m_regionName = hdn->GetAttributeStringOrDefault("regionname","");
  
  
    
}

void UpdateParentIndicies::RegisterFields( PhysicalDomainT& domain )
{

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);

  objectManager.AddKeylessDataField( FieldInfo::globalIndexField,  "parentGlobalIndex", true, true );

}



/**
 * 
 * 
**/

double UpdateParentIndicies::TimeStep( const realT& time,
                                        const realT& dt,
                                        const int cycleNumber,
                                        PhysicalDomainT& domain,
                                        const sArray1d& namesOfSolverRegions ,
                                        SpatialPartition& partition ,
                                        FractunatorBase* const fractunator )
{  
  realT dt_return = dt;

  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);
  gArray1d& parentIndexes = objectManager.GetFieldData<globalIndex>("parentGlobalIndex");

  if(m_setNames.empty()){
	  for( sArray1d::size_type i =0;  i < objectManager.DataLengths(); ++i){
		  localIndex localPI =  objectManager.GetParentIndex(i);
		  parentIndexes[i] =  objectManager.m_localToGlobalMap[localPI];
	  }
  } else {
    for( sArray1d::size_type i =0; i < m_setNames.size(); ++i){
      lSet& set = objectManager.GetSet(m_setNames[i]);
	  for( lSet::iterator si = set.begin();  si != set.end(); ++si){
		  localIndex localPI =  objectManager.GetParentIndex(*si);
		  parentIndexes[*si] =  objectManager.m_localToGlobalMap[localPI];
	  }

    }
  }
  return dt_return;
}


REGISTER_SOLVER( UpdateParentIndicies )
