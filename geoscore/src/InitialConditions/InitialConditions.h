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
 * @file InitialConditions.h
 * @author walsh24
 * @date Feb 28, 2011
 */

#ifndef InitialConditionFACTORY_H_
#define InitialConditionFACTORY_H_

#include "Common/Common.h"
#include "Common/GPException.h"
#include "Common/typedefs.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "Utilities/StringUtilities.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "ObjectManagers/PhysicalDomainT.h"

#include<map>
#include<string>
#include<vector>


class ProblemManagerT;

/// Base class for all initial conditions
class InitialConditionBase
{
public:
  InitialConditionBase( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  InitialConditionBase( const PhysicalDomainT::ObjectDataStructureKeys& objectType, const std::string& regionName,
                        const std::string& fieldName,  const  std::vector<std::string>& setNames,  const FieldType& fieldType);
  virtual ~InitialConditionBase();
  
  virtual void RegisterFields( PhysicalDomainT& domain ) = 0;
  
  virtual void Apply( PhysicalDomainT& domain )=0;

  virtual void ReadXML( TICPP::HierarchicalDataNode* hdn );

  bool IsThisInitialStress();

protected:
  PhysicalDomainT::ObjectDataStructureKeys objectType_;
  std::string regionName_;
  std::string fieldName_; 
  std::vector<std::string> setNames_; 
  FieldType fieldType_;
};


/// Set initial conditions for field based on data in an ascii file
class ReadInitialConditionFromFile: public InitialConditionBase
{
public:
  ReadInitialConditionFromFile( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~ReadInitialConditionFromFile(){};
  
  virtual void RegisterFields( PhysicalDomainT& domain );
  
  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode* hdn );

  static const char* InitialConditionName(){return "ReadInitialConditionFromFile";};

protected:
  std::string filename_;
  bool isIndexedFile_;
};

/// Set initial condition for a field or subset thereof to a constant value
class ConstantInitialCondition: public InitialConditionBase
{
public:
  ConstantInitialCondition( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~ConstantInitialCondition(){};
  
  virtual void RegisterFields( PhysicalDomainT& domain );
  
  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode*  hdn );

  static const char* InitialConditionName(){return "ConstantInitialCondition";};

protected:
  std::string valueStr_;
};

/// Set intial condition for a field or subset thereof based on a table
class InitialConditionTable: public InitialConditionBase
{
public:
  InitialConditionTable( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~InitialConditionTable(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode*  hdn );

  static const char* InitialConditionName(){return "InitialConditionTable";};

protected:
  template <class TABLE, class ARRAY, class ARRAY2>
  inline static bool SetField(const TABLE* ptable, int component, const ARRAY* pos, const NodeManagerT& nodeManager, ElementRegionT& region, const Array1dT<lSet>& sets, ARRAY2* pfield)
  {
    if((!pfield) || (!ptable))
      return false;
    const TABLE& table = *ptable;
    ARRAY2& field = *pfield;
    const bool isFE = pos == 0;
    if(isFE)
    {
      if(sets.size() > 0)
      {
        for (Array1dT<lSet>::const_iterator iset = sets.begin(); iset != sets.end(); ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin(); si != set.end(); ++si)
          {
            R1Tensor cpos(region.GetElementCenter(*si, nodeManager));
            field[*si][component] = table.Lookup(cpos);
          }
        }
      }
      else
      {
        for (localIndex i = 0; i < region.DataLengths(); ++i)
        {
          R1Tensor cpos(region.GetElementCenter(i, nodeManager));
          field[i][component] = table.Lookup(cpos);
        }
      }
    }
    else
    {
      if(sets.size() > 0)
      {
        for (Array1dT<lSet>::const_iterator iset = sets.begin(); iset != sets.end(); ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin(); si != set.end(); ++si)
            field[*si][component] = table.Lookup((*pos)[*si]);
        }
      }
      else
      {
        localIndex i = 0;
        for (typename ARRAY::const_iterator it = pos->begin(); it != pos->end(); ++it, ++i)
          field[i][component] = table.Lookup(*it);
      }
    }
    return true;
  }

  template <class TABLE, class ARRAY, class ARRAY2>
  inline static bool SetField(const TABLE* ptable, const ARRAY* pos, const NodeManagerT& nodeManager, ElementRegionT& region, const Array1dT<lSet>& sets, ARRAY2* pfield)
  {
    if((!pfield) || (!ptable))
      return false;
    const TABLE& table = *ptable;
    ARRAY2& field = *pfield;
    const bool isFE = pos == 0;
    if(isFE)
    {
      if(sets.size() > 0)
      {
        for (Array1dT<lSet>::const_iterator iset = sets.begin(); iset != sets.end(); ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin(); si != set.end(); ++si)
          {
            R1Tensor cpos(region.GetElementCenter(*si, nodeManager));
            field[*si] = table.Lookup(cpos);
          }
        }
      }
      else
      {
        for (localIndex i = 0; i < region.DataLengths(); ++i)
        {
          R1Tensor cpos(region.GetElementCenter(i, nodeManager));
          field[i] = table.Lookup(cpos);
        }
      }
    }
    else
    {
      if(sets.size() > 0)
      {
        for (Array1dT<lSet>::const_iterator iset = sets.begin(); iset != sets.end(); ++iset)
        {
          const lSet& set = *iset;
          for (lSet::const_iterator si = set.begin(); si != set.end(); ++si)
            field[*si] = table.Lookup((*pos)[*si]);
        }
      }
      else
      {
        localIndex i = 0;
        for (typename ARRAY::const_iterator it = pos->begin(); it != pos->end(); ++it, ++i)
          field[i] = table.Lookup(*it);
      }
    }
    return true;
  }

  std::string m_tableName;
  sArray1d m_variableNames;
  Array1dT<FieldType> m_variableTypes;
  int m_component;
};

/// Set intial condition for a field or subset thereof based on a function of the other field variables
class InitialConditionFunction: public InitialConditionBase
{
public:
  InitialConditionFunction( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~InitialConditionFunction(){};
  
  virtual void RegisterFields( PhysicalDomainT& domain );
  
  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode* hdn );

  static const char* InitialConditionName(){return "InitialConditionFunction";};

protected:
  std::string functionName_;
  sArray1d variableNames_;
  Array1dT<FieldType> variableTypes_;
  int component_;

};

///// Set intial condition for a field or subset thereof based on a 3D table
//class InitialConditionFrom3DTable: public InitialConditionBase
//{
//public:
//  InitialConditionFrom3DTable( const TICPP::HierarchicalDataNode* const InitialConditionNode, const ProblemManagerT* const problemManager);
//  virtual ~InitialConditionFrom3DTable(){};
//
//  virtual void RegisterFields( PhysicalDomainT& domain );
//
//  virtual void Apply( PhysicalDomainT& domain );
//
//  virtual void ReadXML( const TICPP::HierarchicalDataNode* const hdn );
//
//  static const char* InitialConditionName(){return "InitialConditionFrom3DTable";};
//
//protected:
//  std::string m_tableName;
//};
//


/// Calculate the center of each face from the node positions
class CalculateFaceCenters: public InitialConditionBase
{
public:
  CalculateFaceCenters( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~CalculateFaceCenters(){};
  
  virtual void RegisterFields( PhysicalDomainT& domain );
  
  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML(  TICPP::HierarchicalDataNode* hdn  ){};

  static const char* InitialConditionName(){return "CalculateFaceCenters";};

};

/// Calculate the center of each face from the node positions
class CalculateElementCenters: public InitialConditionBase
{
public:
  CalculateElementCenters( TICPP::HierarchicalDataNode* InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~CalculateElementCenters(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML(  TICPP::HierarchicalDataNode* hdn  ){};

  static const char* InitialConditionName(){return "CalculateElementCenters";};

};

/// Calculate the normal of each face from the common plane contact
class CalculateFaceNormals: public InitialConditionBase
{
public:
  CalculateFaceNormals(  TICPP::HierarchicalDataNode*  InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~CalculateFaceNormals(){};
  
  virtual void RegisterFields( PhysicalDomainT& domain );
  
  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML(  TICPP::HierarchicalDataNode*  hdn  ){};

  static const char* InitialConditionName(){return "CalculateFaceNormals";};

};


/// Calculate an aperture of each face from the common plane contact
class CalculateAperture: public InitialConditionBase
{
public:
  CalculateAperture(  TICPP::HierarchicalDataNode*  InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~CalculateAperture(){};
  
  virtual void RegisterFields( PhysicalDomainT& domain );
  
  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode*  hdn  ){};

  static const char* InitialConditionName(){return "CalculateAperture";};

};

/// Join sets
class LinkFractureFaces: public InitialConditionBase
{
public:
  LinkFractureFaces(  TICPP::HierarchicalDataNode*  InitialConditionNode, const ProblemManagerT* const problemManager);
  virtual ~LinkFractureFaces(){};

  virtual void RegisterFields( PhysicalDomainT& domain );

  virtual void Apply( PhysicalDomainT& domain );

  virtual void ReadXML( TICPP::HierarchicalDataNode*  hdn );

  static const char* InitialConditionName(){return "JoinFaces";};

protected:

};

//////////////////////////

// Initial Condition Factory
//
// Consists of the following parts:
//   * The function to generate new pointers: "newInitialCondition"
//   * A base class to derive the functions to generate InitialCondition pointers: "InitialConditionInitializer"
//   * A String-to-InitialCondition-Intializer map hidden behind the getInitialConditionCatalogue function
//   * A template to create InitialCondition initializers: "InitialConditionRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_InitialCondition"
// 
// Most initial conditions will only need to use one or two of the parts:
//   * To register a new InitialCondition in the factory: REGISTER_InitialCondition( InitialConditionClassName )
//   * To load a InitialCondition pointer from the factory:       InitialConditionBase* aInitialConditionPtr = newInitialCondition(InitialConditionString, args );

/// The InitialCondition Factory.
InitialConditionBase* newInitialCondition(const std::string& InitialConditionName, TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);

/// Base class to generate new InitialCondition pointers
class InitialConditionInitializer{
  public:
    virtual InitialConditionBase* initializeInitialCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm) = 0;
    virtual ~InitialConditionInitializer() = 0;
};
inline InitialConditionInitializer::~InitialConditionInitializer() { }

/// Interface to the InitialCondition name -> InitialCondition initializer map
std::map<std::string, InitialConditionInitializer*> & getInitialConditionCatalogue();

/// Return a list of supported InitialCondition names
void getInitialConditionNames( std::vector<std::string>& nameList);

/// Template for creating classes derived from InitialConditionInitializer
template< class InitialConditionType >
class InitialConditionRegistrator : public InitialConditionInitializer{

  public:
    InitialConditionRegistrator(void){
      std::string InitialConditionName = std::string(InitialConditionType::InitialConditionName() );
      getInitialConditionCatalogue() [InitialConditionName] = this;
    }; 

    InitialConditionBase* initializeInitialCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm  ){ 
      return new InitialConditionType( hdn, pm );
    };
};

/// Compiler directive to simplify autoregistration
#define REGISTER_InitialCondition( ClassName ) namespace{ InitialConditionRegistrator<ClassName> reg_##ClassName; }

#endif
