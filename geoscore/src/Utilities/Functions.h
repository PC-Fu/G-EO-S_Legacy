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
#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

/**
 * @file Functions.h
 * @author walsh24
 */

#include "Common/Common.h"
#include "Common/GPException.h"
#include "Common/typedefs.h"
#include "DataStructures/Tables/Table.h"
#include "DataStructures/Tables/TableTypes.h"

#include "fparser.hh"

#include <string>

class ProblemManagerT;
namespace TICPP
{
  class HierarchicalDataNode;
}

/// Base class for all user-defined functions
class Function
{

public:
  Function(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Function(const std::string& name);
  virtual ~Function()
  {
  }

  virtual void ReadXML(TICPP::HierarchicalDataNode* hdn)=0;
  virtual realT operator()(const realT& x)
  {
    return 0.0;
  }

  std::string Name()
  {
    return m_name;
  }

  /// returns name of specific function instance
  /// derived classes must define "static const char* FunctionName()"
  /// which returns name of derived class.

private:
  std::string m_name;
};

/// Class representing a constant function
class ConstantFunction: public Function
{

public:
  ConstantFunction(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  ~ConstantFunction()
  {
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    return m_value;
  }

  static const char* FunctionName()
  {
    return "ConstantFunction";
  }

private:
  realT m_value;
};

/// Class representing a polynomial function
class PolynomialFunction: public Function
{

public:
  PolynomialFunction(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  ~PolynomialFunction()
  {
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x);

  static const char* FunctionName()
  {
    return "PolynomialFunction";
  }

private:
  rArray1d coeffs;
};

/// Class representing a uniform random distribution
class UniformRandomDistribution: public Function
{

public:
  UniformRandomDistribution(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  ~UniformRandomDistribution()
  {
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x);

  static const char* FunctionName()
  {
    return "UniformRandomDistribution";
  }

private:
  realT df;
  realT min;
  bool uniqueAcrossProcessors;
  int rank;
};

/// Function class for user-defined functions
class SymbolicFunction: public Function
{

public:
  SymbolicFunction(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  ~SymbolicFunction()
  {
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    return m_fParser.Eval(&x);
  }

  static const char* FunctionName()
  {
    return "SymbolicFunction";
  }

private:
  FunctionParser m_fParser;
};

/// Function wrapper for 4D tables
class Lookup4DTable: public Function
{

public:
  Lookup4DTable(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Lookup4DTable(const std::string& name, const Table4D* tablePtr) :
      Function(name),
      m_tablePtr(tablePtr)
  {/** Empty **/
  }

  ~Lookup4DTable()
  {/** Empty **/
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    rArray1d xx(3);
    xx[0] = (&x)[0];
    xx[1] = (&x)[1];
    xx[2] = (&x)[2];
    xx[3] = (&x)[3];
    return m_tablePtr->Lookup(xx);
  }

  static const char* FunctionName()
  {
    return "Lookup4DTable";
  }

private:
  const Table4D* m_tablePtr;
};

/// Function wrapper for 3D tables
class Lookup3DTable: public Function
{

public:
  Lookup3DTable(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Lookup3DTable(const std::string& name, const Table3D* tablePtr) :
      Function(name),
      m_tablePtr(tablePtr)
  {/** Empty **/
  }
  ;
  ~Lookup3DTable()
  {/** Empty **/
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    rArray1d xx(3);
    xx[0] = (&x)[0];
    xx[1] = (&x)[1];
    xx[2] = (&x)[2];
    return m_tablePtr->Lookup(xx);
  }

  static const char* FunctionName()
  {
    return "Lookup3DTable";
  }

private:
  const Table3D* m_tablePtr;
};

/// Function wrapper for 2D tables
class Lookup2DTable: public Function
{

public:
  Lookup2DTable(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Lookup2DTable(const std::string& name, const Table2D* tablePtr) :
      Function(name),
      m_tablePtr(tablePtr)
  {/** Empty **/
  }

  ~Lookup2DTable()
  {/** Empty **/
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    rArray1d xx(2);
    xx[0] = (&x)[0];
    xx[1] = (&x)[1];
    return m_tablePtr->Lookup(xx);
  }

  static const char* FunctionName()
  {
    return "Lookup2DTable";
  }

private:
  const Table2D* m_tablePtr;
};

/// Function wrapper for 1D tables
class Lookup1DTable: public Function
{

public:
  Lookup1DTable(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm);
  Lookup1DTable(const std::string& name, const Table1D* tablePtr) :
      Function(name),
      m_tablePtr(tablePtr)
  {/** Empty **/
  }

  ~Lookup1DTable()
  {/** Empty **/
  }

  void ReadXML(TICPP::HierarchicalDataNode* hdn);
  realT operator()(const realT& x)
  {
    rArray1d xx(1);
    xx[0] = x;
    return m_tablePtr->Lookup(xx);
  }

  static const char* FunctionName()
  {
    return "Lookup1DTable";
  }
private:
  const Table1D* m_tablePtr;
};

/// Evaluates a mathematical expression in a string with the warp function parser.
/// Should be used for one-time evaluation of real number expressions only.
/// If a user-defined function must be evaluated multiple times use the 
/// SymbolicFunction class instead. 
realT EvaluateStringFunction(const std::string& fString);

//////////////////////////
//
// Function Factory
//
// Consists of the following parts:
//   * The function to generate new pointers: "newFunction"
//   * A base class to derive the functions to generate Function pointers: "FunctionInitializer"
//   * A String-to-Function-Intializer map hidden behind the getFunctionCatalogue function
//   * A template to create Function initializers: "FunctionRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_Function"
// 
// Most initial conditions will only need to use one or two of the parts:
//   * To register a new Function in the factory: REGISTER_Function( FunctionClassName )
//   * To load a Function pointer from the factory:       Function* aFunctionPtr = newFunction(FunctionString, args );

/// The Function Factory.
Function* newFunction(const std::string& FunctionName, TICPP::HierarchicalDataNode* hdn,
                      const ProblemManagerT* const pm);

/// Base class to generate new Function pointers
class FunctionInitializer
{
public:
  virtual Function* initializeFunction(TICPP::HierarchicalDataNode* hdn,
                                       const ProblemManagerT* const pm) = 0;
  virtual ~FunctionInitializer() {}
};

/// Interface to the Function name -> Function initializer map
std::map<std::string, FunctionInitializer*> & getFunctionCatalogue();

/// Return a list of supported Function names
void getFunctionNames(std::vector<std::string>& nameList);

/// Template for creating classes derived from FunctionInitializer
template<class FunctionType>
class FunctionRegistrator: public FunctionInitializer
{

public:
  FunctionRegistrator(void)
  {
    std::string FunctionName = std::string(FunctionType::FunctionName());
    getFunctionCatalogue()[FunctionName] = this;
  }

  Function* initializeFunction(TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm)
  {
    return new FunctionType(hdn, pm);
  }
};

/// Compiler directive to simplify function autoregistration
#define REGISTER_Function( ClassName ) namespace{ FunctionRegistrator<ClassName> reg_##ClassName; }
#endif
