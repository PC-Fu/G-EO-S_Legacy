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
 * @file SolverFactory.h
 * @author walsh24
 * @date Feb 28, 2011
 */

#ifndef SOLVERFACTORY_H_
#define SOLVERFACTORY_H_

#include "SolverBase.h"
#include "Utilities/StringUtilities.h"

#include<map>
#include<string>
#include<vector>

//////////////////////////

// Solver Factory
//
// Consists of the following parts:
//   * The function to generate new solver pointers: "newSolver"
//   * A base class to derive the functions to generate solver pointers: "SolverInitializer"
//   * A String-to-Solver-Intializer map hidden behind the getSolverCatalogue function
//   * A template to create solver initializers: "SolverRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_SOLVER"
// 
// Most solvers will only need to use one or two of the parts:
//   * To register a new solver in the factory: REGISTER_SOLVER( SolverClassName )
//   * To load a solver pointer from the factory:       SolverBase* aSolverPtr = SolverFactory::NewSolver(solverString, args );

/// Base class to generate new Solver pointers
class SolverInitializer
{
public:
  virtual SolverBase* InitializeSolver(TICPP::HierarchicalDataNode* const hdn, ProblemManagerT* const pm) = 0;

  virtual ~SolverInitializer()
  {
  }
};

typedef std::map<std::string, SolverInitializer*> SolverCatalogueType;

class SolverFactory
{
public:
  /// The Solver Factory.
  static SolverBase* NewSolver(const std::string& solverName,
                               TICPP::HierarchicalDataNode* const hdn,
                               ProblemManagerT* const pm);

  /// Interface to the Solver name -> Solver initializer map
  static SolverCatalogueType& GetSolverCatalogue();

  /// Return a list of supported solver names
  static void GetSolverNames(std::vector<std::string>& nameList);
};

/// Template for creating classes derived from SolverInitializer
template<class SolverType>
class SolverRegistrator: public SolverInitializer
{
public:
  SolverRegistrator(void)
  {
    const std::string solverName(SolverType::SolverName());
    SolverFactory::GetSolverCatalogue()[solverName] = this;
  }

  SolverBase* InitializeSolver(TICPP::HierarchicalDataNode* const hdn, ProblemManagerT* const pm)
  {
    if(!hdn)
      throw GPException("Need to specify a valid HierarchicalDataNode to InitializeSolver");

    const std::string name(hdn->GetAttributeString("name"));

    if(name.empty())
      throw GPException("Need to specify a valid name to InitializeSolver");

    SolverBase* ret = new SolverType(name, pm);
    ret->ReadXML(hdn);
    return ret;
  }
};

/// Compiler directive to simplify autoregistration
#define REGISTER_SOLVER( ClassName ) namespace{ SolverRegistrator<ClassName> reg_##ClassName; }

#endif
