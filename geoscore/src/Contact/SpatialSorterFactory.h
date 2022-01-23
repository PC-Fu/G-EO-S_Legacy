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
 * @file SpatialSorterFactory.h
 * @author Scott Johnson
 * @date Oct 6, 2013
 */

#ifndef SPATIALSORTERFACTORY_H_
#define SPATIALSORTERFACTORY_H_

#include "SpatialSorterBase.h"
#include "Utilities/StringUtilities.h"

#include <map>
#include <string>
#include <vector>

//////////////////////////

// SpatialSorter Factory
//
// Consists of the following parts:
//   * The function to generate new spatialSorter pointers: "NewSpatialSorter"
//   * A base class to derive the functions to generate spatialSorter pointers: "SpatialSorterInitializer"
//   * A String-to-SpatialSorter-Intializer map hidden behind the GetSpatialSorterCatalogue function
//   * A template to create spatialSorter initializers: "SpatialSorterRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_SPATIALSORTER"
// 
// Most spatialSorters will only need to use one or two of the parts:
//   * To register a new spatialSorter in the factory: REGISTER_SPATIALSORTER( SpatialSorterClassName )
//   * To load a spatialSorter pointer from the factory:       SpatialSorterBase* aSpatialSorterPtr = SpatialSorterFactory::NewSpatialSorter(spatialSorterString, args );

/// Base class to generate new SpatialSorter pointers

namespace SpatialSorting
{
  class SpatialSorterInitializer
  {
  public:
    virtual SpatialSorterBase* InitializeSpatialSorter() = 0;

    virtual ~SpatialSorterInitializer()
    {
    }
  };

  typedef std::map<std::string, SpatialSorterInitializer*> SpatialSorterCatalogueType;

  class SpatialSorterFactory
  {
  public:
    /// The SpatialSorter Factory.
    static SpatialSorterBase* NewSpatialSorter(const std::string& spatialSorterName);

    /// Interface to the SpatialSorter name -> SpatialSorter initializer map
    static SpatialSorterCatalogueType& GetSpatialSorterCatalogue();

    /// Return a list of supported spatialSorter names
    static void GetSpatialSorterNames(std::vector<std::string>& nameList);
  };

  /// Template for creating classes derived from SpatialSorterInitializer
  template<class SpatialSorterType>
  class SpatialSorterRegistrator: public SpatialSorterInitializer
  {
  public:
    SpatialSorterRegistrator(void)
    {
      SpatialSorterFactory::GetSpatialSorterCatalogue()[SpatialSorterType::SpatialSorterName()] = this;
    }

    SpatialSorterBase* InitializeSpatialSorter()
    {
      return new SpatialSorterType();
    }
  };
}

/// Compiler directive to simplify autoregistration
#define REGISTER_SPATIALSORTER( ClassName ) namespace SpatialSorting{ SpatialSorterRegistrator<ClassName> reg_##ClassName; }

#endif
