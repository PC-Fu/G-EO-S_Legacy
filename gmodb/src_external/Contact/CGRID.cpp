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

/*
 * CGRID.cpp
 *
 *  Created on: Nov 15, 2011
 *      Author: scottjohnson
 */

#include "CGRID.h"
#include "Utilities/Utilities.h"

#include "Contact/SpatialSorterFactory.h"

namespace SpatialSorting
{

  CGRID::CGRID() :
      swapBegin(),
      swapEnd(),
      binIndices(),
      swap(),
      swapIndices()
  {
    this->SetMask();
  }

  CGRID::~CGRID()
  {
    // TODO Auto-generated destructor stub
  }

  /**
   * @brief Set the minima and bin size
   * @author Scott Johnson
   * @param[in] radii Radii of the potentials
   * @param[in] centers Centers of the potentials
   */
  bool CGRID::SetLimits(const rArray1d& radii, const Array1dT<R1Tensor>& centers,
                        const realT binFactor, const int* const excludeFromContact)
  {
    //get extrema and bin size (as the minimum of the element size)
    R1Tensor upper(-std::numeric_limits<realT>::max());
    this->lower = upper;
    this->lower *= -1.0;
    this->binSize = std::numeric_limits<realT>::max();

    if (excludeFromContact != 0)
    {
      for (localIndex i = 0; i < radii.size(); i++)
      {
        if (excludeFromContact[i] > 0)
          continue;

        R1Tensor tmp = centers[i];
        tmp -= radii[i];
        this->lower.SetMin(tmp);
        upper.SetMax(tmp);
        tmp = 2 * radii[i];
        this->binSize.SetMin(tmp);
      }
    }
    else
    {
      for (localIndex i = 0; i < radii.size(); i++)
      {
        R1Tensor tmp = centers[i];
        tmp -= radii[i];
        this->lower.SetMin(tmp);
        upper.SetMax(tmp);
        tmp = 2 * radii[i];
        this->binSize.SetMin(tmp);
      }
    }

    if (isEqual(this->binSize(0), std::numeric_limits<realT>::max()) || isEqual(
        this->binSize(1), std::numeric_limits<realT>::max()) || isEqual(
        this->binSize(2), std::numeric_limits<realT>::max()))
      return false;
    this->binSize *= binFactor;

    //calculate number of bins
    upper -= this->lower;
    for (unsigned int i = 0; i < nsdof; i++)
    {
      realT tt = upper(i) / this->binSize(i);
      this->nbins[i] = ((localIndex) tt) + 1;
    }
    return true;
  }

  bool CGRID::Update(const rArray1d& radii ,
                     const Array1dT<R1Tensor>& x ,
                     const lSet& toResort ,
                     Array1dT<lArray1d>& neighborList ,
                     Array1dT<lSet>& neighborListInverse ,
                     const int* const excludeFromSorting )
  {
    return false;
  }

  /**
   * @brief Sort the potentials into the given neighbor list
   * @author Scott Johnson
   * @param[in] radii List of potential radii
   * @param[in] center List of potential centers
   * @param[out] neighborList Neighbor list
   * @return return
   */
  bool CGRID::Sort(const rArray1d& radii ,
                   const Array1dT<R1Tensor>& x ,
                   Array1dT<lArray1d>& neighborList ,
                   Array1dT<lSet>& neighborListInverse ,
                   const int* const excludeFromSorting )
  {
//    Bin(radii, centers, neighborList);
//
//    //clear neighbor list
//    //for each dimension update sweep list
//    //go through each dimension in neighbors and do the neighbor comparison where the target's list is augmented by the sweep list
//
//
//    //for each bin, you create neighbor connectivity as an 8-list bin with mapping given by pairs
//    //
//
//    this->neighbors[bin[0]][bin[1]][bin[2]] = Array1dT<lArray1d*>();
//    for(int i = 0; i < this->pairs.size(); i++)
//    {
//      lArray1d n1(this->pairs[i].first);
//      lArray1d n2(this->pairs[i].second);
//      for(int i = 0; i < nsdof; i++)
//      {
//        n1[i] += bin[i];
//        n2[i] += bin[i];
//      }
//      lArray1d* s1 = BinContents(n1);
//      lArray1d* s2 = BinContents(n2);
//      if(s1 == 0 || s2 == 0)
//        continue;
//      //add s2 to neighbor1
//      this->neighbors[n1[0]][n1[1]][n1[2]].push_back(s1);
//    }
    return false;
  }

}

/// Register spatial sorter in the spatial sorter factory
REGISTER_SPATIALSORTER( CGRID )
