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
/************************************************************
 * @file DESS.h
 * @date Nov 11, 2011
 * @author Scott Johnson
 *
 * @brief DESS spatial sorting class (Perkins, Engineering Computations, 2001)
 ************************************************************/

#ifndef DESS_H_
#define DESS_H_

#include "Contact/SpatialSorterBase.h"

namespace SpatialSorting
{
  class DESSBound
  {
  public:
    bool isLower;
    R1Tensor position;
    localIndex indexInInput;
    localIndex indexAlongAxis[nsdof];

    DESSBound() : isLower(false), position(0.0), indexInInput(0)
    {
    }

    ~DESSBound(){}

    template<int T>
    static bool Compare(const DESSBound* const b0, const DESSBound* const b1)
    {
      return b0->position(T) < b1->position(T);
    }

    inline void Set(const realT radius, const R1Tensor& center, const localIndex i, const bool isLowerLocal)
    {
      position = center;
      position += isLowerLocal ? -radius : radius;
      this->isLower = isLowerLocal;
      indexInInput = i;
      indexAlongAxis[0] = std::numeric_limits<localIndex>::max();
      indexAlongAxis[1] = std::numeric_limits<localIndex>::max();
      indexAlongAxis[2] = std::numeric_limits<localIndex>::max();
    };
  };

  class DESS : public SpatialSorterBase
  {
  public:
    DESS();
    virtual ~DESS();

    static std::string SpatialSorterName() { return "DESS"; }

    virtual bool Sort(const rArray1d& radii, const Array1dT<R1Tensor>& x,
                      Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                      const int* const excludeFromSorting = 0);

    virtual bool Update(const rArray1d& radii, const Array1dT<R1Tensor>& centers, const lSet& toResort,
                        Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                        const int* const excludeFromSorting = 0);

    void WriteGrid();
  private:
    void SortSub(const localIndex itgt, const rArray1d& radii,
                 const Array1dT<R1Tensor>& centers, Array1dT<lArray1d>& neighborList,
                 Array1dT<lSet>& neighborListInverse, Array1dT<lSet>& shapesInRange,
                 const bool addLeft = false);

    void AddBounds(const rArray1d& radii, const Array1dT<R1Tensor>& centers, const localIndex ii);

    void Allocate(const localIndex size);

    void ClearLists();

    void Populate(const rArray1d & radii,
                  const Array1dT<R1Tensor> & centers,
                  const int* const excludeFromContact = 0);

    void Populate(const rArray1d & radii,
                  const Array1dT<R1Tensor> & centers,
                  const lSet& toResort);

    template<int T>
    void SortGridAxis();

    void SortGrid();

    void DownSelectShapesInRange(const localIndex itgt,
                                 const lArray1d& order,
                                 Array1dT<lSet>& shapesInRange,
                                 const bool addLowerIndexCandidates = true);

    void SetShapesInRange(const DESSBound& lower,
                          const DESSBound& upper,
                          const localIndex dimension,
                          lSet& shapesInRangeAxis,
                          const bool addLeft = false);

    void SetShapesInRange(const DESSBound& lower,
                          const DESSBound& upper,
                          const lArray1d& order,
                          Array1dT<lSet>& shapesInRange,
                          const bool addLeft = false);

    static void SetOrder(const DESSBound& lower, const DESSBound& upper, lArray1d& order)
    {
      lArray1d diffs(nsdof);
      diffs[0] = upper.indexAlongAxis[0] - lower.indexAlongAxis[0];
      diffs[1] = upper.indexAlongAxis[1] - lower.indexAlongAxis[1];
      diffs[2] = upper.indexAlongAxis[2] - lower.indexAlongAxis[2];
      if(diffs[0] < diffs[1])
      {
        if(diffs[2] < diffs[0])
        {
          order[0] = 2;
          order[1] = 0;
          order[2] = 1;
        }
        else if(diffs[2] < diffs[1])
        {
          order[0] = 0;
          order[1] = 2;
          order[2] = 1;
        }
        else
        {
          order[0] = 0;
          order[1] = 1;
          order[2] = 2;
        }
      }
      else
      {
        if(diffs[2] < diffs[1])
        {
          order[0] = 2;
          order[1] = 1;
          order[2] = 0;
        }
        else if(diffs[2] < diffs[0])
        {
          order[0] = 1;
          order[1] = 2;
          order[2] = 0;
        }
        else
        {
          order[0] = 1;
          order[1] = 0;
          order[2] = 2;
        }
      }
    };

    Array1dT<SpatialSorting::DESSBound> dessBounds;
    Array1dT<SpatialSorting::DESSBound*> dessGrid[nsdof];
    lArray1d validIndices;
  };
}

#endif /* DESS_H_ */
