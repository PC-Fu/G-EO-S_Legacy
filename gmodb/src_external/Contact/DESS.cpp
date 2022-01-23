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
 * @file DESS.cpp
 * @date Nov 11, 2011
 * @author Scott Johnson
 *
 * @brief DESS spatial sorting class
 ************************************************************/

#include "DESS.h"
#include "Contact/SpatialSorterFactory.h"

namespace SpatialSorting
{

  DESS::DESS() : dessBounds()
  {
  }

  DESS::~DESS()
  {
    // TODO Auto-generated destructor stub
  }

  void DESS::Allocate(const localIndex size)
  {
    //allocate
    this->validIndices.reserve(size);
    const localIndex length = 2*size;
    this->dessBounds.resize(length);
    this->dessGrid[0].reserve(length);
    this->dessGrid[1].reserve(length);
    this->dessGrid[2].reserve(length);
  }

  void DESS::ClearLists()
  {
    this->validIndices.clear();
    this->dessBounds.clear();
    this->dessGrid[0].clear();
    this->dessGrid[1].clear();
    this->dessGrid[2].clear();
  }

  void DESS::AddBounds(const rArray1d & radii,
                      const Array1dT<R1Tensor> & centers,
                      const localIndex ii)
  {
    this->validIndices.push_back(ii);

    //lower
    {
      DESSBound& o = this->dessBounds[2*ii];
      o.Set(radii[ii], centers[ii], ii, true);
      this->dessGrid[0].push_back(&o);
      this->dessGrid[1].push_back(&o);
      this->dessGrid[2].push_back(&o);
    }
    //upper
    {
      DESSBound& o = this->dessBounds[2*ii+1];
      o.Set(radii[ii], centers[ii], ii, false);
      this->dessGrid[0].push_back(&o);
      this->dessGrid[1].push_back(&o);
      this->dessGrid[2].push_back(&o);
    }
  }

  void DESS::Populate(const rArray1d & radii,
                      const Array1dT<R1Tensor> & centers,
                      const int* const excludeFromContact)
  {
    ClearLists();

    //add the bound objects to the grid
    if(excludeFromContact != 0)
    {
      for(localIndex ii = 0; ii < radii.size(); ii++)
      {
        if(excludeFromContact[ii]>0)
          continue;
        AddBounds(radii, centers, ii);
      }
    }
    else
    {
      for(localIndex ii = 0; ii < radii.size(); ii++)
      {
        AddBounds(radii, centers, ii);
      }
    }
  }

  void DESS::Populate(const rArray1d & radii,
                      const Array1dT<R1Tensor> & centers,
                      const lSet& toResort)
  {
    //add the bound objects to the grid
    for(lSet::const_iterator it = toResort.begin(); it != toResort.end(); ++it)
    {
      const localIndex ii = *it;

      //lower
      {
        DESSBound& o = this->dessBounds[2*ii];
        o.Set(radii[ii], centers[ii], ii, true);
      }
      //upper
      {
        DESSBound& o = this->dessBounds[2*ii+1];
        o.Set(radii[ii], centers[ii], ii, false);
      }
    }
  }

  void DESS::WriteGrid()
  {
    for(localIndex dim = 0; dim < nsdof; ++dim)
    {
      std::cout << "**DIMENSION:" << (dim == 0 ? "X" : (dim == 1 ? "Y" : "Z")) << std::endl;
      for(localIndex i = 0; i < this->dessGrid[dim].size(); i++)
        std::cout << "(" << i << ") " << this->dessGrid[dim][i]->indexInInput << (this->dessGrid[dim][i]->isLower ? "L" : "U") << ", ";
      std::cout << std::endl;
    }
  }

  template<int T>
  void DESS::SortGridAxis()
  {
    //sort the axis
    std::sort(this->dessGrid[T].begin(), this->dessGrid[T].end(), DESSBound::Compare<T>);

    //assign the backward map
    for(localIndex ii = 0; ii < this->dessGrid[T].size(); ++ii)
      this->dessGrid[T][ii]->indexAlongAxis[T] = ii;
  }

  void DESS::SortGrid()
  {
    SortGridAxis<0>();
    SortGridAxis<1>();
    SortGridAxis<2>();
    //WriteGrid();
  }

  void DESS::DownSelectShapesInRange(const localIndex itgt,
                                     const lArray1d& order,
                                     Array1dT<lSet>& shapesInRange,
                                     const bool addLowerIndexCandidates)
  {
    //use the pivot to find the restricted list
    lSet& pivot = shapesInRange[order[0]];
    lSet::const_iterator it = pivot.begin(), tmp, candidate;
    while(it != pivot.end())
    {
      //get the candidate's iterator and advance the current iterator
      const localIndex current = *it;
      candidate = it;
      ++it;

      //determine whether to remove the candidate
      bool rm = current < itgt && !addLowerIndexCandidates;
      if(!rm)
      {
        for(localIndex i = 1; i < nsdof; i++)
        {
          tmp = shapesInRange[order[i]].find(current);
          if(tmp == shapesInRange[order[i]].end())
          {
            rm = true;
            break;
          }
        }
      }

      //if the removal flag is set, remove the candidate
      if(rm)
        pivot.erase(candidate);
    }
  }

  void DESS::SetShapesInRange(const DESSBound& lower,
                              const DESSBound& upper,
                              const localIndex dimension,
                              lSet& shapesInRangeAxis,
                              const bool addLeft)
  {
    //clear the slate before starting
    shapesInRangeAxis.clear();

    //get everything between the two indices
    if(addLeft)
    {
      for(localIndex i = lower.indexAlongAxis[dimension] + 1; i < upper.indexAlongAxis[dimension]; ++i)
        shapesInRangeAxis.insert(this->dessGrid[dimension][i]->indexInInput);
    }
    else
    {
      for(localIndex i = lower.indexAlongAxis[dimension] + 1; i < upper.indexAlongAxis[dimension]; ++i)
        if(this->dessGrid[dimension][i]->isLower)
          shapesInRangeAxis.insert(this->dessGrid[dimension][i]->indexInInput);
    }
  }

  void DESS::SetShapesInRange(const DESSBound& lower,
                              const DESSBound& upper,
                              const lArray1d& order,
                              Array1dT<lSet>& shapesInRange,
                              const bool addLeft)
  {
    SetShapesInRange(lower, upper, order[0], shapesInRange[order[0]], addLeft);
    SetShapesInRange(lower, upper, order[1], shapesInRange[order[1]], true);
    SetShapesInRange(lower, upper, order[2], shapesInRange[order[2]], true);

    //handle case where an owner is completely enclosed in an owned body's boundaries along a slave projection
    for(localIndex i = 0; i < validIndices.size(); i++)
    {
      if(validIndices[i] == lower.indexInInput)
        continue;
      const localIndex index = validIndices[i] * 2;
      for(localIndex j = 0; j < nsdof; j++)
        if(dessBounds[index].indexAlongAxis[order[j]] < lower.indexAlongAxis[order[j]] && dessBounds[index+1].indexAlongAxis[order[j]] > upper.indexAlongAxis[order[j]])
          shapesInRange[order[j]].insert(validIndices[i]);
    }
  }

  void DESS::SortSub(const localIndex itgt, const rArray1d& radii,
                     const Array1dT<R1Tensor>& centers, Array1dT<lArray1d>& neighborList,
                     Array1dT<lSet>& neighborListInverse, Array1dT<lSet>& shapesInRange,
                     const bool addLeft)
  {
    const localIndex ilower = 2*itgt;
    DESSBound& lower = this->dessBounds[ilower];
    DESSBound& upper = this->dessBounds[ilower+1];

    //get a decent pivot estimate
    lArray1d order(nsdof);
    order[0] = 0;
    order[1] = 1;
    order[2] = 2;
    //SetOrder(lower, upper, order);

    //set the shapes in range along each axis and sort according to index in the position list
    SetShapesInRange(lower, upper, order, shapesInRange, addLeft);

    //reduce the list of eligible shapes in-place
    DownSelectShapesInRange(itgt, order, shapesInRange);

    //now, the pivot dimension of shapesInRange contains the eligible neighbors
    {
      const lSet& curr = shapesInRange[order[0]];
      for(lSet::const_iterator it = curr.begin(); it != curr.end(); ++it)
        SpatialSorterBase::AddIfClose(itgt, *it,
                                  radii, centers,
                                  neighborList, neighborListInverse);
    }
  }

  bool DESS::Update(const rArray1d& radii, const Array1dT<R1Tensor>& centers, const lSet& toResort,
                    Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                    const int* const excludeFromSorting)
  {
    //first, let's remove neighbor list entries involving those that need to be resorted
    SpatialSorterBase::Remove(toResort, neighborList, neighborListInverse); //remove the moved from the neighbor lists

    Populate(radii, centers, toResort);

    //second, let's re-sort the lists ... probably overkill, but I'm lazy tonight (fixme)
    SortGrid();

    //third, we only care about re-establishing neighbors for those to be re-sorted
    Array1dT<lSet> shapesInRange(nsdof);
    for (lSet::const_iterator iter = toResort.begin(); iter != toResort.end();
        ++iter)
    {
      //get the lower and upper bounds of the current target
      SortSub(*iter, radii, centers,
              neighborList, neighborListInverse,
              shapesInRange, true);
    }

    SpatialSorterBase::RemoveDuplicates(neighborList);
    return true;
  }

  bool DESS::Sort(const rArray1d & radii,
                  const Array1dT<R1Tensor> & centers,
                  Array1dT<lArray1d> & neighborList,
                  Array1dT<lSet>& neighborListInverse,
                  const int* const excludeFromContact)
  {
    //clear the neighbor list
    for(localIndex ii = 0; ii < neighborList.size(); ii++)
    {
      neighborList[ii].clear();
      neighborListInverse[ii].clear();
    }

    //allocate the internal data structures
    Allocate(radii.size());

    //populate those data structures
    Populate(radii, centers, excludeFromContact);

    //sort the axes' lists; i.e., dessGrid[i] is now sorted for all i
    SortGrid();

    Array1dT<lSet> shapesInRange(nsdof);
    for(localIndex ii = 0; ii < validIndices.size(); ++ii)
    {
      //get the lower and upper bounds of the current target
      SortSub(validIndices[ii], radii, centers,
              neighborList, neighborListInverse, shapesInRange);
    }

    SpatialSorterBase::RemoveDuplicates(neighborList);
    return false;
  }
}

/// Register spatial sorter in the spatial sorter factory
REGISTER_SPATIALSORTER( DESS )
