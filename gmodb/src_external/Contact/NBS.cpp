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
 * @file NBS.cpp
 * @date Nov 11, 2011
 * @author Scott Johnson
 *
 * @brief NBS spatial sorting class (see Munjiza IJNME 1998)
 ************************************************************/

#include "NBS.h"
#include "Utilities/Utilities.h"
#include "Contact/SpatialSorterFactory.h"

namespace SpatialSorting
{

  NBS::NBS() : bins(), lists(), mask(), nbins(), lower(std::numeric_limits<realT>::max()), binSize(0.0)
  {
  }

  /**
   * @brief Called once at initialization to set the mask template
   * @description Note that the first 7 entries involve the target bin
   * @author Scott Johnson
   */
  void NBS::SetMask()
  {
    this->mask.Allocate(13, 2, nsdof);
    localIndex i = 0;
    // (1) 0,0,0  1,0,0
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 1; this->mask(i,1,1) = 0; this->mask(i,1,2) = 0;
    ++i;
    // (2) 0,0,0  1,1,0
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 1; this->mask(i,1,1) = 1; this->mask(i,1,2) = 0;
    ++i;
    // (3) 0,0,0  0,1,0
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 0; this->mask(i,1,1) = 1; this->mask(i,1,2) = 0;
    ++i;
    // (4) 0,0,0  0,0,1
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 0; this->mask(i,1,1) = 0; this->mask(i,1,2) = 1;
    ++i;
    // (5) 0,0,0  0,1,1
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 0; this->mask(i,1,1) = 1; this->mask(i,1,2) = 1;
    ++i;
    // (6) 0,0,0  1,1,1
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 1; this->mask(i,1,1) = 1; this->mask(i,1,2) = 1;
    ++i;
    // (7) 0,0,0  1,0,1
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 1; this->mask(i,1,1) = 0; this->mask(i,1,2) = 1;
    ++i;
    // (8) 0,1,0  1,0,0
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 1; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 1; this->mask(i,1,1) = 0; this->mask(i,1,2) = 0;
    ++i;
    // (9) 0,1,0  0,0,1
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 1; this->mask(i,0,2) = 0;
    this->mask(i,1,0) = 0; this->mask(i,1,1) = 0; this->mask(i,1,2) = 1;
    ++i;
    // (10) 1,0,1  0,1,0
    this->mask(i,0,0) = 1; this->mask(i,0,1) = 0; this->mask(i,0,2) = 1;
    this->mask(i,1,0) = 0; this->mask(i,1,1) = 1; this->mask(i,1,2) = 0;
    ++i;
    // (11) 0,0,1  1,1,0
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 1;
    this->mask(i,1,0) = 1; this->mask(i,1,1) = 1; this->mask(i,1,2) = 0;
    ++i;
    // (12) 0,0,1  1,0,0
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 0; this->mask(i,0,2) = 1;
    this->mask(i,1,0) = 1; this->mask(i,1,1) = 0; this->mask(i,1,2) = 0;
    ++i;
    // (13) 0,1,1  1,0,0
    this->mask(i,0,0) = 0; this->mask(i,0,1) = 1; this->mask(i,0,2) = 1;
    this->mask(i,1,0) = 1; this->mask(i,1,1) = 0; this->mask(i,1,2) = 0;
    ++i;
  }

  NBS::~NBS()
  {
  }

  /**
   * @brief Get the bin that the potential belongs with
   * @author Scott Johnson
   * @param[in] radii Radii of the potentials
   * @param[in] centers Centers of the potentials
   */
  bool NBS::SetLimits(const rArray1d& radii,
                      const Array1dT<R1Tensor>& centers,
                      const realT binFactor,
                      const int* const excludeFromContact)
  {
    //get extrema and bin size (as the minimum of the element size)
    R1Tensor upper(-std::numeric_limits<realT>::max());
    this->lower = upper;
    this->lower *= -1.0;
    this->binSize = 0;

    if(excludeFromContact != 0)
    {
      for(localIndex  i = 0; i < radii.size(); i++)
      {
        if(excludeFromContact[i]>0)
          continue;

        R1Tensor tmp = centers[i];
        tmp -= radii[i];
        this->lower.SetMin(tmp);
        upper.SetMax(tmp);
        tmp = 2*radii[i];
        this->binSize.SetMax(tmp);
      }
    }
    else
    {
      for(localIndex  i = 0; i < radii.size(); i++)
      {
        R1Tensor tmp = centers[i];
        tmp -= radii[i];
        this->lower.SetMin(tmp);
        upper.SetMax(tmp);
        tmp = 2*radii[i];
        this->binSize.SetMax(tmp);
      }
    }


    if( isZero(Dot(this->binSize, this->binSize)) )
      return false;
    this->binSize *= binFactor;

    //calculate number of bins
    upper -= this->lower;
    for( unsigned int i = 0; i < nsdof; i++)
    {
      realT tt = upper(i) / this->binSize(i);
      this->nbins[i] = ((localIndex)tt) + 1;
    }
    return true;
  }

  /**
   * @brief Bin indices of the potential
   * @author Scott Johnson
   * @param[in] radius Radius of the potential
   * @param[in] center Center of the potential
   * @param[out] bin Bin indices of the potential
   */
  void NBS::Bin(const realT radius, const R1Tensor& center, lArray1d& bin) const
  {
    R1Tensor tmp = center;
    tmp -= this->lower;
    tmp -= radius;
    for( unsigned int i = 0; i < nsdof; i++)
    {
      realT tt = tmp(i)/this->binSize(i);
      bin[i] = (localIndex)tt;
    }
  }

  /**
   * @brief Bin all of the potentials in the list
   * @author Scott Johnson
   * @param[in] radii Radii list
   * @param[in] centers Centers list
   */
  void NBS::Bin(const rArray1d& radii,
                const Array1dT<R1Tensor>& centers,
                Array1dT<lArray1d>& neighborList,
                const realT binFactor,
                const int* const excludeFromContact)
  {
    bins.clear();
    lists.clear();

    SetLimits(radii, centers, binFactor, excludeFromContact);
    lArray1d bin(nsdof, 0);

    if(excludeFromContact != 0)
    {
      for(localIndex i = 0; i < radii.size(); i++)
      {
        neighborList[i].clear();
        if(excludeFromContact[i]>0)
          continue;

        Bin(radii[i], centers[i], bin);
        AddToBin(i, bin);
      }
    }
    else
    {
      for(localIndex i = 0; i < radii.size(); i++)
      {
        neighborList[i].clear();
        Bin(radii[i], centers[i], bin);
        AddToBin(i, bin);
      }
    }
  }

  /**
   * @brief Add the given list index to the bin
   * @author Scott Johnson
   * @param[in] index Index to add to bin
   * @param[in] bin Indices of bin
   */
  void NBS::AddToBin(const localIndex index, lArray1d& bin)
  {
    lArray1d* cnts = BinContents(bin);
    if(cnts == 0) {
      lArray1d* tmp = CreateBin(bin);
      tmp->push_back(index);
    } else {
      cnts->push_back(index);
    }
  }

  /**
   * @brief Return the contents (list indices) of the requested bin
   * @author Scott Johnson
   * @param[in] bin bin indices
   * @return List of list indices associated with the bin
   */
  lArray1d* NBS::BinContents(lArray1d& bin)
  {
    //Dimension 0
    std::map<localIndex, std::map<localIndex, std::map<localIndex, lArray1d> > >::iterator iter0 = this->bins.find(bin[0]);
    if(iter0 == this->bins.end())
      return 0;

    //Dimension 1
    std::map<localIndex, std::map<localIndex, lArray1d> >::iterator iter1 = iter0->second.find(bin[1]);
    if(iter1 == iter0->second.end())
      return 0;

    //Dimension 2
    std::map<localIndex, lArray1d>::iterator iter2 = iter1->second.find(bin[2]);
    return iter2 != iter1->second.end() ? &iter2->second : 0;
  }

  /**
   * @brief Create the requested bin
   * @author Scott Johnson
   * @param[in] bin Bin to request
   * @return Pointer to list of indices
   */
  lArray1d* NBS::CreateBin(lArray1d& bin)
  {
    this->bins[bin[0]][bin[1]][bin[2]] = lArray1d();
    //add links to the other 7 lists here
    SetMask(bin);
    return &this->bins[bin[0]][bin[1]][bin[2]];
  }

  /**
   * @brief Connect the bin with the neighbors
   * @author Scott Johnson
   * @param[in] bin Bin to connect in to neighbors
   */
  void NBS::SetMask(lArray1d& bin  )
  {
    lArray1d tmp0(nsdof, 0);
    lArray1d tmp1(nsdof, 0);
    lArray1d tmpTgt(nsdof, 0);
//    lArray1d* s0, *s1, *st;

//    for(localIndex i = 0; i < mask.Dimension(0); i++)
//    {
//      //me as first entry
//      for(int j = 0; j < nsdof; j++)
//      {
//        tmpTgt[j] = bin[j] - mask(i,0,j);
//        tmp0[j] = bin[j];
//        tmp1[j] = tmpTgt[j] + mask(i,1,j);
//      }
//      st = this->BinContents(tmpTgt);
//      s0 = this->BinContents(tmp0);
//      s1 = this->BinContents(tmp1);
//      if(s0 != 0 && s1 != 0)
//      {
//        lists[bin[0]][bin[1]][bin[2]].push_back(Array1dT<Array1dT<lArray1d*> >());
//        Array1dT<Array1dT<lArray1d*> >& pairs = lists[bin[0]][bin[1]][bin[2]].back();
//        pairs.push_back(Array1dT<lArray1d*>());
//        pairs.back().resize(2,0);
//        pairs.back()[0] = s0;
//        pairs.back()[1] = s1;
//      }
//
//      //me as second entry
//      for(int j = 0; j < nsdof; j++)
//      {
//        tmpTgt[j] = bin[j] - mask(i,1,j);
//        tmp1[j] = bin[j];
//        tmp0[j] = tmpTgt[j] + mask(i,0,j);
//      }
//      st = this->BinContents(tmpTgt);
//      s0 = this->BinContents(tmp0);
//      s1 = this->BinContents(tmp1);
//      if(s0 != 0 && s1 != 0)
//      {
//        lists[bin[0]][bin[1]][bin[2]].push_back();
//        Array1dT<lArray1d*>& pair = lists[bin[0]][bin[1]][bin[2]].back();
//        pairs.resize(2,0);
//        pair[0] = s0;
//        pair[1] = s1;
//      }
//    }
  }

  bool NBS::Update(const rArray1d& radii, const Array1dT<R1Tensor>& x, const lSet& toResort,
                    Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                    const int* const excludeFromSorting)
  {
    return Sort(radii, x, neighborList, neighborListInverse, excludeFromSorting);
  }

  /**
   * @brief Sort the potentials into the given neighbor list
   * @author Scott Johnson
   * @param[in] radii List of potential radii
   * @param[in] center List of potential centers
   * @param[out] neighborList Neighbor list
   * @return return
   */
  bool NBS::Sort(const rArray1d& radii,
                 const Array1dT<R1Tensor>& centers,
                 Array1dT<lArray1d>& neighborList,
                 Array1dT<lSet>& neighborListInverse,
                 const int* const excludeFromSorting)
  {
    const realT binFactor = 1.1;

    Bin(radii, centers, neighborList,  binFactor, excludeFromSorting);

    std::map<localIndex, std::map<localIndex, std::map<localIndex, Array1dT<Array1dT<lArray1d*> > > > >::iterator it0;
    std::map<localIndex, std::map<localIndex, Array1dT<Array1dT<lArray1d*> > > > ::iterator it1;
    std::map<localIndex, Array1dT<Array1dT<lArray1d*> > >::iterator it2;
    Array1dT<Array1dT<lArray1d*> >::iterator it3;
    for(it0 = lists.begin(); it0 != lists.end(); ++it0)
    {
      for(it1 = it0->second.begin(); it1 != it0->second.end(); ++it1)
      {
        for(it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
        {
          for(it3 = it2->second.begin(); it3 != it2->second.end(); ++it3)
          {
            for( localIndex ii = 0; ii < it3->size(); ii++)
            {
              for( localIndex i = 0; i < (*it3)[ii][0].size(); i++) {
                localIndex ia = (*it3)[ii][0][i];
                for( localIndex j = 0; j < (*it3)[ii][1].size(); j++) {
                  const localIndex ib = (*it3)[ii][1][j];
                  SpatialSorterBase::AddIfClose(ia, ib, radii, centers, neighborList, neighborListInverse);
                }
              }
            }
          }
        }
      }
    }
    return false;
  }
}

/// Register spatial sorter in the spatial sorter factory
REGISTER_SPATIALSORTER( NBS )
