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
 * @file NBS.h
 * @date Nov 11, 2011
 * @author Scott Johnson
 *
 * @brief NBS spatial sorting class (see Munjiza IJNME 1998)
 ************************************************************/

#ifndef NBS_H_
#define NBS_H_

#include "Contact/SpatialSorterBase.h"

namespace SpatialSorting
{
  class NBS: public SpatialSorterBase
  {
  public:
    NBS();
    ~NBS();

    static std::string SpatialSorterName() { return "NBS"; }

    void SetMask();

    virtual bool Sort(const rArray1d& radii, const Array1dT<R1Tensor>& x,
                      Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                      const int* const excludeFromSorting = 0);

    virtual bool Update(const rArray1d& radii, const Array1dT<R1Tensor>& x, const lSet& toResort,
                        Array1dT<lArray1d>& neighborList, Array1dT<lSet>& neighborListInverse,
                        const int* const excludeFromSorting = 0);
  protected:
    virtual bool SetLimits(const rArray1d& radii, const Array1dT<R1Tensor>& centers,
                           const realT binFactor, const int* const excludeFromContact = 0);

    void Bin(const realT radius, const R1Tensor& center, lArray1d& bin) const;

    void Bin(const rArray1d& radii, const Array1dT<R1Tensor>& centers,
             Array1dT<lArray1d>& neighborList, const realT binFactor,
             const int* const excludeFromContact = 0);

    void AddToBin(const localIndex index, lArray1d& bin);

    lArray1d* BinContents(lArray1d& bin);

    lArray1d* CreateBin(lArray1d& bin);

    void SetMask(lArray1d& bin);

    ///List of indices associated with each bin
    std::map<localIndex, std::map<localIndex, std::map<localIndex, lArray1d> > > bins;

    ///Set of list (pointer) pairs associated with each bin
    std::map<localIndex,
        std::map<localIndex, std::map<localIndex, Array1dT<Array1dT<lArray1d*> > > > > lists;

    ///Mask template associated with any bin (constant)
    Array3dT<localIndex> mask;

    ///Number of bins in each Cartesian direction
    lArray1d nbins;

    ///Lowest point in the domain and size of the bin in each Cartesian direction
    R1Tensor lower, binSize;
  };
}

#endif /* NBS_H_ */
