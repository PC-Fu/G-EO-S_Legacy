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
/******************************************************************************
 * R2SymTensorT.h - Rank 2 Symmetric Tensor Class
 *
 * created: RRS (01/18/2002)
 *****************************************************************************/

#ifndef _R4_SYM_TENSOR_T_H_
#define _R4_SYM_TENSOR_T_H_

#include "TensorBaseT.h"

template<int N>
  struct SymSize {
    enum {
      value = N + SymSize<N - 1>::value
    };
  };

template<>
  struct SymSize<1> {
    enum {
      value = 1
    };
  };

template<int T_dim>
  class R2SymTensorT;

//*****************************************************************************
//***** R4SymTensorT Declaration **********************************************
//*****************************************************************************
template<int T_dim>
  class R4SymTensorT : public TensorBaseT<SymSize<T_dim>::value> {
  public:
    //**** CONSTRUCTORS AND DESTRUCTORS *****************************************
    R4SymTensorT(void);
    ~R4SymTensorT(void);
    R4SymTensorT(const R4SymTensorT<T_dim>& rhs);

    //***** ASSIGNMENT OPERATORS **************************************************
    R4SymTensorT<T_dim>&
    operator=(const int& rhs);
    R4SymTensorT<T_dim>&
    operator=(const realT& rhs);
#if __LONG_real
    R4SymTensorT<T_dim>& operator=( const realT& rhs )
    { operator=(static_cast<realT>(rhs));
      return *this;}
#endif

    R4SymTensorT<T_dim>&
    operator=(const R4SymTensorT<T_dim>& rhs);

    //***** ACCESS OPERATORS ****************************************************
    inline realT
    operator()(const int i, const int j) const;
    inline realT&
    operator()(const int i, const int j);

    //***** MULTIPLICATION OPERATIONS *******************************************


    //****** TENSOR OPERATIONS **************************************************


    friend class R2TensorT<T_dim> ;

  private:
    R4SymTensorT(R4SymTensorT<T_dim>&);

  };


template<int T_dim>
  R4SymTensorT<T_dim>::R4SymTensorT(void) :
    TensorBaseT<SymSize<T_dim>::value> ()
  {
  }

template<int T_dim>
  R4SymTensorT<T_dim>::~R4SymTensorT(void)
  {
  }

template<int T_dim>
  R4SymTensorT<T_dim>::R4SymTensorT(const R4SymTensorT<T_dim>& rhs) :
    TensorBaseT<SymSize<T_dim>::value> ()
  {
    TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
  }

//***** ASSIGNMENT OPERATORS **************************************************

// Assigns all components to an integer
template<int T_dim>
  R4SymTensorT<T_dim>&
  R4SymTensorT<T_dim>::operator=(const int& rhs)
  {
    TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
    return *this;
  }

// Assigns all components to a realT
template<int T_dim>
  R4SymTensorT<T_dim>&
  R4SymTensorT<T_dim>::operator=(const realT& rhs)
  {
    TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
    return *this;
  }

// Assigns all components to another TensorBaseT's (Copy Constructor)
template<int T_dim>
  R4SymTensorT<T_dim>&
  R4SymTensorT<T_dim>::operator=(const R4SymTensorT<T_dim>& rhs)
  {
    TensorBaseT<SymSize<T_dim>::value>::operator=(rhs);
    return *this;
  }

template<int T_dim>
  realT
  R4SymTensorT<T_dim>::operator()(const int i, const int j) const
  {
    int index = 0;
    int i_sym = i;
    int j_sym = j;

    if (j > i)
    {
      i_sym = j;
      j_sym = i;
    }

    for (int k = 1; k < i_sym; ++k)
      index += k;
    index += (j_sym - 1);

    return this->t_data[index];
  }

template<int T_dim>
  inline realT&
  R4SymTensorT<T_dim>::operator()(const int i, const int j)
  {
    int index = 0;
    int i_sym = i;
    int j_sym = j;

    if (j > i)
    {
      i_sym = j;
      j_sym = i;
    }

    for (int k = 1; k < i_sym; ++k)
      index += k;
    index += (j_sym - 1);

    return this->t_data[index];
  }

#endif

