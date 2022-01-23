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

 
#ifndef _KINEMATICS_H_
#define _KINEMATICS_H_

#include "Common/Common.h"
#include "assert.h"

//*****************************************************************************
//***** DECLARATIONS **********************************************************
//*****************************************************************************
void IncrementalKinematics( const R2TensorT<nsdof>& A ,
                            R2SymTensorT<nsdof>& Dadt ,
                            R2TensorT<nsdof>& Rhat );

void IncrementalRotation( const R2TensorT<nsdof>& A ,
                          R2TensorT<nsdof>& Rot );

inline void CalculateGradient( R2TensorT<nsdof>& Gradient ,
                        const int* bConnectivity,
                        const Array1dT<R1TensorT<nsdof> >& disp,
                        const Array1dT<R1TensorT<nsdof> >& dNdX )

{
  Gradient = 0.0;
  for( int a=0 ; a<8 ; ++a )
    Gradient.plus_dyadic_ab( disp(bConnectivity[a]) , dNdX(a));

}

inline void CalculateGradient( R2Tensor& Gradient ,
                        const Array1dT<R1Tensor >& disp,
                        const Array1dT<R1Tensor >& dNdX )

{
  //Gradient = 0.0;
  //for( int a=1 ; a<=8 ; ++a )
  //Gradient.plus_dyadic_ab( disp(bConnectivity[a-1]) , dNdX(a));

  assert( disp.size() == dNdX.size() );

  Gradient.dyadic_ab( disp(0) , dNdX(0) );
  for( Array1dT<R1Tensor >::size_type a=1 ; a<disp.size() ; ++a )
  {
    Gradient.plus_dyadic_ab( disp(a) , dNdX(a) );
  }
}

inline void CalculateGradient( R2Tensor& Gradient ,
                        const Array1dT<R1Tensor >& disp,
                        const R1Tensor* const dNdX )

{

  Gradient.dyadic_ab( disp(0) , dNdX[0] );
  for( Array1dT<R1Tensor >::size_type a=1 ; a<disp.size() ; ++a )
  {
    Gradient.plus_dyadic_ab( disp(a) , dNdX[a] );
  }
}



void CalculatePhantomGradient( R2TensorT<nsdof>& Gradient ,
                               const int* bConnectivity,
                               const Array1dT<R1TensorT<nsdof> >& disp,
                               const Array2dT<R1TensorT<nsdof> >& dNdX );


//*****************************************************************************

                                            


#endif
