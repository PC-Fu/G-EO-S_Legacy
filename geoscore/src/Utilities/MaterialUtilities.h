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
/*
 * MaterialUtilities.h
 *
 *  Created on: Nov 27, 2013
 *      Author: walsh24
 */

#ifndef MATERIAL_UTILITIES_H_
#define MATERIAL_UTILITIES_H_

#include "Common/typedefs.h"
#include <sys/resource.h>
#include "Common/GPException.h"

#include <map>

/*
 * FillLinearElasticModuli
 *
 *
 * @param K bulk Modulus
 * @param G Shear Modulus
 * @param E Youngs Modulus
 * @param nu Poissions ratio
 * @param lame Lames constant
 * @param M P-wave modulus
 * @return error code 1 or 2 if failed
 */
inline
int FillLinearElasticModuli(realT& K, realT& G, realT& E, realT& nu,realT& lame, realT& M, const bool forceRecalculate = false){

  if(K <= 0){
    if(G>0 && E >0){
      K = E*G/(3.0*(3.0*G-E));
    } else if (G>0  &&  nu >0){
        K = 2.0*G*(1.0+nu)/(3.0*(1.0-2.0*nu));
    } else if (G>0  &&  lame >0){
        K = lame + 2.0*G/3.0;
    } else if (G>0  &&  M >0){
       K = M- 4.0*G/3.0;
    } else if(E>0  &&  nu >0){
      K = E/(3.0*(1.0-2.0*nu));
    } else if (lame>0  &&  nu >0){
      K = lame*(1.0+nu)/(3.0*nu);
    } else {
      return 1;
    }
  }


  if(G<=0){
    if(E>0){
      G = 3.0*K*E/(9.0*K-E);
    } else if(lame > 0){
      G = 3.0*(K-lame)/2.0;
    } else if(nu > 0){
      G = 3.0*K*(1.0-2*nu)/(2.0*(1+nu));
    } else if(M > 0){
      G = 0.75*(M-K);
    } else {
      return 2;
    }
  }


  //K & G are set - fill in the others (but only if not already set to avoid roundoff)
  if(forceRecalculate)
    E = nu = lame = M = 0;

  if(E<=0)
    E = 9*K*G/(3.0*K+G);
  if(nu<=0)
    nu = (3*K-2*G)/(2.0*(3.0*K+G));
  if(lame<=0)
    lame = K - 2.0*G/3.0;
  if(M<=0)
    M = K + 4.0*G/3.0;

  return 0;
}



#endif /* UTILITIES_H_ */
