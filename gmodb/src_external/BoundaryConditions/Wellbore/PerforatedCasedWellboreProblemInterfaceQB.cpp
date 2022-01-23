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
/**
 * @file PerforatedCasedWellboreProblemInterfaceQB.cpp
 * @author johnson346
 * @date July 4, 2014
 */

#include "PerforatedCasedWellboreProblemInterfaceQB.h"

bool PerforatedCasedWellboreProblemInterfaceQB::computeF(const Epetra_Vector& x,
                                                         Epetra_Vector& f,
                                                         NOX::Epetra::Interface::Required::FillType F)
{
  //
  //Solve the system: assume incompressible flow in a simple well-bore, where perforations are clustered
  //
  //The system is 3*N-1 unknowns (N=number of perforation clusters): IN ORDER
  //  N fluxes through the casing into the initiated fractures
  //  N-1 fluxes along the well-bore between clusters
  //  N pressures in the well-bore at the points where the clusters are placed
  //
  //The equations can be formulated as:
  // (1: i=1..N)     q_f,i + q_w,i - q_w,i-1
  // (2: i=1..N)     p_w,i - p_f,i - alpha_i * |q_f,i| * q_f,i
  // (3: 1..N-1) p_w,i+1 - p_w,i - beta_i * |q_w,i| * q_w,i
  //
  // x = {q_f:i=1..N,
  //      q_w:i=1..N-1,
  //      p_w:i=1..N}
  //
  const size_t nclusters = m_fracturePressures.size();
  for (localIndex i = 0; i < nclusters; i++)
  {
    const realT pwi = GetPw(i, nclusters, x);
    const realT pfi = m_fracturePressures[i];
    const realT qfi = GetQf(i, nclusters, x);
    const realT qwi = i < (nclusters-1) ? GetQw(i, nclusters, x) : 0;
    const realT qwim1 = i > 0 ? GetQw(i - 1, nclusters, x) : m_qwellbore;

#if 0
    //DEBUG
    std::cout << "pw[" << i << "]=" << pwi << " pf[" << i << "]=" << pfi << " qf[" << i << "]=" << qfi << " qw[" << i << "]=" << qwi;
    if(i < 1)
      std::cout << " qw[-1";
    else
      std::cout << " qw[" << (i-1);
    std::cout << "]=" << qwim1 << std::endl;
#endif

    //first N eqns
    // (1: i=1..N)     q_f,i + q_w,i - q_w,i-1
    f[i] = qfi + qwi - qwim1;

    //second N eqns
    // (2: i=1..N)     p_w,i - p_f,i - alpha_i * |q_f,i| * q_f,i
    f[i + nclusters] = pwi - pfi - m_alpha * fabs(qfi) * qfi;
    f[i + nclusters] /= m_alpha;

    //third N-1 eqns
    if (i < (nclusters-1))
    {
      // (3: 1..N-1) p_w,i+1 - p_w,i - beta_i * |q_w,i| * q_w,i
      const realT pwip1 = GetPw(i + 1, nclusters, x);
      f[i + 2 * nclusters] = pwip1 - pwi - m_beta[i] * fabs(qwi) * qwi;
      f[i + 2 * nclusters] /= m_alpha;
    }
  }

#if 0
  //DEBUG
  for (localIndex i = 0; i < (3*nclusters-1); i++)
    std::cout << " f[" << i << "] = " << f[i] << std::endl;
  throw GPException("DEBUG");
#endif

  return true;
}

void PerforatedCasedWellboreProblemInterfaceQB::Jacobian(const localIndex ieqn, const Epetra_Vector& x, rArray1d& values) const
{
  //The equations can be formulated as:
  // (1: i=1..N)     q_f,i + q_w,i - q_w,i-1
  // (2: i=1..N)     p_w,i - p_f,i - alpha_i * |q_f,i| * q_f,i
  // (3: 1..N-1) p_w,i+1 - p_w,i - beta_i * |q_w,i| * q_w,i
  //
  // x = {q_f:i=1..N,
  //      q_w:i=1..N-1,
  //      p_w:i=1..N}
  //
  //NOTE: d(x|x|)/dx = x*x/|x| + |x|
  const size_t nclusters = m_fracturePressures.size();

  const localIndex oqw = nclusters, opw = 2 * nclusters - 1;
  values = 0.0;

  if (ieqn < nclusters)
  {
    //////////////////////////////////////////
    //first, i<N eqns, N
    // (1: i=1..N)     q_f,i + q_w,i - q_w,i-1
    values[ieqn] = 1.0;       //dfi/dq_f,i
    if(ieqn < (nclusters - 1))
      values[ieqn + oqw] = 1.0; //dfi/dq_w,i
    if (ieqn > 0)
      values[ieqn + oqw - 1] = -1.0; //df1/dq_w,i-1
  }
  else
  {
    if (ieqn < 2*nclusters)
    {
      //////////////////////////////////////////
      //second set: i=[N,2N-1], N
      // (2: i=1..N)     p_w,i - p_f,i - alpha_i * |q_f,i| * q_f,i
      const localIndex iq = ieqn - nclusters;
      const realT qfi = GetQf(iq, nclusters, x);
      const realT adqabsq = PerforatedCasedWellboreProblemInterface::dxabsxdx(qfi);
      values[iq + opw] = 1.0 / m_alpha; //df2/dp_w,i
      values[iq] = -adqabsq;  //df2/dq_f,i
    }
    else
    {
      //////////////////////////////////////////
      //third set: i=[2N,3*N-2], N-1
      const localIndex iq = ieqn - 2 * nclusters;
      const realT qwi = GetQw(iq, nclusters, x);
      // (3: 1..N-1) p_w,i+1 - p_w,i - beta_i * |q_w,i| * q_w,i
      values[iq + oqw] = -m_beta[iq] * PerforatedCasedWellboreProblemInterface::dxabsxdx(qwi) / m_alpha; //df3/dq_w,i
      values[iq + opw] = -1.0/m_alpha; //df2/dp_w,i
      values[iq + opw + 1] = 1.0/m_alpha; //df2/dp_w,i+1
    }
  }

#if 0
    //DEBUG
    std::cout << "J[" << ieqn << "] = [";
    for(rArray1d::const_iterator it = values.begin(); it != values.end(); ++it)
      std::cout << " " << *it;
    std::cout << " ]" << std::endl;
#endif
}
