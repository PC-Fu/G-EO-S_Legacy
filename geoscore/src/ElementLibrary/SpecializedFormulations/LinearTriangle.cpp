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
 * @file LinearTriangle.cpp
 * @author Fu, Pengcheng
 * @date July 3, 2012
 */

#include "LinearTriangle.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


LinearTriangle::LinearTriangle():
FiniteElement<2>(1,3,0)
{
  m_nodeOrdering.resize(3);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 2;

}

LinearTriangle::~LinearTriangle()
{
  // TODO Auto-generated destructor stub
}


/**
 * Reinitialize the finite element basis on a particular element.
 * We use the coordinates of the support points in real space to
 * construct the forward mapping from the parent coordinate system.  The
 * support points are assumed to follow a lexicographic ordering:
 * On the parent element, we loop over the x-coordinate fastest,
 * the y, then z (depending on the desired spatial dimension of the
 * element).
 */

void LinearTriangle::reinit(const std::vector<R1TensorT<3> > &mapped_support_points)
{

  assert(mapped_support_points.size() == n_dofs);

  //See Chapter 15 of Int. FEM of U Colorado by Carlos Felippa for detailed formulation.
  //http://www.colorado.edu/engineering/cas/courses.d/IFEM.d/
  //Accessed in July 2012

  const std::vector<R1TensorT<3> >& X = mapped_support_points;

  realT V;
  const realT half = 1.0 / 2.0;



  V = (X[1][0] * X[2][1] - X[2][0] * X[1][1]) + (X[2][0] * X[0][1] - X[0][0] * X[2][1]) + (X[0][0] * X[1][1] - X[1][0] * X[0][1]);
  V *= half;
  data[0].jacobian_determinant = V ;

  data[0].mapped_gradients[0](0) = (X[1][1] - X[2][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[0](1) = (X[2][0] - X[1][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[1](0) = (X[2][1] - X[0][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[1](1) = (X[0][0] - X[2][0]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[2](0) = (X[0][1] - X[1][1]) * half / data[0].jacobian_determinant;
  data[0].mapped_gradients[2](1) = (X[1][0] - X[0][0]) * half / data[0].jacobian_determinant;



}

