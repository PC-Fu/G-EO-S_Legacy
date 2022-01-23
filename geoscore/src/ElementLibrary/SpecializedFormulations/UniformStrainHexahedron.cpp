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
 * @file UniformStrainHexahedron.cpp
 * @author settgast1
 * @date Jun 6, 2011
 */

#include "UniformStrainHexahedron.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


static void CalculateShapeFunctionDerivative( const realT y[8] ,
                                  const realT z[8] ,
                                  realT b[8] ) ;


static void CalculateFBHourGlassModes( const Array1dT<R1Tensor>& xpos,
                                       const Array1dT<R1Tensor>& dNdx,
                                       realT hourgamma[4][8] );

static void
CalcFBHourForce( const Array1dT<R1Tensor>& vel,
                 const realT hourgamma[4][8],
                 const realT& dampcoef,
                 const realT& stiffcoef,
                 const realT& rho,
                 const realT& modulus,
                 const realT& volume,
                 const Array1dT<R1Tensor>& dNdx,
                 const realT& dt,
                 Array1dT<R1Tensor>& Qstiffness,
                 Array1dT<R1Tensor>& hgforce);


UniformStrainHexahedron::UniformStrainHexahedron():
FiniteElement<3>(1,8,4)
{
  m_nodeOrdering.resize(8);

  m_nodeOrdering[0] = 0;
  m_nodeOrdering[1] = 1;
  m_nodeOrdering[2] = 3;
  m_nodeOrdering[3] = 2;
  m_nodeOrdering[4] = 4;
  m_nodeOrdering[5] = 5;
  m_nodeOrdering[6] = 7;
  m_nodeOrdering[7] = 6;


}

UniformStrainHexahedron::~UniformStrainHexahedron()
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

void UniformStrainHexahedron::reinit(const std::vector<R1TensorT<3> > &mapped_support_points)
{
  assert(mapped_support_points.size() == n_dofs);

  const unsigned int q = 0;



  const realT x[8] = {  mapped_support_points[m_nodeOrdering[0]](0),
                        mapped_support_points[m_nodeOrdering[1]](0),
                        mapped_support_points[m_nodeOrdering[2]](0),
                        mapped_support_points[m_nodeOrdering[3]](0),
                        mapped_support_points[m_nodeOrdering[4]](0),
                        mapped_support_points[m_nodeOrdering[5]](0),
                        mapped_support_points[m_nodeOrdering[6]](0),
                        mapped_support_points[m_nodeOrdering[7]](0) };

  const realT y[8] = {  mapped_support_points[m_nodeOrdering[0]](1),
                        mapped_support_points[m_nodeOrdering[1]](1),
                        mapped_support_points[m_nodeOrdering[2]](1),
                        mapped_support_points[m_nodeOrdering[3]](1),
                        mapped_support_points[m_nodeOrdering[4]](1),
                        mapped_support_points[m_nodeOrdering[5]](1),
                        mapped_support_points[m_nodeOrdering[6]](1),
                        mapped_support_points[m_nodeOrdering[7]](1) };

  const realT z[8] = {  mapped_support_points[m_nodeOrdering[0]](2),
                        mapped_support_points[m_nodeOrdering[1]](2),
                        mapped_support_points[m_nodeOrdering[2]](2),
                        mapped_support_points[m_nodeOrdering[3]](2),
                        mapped_support_points[m_nodeOrdering[4]](2),
                        mapped_support_points[m_nodeOrdering[5]](2),
                        mapped_support_points[m_nodeOrdering[6]](2),
                        mapped_support_points[m_nodeOrdering[7]](2) };



  realT b[3][8];

  CalculateShapeFunctionDerivative(  y , z , b[0] ) ;
  CalculateShapeFunctionDerivative(  z , x , b[1] ) ;
  CalculateShapeFunctionDerivative(  x , y , b[2] ) ;

  /*
  for( int i=0 ; i<3 ; ++i )
  { 
    std::cout<<"i, b: "<<i;
    for( int a=0 ; a<8 ; ++a )
    {
      std::cout<<", "<<b[i][a];
    }
    std::cout<<std::endl;
  }
  */

  data[q].jacobian_determinant = 0.0;
  for( int a=0 ; a<8 ; ++a )
  {
    data[q].jacobian_determinant += x[a]*b[0][a] + y[a]*b[1][a] + z[a]*b[2][a] ;
  }
  data[q].jacobian_determinant /= 3.0;

  //std::cout<<"data[q].jacobian_determinant = "<<data[q].jacobian_determinant<<std::endl;


  for( int a=0 ; a<8 ; ++a )
  {
    data[q].mapped_gradients[a](0) = b[0][m_nodeOrdering[a]] / data[q].jacobian_determinant;
    data[q].mapped_gradients[a](1) = b[1][m_nodeOrdering[a]] / data[q].jacobian_determinant;
    data[q].mapped_gradients[a](2) = b[2][m_nodeOrdering[a]] / data[q].jacobian_determinant;

    //std::cout<<"data[q].mapped_gradients[a] = "<<data[q].mapped_gradients[a]<<std::endl;
  }
  //std::cout<<"data[q].jacobian_determinant = "<<data[q].jacobian_determinant<<std::endl;

}

void UniformStrainHexahedron::zero_energy_mode_control( const Array1dT<R1Tensor>& dNdx,
                                                        const realT& volume,
                                                        const Array1dT<R1Tensor>& x,
                                                        const Array1dT<R1Tensor>& vel,
                                                        const realT& dampcoef,
                                                        const realT& stiffcoef,
                                                        const realT& rho,
                                                        const realT& modulus,
                                                        const realT& dt,
                                                        Array1dT<R1Tensor>& Qstiffness,
                                                        Array1dT<R1Tensor>& force )
{
  realT hourgamma[4][8];

  CalculateFBHourGlassModes( x, dNdx, hourgamma );


  CalcFBHourForce( vel, hourgamma, dampcoef, stiffcoef, rho, modulus, volume, dNdx, dt, Qstiffness, force );


}




void CalculateShapeFunctionDerivative( const realT y[8] ,
                                       const realT z[8] ,
                                       realT b[8] )
{
  const realT y0 = y[0];
  const realT y1 = y[1];
  const realT y2 = y[2];
  const realT y3 = y[3];
  const realT y4 = y[4];
  const realT y5 = y[5];
  const realT y6 = y[6];
  const realT y7 = y[7];

  const realT z0 = z[0];
  const realT z1 = z[1];
  const realT z2 = z[2];
  const realT z3 = z[3];
  const realT z4 = z[4];
  const realT z5 = z[5];
  const realT z6 = z[6];
  const realT z7 = z[7];
  const realT twelfth = 1.0/12.0;

  b[0] = ( y1*((z5-z2)-(z3-z4))
          +y2*(z1-z3)
          +y3*((z2-z7)-(z4-z1))
          +y4*((z7-z5)-(z1-z3))
          +y5*(z4-z1)
          +y7*(z3-z4) )*twelfth;

  b[1] = ( y2*((z6-z3)-(z0-z5))
          +y3*(z2-z0)
          +y0*((z3-z4)-(z5-z2))
          +y5*((z4-z6)-(z2-z0))
          +y6*(z5-z2)
          +y4*(z0-z5) )*twelfth;

  b[2] = ( y3*((z7-z0)-(z1-z6))
          +y0*(z3-z1)
          +y1*((z0-z5)-(z6-z3))
          +y6*((z5-z7)-(z3-z1))
          +y7*(z6-z3)
          +y5*(z1-z6))*twelfth;

  b[3] = ( y0*((z4-z1)-(z2-z7))
          +y1*(z0-z2)
          +y2*((z1-z6)-(z7-z0))
          +y7*((z6-z4)-(z0-z2))
          +y4*(z7-z0)
          +y6*(z2-z7))*twelfth;

  b[4] = ( y7*((z3-z6)-(z5-z0))
          +y6*(z7-z5)
          +y5*((z6-z1)-(z0-z7))
          +y0*((z1-z3)-(z7-z5))
          +y3*(z0-z7)
          +y1*(z5-z0))*twelfth;

  b[5] = ( y4*((z0-z7)-(z6-z1))
          +y7*(z4-z6)
          +y6*((z7-z2)-(z1-z4))
          +y1*((z2-z0)-(z4-z6))
          +y0*(z1-z4)
          +y2*(z6-z1))*twelfth;

  b[6] = ( y5*((z1-z4)-(z7-z2))
          +y4*(z5-z7)
          +y7*((z4-z3)-(z2-z5))
          +y2*((z3-z1)-(z5-z7))
          +y1*(z2-z5)
          +y3*(z7-z2))*twelfth;

  b[7] = ( y6*((z2-z5)-(z4-z3))
          +y5*(z6-z4)
          +y4*((z5-z0)-(z3-z6))
          +y3*((z0-z2)-(z6-z4))
          +y2*(z3-z6)
          +y0*(z4-z3))*twelfth;
}


void CalculateFBHourGlassModes( const Array1dT<R1Tensor>& xpos,
                                const Array1dT<R1Tensor>& dNdx,
                                realT hourgamma[4][8] )
{
/*
  const realT gamma[4][8] =  { { 1,  1, -1, -1, -1, -1, 1,  1},
                               { 1, -1, -1,  1, -1,  1, 1, -1},
                               { 1, -1,  1, -1,  1, -1, 1, -1},
                               {-1,  1, -1,  1,  1, -1, 1, -1} };
*/
  const realT gamma[4][8] =  { { 1,  1, -1, -1, -1, -1,  1,  1},
                               { 1, -1,  1, -1, -1,  1, -1,  1},
                               { 1, -1, -1,  1,  1, -1, -1,  1},
                               {-1,  1,  1, -1,  1, -1, -1,  1} };


/*
  Array1dT<R1Tensor> b(8);

  b(0) = dNdx(0);
  b(1) = dNdx(1);
  b(2) = dNdx(3);
  b(3) = dNdx(2);
  b(4) = dNdx(4);
  b(5) = dNdx(5);
  b(6) = dNdx(7);
  b(7) = dNdx(6);*/

  // compute the hourglass modes
  for( int mode=0 ; mode<4 ; ++mode )
  {
    R1Tensor temp;
    R1Tensor hourmod;

    for( int a=0 ; a<8 ; ++a )
    {
      temp = xpos(a);
      temp *= gamma[mode][a];
      hourmod += temp;
    }

    for( int a=0 ; a<8 ; ++a )
    {
      hourgamma[mode][a] = gamma[mode][a] - Dot( dNdx(a), hourmod);
    }
  }
}


void
CalcFBHourForce( const Array1dT<R1Tensor>& vel,
                 const realT hourgamma[4][8],
                 const realT& dampcoef,
                 const realT& stiffcoef,
                 const realT& rho,
                 const realT& modulus,
                 const realT& volume,
                 const Array1dT<R1Tensor>& dNdx,
                 const realT& dt,
                 Array1dT<R1Tensor>& Qstiffness,
                 Array1dT<R1Tensor>& hgforce )
{


  R1Tensor q[4], Q;
  R1Tensor temp;

  const realT BB = Dot(dNdx[0],dNdx[0]) + Dot(dNdx[1],dNdx[1]) + Dot(dNdx[2],dNdx[2]) + Dot(dNdx[3],dNdx[3]) +
                   Dot(dNdx[4],dNdx[4]) + Dot(dNdx[5],dNdx[5]) + Dot(dNdx[6],dNdx[6]) + Dot(dNdx[7],dNdx[7]) ;

  const realT Cdamp  = dampcoef * sqrt( rho*modulus*BB / 6.0 ) * volume;

  const realT Cstiff = stiffcoef * modulus * BB * volume / 3.0 ;

  for( int mode=0 ; mode<4 ; ++mode )
  {
    q[mode] = 0.0;
    for( int a=0 ; a<8 ; ++a )
    {
      temp = vel(a);
      temp *= hourgamma[mode][a];
      q[mode] += temp;
    }

    q[mode] *= 1.0/sqrt(8.0);

    temp = q[mode];
    temp *= Cstiff * dt;
    Qstiffness[mode] += temp;

    q[mode] *= Cdamp;
    q[mode] += Qstiffness[mode];
    q[mode] *= 1.0/sqrt(8.0);
  }

  for( int a=0 ; a<8 ; ++a )
  {
    hgforce(a) = 0.0;
    for( int mode=0 ; mode<4 ; ++mode )
    {
      temp  = q[mode];
      temp *= hourgamma[mode][a];
      hgforce[a] -= temp;
    }
  }
}


