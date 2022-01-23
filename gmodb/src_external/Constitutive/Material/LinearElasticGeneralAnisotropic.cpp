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


/*
 * LinearElasticGeneralAnisotropic.cpp
 *
 *  Created on: Sun Jun 29 15:11:06 PDT 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"

#include "LinearElasticGeneralAnisotropic.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "Constitutive/Material/MaterialFactory.h"
#include <typeinfo>
#include <assert.h>

LinearElasticGeneralAnisotropic::LinearElasticGeneralAnisotropic( ):
MaterialBase( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

LinearElasticGeneralAnisotropic::~LinearElasticGeneralAnisotropic()
{

}

void LinearElasticGeneralAnisotropic::StrainDrivenUpdateMember( const localIndex index0,
                                                                const localIndex index1,
                                                                const R2SymTensorT<3>& Ddt,
                                                                const R2TensorT < 3 >& L,
                                                                const R2Tensor& Rot,
                                                                const realT dt)
{
  const localIndex paramIndex = m_parameterData.size() > 1 ? index0 : 0;
  const LinearElasticGeneralAnisotropicParameterData& matParams = m_parameterData[paramIndex];
  LinearElasticGeneralAnisotropicStateData& matState = m_stateData(index0, index1);

  realT& pressure = matState.pressure;
  R2SymTensorT<3>& devStress = matState.devStress;


  // use devstress to store stress
  devStress.PlusIdentity( pressure );
  realT* const stress = devStress.Data();
  const realT* D = Ddt.Data();

  stress[0] += matParams.c11*D[0] + matParams.c16*D[1] + matParams.c12*D[2] + matParams.c15*D[3] + matParams.c14*D[4] + matParams.c13*D[5];
  stress[2] += matParams.c12*D[0] + matParams.c26*D[1] + matParams.c22*D[2] + matParams.c25*D[3] + matParams.c24*D[4] + matParams.c23*D[5];
  stress[5] += matParams.c13*D[0] + matParams.c36*D[1] + matParams.c23*D[2] + matParams.c35*D[3] + matParams.c34*D[4] + matParams.c33*D[5];
  stress[4] += matParams.c14*D[0] + matParams.c46*D[1] + matParams.c24*D[2] + matParams.c45*D[3] + matParams.c44*D[4] + matParams.c34*D[5];
  stress[3] += matParams.c15*D[0] + matParams.c56*D[1] + matParams.c25*D[2] + matParams.c55*D[3] + matParams.c45*D[4] + matParams.c35*D[5];
  stress[1] += matParams.c16*D[0] + matParams.c66*D[1] + matParams.c26*D[2] + matParams.c56*D[3] + matParams.c46*D[4] + matParams.c36*D[5];

  pressure = devStress.Trace() / 3.0 ;
  devStress.PlusIdentity( -pressure );
  matState.RotateState(Rot);
  return;
}

void LinearElasticGeneralAnisotropic::StrainDrivenUpdateMember( const localIndex index0,
                                                                const localIndex index1,
                                                                const R2SymTensorT<3>& Ddt,
                                                                const R2TensorT < 3 >& L,
                                                                const R2Tensor& Rot,
                                                                const realT& volume_n,
                                                                const realT& volume_np1,
                                                                const realT dt)
{
  LinearElasticGeneralAnisotropicStateData& matState = m_stateData(index0, index1);

  realT& pressure = matState.pressure;
  R2SymTensorT<3>& devStress = matState.devStress;

  const realT trDdt = Ddt.Trace();

  realT StressPowerIncrement = 0.5 * (Dot(devStress, Ddt) + pressure * trDdt) * volume_n;
  //  realT strainEnergyIncrement = - 0.5 * ( devStress.Inner() / (2*G) + pow(pressure,2)/K ) * volume_n;

  StrainDrivenUpdateMember( index0, index1, Ddt, L, Rot, dt);

  StressPowerIncrement += 0.5 * (Dot(devStress, Ddt) + pressure * trDdt) * volume_np1;
  //  strainEnergyIncrement += 0.5 * ( devStress.Inner() / (2*G) + pow(pressure,2)/K ) * volume_np1;

  matState.ElasticStrainEnergy += StressPowerIncrement; //* 0.5 * ( volume_n + volume_np1);
  matState.StressPower += StressPowerIncrement; //* 0.5 * ( volume_n + volume_np1 );

  matState.RotateState(Rot);
  return;
}


void LinearElasticGeneralAnisotropic::MeanPressureDevStress( const localIndex index, realT& pressure,
                                                      R2SymTensor& devStress) const
{
  MeanPressureDevStressFromDerived<LinearElasticGeneralAnisotropic>(index, pressure, devStress);
}


/// Register class in the class factory
REGISTER_MATERIAL( LinearElasticGeneralAnisotropic )
