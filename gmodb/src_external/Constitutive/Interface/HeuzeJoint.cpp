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
 * HeuzeJoint.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"

#include "HeuzeJoint.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "Constitutive/Interface/InterfaceFactory.h"
#include <typeinfo>
#include <assert.h>

HeuzeJoint::HeuzeJoint( ):
JointIntermediate( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

HeuzeJoint::~HeuzeJoint()
{

}

void
HeuzeJoint::Initialize( const localIndex index,
                        const realT stressNormal,
                        const realT stressShear)
{
  //
  //NOTE: BE AWARE OF THE ASSUMPTIONS OF THIS PROCEDURE
  //
  //1) NORMAL DISPLACEMENT OCCURS BEFORE ANY SHEAR DISPLACEMENT
  //2) THE JOINT IS EQUILIBRIUM AT INIT ... THIS IS ENFORCED
  //3) THE JOINT HAS ONLY ELASTIC DILATION ... THIS IS ENFORCED

  HeuzeJointParameterData& matParams = *this->ParameterData(index);
  HeuzeJointStateData& matState = *this->StateData(index, 0);

  matState.stress = stressNormal;
  matState.stressShear = stressShear;

  //set shear strain
  matState.xs = -matState.stressShear / matParams.kshear;

  //set strength parameters, if necessary to keep it in equilibrium at initialization
  const realT xsmag = fabs(matState.xs);
  {
    if(xsmag >= matParams.xsResidual)
      matParams.xsResidual = xsmag * (1.0 + 1e-6);
    matState.mu = matParams.mu0;
    UpdateFriction(matParams, matState);
    const realT strength = ShearStrength( matParams, matState);
    if(strength < fabs(matState.stressShear))
    {
      const realT mu_new = fabs(matState.stressShear) / xsmag;
      matParams.muResidual += mu_new - matState.mu;
      matState.mu = mu_new;
    }
  }

  //now that we know sig_n, sig_s, and x_s, we need x_n
  //for now, let's force elastic dilation ...
  if(stressNormal > matParams.stressDilationLimit)
  {
    matParams.stressDilationLimit = stressNormal;
  }

  //1) Dilation elastic
  {
    const realT f1 = 0.4 * matParams.dilationCoefficient0 * xsmag * xsmag * sqrt(xsmag);
    const realT f3 = 2.0 / 3.0;

    //1a) Is this normal elastic? ... assume normal before
    realT xnelastic = 0;
    {
      //solve the cubic equation for x^1/2 -> elastic
      const realT f0 = -stressNormal / matParams.kCoefficientElastic;
      //solve for x^1/2 -> plastic

      const realT params[] = {f0, f1, f3};
      const realT xnelastic_sqrt = FindRoots::FindRoot(ZerosFunctionElasticElastic,
                                                       params,
                                                       0, 2 * matParams.normalApproachYield,
                                                       1e-6);
      xnelastic = xnelastic_sqrt * xnelastic_sqrt;
    }

    //check ...
    if(xnelastic <= matParams.normalApproachYield)
    {
      //Yes, it is still elastic
      matState.normalApproach = xnelastic;
    }
    else
    {
      //No, it is not elastic, so let's get the plastic part
      matState.normalApproach = matParams.kplasticLoad * f1;
      matState.normalApproach += f3 * matParams.normalApproachYield * sqrt(matParams.normalApproachYield);
      matState.normalApproach /= matParams.kplasticLoad;
      matState.normalApproach += matParams.normalApproachYield;
    }
  }
//  else
//  {
//    //So, now we know we may have a plastic dilation ... more complex
//
//  }
}
void
HeuzeJoint::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                   InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const HeuzeJointParameterData& matParams = static_cast < const HeuzeJointParameterData& > ( matParamsBase );
  HeuzeJointStateData& matState = static_cast < HeuzeJointStateData& > ( matStateBase );

  const realT xsmag = fabs(matState.xs);
  if(xsmag < matParams.xsResidual)
  {
    //not yet failed ...
    const realT slope = xsmag / (matParams.xsResidual-xsmag);
    const realT phi1 = (1.0 - slope) * (matState.mu - matParams.muResidual);
    matState.mu = matParams.muResidual + ( (phi1 < 0.0) ? 0.0 : phi1 );
  }
  else
  {
    //just failed ...
    matState.ifail = 1;
    matState.mu = matParams.muResidual;
  }
}

realT
HeuzeJoint::DilationalStressIncrement( const InterfaceBaseParameterData& matParamsBase,
                                       InterfaceBaseStateData& matStateBase,
                                       const realT dxs) const
{
  HeuzeJointStateData& matState = static_cast<HeuzeJointStateData&> (matStateBase);
  const HeuzeJointParameterData& matParams = static_cast<const HeuzeJointParameterData&> (matParamsBase);

  const realT xsmag = fabs(matState.xs);
  const realT dstress = matState.kcurrent * matParams.dilationCoefficient0 * sqrt(xsmag) * dxs;
  const realT stress0 = matState.stress + matState.stressDilation;
  const realT stress1 = stress0 + dstress;

  return stress1 < matParams.stressDilationLimit ? dstress :
      ((matParams.stressDilationLimit - stress0) > 0 ?
          matParams.stressDilationLimit - stress0 :
          0.0);
}

realT
HeuzeJoint::NormalApproachAtInitialization( const realT stressNormal,
                                            const realT knCoefficientElastic,
                                            const realT knPlasticLoading,
                                            const realT knPlasticUnloading,
                                            realT& kcurrent,
                                            realT& normalGap,
                                            const realT normalApproachNormalYield) const
{
  //Let's make the assumption that the loading monotonically increasing in time before initialization
  //otherwise, must assume that deformation has been elastic in the past ... if not, we're indeterminate
  if(stressNormal <= 0.0)
    return 0.0;

  realT normalApproach = pow(stressNormal / knCoefficientElastic, 2.0/3.0);
  kcurrent = knCoefficientElastic * sqrt(normalApproach);

  if(normalApproach <= normalApproachNormalYield)
  {
    return normalApproach;
  }

  if (knPlasticLoading > knPlasticUnloading)
    throw GPException(
        "Cannot have a kplasticLoad > kplasticUnload without injecting energy!");
  else if (isZero(knPlasticUnloading))
    throw GPException("Cannot have a zero slope plastic unloading curve");

  //get the solution
  if (knPlasticLoading > 0.0)
  {
    normalApproach = stressNormal;
    normalApproach -= knCoefficientElastic * sqrt(normalApproachNormalYield) * normalApproachNormalYield;
    normalApproach /= knPlasticLoading;
    normalApproach += normalApproachNormalYield;
  } //else, normal approach was already calculated correctly!

  const realT fpu = knPlasticUnloading * (normalApproach - normalGap);
  if (fpu > stressNormal)
    normalGap = normalApproach - (stressNormal / knPlasticUnloading);

  kcurrent = knPlasticUnloading;
  return normalApproach;
}

realT
HeuzeJoint::ZerosFunctionElasticElastic(const realT x, const realT* params )
{
  return params[2] * x * x * x + params[1]  * x + params[0];
}

/// Register class in the class factory
REGISTER_INTERFACE( HeuzeJoint )
