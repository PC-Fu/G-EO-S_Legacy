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
 * JointIntermediate.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"

#include "JointIntermediate.h"
#include <typeinfo>
#include <assert.h>

JointIntermediate::JointIntermediate( const int paramSize, const int stateSize ):
InterfaceBase( paramSize, stateSize )
{
  // TODO Auto-generated constructor stub

}

JointIntermediate::~JointIntermediate()
{

}

void
JointIntermediate::Initialize( const localIndex index,
                               const realT stressNormal,
                               const realT stressShear)
{
  JointIntermediateParameterData& matParams = *this->ParameterData(index);
  JointIntermediateStateData& matState = *this->StateData(index, 0);

  //Normal stress
  matState.stress = stressNormal;

  //Shear stress
  UpdateFriction(matParams, matState);
  matState.stressShear = stressShear;
  const realT strength = ShearStrength( matParams, matState);
  if(strength < fabs(matState.stressShear))
  {
    //At initialization, the system is necessarily at equilibrium, so we
    //need to set the joint at critical stress ... by increasing strength
    matState.mu = matParams.mu0 = fabs(matState.stressShear) / matParams.kshear;
  }
  matState.xs = -matState.stressShear / matParams.kshear;

  //Normal stress
  matState.normalApproach = NormalApproachAtInitialization(stressNormal,
                                                           matParams.kCoefficientElastic,
                                                           matParams.kplasticLoad,
                                                           matParams.kplasticUnload,
                                                           matState.kcurrent,
                                                           matState.gap,
                                                           matParams.normalApproachYield);
}

void
JointIntermediate::StrainDrivenUpdate( const localIndex index )
{
  const JointIntermediateParameterData& matParams = *this->ParameterData(index);
  JointIntermediateStateData& matState = *this->StateData(index, 0);

  //(1) evolve the normal stiffness state
  matState.kcurrent = NormalStiffnessRock(matParams.kCoefficientElastic, matParams.kplasticLoad, matParams.kplasticUnload,
            matState.normalApproach, matState.gap, matParams.normalApproachYield);
  //const realT normalApproachNext = matState.normalApproach + matState.dxndt * matState.dt;
  //const realT kNext = normalApproachNext > matState.normalApproach ? NormalStiffnessRock(matParams.kCoefficientElastic, matParams.kplasticLoad, matParams.kplasticUnload,
  //                                                                                   normalApproachNext, matState.gap, matParams.normalApproachYield) : matState.kcurrent;

  //(2) evolve the normal stress
  matState.stress = matState.kcurrent * (matState.normalApproach - matState.gap);

  //(3) evolve the normal stress due to dilation
  const realT dxs = matState.dxsdt * matState.dt;
  matState.stressDilation += DilationalStressIncrement( matParams, matState, dxs);
  matState.stress += matState.stressDilation;

  //(4) evolve the friction
  UpdateFriction(matParams, matState);

  //(5) evolve the shear stress
  matState.stressShear += -matParams.kshear * dxs;
  matState.xs += dxs;

  //(6) determine whether strength exceeded and return to failure surface if necessary
  ThresholdToFailureSurface(matParams, matState, dxs);

  matState.stressShearVector.Normalize();
  matState.stressShearVector *= matState.stressShear;
  return;
}

realT
JointIntermediate::StiffnessProjected(const localIndex index)
{
  const JointIntermediateParameterData& matParams = *this->ParameterData(index);
  JointIntermediateStateData& matState = *this->StateData(index, 0);
  const realT normalApproachEstimate = matState.normalApproach +
      (matState.dxndt > 0 ? matState.dxndt * matState.dt : 0.0);
  return NormalStiffnessRock(matParams.kCoefficientElastic, matParams.kplasticLoad, matParams.kplasticUnload,
           normalApproachEstimate, matState.gap, matParams.normalApproachYield);
}
void
JointIntermediate::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                   InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const JointIntermediateParameterData& matParams = static_cast < const JointIntermediateParameterData& > ( matParamsBase );
  JointIntermediateStateData& matState = static_cast < JointIntermediateStateData& > ( matStateBase );

  matState.mu = matParams.mu0;
}

realT
JointIntermediate::ShearStrength(const InterfaceBaseParameterData& ,
                                 InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  JointIntermediateStateData& matState = static_cast < JointIntermediateStateData& > ( matStateBase );

  return matState.stress * matState.mu;
}

realT
JointIntermediate::DilationalStressIncrement( const InterfaceBaseParameterData&,
                                              InterfaceBaseStateData& ,
                                              const realT ) const
{
  return 0.0;
}

realT
JointIntermediate::NormalApproachAtInitialization( const realT stressNormal,
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
    return normalApproach;

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
JointIntermediate::NormalStiffnessRock( const realT knCoefficientElastic,
                                        const realT knPlasticLoading,
                                        const realT knPlasticUnloading,
                                        const realT normalApproach,
                                        realT& normalGap,
                                        const realT normalApproachNormalYield) const
{
  if (normalApproach < 0.0)
    return 0.0;

  const bool noGap = isZero(normalGap);

  //if elastic, give the elastic stiffness, otherwise the plastic unloading stiffness
  if (normalApproach <= normalApproachNormalYield && noGap)
    return knCoefficientElastic * sqrt(normalApproach);

  //adjust the gap appropriately
  if (normalApproach > normalApproachNormalYield)
  {
    //handle case of simple linear contact
    if (isZero(knPlasticLoading - knPlasticUnloading) && noGap)
      return knPlasticLoading;
    //handle case of energy-injecting contact
    else if (knPlasticLoading > knPlasticUnloading)
      throw GPException(
          "Cannot have a kplasticLoad > kplasticUnload without injecting energy!");
    else if (isZero(knPlasticUnloading))
      throw GPException("Cannot have a zero slope plastic unloading curve");

    const realT fpu = knPlasticUnloading * (normalApproach - normalGap);

    realT fpl = 0.0;
    if (knPlasticLoading > 0.0)
    {
      fpl = knPlasticLoading * (normalApproach - normalApproachNormalYield);
      fpl += knCoefficientElastic * sqrt(normalApproachNormalYield) * normalApproachNormalYield;
    }
    else
    {
      fpl = knCoefficientElastic * sqrt(normalApproach);
      if (fpl > knPlasticUnloading)
        fpl = knPlasticUnloading;
      fpl *= normalApproach;
    }

    if (fpu > fpl)
    {
      normalGap = normalApproach - (fpl / knPlasticUnloading);
    }
  }
  return knPlasticUnloading;
}

realT
JointIntermediate::ShearStrength(const realT stress,
                                 const realT normalStressAtDilationLimit,
                                 const realT tanFrictionCoefficientInitial,
                                 const realT tanFrictionCoefficientResidual,
                                 const realT cohesion,
                                 const int ifail) const
{
  //note: this is the effect of Heuze and Itasca models
  //don't accept if the friction coefficient too close to 90-degrees
  if(tanFrictionCoefficientInitial > 1.0e5)
    return 0.0;

  // NOTE: Otis Walton (Private communication to Joe Morris: 22.FEB.2006) points out that
  // the shear strength should be determined by the SUM of elastic and viscous terms
  // The issue is that at high rate, viscous stresses will be a significant fraction of
  // the total stress and you can no longer argue that viscosity is simply a numerical
  // convenience
  if(isZero(tanFrictionCoefficientResidual - tanFrictionCoefficientInitial))
    return tanFrictionCoefficientInitial * stress + cohesion;

  realT shearStrength = 0.0;
  if (stress > 0.0)
  {
    if(ifail == 1)
    {
      shearStrength += tanFrictionCoefficientResidual * stress;
    }
    else
    {
      // We must figure where on the shear strength
      // envelope we are
      if (stress < normalStressAtDilationLimit)
      {
        // We are below critical normal stress
        // Use initial friction angle
        shearStrength += tanFrictionCoefficientInitial * stress;
      }
      else
      {
        // Initial portion
        shearStrength += tanFrictionCoefficientInitial * normalStressAtDilationLimit;
        // Portion past critical
        shearStrength += tanFrictionCoefficientResidual * (stress - normalStressAtDilationLimit);
      }
    }
  }
  //Add in the contribution due to cohesion
  shearStrength += cohesion;
  return shearStrength;
}
