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
 * Linearized.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"

#include "Linearized.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "Constitutive/Interface/InterfaceFactory.h"
#include <typeinfo>
#include <assert.h>

Linearized::Linearized( ):
HertzianIntermediate( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

Linearized::~Linearized()
{

}

void
LinearizedStateData::Update(const realT curvature1,
                             const realT curvature2)
{
  HertzianIntermediateStateData::Update(curvature1, curvature2);

  //calculate the yield force
  realT fy = 1.0 / 6.0;
  {
    //Hertz solution is F = pi^3*R^2/(6E^2) * p0^3
    //under von Mises and Tresca yield criteria the
    //yield occurs when the deviatoric stress reaches
    //the yield stress, which occurs where p0 = 1.60*y0
    realT tmp = 1.6 * yield * 3.141592653589793;
    tmp *= tmp * tmp;
    fy *= tmp;
    fy /= youngs * youngs;
    fy *= radius * radius;
  }

  //calculate the critical overlap (onset of plastic deformation)
  ac = pow(fy / hertzCf, 2.0 / 3.0);

  //calculate the elastic energy of the Hertzian contact up to plastic yield
  //REAL tmpS = 1.5*hertzCf*sqrt(ac);
  //f(a) = hcf*(ac^1.5) -> ep = integral(f(a)*da) = (2/5)hcf*(ac^2.5)
  ecrit = 0.4 * hertzCf * ac * ac * sqrt(ac);

  //differentiate Hertzian from linearized
  linearStiffnessElastic = 1.5 * fy / ac;
  if (false)
  {
    fcrit = fy;
  }
  else
  {
    //assume that elastic energy in both the Hertzian and linear is the same
    //f(a) = k*ac -> ep = integral(ka*da) = 0.5*k*ac^2
    ac = 4.0 * ecrit / (3.0 * fy);
    fcrit = linearStiffnessElastic * ac;
  }
}

void
LinearizedStateData::Initialize(const realT curvature1, const realT curvature2,
                                 const realT poissons1, const realT poissons2,
                                 const realT youngs1, const realT youngs2,
                                 const realT mass1, const realT mass2,
                                 const realT rest1, const realT rest2,
                                 const realT yield1, const realT yield2,
                                 const realT velHalf1, const realT velHalf2,
                                 const realT surfaceEnergy1, const realT surfaceEnergy2,
                                 const realT cement1, const realT cement2)
{
  youngs = EffectiveYoungsModulus(poissons1, poissons2, youngs1, youngs2);
  mass = EffectiveMass(mass1, mass2);
  coefficientOfRestitution2 = EffectiveCoefficientOfRestitutionSquared(youngs, rest1, rest2,
                                                                       poissons1, poissons2,
                                                                       youngs1, youngs2);
  yield = EffectiveYieldStrength(yield1, yield2);
  halfRestVel = HalfCoefficientOfRestitutionVelocity(youngs, velHalf1, velHalf2,
                                                     poissons1, poissons2,
                                                     youngs1, youngs2);
  //set cohesive parameters
  dupre = 0.0;
  if(cement1 > 0 && cement2 > 0)
    isCohesive = 2;
  else
  {
    isCohesive = (surfaceEnergy1 > 0 && surfaceEnergy2 > 0) ? 1 : 0;
    if(isCohesive > 0)
    {
      dupre = surfaceEnergy1 + surfaceEnergy2
          - 0.0;
      //0.0 should be replaced with interface energy
      //the following assumes a JKR model
      {
        //The following is the direct calculation using JKR
        //  REAL w2K2_two_thirds = dupre*dupre/(Eeff*Eeff);
        //  w2K2_two_thirds = pow(w2K2_two_thirds,2.0/3.0);
        //  dcoh0 = w2K2_two_thirds * pow(reff,2.0/3.0) * (7.08273 - 22.2956*w2K2_two_thirds);
        //...however, we would like to linearize the relationship while still dissipating the correct amount of energy
        //looking at the functional form of the JKR model, the energy in the hysteretic part of the force displacement
        //curve is ...
        dcoh0 = 0.0;
      }
    }
  }
  pullOffForce = PullOffForce(cement1, cement2);

  //set mechanical properties of the contact
  Update(curvature1, curvature2);
}

realT
LinearizedStateData::PullOffForce(const realT cement1,
                                  const realT cement2) const
{
  realT pof = 0.0;
  if (isCohesive > 0)
  {
    if (cement1 > 0 && cement2 > 0)
      pof = cement1 < cement2 ? -cement1 : -cement2;
    else
      pof = -1.5 * dupre * 3.1415926535897933 * radius;
  }
  return pof;
}


realT
LinearizedStateData::EffectiveCoefficientOfRestitutionSquared(const realT Eeff,
                                                              const realT rest1, const realT rest2,
                                                              const realT poissons1, const realT poissons2,
                                                              const realT youngs1, const realT youngs2) const
{
  realT ret = Eeff;
  ret *= rest1 * rest1 * (1 - poissons1 * poissons1) / youngs1
      + rest2 * rest2 * (1 - poissons2 * poissons2) / youngs2;
  return ret;
}

realT
LinearizedStateData::EffectiveYieldStrength(const realT yield1,
                                            const realT yield2) const
{
  return 2.0 / (1.0 / yield1 + 1.0 / yield2);
  //note that for similar yield strengths this would yield the same yield strength
}

realT
LinearizedStateData::HalfCoefficientOfRestitutionVelocity(const realT Eeff,
                                                          const realT velHalf1, const realT velHalf2,
                                                          const realT poissons1, const realT poissons2,
                                                          const realT youngs1, const realT youngs2) const
{
  return Eeff
      * (velHalf1 * velHalf1 * (1 - poissons1 * poissons1) / youngs1
          + velHalf2 * velHalf2 * (1 - poissons2 * poissons2) / youngs2);
}


realT
LinearizedStateData::min_dt(const realT stiffness,
                            const realT tmass) const
{
  realT keq = stiffness; //(1+1.5*rattc)*stiffness;
  realT dt_ = 0.0484 * tmass / keq; //denominator is k2 + shear_stiffness
  return dt_;
}


realT
LinearizedStateData::ThetaK(const realT Tk_old,
                            const R1Tensor& shearDirection,
                            const R1Tensor& tangentialForces,
                            const realT mu_,
                            const realT fn_mag,
                            const realT dfn_mag,
                            int& loading_) const
{
  int k = 0;
  realT Tk_ = Tk_old;
  if (Tk_old > 0.0)
  {
    //if the maximum tangential loading reached is non-zero,
    //this is either unloading or reloading; for unloading the
    //new shear direction will be in the opposite direction as the new one
    //for reloading
    loading_ = Dot(shearDirection, shearDirection) < 0 ? -loading_ : loading_;
    k = loading_ > 0 ? 2 : 1;
  }
  else
  {
    if (Dot(shearDirection, shearDirection) < 0)
    {
      loading_ = -loading_;
      Tk_ = tangentialForces.L2_Norm();
      k = 2;
    }
  }

  //now, find the ThetaK term
  realT tmp = tangentialForces.L2_Norm();
  switch (k)
  {
    case 0:
      tmp += mu_ * dfn_mag / fn_mag;
      break;
    case 1:
      tmp -= Tk_;
      tmp *= -1.0;
      tmp += 2.0 * mu_ * dfn_mag;
      tmp /= 2.0 * mu_;
      break;
    default:
      tmp -= Tk_;
      tmp += 2.0 * mu_ * dfn_mag;
      tmp /= 2.0 * mu_;
      break;
  }
  tmp = pow(1 - tmp, (1.0 / 3.0));
  return tmp;
}
realT
Linearized::NormalStiffness(const InterfaceBaseParameterData& ,
                            InterfaceBaseStateData& matStateBase,
                            const realT normalApproach,
                            const bool setForces) const
{
  LinearizedStateData& matState = static_cast < LinearizedStateData& > ( matStateBase );

  const realT youngs2 = matState.youngs * matState.youngs;

  realT fnmag_coh = 0.0;
  realT fnmag_mat = 0.0;

  //------COHESION------
  //if < a0, then we are (potentially) cohesive
  if (normalApproach < matState.a0)
  {
    //this is now out of contact with the other body ... unless it is cohesive
    realT coh_stiffness = std::numeric_limits<realT>::min();
    if (matState.isCohesive > 0)
    {
      //cohesive stiffness
      coh_stiffness = matState.at <= matState.ac ?
          matState.linearStiffnessElastic :
          PlasticUnloadingSlope(matState.youngs, youngs2,
                                matState.mass, matState.vark2,
                                matState.linearStiffnessElastic,
                                matState.at, matState.ac, matState.halfRestVel);
      //get the cohesion distance
      if (setForces)
      {
        realT dcoh0 = matState.pullOffForce / (-1.0 * coh_stiffness);
        if(normalApproach >= (matState.a0 - dcoh0))
        {
          fnmag_coh = coh_stiffness * (matState.a0 - dcoh0 - normalApproach);
          matState.ElasticStrainEnergy += 0.5 * fnmag_coh * fnmag_coh / coh_stiffness;
          matState.stress += fnmag_coh;
        }
      }
    }
    return coh_stiffness;
  }

  //------GRAIN LOADING/UNLOADING------
  //otherwise, contact is in a different regime ... note dcoh0 = 0 for !is_cohesive
  realT at; //at >= ac ... if > ac, then it is yielding
  if(setForces)
  {
    matState.at = matState.at < matState.ac ? matState.ac : matState.at;
    at = matState.at;
  }
  else
  {
    at = matState.at < matState.ac ? matState.ac : matState.at;
  }

  //previous update stored in the "store" state
  if (normalApproach > at)
  {
    //------POST-YIELD LOADING------
    at = normalApproach;
    const realT k2 = PlasticUnloadingSlope(matState.youngs, youngs2,
                                           matState.mass, matState.vark2,
                                           matState.linearStiffnessElastic,
                                           at, matState.ac, matState.halfRestVel);
    if(setForces)
    {
      matState.at = at;
      if ( !isEqual( k2, matState.linearStiffnessElastic ) )
        matState.a0 = (k2 - matState.linearStiffnessElastic) * at / k2;

      fnmag_mat = k2 * (normalApproach - matState.a0);
      matState.ElasticStrainEnergy += 0.5 * fnmag_mat * fnmag_mat / k2;
      //why was I adding in a cohesive force during post-yield loading????
      //fnmag_coh = matState.pullOffForce; //just a constant value here ... could alter this in the future
      matState.stress += fnmag_mat;//+fnmag_coh
    }
    return k2;
  }
  else if (at > matState.ac)
  {
    //------POST-YIELD UNLOADING------
    const realT k2 = PlasticUnloadingSlope(matState.youngs, youngs2,
                                           matState.mass, matState.vark2,
                                           matState.linearStiffnessElastic,
                                           at, matState.ac,
                                           matState.halfRestVel);
    if (setForces)
    {
      fnmag_mat = k2 * (normalApproach - matState.a0);
      matState.ElasticStrainEnergy += 0.5 * fnmag_mat * fnmag_mat / k2;
      matState.stress += fnmag_mat;
    }

    return k2;
  }
  else
  {
    //------ELASTIC LOADING/UNLOADING------
    if (setForces)
    {
      //still elastic
      fnmag_mat = matState.linearStiffnessElastic * normalApproach;
      matState.ElasticStrainEnergy += 0.5 * fnmag_mat * fnmag_mat / matState.linearStiffnessElastic;
      matState.stress += fnmag_mat;
    }
    return matState.linearStiffnessElastic;
  }
}

realT
Linearized::PlasticUnloadingSlope(const realT Eeff, const realT Eeff2, const realT tmass,
                                  const int vark2, const realT linearStiffnessElastic,
                                  const realT at, const realT ac, const realT velHalf) const
{
  if (vark2!=0)
  {
    //Get variable latching term
    const realT slope = 3 * sqrt(linearStiffnessElastic / tmass) / velHalf;

    //Calculate k2 with saturation at e = 0.1
    realT k2 = at > ac ? linearStiffnessElastic + slope * linearStiffnessElastic * (at - ac) : linearStiffnessElastic;

    if (k2 > 100.0 * linearStiffnessElastic) //saturate at 100*K1, e = sqrt(K1/K2) = 0.1
      k2 = 100.0 * linearStiffnessElastic;
    return k2;
  }
  else
  {
    return linearStiffnessElastic / Eeff2;
  }
}

/// Register class in the class factory
REGISTER_INTERFACE( Linearized )
