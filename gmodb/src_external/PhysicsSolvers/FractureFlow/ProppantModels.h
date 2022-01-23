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
 * @file ProppantModels.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef PROPPANTMODELS_H_
#define PROPPANTMODELS_H_

#include "Common/typedefs.h"
#include "IO/ticpp/HierarchicalDataNode.h"

enum EffectiveViscosityModelType
{
  powerlaw = 0,
  einstein = 1,
  mooney = 2,
  maronpierce = 3,
  kriegerdougherty = 4,
  numEffectiveViscosityModels = 5
};

namespace EffectiveViscosityStrings{
  extern const std::string Powerlaw;
  extern const std::string Einstein;
  extern const std::string Mooney;
  extern const std::string MaronPierce;
  extern const std::string KriegerDougherty;
}


// override the fromString template for EffectiveViscosityModelType
template <>
inline EffectiveViscosityModelType fromString<EffectiveViscosityModelType>(std::string theString){
  if(theString==EffectiveViscosityStrings::Powerlaw){
  	return powerlaw;
  } else if(theString==EffectiveViscosityStrings::Einstein){
	return einstein;
  } else if(theString==EffectiveViscosityStrings::Mooney){
	return mooney;
  } else if(theString==EffectiveViscosityStrings::MaronPierce){
	return maronpierce;
  } else if(theString==EffectiveViscosityStrings::KriegerDougherty){
	return kriegerdougherty;
  } else {
    throw GPException("Error unrecognized effective viscosity type: " + theString +".");
  }

  // should never get here
  throw GPException("Error unrecognized effective viscosity type: " + theString +".");
  return EffectiveViscosityModelType::numEffectiveViscosityModels;
}

//class SlurryBase;
class ProppantBase;

class SlurryBase {

  public:
      SlurryBase();
      ~SlurryBase(){};
      virtual void ReadXML( TICPP::HierarchicalDataNode* hdn, ProppantBase& proppantData) ;
      virtual realT GetViscosity(realT vf);

      virtual realT HinderedSettlingSpeed(realT phi){
    	   return m_singleParticleSettlingSpeed*exp(m_hinderedSettlingCoefficient*phi);
      };
  private:
      realT m_beta;  /// effective viscosity power  mu = muo (1-vf/max_vf)^-beta
      realT m_muo;   /// fluid viscosity
      realT m_maxVf; /// max volume fraction

      realT m_singleParticleSettlingSpeed;
      realT m_hinderedSettlingCoefficient;

      EffectiveViscosityModelType m_effectiveViscosityModelType;

};

class ProppantBase
    {

  public:
      ProppantBase();
      ~ProppantBase(){};
      virtual void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;

      realT BuoyancyDeltaP(realT vf,realT dz){
    	  return -m_deltaRhoG*vf*dz;
      };

      realT PackPermeability(realT vf, realT mu);

      SlurryBase m_slurryModel;

      realT m_density;
      realT m_fluidDensity;
      realT m_rhoRatio;
      realT m_fluidViscosity;
      realT m_diameter;
      realT m_sphericity;
      realT m_gravity;
      realT m_maxVf;  // maximum volume fraction

  private:
      realT m_deltaRhoG; // (delta rho_p-delta rho_l)*g
      realT m_rholG;  // rho_l*g
      realT m_sphericityDiameterSqrdOn180; // squared sphericity*diameter/(180)

  friend class SlurryBase;
};

/// Base for erosion models
class ScourErosionBase
{
public:
	ScourErosionBase();
	~ScourErosionBase(){};
	virtual void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;

	// erosion rate for one surface of the fracture (will be multiplied by two later)
	virtual realT ErosionRate(realT tau, realT conc){  // erosion rate as a function of shear stress
	   realT rv(0.0);
	   if(tau > m_tau_o) rv =  m_K*std::pow(tau-m_tau_o,m_beta)/m_density;
       return rv;
	}

	 /*
	virtual realT VelocityBasedErosionRate(realT velocityMag, realT conc){  // erosion rate as a function of fluid velocity
		return m_K*conc*std::pow(velocityMag,m_beta);
		//K*c*v^beta
	}
	*/
private:
	realT m_K;  // erosion coefficent
	realT m_beta;  // powerlaw behavior
	realT m_tau_o; // shear stress threshold
    realT m_density; // density of solid
};

/// Base for erosion models due to flow from the matrix
class MatrixFluxErosionBase
{
public:
	MatrixFluxErosionBase();
	~MatrixFluxErosionBase(){};
	virtual void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;

	virtual realT ErosionRate(realT fluidFlux, realT phi); // Erosion rate as a function of fluid flux and matrix porosity


private:
	realT m_K;
	realT m_beta;
	realT m_lambda;

};

/// Near wall lift force models
class InertialShearLiftBase
{
public:
	InertialShearLiftBase();
	~InertialShearLiftBase(){};
	virtual void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;

	// F+ = K d+^\beta
	// F+ = F/\rho \nu^2
	// d+ = d u/ \nu
	realT LiftForce(realT shearVelocity){
		realT dPlus = m_d*shearVelocity/m_nu;
		return m_rho*m_nu*m_nu*m_K*pow(dPlus,m_beta);
	};

private:
	realT m_K;  //  coefficient
	realT m_beta; // power
	realT m_rho; // fluid density
	realT m_d; // particle diameter
    realT m_nu; // kinematic viscosity
};

// Proppant Forces new

class ProppantForceBase
{
public:
	ProppantForceBase();
	~ProppantForceBase(){};
	virtual void ReadXML( TICPP::HierarchicalDataNode* hdn ) {};

	virtual R1Tensor GetForce(R1Tensor& proppantVelocity, R1Tensor& FluidVelocity, R1Tensor& FaceNormal) = 0;
};

#endif /* PROPPANTMODELS_H_ */
