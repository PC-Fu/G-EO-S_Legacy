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
 * @file ProppantModels.cpp
 * @author walsh24
 * @date June 1, 2011
 */
#include "ProppantModels.h"

namespace EffectiveViscosityStrings{
  const std::string Powerlaw = "Powerlaw";
  const std::string Einstein = "Einsten";
  const std::string Mooney = "Mooney";
  const std::string MaronPierce = "MaronPierce";
  const std::string KriegerDougherty = "KriegerDougherty";
}


namespace{
    const realT TINY = 1e-64;

    // dp = particle diameter
    // rho_l = liquid density
    // v = velocity
    // mu = viscosity
    realT ParticleReynoldsNumber(realT dp, realT rho_l, realT v, realT mu){
      realT Rep = dp * rho_l*fabs(v)/(mu+TINY);
	  return Rep;
    }

    realT SettlingVelocityFromDragCoefficient(realT Cd,realT rho_s,realT rho_l,realT g, realT d){
      realT vs = sqrt( 4.0/3.0 * fabs(rho_s-rho_l)*g*d/((Cd+TINY)*rho_l) );
      return vs;
    }

	// Schiller Nauman Drag Coefficient
	realT SchillerNaumanDragCoefficient(realT Rep){
		realT Cd = 0.44;
		if(Rep < 1000.0 ) Cd= (24.0/Rep)*(1.0 + 0.15*pow(Rep,0.687) );
		return Cd;
	}

	//Ishii Mishima 1984 DragCoefficient (for Re below 1000)
	realT IshiiMishimaDragCoefficient(realT Rep){
		realT Cd = 0.44; // Schiller Nauman for Re above 1000
		if(Rep < 1000.0 ) Cd= (24.0/Rep)*(1.0 + 0.1*pow(Rep,0.75) );

		// above 1000 should be
		// Cd = 0.45*(1+17.67*pow(f,6/7) )/18.67*f
		// where f = sqrt(1 - ap) mu_c/mu_m
		// ap = particle volume fraction, mu_m = mixed fluid viscosity
		return Cd;
	}

	//Morrison Drag Coefficient (~standard drag coefficient) "Data Correlation for Drag Coefficient for Sphere" 2010
	realT MorrisonDragCoefficient(realT Rep){
		realT Cd = 24/Rep;
		if(Rep > 2.0 ) Cd +=
				   2.6*(Rep/5.0)/( 1+pow(Rep/5.0,1.52)) +
				   0.411*pow(Rep/263000,-7.94)/( 1+pow(Rep/263000,-8.0) ) +
				   pow(Rep,0.80)/ 461000;
		return Cd;
	}

}


ProppantBase::ProppantBase(){
/*Empty*/
}

void ProppantBase::ReadXML( TICPP::HierarchicalDataNode* proppantNode ){

  m_density    = proppantNode->GetAttributeOrDefault("density","2.5 g/cc");
  m_fluidDensity    = proppantNode->GetAttributeOrDefault("fluidDensity","1.0 g/cc");
  m_fluidViscosity   = proppantNode->GetAttributeOrDefault("fluidViscosity","1.0e-3 Pa.s");
  m_diameter    = proppantNode->GetAttributeOrDefault("diameter","200 um");
  m_sphericity    = proppantNode->GetAttributeOrDefault("sphericity","1");
  m_gravity   = proppantNode->GetAttributeOrDefault("gravity","9.8 m/s^2");
  m_maxVf   = proppantNode->GetAttributeOrDefault("maxPackingFraction","0.62");

  m_deltaRhoG = m_gravity*(m_density-m_fluidDensity);
  m_rholG = m_gravity*m_fluidDensity;
  m_sphericityDiameterSqrdOn180 = pow(m_sphericity*m_diameter,2)/180.0;
  m_rhoRatio = m_fluidDensity/m_density;

  TICPP::HierarchicalDataNode* slurryNode = proppantNode->GetChild("SlurryModel");
  if(!slurryNode){
    throw GPException("Error ProppantData: SlurryModel not provided.");
  } else {
    m_slurryModel.ReadXML(slurryNode,*this);
  }
}


realT ProppantBase::PackPermeability(realT vf, realT mu){
	 return m_sphericityDiameterSqrdOn180*(vf*vf*vf)/(mu*(1-vf*vf));
	 // Kozeny Carmen permeability
}


//////////////////////////


SlurryBase::SlurryBase(){

}

realT SlurryBase::GetViscosity(realT vf){
	// mu = muo(1-c/c_max)^-\beta
	//  return (vf < m_maxVf)? m_muo*pow(1.0 -  vf/m_maxVf,-m_beta) : 1e10;
	realT rv = 0.0;
	switch (m_effectiveViscosityModelType) {
	case powerlaw:
		rv = (vf < m_maxVf)? m_muo*pow(1.0 -  vf/m_maxVf,-m_beta) : 1e10;
		break;
	case einstein:
		rv = m_muo*(1.0 +2.5* vf);
		break;
	case mooney:
		rv = (vf < m_maxVf)? m_muo*exp(2.5*vf/(1.0 -  vf/m_maxVf)): 1e10;
		break;
	case maronpierce:
	    {
		realT xx = (1.0 -  vf/m_maxVf);
		rv = (vf < m_maxVf)? m_muo/(xx*xx): 1e10;
		break;
   	    }
	case kriegerdougherty:
		rv = (vf < m_maxVf)?m_muo*pow(1.0 -  vf/m_maxVf,-2.5*m_maxVf): 1e10 ;
		break;
	default :
		throw GPException("Error SlurryBase: Unsupported effective viscosity model.");
	}

	return rv;

}

void SlurryBase::ReadXML( TICPP::HierarchicalDataNode* slurryNode, ProppantBase& proppantData){
    std::string singleParticleModelType = slurryNode->GetAttributeStringOrDefault("singleParticleSettlingModel","stokes");
    m_hinderedSettlingCoefficient = slurryNode->GetAttributeOrDefault("hinderedSettlingCoefficient","-5.9");
    std::string effectiveViscosityModelType = slurryNode->GetAttributeStringOrDefault("effectiveViscosityModel","powerlaw");
    m_effectiveViscosityModelType = fromString<EffectiveViscosityModelType>(effectiveViscosityModelType);

    // mu = muo(1-c/c_max)^-\beta
    m_beta = slurryNode->GetAttributeOrDefault("beta","1.5"); m_beta = fabs(m_beta);
    m_muo = proppantData.m_fluidViscosity;
    m_maxVf = proppantData.m_maxVf;

    realT g = proppantData.m_gravity;
    realT mu =  proppantData.m_fluidViscosity;
    realT d = proppantData.m_diameter;
    realT rhos = proppantData.m_density;
    realT rhol = proppantData.m_fluidDensity;

    toLower(singleParticleModelType);
    if( ieq(singleParticleModelType,"stokes") ){
    	  m_singleParticleSettlingSpeed = g*(rhos-rhol)*d*d/(18.0*mu);
    } else if( ieq(singleParticleModelType,"newton") ) {
      	  m_singleParticleSettlingSpeed = 1.74*sqrt(d)*sqrt(g*(rhos-rhol)/rhol);
    } else if( ieq(singleParticleModelType,"allen") ) {
      	  m_singleParticleSettlingSpeed = 0.2*pow(d,1.18)*pow( (g*(rhos-rhol)/rhol),0.72)/pow(mu/rhol,0.45);
    } else if( ieq(singleParticleModelType,"SchillerNauman") ) {
    	  m_singleParticleSettlingSpeed = g*(rhos-rhol)*d*d/(18.0*mu); // use stokes as first guess
    	  // Update settling speed estimate using picard iteration
    	  realT tol = 1e-8;
          realT dvs = m_singleParticleSettlingSpeed;
          while(dvs > tol*m_singleParticleSettlingSpeed){
        	  realT vold = m_singleParticleSettlingSpeed;
        	  realT Rep = ParticleReynoldsNumber(d, rhol, vold, mu);
        	  realT Cd = SchillerNaumanDragCoefficient(Rep);
        	  m_singleParticleSettlingSpeed = SettlingVelocityFromDragCoefficient(Cd,rhos,rhol,g,d);
              dvs = fabs(m_singleParticleSettlingSpeed- vold);
              std::cout << m_singleParticleSettlingSpeed << " " << dvs << std::endl;
          }

    } else if( ieq(singleParticleModelType,"IshiiMishima") ) {
  	  m_singleParticleSettlingSpeed = g*(rhos-rhol)*d*d/18.0*mu; // use stokes as first guess
  	  // Update settling speed estimate using picard iteration
  	  realT tol = 1e-8;
        realT dvs = m_singleParticleSettlingSpeed;
        while(dvs > tol*m_singleParticleSettlingSpeed){
      	  realT vold = m_singleParticleSettlingSpeed;
      	  realT Rep = ParticleReynoldsNumber(d, rhol, vold, mu);
      	  realT Cd = IshiiMishimaDragCoefficient(Rep);
    	  m_singleParticleSettlingSpeed = SettlingVelocityFromDragCoefficient(Cd,rhos,rhol,g,d);
            dvs = fabs(m_singleParticleSettlingSpeed- vold);
            std::cout << m_singleParticleSettlingSpeed << " " << dvs << std::endl;
        }

    } else if( ieq(singleParticleModelType,"Standard") ) {
    	  m_singleParticleSettlingSpeed = g*(rhos-rhol)*d*d/18.0*mu; // use stokes as first guess
    	  // Update settling speed estimate using picard iteration
    	  realT tol = 1e-8;
          realT dvs = m_singleParticleSettlingSpeed;
          while(dvs > tol*m_singleParticleSettlingSpeed){
        	  realT vold = m_singleParticleSettlingSpeed;
        	  realT Rep = ParticleReynoldsNumber(d, rhol, vold, mu);
          	  realT Cd = MorrisonDragCoefficient(Rep);
        	  m_singleParticleSettlingSpeed = SettlingVelocityFromDragCoefficient(Cd,rhos,rhol,g,d);
              dvs = fabs(m_singleParticleSettlingSpeed- vold);
              std::cout << m_singleParticleSettlingSpeed << " " << dvs << std::endl;
          }

    } else if( ieq(singleParticleModelType,"constant") ) {
    	  m_singleParticleSettlingSpeed = slurryNode->GetAttributeValue<realT>(std::string("singleParticleSettlingSpeed"));
    } else {
          throw GPException("Error SlurryBase: Unsupported singleParticleSettlingModel:" +singleParticleModelType +".");
    }
    std::cout << "Single particle terminal velocity (" << singleParticleModelType << "): " << m_singleParticleSettlingSpeed << std::endl;
}



//////////////////////////

// sand erosion from scour
ScourErosionBase::ScourErosionBase(){
/* empty */
}


void ScourErosionBase::ReadXML( TICPP::HierarchicalDataNode* hdn ){
    m_K = hdn->GetAttributeOrDefault("K","0.01 s/m");
    m_beta = hdn->GetAttributeOrDefault("beta","1.0");
    m_tau_o = hdn->GetAttributeOrDefault("tau_o","10 Pa");
    m_density = hdn->GetAttributeOrDefault("rho","2700.0 kg/m^3");

}

//////////////////////////

// simple model of sand erosion from porous medium flow
MatrixFluxErosionBase::MatrixFluxErosionBase(){
/* empty */
}


void MatrixFluxErosionBase::ReadXML( TICPP::HierarchicalDataNode* hdn ){
    m_lambda = hdn->GetAttributeOrDefault("lambda","0.005m^-1");
    m_beta = hdn->GetAttributeOrDefault("beta","1");
}

realT MatrixFluxErosionBase::ErosionRate(realT fluidFlux, realT phi){
	realT rv(0.0);

	if(fluidFlux > 0.0) rv =  m_lambda*(1-phi)*std::pow(fluidFlux,m_beta);

	return rv;
}


//////////////////////////


InertialShearLiftBase::InertialShearLiftBase(){
/* empty */
}


void InertialShearLiftBase::ReadXML( TICPP::HierarchicalDataNode* hdn ){

	// Mollinger and Nieuwstadt (1996)
    m_K = hdn->GetAttributeOrDefault("K","15.57");
    m_beta = hdn->GetAttributeOrDefault("beta","1.87");

    m_rho= hdn->GetAttributeOrDefault("rho","1e-6");
    m_d = hdn->GetAttributeOrDefault("diameter","200 um");
    realT mu = hdn->GetAttributeOrDefault("mu","1e-3 Pa.s");
    m_nu = mu/m_rho;

}





