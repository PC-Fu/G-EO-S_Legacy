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
 * @file ParallelPlateFlowSolverBase.cpp
 * @author settgast1
 * @date Aug 3, 2011
 */

#include "Common/Common.h"
#include "PhysicsSolverStrings.h"

#include "ParallelPlateFlowSolverBase.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"


using namespace PS_STR;

ParallelPlateFlowSolverBase::ParallelPlateFlowSolverBase( const std::string& name,
                                                          ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

ParallelPlateFlowSolverBase::~ParallelPlateFlowSolverBase()
{}


void ParallelPlateFlowSolverBase::ReadXML( TICPP::HierarchicalDataNode* const hdn  )
{
  SolverBase::ReadXML(hdn);

  // Permeability
  m_mu = hdn->GetAttributeOrDefault("mu","1.0e-3 N.s/m^2");
  m_SHP_FCT = hdn->GetAttributeOrDefault<realT>("shapefactor",1.0);
  m_min_aperture = hdn->GetAttributeOrDefault(MinimumApertureStr,"1 um");
  m_max_aperture = hdn->GetAttributeOrDefault(MaximumApertureStr,"4 mm");
  m_zeroApertureVolume = hdn->GetAttributeOrDefault("ZeroApertureVolume","1 mm");

  R1Tensor zeroVector;
  zeroVector *= 0.0;
  m_gravityVector = hdn->GetAttributeOrDefault<R1Tensor>("gravityVector", zeroVector);

  // Fluid EOS
  /*
  m_bulk_modulus = hdn->GetAttributeOrDefault(BulkModulusStr,"2.0e9 Pa");
  m_rho_o = hdn->GetAttributeOrDefault("rho_o","1 kg/L");
  m_press_o = hdn->GetAttributeOrDefault("press_o","0.0 Pa");
  m_pressureCap = hdn->GetAttributeOrDefault("pressurecap","1.0e8 Pa");
  if (m_pressureCap > m_bulk_modulus)
  {
    throw GPException("The pressure cap is higher than the bulk modulus!");
  }
  */

  TICPP::HierarchicalDataNode* fluidEOShdn = hdn->GetChild("FluidEOS");
  if(fluidEOShdn){
	std::string fluidEOSname = fluidEOShdn->GetAttributeString("name");
    m_fluidEOS = PPFS::newFluidEOS(fluidEOSname,fluidEOShdn);
  } else {
	  // build pressure cap model
	  m_bulk_modulus = hdn->GetAttributeOrDefault(BulkModulusStr,"2.0e9 Pa");
	  m_rho_o = hdn->GetAttributeOrDefault("rho_o","1 kg/L");
	  m_press_o = hdn->GetAttributeOrDefault("press_o","0.0 Pa");
	  m_pressureCap = hdn->GetAttributeOrDefault("pressurecap","1.0e8 Pa");
	  if (m_pressureCap > m_bulk_modulus)
	  {
	    throw GPException("The pressure cap is higher than the bulk modulus!");
	  }
	  m_fluidEOS = new PPFS::PressureCapEOS(m_rho_o, m_bulk_modulus, m_pressureCap);

  }


  // Viscosity model

  m_usePowerlawFluid = false;
  TICPP::HierarchicalDataNode* powerLawFluidNode = hdn->GetChild("PowerlawFluid");
  if(powerLawFluidNode){
	  //m_powerlawFluid.ReadXML(powerLawFluidNode);
	  m_usePowerlawFluid = true;
	  m_fluidModelPtr = new PowerlawFluidModel();
	  m_fluidModelPtr->ReadXML(powerLawFluidNode);
  }

  TICPP::HierarchicalDataNode* hbFluidNode = hdn->GetChild("HerschelBulkleyFluid");
  if(hbFluidNode){
	  m_usePowerlawFluid = true;
	  m_fluidModelPtr = new HerschelBulkleyParallelPlateFluidModel();
	  m_fluidModelPtr->ReadXML(hbFluidNode);

  }
  if(!m_usePowerlawFluid){
    m_fluidModelPtr = new PowerlawFluidModel();
  }

  m_updateFaceArea = hdn->GetAttributeOrDefault("updateFaceArea","0");


  // Timestep
  m_courant = hdn->GetAttributeOrDefault<realT>("courant",0.5);

  // Faceset
  m_flowFaceSetName = hdn->GetAttributeString("flowFaceSet");
  if(m_flowFaceSetName.empty()) m_flowFaceSetName = hdn->GetAttributeString("faceset");
}


void ParallelPlateFlowSolverBase::RegisterFields( FaceManagerT& faceManager, EdgeManagerT& edgeManager )
{

  faceManager.AddKeylessDataField<realT>(ApertureStr,true,true);
  faceManager.AddKeylessDataField<R1Tensor>(FaceCenterStr,true,true);
  faceManager.AddKeyedDataField<FieldInfo::pressure>();
  faceManager.AddKeyedDataField<FieldInfo::density>();
  faceManager.AddKeyedDataField<FieldInfo::mass>();
  faceManager.AddKeylessDataField<realT>( "appliedPressure", true, true );
  faceManager.AddKeylessDataField<int>( "flowFaceType", true, true );
  faceManager.AddKeylessDataField<realT>( "faceArea", true, false );

  edgeManager.AddKeylessDataField<int>( "flowEdgeType", true, true );
  edgeManager.AddMap< UnorderedVariableOneToManyRelation >( "edgeToFlowFaces");
  edgeManager.AddKeylessDataField<realT>(PermeabilityStr,true,true);
  edgeManager.AddKeylessDataField<realT>(PressureStr,true,true);
  edgeManager.AddKeylessDataField<realT>(VolumetricFluxStr,true,true);
  edgeManager.AddKeylessDataField<realT>("length",true,true);
  edgeManager.AddKeylessDataField<R1Tensor>("center",true,true);


}




/**
 *  Power law fluid model
 */

void PowerlawFluidModel::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
	m_M = hdn->GetAttributeOrDefault("ConsistencyIndex","1.0e-3 N.s/m^2");
	m_n = hdn->GetAttributeOrDefault("FluidBehaviorIndex","1");
	m_phiM = PPFS::CalculateModifiedConsistencyIndex(m_M,m_n);
}

realT PowerlawFluidModel::CalculatePermeability(const realT la, const realT lb,
		                                        const realT apa, const realT apb,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){
		return PPFS::CalculatePermeability_PowerLawFluid(la,lb,apa,apb, w, qMag, m_phiM, m_n, SHP_FCT);
}

// one sided permeability
realT PowerlawFluidModel::CalculatePermeability(const realT l, const realT ap,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){
		return PPFS::CalculatePermeability_PowerLawFluid(l,ap, w, qMag, m_phiM, m_n, SHP_FCT);
}



/**
 *  Herschel Bulkley Parallel Plate Fluid Model
 */

void HerschelBulkleyParallelPlateFluidModel::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
	PowerlawFluidModel::ReadXML(hdn);
	m_tau_y = hdn->GetAttributeOrDefault("YieldStress","0.0");
	m_phi = m_phiM/m_M;

	/* build lookup table */
    rArray1d xs;
    rArray1d values;
    xs.resize(16);
    values.resize(16);
    for(int i =0 ; i < 16; i++){
     realT x = -3 + 0.2*i;
     xs[i] = x;
     realT kappa = pow(10,x);
     values[i] = CalculateZp(kappa,m_n);
    }
	m_zp_kappaLookup.SetAxisValues(0,xs);
	m_zp_kappaLookup.SetValues(values);

    for(int i =0 ; i < 16; i++){
    	std::cout <<  "k: " <<  pow(10,xs[i])
    			  << " v: " <<  values[i]
                  << "\n";
    }


}

// analytical solution relating dimensionless parameter V to zp (dimensionless plug thickness)
realT HerschelBulkleyParallelPlateFluidModel::analyticalVFunc(realT zp, realT n){
	realT V = pow(1-zp,n+1)*pow(n/(n+1)*zp + 1,n);
	return V;
}
// derivative of analyticalVFunc wrt zp
realT HerschelBulkleyParallelPlateFluidModel::dVdz(realT zp, realT n){

	  realT deriv = -(n+1)*pow(1-zp,n) *pow(n/(n+1)*zp + 1,n) +
	          pow(1-zp,n+1)*(n*n/(n+1))*pow(n/(n+1)*zp + 1,n-1);
	  return deriv;
}

// calculate Zp as a function of kappa = Kq^n
realT HerschelBulkleyParallelPlateFluidModel::CalculateZp(realT Kqn, realT n){

  realT zp = 0.5;
  realT dzp = 1;
  int numIters = 0;
  realT tol = 1e-8;
  // newton's method
  while(fabs(dzp) > tol*zp && numIters < 500){
    realT f = Kqn*zp - analyticalVFunc(zp,n);
	realT dfdz = Kqn - dVdz(zp,n);
	dzp = -f/dfdz;
	zp += dzp;
	numIters += 1;
  }
  if(numIters == 500)
    throw GPException("HerschelBulkleyParallelPlateFluidModel: Zp calculation failed to converge");
  return zp;
}

realT HerschelBulkleyParallelPlateFluidModel::CalculatePermeability(const realT la, const realT lb,
		                                        const realT apa, const realT apb,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){

        realT zp(0.0);
        if(m_tau_y > 0.0){
	      realT kappa = pow( (2*m_n+1)/(2*m_n),m_n)*pow( 4*qMag/w,m_n)* m_M/m_tau_y;
          if(kappa > 1000.0){
        	zp = 1.0/kappa;
          } else if(kappa < 1e-3){
        	realT denom = pow(1 + m_n/(m_n+1),m_n);
        	zp = 1 - pow(kappa/denom , 1.0/(m_n+1));
          } else {
    	    realT xx[1];
    	    xx[0] = log10(kappa);
        	zp = m_zp_kappaLookup.Lookup(xx);
          }
        }
		return PPFS::CalculatePermeability_HerschelBulkleyFluid(la,lb,apa,apb, w, qMag, m_phiM, m_n, zp, SHP_FCT);
}

// one sided permeability
realT HerschelBulkleyParallelPlateFluidModel::CalculatePermeability(const realT l, const realT ap,
                                                const realT w, const realT qMag,
                                                const realT SHP_FCT){

    realT zp(0.0);
    if(m_tau_y > 0.0){
      realT kappa = pow( (2*m_n+1)/(2*m_n),m_n)*pow( 4*qMag/w,m_n)* m_M/m_tau_y;
      if(kappa > 1000.0){
    	zp = 1.0/kappa;
      } else if(kappa < 1e-3){
     	realT denom = pow(1 + m_n/(m_n+1),m_n);
    	zp = 1 - pow(kappa/denom , 1.0/(m_n+1));
      } else {
	    realT xx[1];
	    xx[0] = log10(kappa);
    	zp = m_zp_kappaLookup.Lookup(xx);
      }
    }
	return PPFS::CalculatePermeability_HerschelBulkleyFluid(l,ap, w, qMag, m_phiM, m_n,zp, SHP_FCT);
}



//// Fluid EOS
///////////////

void PPFS::FluidEOSBase::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
	m_rho_o = hdn->GetAttributeOrDefault("rho_o","1 kg/L");
}


//// Linear EOS

PPFS::LinearEOS::LinearEOS(TICPP::HierarchicalDataNode* hdn):
FluidEOSBase(hdn),
m_K_o(){
	  ReadXML(hdn);
}

void PPFS::LinearEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
	FluidEOSBase::ReadXML(hdn);
	m_K_o = hdn->GetAttributeOrDefault(BulkModulusStr,"2.0e9 Pa");
}

REGISTER_FLUID_EOS( LinearEOS )



//// Pressure cap

PPFS::PressureCapEOS::PressureCapEOS(TICPP::HierarchicalDataNode* hdn):
		LinearEOS(hdn),
        m_pressureCap(){
	  ReadXML(hdn);

}

void PPFS::PressureCapEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
	LinearEOS::ReadXML(hdn);
	m_pressureCap = hdn->GetAttributeOrDefault(BulkModulusStr,"1.0e8 Pa");

	if (m_pressureCap > m_K_o)
	{
	    throw GPException("PressureCapEOS: The pressure cap is higher than the bulk modulus!");
	}
}
REGISTER_FLUID_EOS( PressureCapEOS )

//// Adiabatic eos

PPFS::AdiabaticEOS::AdiabaticEOS(TICPP::HierarchicalDataNode* hdn):
		FluidEOSBase(hdn),
		m_P_o(),
		m_P_ref(),
		m_gamma(){
	  ReadXML(hdn);

}

void PPFS::AdiabaticEOS::ReadXML( TICPP::HierarchicalDataNode* const hdn  ){
	FluidEOSBase::ReadXML(hdn);
	m_P_o = hdn->GetAttributeOrDefault("P_o","1e5 Pa");
	m_P_ref = hdn->GetAttributeOrDefault("P_ref","0");
	m_rho_o = hdn->GetAttributeOrDefault("rho_o","1.204 kg/m^3");  // reset default rho_o to reflect gas
	m_gamma = hdn->GetAttributeOrDefault("gamma","7.0/5.0"); // diatomic gas
}

REGISTER_FLUID_EOS( AdiabaticEOS )


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////////
//
// Fluid EOS factory

typedef std::map<std::string, PPFS::FluidEOSInitializer*> FluidEOSCatalogueType;

FluidEOSCatalogueType & PPFS::getFluidEOSCatalogue(){
  static FluidEOSCatalogueType theCatalogue ;
  return theCatalogue;
}

void PPFS::getFluidEOSNames( std::vector<std::string>& nameList){

  using namespace PPFS;
  for(FluidEOSCatalogueType::const_iterator it = getFluidEOSCatalogue().begin();
      it != getFluidEOSCatalogue().end(); ++it){
        nameList.push_back(it->first);
  }
}

PPFS::FluidEOSBase* PPFS::newFluidEOS(const std::string& FluidEOSName , TICPP::HierarchicalDataNode* hdn)
{
  using namespace PPFS;

  FluidEOSInitializer* FluidEOSInitializer = getFluidEOSCatalogue()[FluidEOSName];
  FluidEOSBase *theNewFluidEOS = NULL;

  if(!FluidEOSInitializer)
      throw GPException("Could not create unrecognized FluidEOS: "+ FluidEOSName);

  theNewFluidEOS = FluidEOSInitializer->initializeFluidEOS( hdn );

  return theNewFluidEOS;
}

