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
 * @file ParallelPlateFlowSolverBase.h
 * @author settgast1
 * @date Aug 3, 2011
 */

#ifndef PARALLELPLATEFLOWSOLVERBASE_H_
#define PARALLELPLATEFLOWSOLVERBASE_H_

#include "SolverBase.h"
#include "DataStructures/Tables/Table.h"

namespace PPFS
{
class FluidEOSBase;
}

class ParallelPlateFluidModelBase{

    public:
	ParallelPlateFluidModelBase(){/** empty **/};
    virtual ~ParallelPlateFluidModelBase(){/** empty **/};

	virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn){};
	virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                              const realT apb,const realT w, const realT qMag, const realT SHP_FCT) = 0;
	virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                  const realT qMag, const realT SHP_FCT) = 0; // one sided

};





class PowerlawFluidModel: public ParallelPlateFluidModelBase{

  public:
	PowerlawFluidModel():
		ParallelPlateFluidModelBase(), m_n(),m_M(),m_phiM(){/** empty **/};
	~PowerlawFluidModel(){/** empty **/};

	virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);
	virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                                const realT apb,const realT w, const realT qMag, const realT SHP_FCT);
	virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                    const realT qMag, const realT SHP_FCT); // one sided

	realT m_n; // fluid behavior index
	realT m_M; // consistency index
	realT m_phiM;// modified consistency index (phi accounts for flow between parallel plates)
};

/*
 * The Herschel Bulkley model describes a non-Newtonian fluid with a Bingham law-like yield stress
 * accompanied by a power-law fluid like growth in stress as a function of shear strain rate:
 *
 * \tau = \tau_y + k |du/dr|^n   for |\tau| > \tau_y
 * du/dr = 0                     for |\tau| < \tau_y
 *
 * The parallel plate model implemented here follows the approximate solution provided in
 * Wang and Gordaninejad 1999, Flow analysis of field controllable, electro- and magneto-rheological fluids
 * using Herschel-Bulkley model.
 *
 */
class HerschelBulkleyParallelPlateFluidModel: public PowerlawFluidModel{

    public:
	  HerschelBulkleyParallelPlateFluidModel():
	  	  PowerlawFluidModel(),m_tau_y(0){/** empty **/};
	  ~HerschelBulkleyParallelPlateFluidModel(){/** empty **/};

	virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);
	virtual realT CalculatePermeability(const realT la, const realT lb, const realT apa,
                              const realT apb,const realT w, const realT qMag, const realT SHP_FCT);
	virtual realT CalculatePermeability(const realT l, const realT ap,const realT w,
                                  const realT qMag, const realT SHP_FCT); // one sided

	realT m_tau_y; // fluid yield stress
	realT m_phi; // modified consistency index factor

    private:
	Table<1, realT> m_zp_kappaLookup;

	realT CalculateZp(realT kappa, realT n);
	/*
	 lookup table relating the dimensionless plug thickness (zp) to the dimensionless variable
	 kappa = K q^n where
	 K = ((2n+1)/n)^n 2^n/h^(n+1) * M/\tau_y = phi/2* 1/h^(n+1) *M/\tau_y
	 */
	// functions used in Newton's method
	realT analyticalVFunc(realT zp, realT n);
	realT dVdz(realT zp, realT n);

};


class ParallelPlateFlowSolverBase: public SolverBase
{
public:

  ParallelPlateFlowSolverBase( const std::string& name,
                               ProblemManagerT* const pm );

  virtual ~ParallelPlateFlowSolverBase();


  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn) ;
  virtual void RegisterFields( PhysicalDomainT& domain ) = 0 ;
  virtual void Initialize( PhysicalDomainT& domain, SpatialPartition& partition) = 0;


  realT BoundedAperture( const realT& aper )
  {
    realT rval = aper;
    const int b = 1;

    if( aper< b * m_min_aperture )
    {
      const realT pi = 3.14159265358979323846;
      rval = m_min_aperture * ( 2.0/pi * atan( pi/2.0 * (aper/m_min_aperture - b) ) + b );
    }

    return rval;
  }

  realT BoundedApertureDerivative( const realT& aper )
  {
    realT rval = 1;
    const int b = 1;

    if( aper< b * m_min_aperture )
    {
      const realT pi = 3.14159265358979323846;
      rval = 1.0 / (1.0 + pow( pi / 2.0 * ( b - aper / m_min_aperture ) , 2 ) );
    }

    return rval;
  }



  virtual void UpdateEOS( const realT time, const realT dt, PhysicalDomainT& domain)
  {
    throw GPException("ParallelPlateFlowSolverBase::UpdateEOS() not overridden");
  }


  // functions for implicit
  virtual realT Assemble ( PhysicalDomainT& domain,
                           Epetra_System& epetraSystem,
                           const realT& time,
                           const realT& dt )
  {
    throw GPException("ParallelPlateFlowSolverBase::Assemble() not overridden");
    return 0;
  }


  virtual void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition, const realT& time, const realT& dt)
  {
    throw GPException("ParallelPlateFlowSolverBase::Assemble() not overridden");
  }

  virtual void RegisterTemporaryFields( PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::RegisterTemporaryFields() not overridden");
  }

  virtual void DeregisterTemporaryFields( PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::DeregisterTemporaryFields() not overridden");
  }

  virtual void FillTemporaryFields( PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::FillTemporaryFields() not overridden");
  }

  virtual void OverwriteFieldsWithTemporaryFields( PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::OverwriteFieldsWithTemporaryFields() not overridden");
  }

  virtual void SetNumRowsAndTrilinosIndices( PhysicalDomainT& domain,
                                             SpatialPartition& partition,
                                             int& numLocalRows,
                                             int& numGlobalRows )
  { throw GPException("ParallelPlateFlowSolverBase::SetNumRowsAndTrilinosIndices() not overridden"); }

  // functions for explicit
  virtual void UpdateFlux( const realT time,
                           const realT dt,
                           PhysicalDomainT& domain,SpatialPartition& partition)
  {
    throw GPException("ParallelPlateFlowSolverBase::UpdateFlux() not overridden");
  }

  // functions for proppant
  virtual void GenerateSlurryParameters( PhysicalDomainT& domain)
  {
    throw GPException("ParallelPlateFlowSolverBase::GenerateSlurryParameters() not overridden");
  }


  virtual void CalculateMassRate(PhysicalDomainT& domain, SpatialPartition& partition,realT time,realT dt ){
	    throw GPException("ParallelPlateFlowSolverBase::CalculateMassRate() not overridden");
  }


  // boundary conditions
  virtual void CalculateAndApplyMassFlux( const realT dt, PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::CalculateAndApplyMassFlux() not overridden");
  }


  virtual void CalculateCarterLeakOff( const realT time, const realT dt, PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::CalculateCarterLeakOff() not overridden");
  }

  virtual void CalculateMatrixFlowLeakOff( const realT time, const realT dt, PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::CalculateMatrixFlowLeakOff() not overridden");
  }

  virtual void ApplyFluxBoundaryCondition( const realT time, const realT dt, const int cycleNumber, const int rank, PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::ApplyFluxBoundaryCondition() not overridden");
  }

  virtual void FlowControlledBoundaryCondition( PhysicalDomainT& domain, ObjectDataStructureBaseT& object ,
          BoundaryConditionBase* bc ,
          const lSet& set,
          realT time )
  {
    throw GPException("ParallelPlateFlowSolverBase::FlowControlledBoundaryCondition() not overridden");
  }


  virtual void GenerateParallelPlateGeometricQuantities( PhysicalDomainT& domain,
          realT time,
          realT dt )
  {
    throw GPException("ParallelPlateFlowSolverBase::GenerateParallelPlateGeometricQuantities() not overridden");
  }

  virtual void OverwriteOldGeometricQuantities( PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::OverwriteOldGeometricQuantities() not overridden");
  }



  virtual void DefineFlowSets( PhysicalDomainT& domain )
  {
    throw GPException("ParallelPlateFlowSolverBase::DefineFlowSets() not overridden");
  }


  virtual void CalculateNodalPressure ( PhysicalDomainT& domain, SpatialPartition& partition)
  {
    throw GPException("ParallelPlateFlowSolverBase::CalculateNodalPressure() not overridden");
  }





  Array1dT< rArray1d > m_dwdu;
  Array1dT< iArray1d > m_dwdu_dof;
  rArray1d m_dwdw;
  int m_boundPhysicalAperture;
  std::map<localIndex,lArray1d> m_edgesToFaces;

  realT m_leakoffCoef;
  int m_overLeakCompensation;
  // If =1, we allow flow cells to temporarily have a negative fluid mass.
  R1Tensor m_gravityVector;


protected:

  // permeability calculation
  realT m_SHP_FCT;
  realT m_mu;
  realT m_min_aperture;
  realT m_max_aperture;
  realT m_zeroApertureVolume;

  // fluid eos
  realT m_bulk_modulus;
  realT m_rho_o; // fluid density @ pressure=0 (for eos)
  realT m_press_o;

  PPFS::FluidEOSBase * m_fluidEOS;

  int m_updateFaceArea;

  // powerlaw fluid
  bool m_usePowerlawFluid;  // ultimately want a generic fluid model here
  //PowerlawFluidModel m_powerlawFluid;
  ParallelPlateFluidModelBase* m_fluidModelPtr;

  // This is the maximum pressure allowed in the flow solver.  The purpose of having this cap is to prevent sudden pressure spikes.
  // If the pressure value according to EOS is lower than half of this limiting value, the EOS will be directly used.  Otherwise, some special treatment will be adopted.
  // Its value should be significantly lower than the bulk modulus.
  // It should be at least twice of the minimum principal stress in the rock matrix.
  realT m_pressureCap;


  std::string m_flowFaceSetName;


  realT m_dT;
  realT m_farFieldPorePressure;
  int m_pressureDependentLeakoff;


  void RegisterFields( FaceManagerT& faceManager, EdgeManagerT& edgeManager );


private:
  ParallelPlateFlowSolverBase();
};



namespace PPFS
{
  // Equation of state for pressure
  inline realT P_EOS(realT rho, realT K_o, realT rho_o)
  {
    return K_o*(rho - rho_o)/rho_o;
  }

  inline realT Inverse_EOS(realT P, realT K_o, realT rho_o)
  {
    return (P + K_o) * rho_o / K_o;
  }

  // With pressure cap.
  inline realT P_EOS(realT rho, realT K_o, realT rho_o, realT pressureCap)
  {
    realT P;
    if (rho < rho_o)
    {
      P = 0.0;
    }
    else
    {
      P = K_o*(rho - rho_o)/rho;
      if ( P > 0.5 * pressureCap)
      {
        P = 0.5 * pressureCap + ( P - 0.5 * pressureCap) / (K_o - 0.5 * pressureCap) * 0.5 * pressureCap;
      }

    }
    return (P);
  }

  inline realT Inverse_EOS(realT P, realT K_o, realT rho_o, realT pressureCap)
  {
    realT rho;
    if (P <= pressureCap * 0.5)
    {
      rho = K_o * rho_o / (K_o - P);
    }
    else
    {
      rho = K_o * rho_o / (K_o - 0.5*pressureCap - (2.0*P/pressureCap - 1.0) * (K_o - 0.5*pressureCap));
    }
    return(rho);
  }


  // Equation for d Pressure/ d rho
  inline realT dPdRho_EOS(realT K_o, realT rho_o)
  {
    return K_o/rho_o;
  }

  // with pressure cap
  inline realT dPdRho_EOS(realT rho, realT K_o, realT rho_o, realT pressureCap)
  {
	 realT P = K_o*(rho - rho_o)/rho;
	 realT dPdRho(0);
     if ( P < 0.5 * pressureCap)
     {
	   dPdRho = K_o*rho_o/(rho*rho);
     } else {
		 dPdRho *= 0.5 * pressureCap/ (K_o - 0.5 * pressureCap);
	 }

	 return dPdRho;
  }

  // Equation for K: K = rho* d Pressure/ d rho
  inline realT BulkModulus_EOS(realT rho,realT K_o, realT rho_o)
  {
    return rho*dPdRho_EOS(K_o, rho_o);
  }

  // Equation of state for density
  inline realT rho_EOS(realT P,realT K_o, realT rho_o)
  {
    return (1.0 + P/K_o)*rho_o;
  }

  // with pressure cap
  inline realT rho_EOS(realT P, realT K_o, realT rho_o, realT pressureCap)
  {
	    realT rho;
	    if (P <= pressureCap * 0.5)
	    {
	      rho = K_o * rho_o / (K_o - P);
	    }
	    else
	    {
	      rho = K_o * rho_o / (K_o - 0.5*pressureCap - (2.0*P/pressureCap - 1.0) * (K_o - 0.5*pressureCap));
	    }
	    return(rho);
  }
  /**
   *
   * @param la length of the path between centers of two flow elements projected on normal attached to element a
   * @param lb length of the path between centers of two flow elements projected on normal attached to element b
   * @param apa aperture of the first flow element
   * @param apb aperture of the second flow element
   * @param w width of the flow path
   *
   * The "Permeability" returned is scaled by the cell face area
   * and divided by the product of the distance between the cells and kinematic viscosity.
   *
   *   k' = k*(h*w)/(mu*L)
   *
   * so that the mass flux from cell a to cell b is given by
   *
   *   q_Mass = k' * rho * (Pa-Pb)
   *
   */
  inline realT CalculatePermeability(const realT la,
                                     const realT lb,
                                     const realT apa,
                                     const realT apb,
                                     const realT w,
                                     const realT mu,
                                     const realT SHP_FCT)
  {
    if (apa <= 0 || apb <= 0)
      return 0.0;


    const realT ka = apa*apa*apa;
    const realT kb = apb*apb*apb;

    realT permeability = ka * kb * w/(12.0 * mu *(ka*lb + kb*la)); // harmonic mean

    permeability *= SHP_FCT;
    return permeability;
  }

  /**
   * One sided permeability calculation
   *
   * @author walsh24
   *
   * @param la length of the path between centers of two flow elements projected on normal attached to element a\
   * @param h aperture of the flow element
   * @param w width of the flow path
   *
   * NB formulation differs from that in Johnson Morris 2009
   * which uses distance between adjacent elements (2*la) rather than element and boundary (la).
   *
   * The "Permeability" returned is scaled by the cell face area
   * and divided by the product of the distance between the cells and kinematic viscosity.
   *
   *   k' = k*(h*w)/(mu*L)
   *
   * so that the mass flux from cell a to cell b is given by
   *
   *   q_Mass = k' * rho * (Pa-Pb)
   *
   */
  inline realT CalculatePermeability(const realT la,
                                     const realT h,
                                     const realT w,
                                     const realT mu,
                                     const realT SHP_FCT)
  {
    if (h <= 0) return 0.0;

    realT permeability = h*h*h *w / ( 12.0 * mu * la );

    permeability *= SHP_FCT;
    return permeability;
  }


  inline realT CalculatePRhoGravity(const R1Tensor xa,
                                     const R1Tensor xb,
                                     const realT rho_a,
                                     const realT rho_b,
                                     const R1Tensor g)
  {
    R1Tensor vec = xb;
    vec -= xa;
    realT dist = Dot(vec, g);
    realT Prho = dist * (rho_a + rho_b) * 0.5;
    Prho *= (rho_a + rho_b) * 0.5;
    return Prho;
  }

  inline realT CalculatePRhoGravity(const R1Tensor xa,
                                     const R1Tensor xb,
                                     const realT rho_a,
                                     const R1Tensor g)
  {
    R1Tensor vec = xb;
    vec -= xa;
    realT dist = Dot(vec, g);
    realT Prho = dist * rho_a * rho_a;
    return Prho;
  }

  /**
     *
     * @param la length of the path between centers of two flow elements projected on normal attached to element a
     * @param lb length of the path between centers of two flow elements projected on normal attached to element b
     * @param apa aperture of the first flow element
     * @param apb aperture of the second flow element
     * @param w width of the flow path
     * @param phiM = M' = modified consistency index     (units Pa.s^n)
     * @param qMag magnitude of fluid velocity
     * @param n fluid behavior index  (dimensionless)
     *
     * Power Law fluid: \sigma_{xy} = M(2 \dot{\epsilon}_{xy})^{n}
     * n = 1 : Newtonian fluid
     * n < 1 : Shear thinning
     * n > 1 : Shear thickening
     *
     * For a powerlaw fluid between two parallel plates:
     * dp/dx_{i} = - \phi M q_i |q|^{n-1}/h^{2n+1}
     * Where \phi = 2^{n+1}(2*n+1)^n/n^n
     * See Adachi and Detournay,
     * Self-similar solution of a plane-strain fracture driven by a power-law fluid (2002)
     *
     * The "Permeability" returned is scaled by the cell face area
     * and divided by the product of the distance between the cells and kinematic viscosity.
     *
     *   k' = k*(h*w)/(mu*L)
     *
     * so that the mass flux from cell a to cell b is given by
     *
     *   q_Mass = k' * rho * (Pa-Pb)
     *
     */
    inline realT CalculatePermeability_PowerLawFluid(const realT la,
                                                     const realT lb,
                                                     const realT apa,
                                                     const realT apb,
                                                     const realT w,
                                                     const realT qMag,
                                                     const realT phiM,
                                                     const realT n,
                                                     const realT SHP_FCT)
    {
      if (apa <= 0 || apb <= 0)
        return 0.0;


      const realT ka = pow(apa,2*n+1);
      const realT kb = pow(apb,2*n+1);

      realT permeability = ka * kb * w*pow(qMag,1-n)/(phiM *(ka*lb + kb*la)); // harmonic mean

      permeability *= SHP_FCT;
      return permeability;
    }

    // one sided permeability calculation for powerlaw fluid
    inline realT CalculatePermeability_PowerLawFluid(const realT la,
                                                     const realT h,
                                                     const realT w,
                                                     const realT qMag,
                                                     const realT phiM,
                                                     const realT n,
                                                     const realT SHP_FCT)
    {
      if (h <= 0) return 0.0;

      realT permeability = pow(h,2*n+1) * w*pow(qMag,1-n)/(phiM*la);

      permeability *= SHP_FCT;
      return permeability;
    }

    // modified consistency index is required when calculating the power law fluid permeability
    inline realT CalculateModifiedConsistencyIndex(const realT M, const realT n){
    	realT phi =2*pow(2*(2*n+1)/n,n);
    	return phi*M;
    }

    /**
        *
        * @param la length of the path between centers of two flow elements projected on normal attached to element a
        * @param lb length of the path between centers of two flow elements projected on normal attached to element b
        * @param apa aperture of the first flow element
        * @param apb aperture of the second flow element
        * @param w width of the flow path
        * @param phiM = M' = modified consistency index     (units Pa.s^n)
        * @param qMag magnitude of fluid velocity
        * @param n fluid behavior index  (dimensionless)
        * @param z_p dimensionless plug width = |2 \tau_y|/(|p_{,i}|h) where h = aperture, \tau_y is the yield stress and |p_{,i}| is the magnitude of the pressure gradient
        *
        * Herschel-Bulkley fluid:
        *    \tau = \tau_y + M(2 \dot{\epsilon}_{xy})^{n}  for \tau > \tau_y
        *    \dot{\epsilon}_{xy} = 0 for \tau < \tau_y
        *
        * See Wang and Gordaninejad 1999 "Flow analysis of Field controllable Electro and Magneto-Rheological fluids using Herschel Bulkley Model",
        * Self-similar solution of a plane-strain fracture driven by a power-law fluid (2002)
        *
        * The "Permeability" returned is scaled by the cell face area
        * and divided by the product of the distance between the cells and kinematic viscosity.
        *
        *   k' = k*(h*w)/(mu*L)
        *
        * so that the mass flux from cell a to cell b is given by
        *
        *   q_Mass = k' * rho * (Pa-Pb)
        *
        */
       inline realT CalculatePermeability_HerschelBulkleyFluid(const realT la,
                                                        const realT lb,
                                                        const realT apa,
                                                        const realT apb,
                                                        const realT w,
                                                        const realT qMag,
                                                        const realT phiM,
                                                        const realT n,
                                                        const realT zp,
                                                        const realT SHP_FCT)
       {
         if (apa <= 0 || apb <= 0 || zp >= 1)
           return 0.0;

         realT permeability =CalculatePermeability_PowerLawFluid(la, lb,apa, apb,w, qMag, phiM, n, SHP_FCT);
         realT denom = pow(1-zp , n+1)*pow(n/(n+1)*zp +1,n);

         permeability /= denom;

         return permeability;
       }

       // one sided permeability calculation for Herschel-Bulkley fluid
       inline realT CalculatePermeability_HerschelBulkleyFluid(const realT la,
                                                        const realT h,
                                                        const realT w,
                                                        const realT qMag,
                                                        const realT phiM,
                                                        const realT n,
                                                        const realT zp,
                                                        const realT SHP_FCT)
       {
         if (h <= 0|| zp >= 1) return 0.0;


         realT permeability =CalculatePermeability_PowerLawFluid(la,h,w, qMag, phiM, n, SHP_FCT);
         realT denom = pow(1-zp , n+1)*pow(n/(n+1)*zp +1,n);

         permeability /= denom;

         return permeability;
       }




       ////////////////////////
       // Equation of state  //
       ////////////////////////

       // base class for fluid equation of state
       class FluidEOSBase{
    	   public:

    	     FluidEOSBase( TICPP::HierarchicalDataNode* hdn){};
    	     FluidEOSBase(realT rho):m_rho_o(rho){/*empty*/};
             virtual ~FluidEOSBase(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             virtual realT pressure(realT rho) = 0;
             virtual realT density(realT P) = 0;
             virtual realT dPdRho(realT rho) = 0;
             virtual realT rho_o(){return m_rho_o;};
    	   protected:
             realT m_rho_o;
       };

       // Linear EOS
       class LinearEOS: public FluidEOSBase{
    	   public:

    	     LinearEOS(TICPP::HierarchicalDataNode* hdn);
    	     LinearEOS(realT rho, realT Ko): FluidEOSBase(rho), m_K_o(Ko){ /* empty */};
             virtual ~LinearEOS(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             virtual realT pressure(realT rho){return m_K_o*(rho - m_rho_o)/m_rho_o;};
             virtual realT density(realT P){return (P + m_K_o) * m_rho_o / m_K_o; };
             virtual realT dPdRho(realT rho){return m_K_o/m_rho_o;};

             static const char* FluidEOSName(){return "LinearEOS";};
    	   protected:
             realT m_K_o;
       };

       // PressureCap EOS
       class PressureCapEOS: public LinearEOS{
    	   public:

    	     PressureCapEOS(TICPP::HierarchicalDataNode* hdn);
    	     PressureCapEOS(realT rho, realT Ko, realT pCap):LinearEOS(rho,Ko), m_pressureCap(pCap){/* empty */};
             virtual ~PressureCapEOS(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             static const char* FluidEOSName(){return "PressureCapEOS";};

             virtual realT pressure(realT rho){
            	    realT P;
            	    if (rho < m_rho_o)
            	    {
            	      P = 0.0;
            	    }
            	    else
            	    {
            	      P = m_K_o*(rho - m_rho_o)/rho;
            	      if ( P > 0.5 * m_pressureCap)
            	      {
            	        P = 0.5 * m_pressureCap + ( P - 0.5 * m_pressureCap) / (m_K_o - 0.5 * m_pressureCap) * 0.5 * m_pressureCap;
            	      }

            	    }
            	    return (P);
             };

             virtual realT density(realT P){
            	    realT rho;
            	    if (P <= m_pressureCap * 0.5)
            	    {
            	      rho = m_K_o * m_rho_o / (m_K_o - P);
            	    }
            	    else
            	    {
            	      rho = m_K_o * m_rho_o / (m_K_o - 0.5*m_pressureCap - (2.0*P/m_pressureCap - 1.0) * (m_K_o - 0.5*m_pressureCap));
            	    }
            	    return(rho);
             };

             virtual realT dPdRho(realT rho){
            	 realT P = m_K_o*(rho - m_rho_o)/rho;
            	 realT dPdRho(0);
            	 if ( P < 0.5 * m_pressureCap)
            	 {
            	     dPdRho = m_K_o*m_rho_o/(rho*rho);
            	 } else {
            	 	 dPdRho *= 0.5 * m_pressureCap/ (m_K_o - 0.5 * m_pressureCap);
            	 }
            	 return dPdRho;
             };

    	   protected:
             realT m_pressureCap;
       };

       // Adiabatic EOS
       class AdiabaticEOS: public FluidEOSBase{
    	   public:

    	     AdiabaticEOS(TICPP::HierarchicalDataNode* hdn);
             virtual ~AdiabaticEOS(){};
             virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn);

             virtual realT pressure(realT rho){
            	 realT rv = (rho > 0)? m_P_o*std::pow(rho/m_rho_o,m_gamma): 0.0;
            	 rv -= m_P_ref;
            	 if(rv < 0.0) rv = 0.0;
                 return rv;
             };
             virtual realT density(realT P){return (P+m_P_ref > 0)? m_rho_o*std::pow( (P+m_P_ref)/m_P_o, 1/m_gamma ): 0.0; };
             virtual realT dPdRho(realT rho){
            	 realT dPdRho = 0;
            	 //if(rho > 0) dPdRho = m_gamma*m_P_o * std::pow(rho/m_rho_o,m_gamma)/rho;
            	 realT p = pressure(rho);
            	 if(p > 0) dPdRho = m_gamma*m_P_o * std::pow(rho/m_rho_o,m_gamma)/rho;
            	 return dPdRho;
             };

             static const char* FluidEOSName(){return "AdiabaticEOS";};
    	   protected:
             realT m_P_o;
             realT m_gamma;
             realT m_P_ref;
       };


//////////////////////////

// Fluid EOS Factory
//
// Consists of the following parts:
//   * The function to generate new pointers: "newFluidEOS"
//   * A base class to derive the functions to generate FluidEOS pointers: "FluidEOSInitializer"
//   * A String-to-FluidEOS-Intializer map hidden behind the getFluidEOSCatalogue function
//   * A template to create FluidEOS initializers: "FluidEOSRegistrator"
//   * A compiler directive to simplify autoregistration: "REGISTER_FLUID_EOS"
//
// Most initial conditions will only need to use one or two of the parts:
//   * To register a new FluidEOS in the factory: REGISTER_FluidEOS( FluidEOSClassName )
//   * To load a FluidEOS pointer from the factory:       FluidEOSBase* aFluidEOSPtr = newFluidEOS(FluidEOSString, args );

/// The FluidEOS Factory.
FluidEOSBase* newFluidEOS(const std::string& EOSName, TICPP::HierarchicalDataNode* hdn);

/// Base class to generate new FluidEOS pointers
class FluidEOSInitializer{
public:
	virtual FluidEOSBase* initializeFluidEOS( TICPP::HierarchicalDataNode* hdn) = 0;
	virtual ~FluidEOSInitializer() = 0;
};
inline FluidEOSInitializer::~FluidEOSInitializer() { }

/// Interface to the FluidEOS name -> FluidEOS initializer map
std::map<std::string, FluidEOSInitializer*> & getFluidEOSCatalogue();

/// Return a list of supported FluidEOS names
void getFluidEOSNames( std::vector<std::string>& nameList);

/// Template for creating classes derived from FluidEOSInitializer
template< class FluidEOSType >
class FluidEOSRegistrator : public FluidEOSInitializer{

public:
	FluidEOSRegistrator(void){
		std::string FluidEOSName = std::string(FluidEOSType::FluidEOSName() );
		getFluidEOSCatalogue() [FluidEOSName] = this;
	};

	FluidEOSBase* initializeFluidEOS( TICPP::HierarchicalDataNode* hdn ){
		return new FluidEOSType( hdn );
	};
};



} // end namespace

/// Compiler directive to simplify autoregistration
#define REGISTER_FLUID_EOS( ClassName ) namespace{ PPFS::FluidEOSRegistrator<PPFS::ClassName> reg_##ClassName; }


#endif /* PARALLELPLATEFLOWSOLVERBASE_H_ */
