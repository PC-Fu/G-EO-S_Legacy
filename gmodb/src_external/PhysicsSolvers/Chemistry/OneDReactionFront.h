//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2013, Lawrence Livermore National Security, LLC.
//
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)     Stuart Walsh(walsh24@llnl.gov)
//  Scott Johnson (johnson346@llnl.gov)        Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)           
//
//  All rights reserved.
//
//  This file is part of GPAC.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file OneDReactionFront.h
 * @author walsh24
 * @date July 27, 2011
 */

#ifndef REACTIONFRONTSOLVER_H_
#define REACTIONFRONTSOLVER_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/ChemistryManager.h"
#include "ObjectManagers/LogManager.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

#include "Common/Common.h"
#include "Utilities/Functions.h"

namespace ODRF{

  /// Base class for non-linear systems of equations
  class NonLinearSystemBase{

    public:
      NonLinearSystemBase(){};
      virtual ~NonLinearSystemBase(){};
      virtual void GetFunction(rArray1d& X,rArray1d& F) = 0;
      virtual rArray2d& GetJacobian(const rArray1d& X) = 0;
      virtual void RescaleSolution(rArray1d& dX, rArray1d& X,rArray1d& F){/** empty  **/};
 
      rArray1d m_F;
      rArray1d m_g;
      rArray1d m_dx;

  };
	
  /// 
  class EquilibriumSystem: public NonLinearSystemBase{
    public:
      EquilibriumSystem():
        NonLinearSystemBase(),
        m_K(),m_p(),m_A(),m_b(){};
      EquilibriumSystem(const rArray2d& K, const rArray1d& p, 
                        const rArray2d& A, const rArray1d& b,
                        const sArray1d& species, realT Tc):
        NonLinearSystemBase(),
        m_K(K),m_p(p),m_A(A),m_b(b),
        m_ionicStrengthCoefficients( GPChemistry::IonicStrengthCoefficients(species) ),
        m_activityFunction(species,Tc),
        m_doScaling(false),
        m_isolatedDofs()
        {
          BuildJacobian();
        };
      ~EquilibriumSystem(){};

      rArray1d& Get_b(void){return m_b;};      

      void Initialize(const rArray2d& K, const rArray1d& p,
                      const rArray2d& A, const rArray1d& b,
                      const sArray1d& species, realT Tc){
        m_K = K;
        m_p = p;
        m_A = A;
        m_b = b;
        m_activityFunction = ExtendedDebyeHuckelActivityFunction(species,Tc);
        m_ionicStrengthCoefficients = GPChemistry::IonicStrengthCoefficients(species);
        BuildJacobian();
      }

      void GetFunction(rArray1d& X,rArray1d& F);
      rArray2d& GetJacobian(const rArray1d& X);
    private:
      rArray2d m_K;
      rArray1d m_p;
      rArray2d m_A;
      rArray1d m_b;

      // Activity function
      rArray1d m_ionicStrengthCoefficients;
      //DaviesActivityFunction m_activityFunction;
      ExtendedDebyeHuckelActivityFunction m_activityFunction;
      
      rArray2d m_gradF;

      bool m_doScaling;
      lArray1d m_isolatedDofs;

      void BuildJacobian(void);
      void CheckActivities(void);
  };
  
  /// 
  class ReactionFrontSystem: public NonLinearSystemBase{
    public:
      ReactionFrontSystem():
        m_K(),m_p(),m_A(),m_b(),m_numFronts(){};
      ReactionFrontSystem(const rArray2d& K, const rArray1d& p, 
                          const rArray2d& Kf, const rArray1d& pf, 
                          const rArray2d& A, const rArray1d& b, int numFronts,
                          const sArray1d& species,  realT Tc):
        m_K(K),m_Kf(Kf),m_p(p),m_pf(pf),m_A(A),m_b(b),
        m_ionicStrengthCoefficients( GPChemistry::IonicStrengthCoefficients(species) ),
        m_activityFunction(species,Tc),
        m_numFronts(numFronts),
        m_doScaling(false),
        m_isolatedDofs()
{BuildJacobian(); };
      ~ReactionFrontSystem(){};
      
      void Initialize(const rArray2d& K, const rArray1d& p,
                      const rArray2d& Kf, const rArray1d& pf,
                      const rArray2d& A, const rArray1d& b, int numFronts,
                      const sArray1d& species,  realT Tc){
        m_K=K;    m_p= p;
        m_Kf =Kf; m_pf = pf;
        m_A =A;   m_b = b;
        m_numFronts = numFronts;
        //m_activityFunction = DaviesActivityFunction(species,Tc);
        m_activityFunction = ExtendedDebyeHuckelActivityFunction(species,Tc);
        m_ionicStrengthCoefficients = GPChemistry::IonicStrengthCoefficients(species);
        BuildJacobian();
      };

      rArray2d& Get_A(void){return m_A;};
      rArray1d& Get_b(void){return m_b;};
      
      void GetFunction(rArray1d& X,rArray1d& F);
      rArray2d& GetJacobian(const rArray1d& X);
      void UpdateJacobian(void);
      void RescaleSolution(rArray1d& dX, rArray1d& X,rArray1d& F);

      realT GetIonicStrengthCoefficient(int i){return m_ionicStrengthCoefficients[i];};
      
    private:
      rArray2d m_K;
      rArray2d m_Kf;
      rArray1d m_p;
      rArray1d m_pf;
      rArray2d m_A;
      rArray1d m_b;

      // Activity function
      rArray1d m_ionicStrengthCoefficients;
      //DaviesActivityFunction m_activityFunction;
      ExtendedDebyeHuckelActivityFunction m_activityFunction;

      unsigned m_numFronts;
      
      rArray2d m_gradF;
      
      bool m_doScaling;

      lArray1d m_isolatedDofs;

      void BuildJacobian(void);
      void CheckActivities(void);
  };
	
}

/// 
class ReactionFrontSolver : public SolverBase
{
public:
  ReactionFrontSolver(  const std::string& name,
                        ProblemManagerT* const pm );
  virtual ~ReactionFrontSolver();
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );


  virtual void SetMaxStableTimeStep( PhysicalDomainT& domain,
                                     const sArray1d& namesOfSolverRegions,
                                     SpatialPartition& partition);


  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "ReactionFrontSolver";};
  
  
  // might make the following private
  int FindReactionFrontEquilibrium(bool haveInitialGuess = true);
  void UpdateInitialGuess(unsigned face);
  void SetRegionLengths();
  void BuildReactionFrontEquations(rArray2d& AA, rArray1d& bb);

private:
   
   ODRF::EquilibriumSystem m_brineSystem;  // Brine equilibrium
   ODRF::ReactionFrontSystem m_rfSystem; // Reaction front equilibrium 
   
   // Data
   struct{
     sArray1d names;
     unsigned size;
     std::map<std::string,int> indices;
     iArray1d valences;
     rArray1d ionicStrengthCoefficients; // 0.5*valence^2
     sArray1d eqConstants;
     rArray1d brineConcentrations;  // Concentrations followed by their natural log
     rArray1d frontConcentrations;  // Concentrations followed by their natural log
   } m_species;

   struct{
     sArray1d names;
     unsigned size;
     std::map<std::string,int> indices;
     std::vector<rArray1d > species;
     rArray1d brineConcentrations;  // Net element concentrations
   } m_elements;

   struct{
     sArray1d names;
     unsigned size;
     std::map<std::string,int> indices;
     rArray1d vfs;
     rArray1d MMs;
     rArray1d rhos;
   } m_solidphases;

   struct{
     std::map< std::string, rArray1d > equations;
     std::map< std::string, rArray1d > log10Kdata;
   } m_equilibriumData;

   struct{
     sArray1d names;
     unsigned size;
     sArray1d eqConstants;
     rArray1d initial_positions;
     rArray1d current_positions;
     rArray1d velocities;
     rArray1d lengths;
     std::vector< std::vector< rArray1d > > front_conditions;
     std::vector< rArray1d > dissolution_rates;
     rArray1d molar_densities;
     realT minimum_width;
   } m_fronts;

   struct{
     sArray1d names;
     std::vector<sArray1d > solidPhases;
     rArray1d phis;
     rArray1d tortuosities;
     rArray1d lengths;
   } m_regions;
   
   // Effective diffusivity;
   realT m_D; // diffusivity
   Function* m_effDiffFuncPtr;
   std::string m_EffectiveDiffusivityFunction;
   realT EffectiveDiffusivity(realT phi,realT tau);
   
   realT ReactionFrontVelocity(realT C,realT C_m,realT C_p,
                               realT dx_m,realT dx_p,
                               realT phi_m,realT phi_p,
                               realT tau_m,realT tau_p,
                               realT MolarDensity);

   // Activity function
   //DaviesActivityFunction m_activityFunction;
   ExtendedDebyeHuckelActivityFunction m_activityFunction;
   
   // Face set
   const lSet* m_faceSet;
   std::string m_faceSetName;
   localIndex m_numFaces;
  
   // Temperature   
   realT m_Tc; // Temperature in Celcius
 
   // equilibrium equations
   rArray2d m_KK; // equilibrium data for aqueous phases
   rArray1d m_pp;
   
   rArray2d m_Kf;  // equilibrium data for front (solid) phases
   rArray1d m_pf;
   
   rArray1d m_X_br; // initial guess for brine 
   rArray1d m_X_s;  // initial guess for solid
 
   // data pointers
   std::vector<Array1dT<realT>* > m_elementConcPtrs;
   std::vector<Array1dT<realT>* > m_elementReactionRatePtrs;
   std::vector<Array1dT<realT>* > m_frontPositionPtrs;
   
   // conversion factors
   realT m_convertToMolPerL;
   realT m_convertFromMolPerL;
   
   // parallel communication
   std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;
   
   // record species concentrations
   bool m_recordBrineSpecies;
   bool m_recordFrontSpecies;
   std::vector<Array1dT<realT>* > m_brineSpeciesConcPtrs; // only used if brine species are recorded
   std::vector<Array1dT<realT>* > m_frontSpeciesConcPtrs; // only used if front species are recorded


   // Cylindrical growth model
   bool m_useCylindricalGrowth;
   realT m_Ro;

   // Aperture flux limiter - to damp oscillatory concentrations
   bool m_useApertureBasedFluxLimiter;
   const Array1dT<realT>* m_FaceAperturesPtr;

   // timestep
   realT m_min_dt;
   realT m_max_dt;

   // verbose
   logLevelT m_verbosity;

};

void NewtonRaphson( ODRF::NonLinearSystemBase& nls, 
                    rArray1d& x, unsigned maxNumItrs, realT tolx, bool& rootFound, unsigned& numItrs, realT& f);

inline void NewtonRaphson( ODRF::NonLinearSystemBase& nls, 
                    rArray1d& x, unsigned maxNumItrs, realT tolx){
  bool rootFound(false); unsigned numItrs(0); realT f(0.0);
  NewtonRaphson( nls, x,  maxNumItrs,  tolx, rootFound, numItrs,  f);

}


#endif /* REACTIONFRONTSOLVER_H_ */
