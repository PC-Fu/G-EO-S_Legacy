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
 * @file ChemistryManager.h
 * @author walsh24
 * @date July 20, 2011
 */

#ifndef CHEMISTRYMANAGER_H_
#define CHEMISTRYMANAGER_H_

#include "Common/typedefs.h"
#include "Utilities/Utilities.h"
#include "DataStructures/Tables/Table.h"

#include "IO/ticpp/HierarchicalDataNode.h"


//#include "ChemistryData.h"


#include <map>
#include <utility>

#if GPAC_MPI
#include <mpi.h>
#endif

#include "DataStructures/Tables/TableTypes.h"

namespace GPChemistry
{

 /*
  * @brief Returns the base 10 log of the reaction rate constant
  * 
  * @param A[5] array of equilibrium constant coefficients, such that
  *    log_10 K = A0 + A1*T_K + A2*T_K^{-1} + A3*log_10(T_K) + A4*T_K^{-2} 
  *    where T_K is the temperature in Kelvin
  * @param TC Temperature in Celsius
  */
  inline
  realT Log10K(realT A[5], realT TC){
	realT TK = TC + 273.15;
	return A[0] + A[1]*TK + A[2]/TK + A[3]*log10(TK) + A[4]/(TK*TK);
  }
  
  inline
  realT Log10K(rArray1d A, realT TC){
	realT TK = TC + 273.15;
	return A[0] + A[1]*TK + A[2]/TK + A[3]*log10(TK) + A[4]/(TK*TK);
  }

	
  struct MasterSpecies {
    std::string element;
    std::string species;
    realT alkalinity;
    realT gfw;  // gram-formula-weight
    std::string gfw_formula; // gram-formula-weight formula
  };

  struct ChemicalReaction {
    std::string name;
    std::string formula;
    realT llnl_gamma;
    realT log_K;   // log_10 of equilibrium constant @ 25C
    realT delta_H;                     
    realT log_K_coeffs[5]; // log_10 K = A0 + A1*T_K + A2*T_K^{-1} + A3*log_10(T_K) + A4*T_K^{-2} 
    
    // Species indicies, valences and stoichiometric coefficients for the reaction
    // -ve stoichiometric coefficient indicates a reactant, +ve indicates a product
    sArray1d species;
    lArray1d speciesIndices; 
    rArray1d valences;
    rArray1d stoich_coeffs;
    
    ChemicalReaction():
      name(), formula(),
      llnl_gamma(0.0),log_K(0.0),delta_H(0.0)
    {
      for(int i =0; i < 5; ++i) log_K_coeffs[i] = 0.0;
    }
    
    
    // returns the equilibrium constant for the reaction
    realT EquilibriumConstant(realT TC)
    {
      realT rv = 0.0;
      if ( !isZero( log_K_coeffs[0]) || !isZero( log_K_coeffs[1] ) ||
           !isZero( log_K_coeffs[2]) || !isZero( log_K_coeffs[3] ) ||
           !isZero( log_K_coeffs[4] ) ){
    	rv = Log10K(log_K_coeffs, TC);
      } else {
      	rv = log_K; // FIXME need to add debye-huckle
      }
      return rv;
    };
  };

      
  /// Return a count of the elements named in a chemical species as a map, eg.
  /// Fe(OH)2 -> {{"Fe",1},{"O",2},{"H",2}}
  void GetElementCountFromFormula(std::string formula, std::map<std::string, realT>& elementCount);
  
  /// Returns the valence of a species string
  /// "H+"-> 1, "Fe+++" -> 3, "SO4-2" -> -2
  realT GetValence(std::string formula);

  inline
  rArray1d GetValences(const sArray1d& formula){
    unsigned n = formula.size();
    rArray1d valences(n);
    for(unsigned i = 0; i < n; ++i){
      valences[i] = GetValence(formula[i]);
    }
    return valences;
  }
  
  void ReadXMLEquilibriumData(TICPP::HierarchicalDataNode* hdn,
                              std::map< std::string, unsigned>& speciesIndicies,
                              std::map< std::string, rArray1d >& equations,
                              std::map< std::string, rArray1d >& log10Kdata);

  // m = molality of the species (mol/Kg)
  // Z = valence of the species
  template <typename ArrayT>
  inline realT IonicStrength(const ArrayT& m,const rArray1d& Z){
    realT sum =0.0;
    for(unsigned i =0; i < Z.size();++i) sum += fabs(m[i])*Z[i]*Z[i]; // fabs as a safety net
    return 0.5*sum;
  }

  inline
  rArray1d IonicStrengthCoefficients(const sArray1d& formula){
    unsigned n = formula.size();
    rArray1d iac(n);
    for(unsigned i = 0; i < n; ++i){
      realT v = GetValence(formula[i]);
      iac[i] = 0.5*v*v;
    }
    return iac;
  }


}


class ChemistryManager
{
public:

  static ChemistryManager& Instance()
  {
      static ChemistryManager theChemistyManager;

    return theChemistyManager;
  }
  void ReadXML(TICPP::HierarchicalDataNode* hdn);

  void ReadPhreeqcDatabase(const std::string& filename);
  
  void ReadSITDatabase(const std::string& filename);

  void BuildReactionList(std::set<std::string>& reactionSpecies,
                         std::vector<GPChemistry::ChemicalReaction>& activeSolutionReactions,
                         std::vector<GPChemistry::ChemicalReaction>& activeSolidPhaseReactions);

  std::vector<GPChemistry::MasterSpecies> m_solutionMasterSpecies;  
  std::vector<GPChemistry::MasterSpecies> m_exchangeMasterSpecies;  
  std::vector<GPChemistry::MasterSpecies> m_surfaceMasterSpecies;
  std::vector<GPChemistry::ChemicalReaction> m_solutionReactions;
  std::vector<GPChemistry::ChemicalReaction> m_phases;
  std::vector<GPChemistry::ChemicalReaction> m_exchangeReactions;
  std::vector<GPChemistry::ChemicalReaction> m_surfaceReactions;

  std::map<std::string, localIndex > m_speciesIndexMap;

  // Lookup tables 
  Table1D m_DebyeHuckelA; /// temperature (C) lookup table for debye_huckel_A parameter
  Table1D m_DebyeHuckelB; /// temperature (C) lookup table for debye_huckel_B parameter
  Table1D m_WATEQBdot; /// temperature (C) lookup table for WATEQ B dot parameter
  std::map<std::string,realT> m_atomicWeights; /// element atomic weights.
  std::map<std::string,realT> m_ionSizeParameters; /// species ion size parameters for use in extended-debye-huckel equation
  std::map<std::pair<std::string,std::string>,realT > m_SITParameters; /// specific ion interaction theory parameters for use in the BSG IAC
  
private:

  ChemistryManager():
    m_solutionMasterSpecies(),
    m_solutionReactions(), m_phases(), m_exchangeReactions(), m_surfaceReactions(),
    m_speciesIndexMap(),
    m_isReducingEnvironment(false)
  {
    ResetLookupTables();
  };
  ~ChemistryManager() {}
  ChemistryManager( const ChemistryManager& );
  ChemistryManager& operator=( const ChemistryManager& );
  
  void RecordReactionSpecies(std::vector<GPChemistry::ChemicalReaction>& reactions);
  void UpdateReactionSpeciesIndicies(std::vector<GPChemistry::ChemicalReaction>& reactions);
  bool CheckIfReactionIsPossible(const GPChemistry::ChemicalReaction& reaction,
                               std::set<std::string>& species,
                               std::vector<GPChemistry::ChemicalReaction>& reactionList,
                               bool recordMasterSpeciesReactions = false);
  void ResetLookupTables();

  int m_rank;
  bool m_isReducingEnvironment;
  
};

/// Activity functions
//////////////////////

namespace{
 const realT ln10 = log(10.0);
}

class ActivityFunction{
  public:
    ActivityFunction(){};
    ~ActivityFunction(){/*empty*/};
    virtual void CalculateActivities(rArray1d& activities,              
                                     const rArray1d& molalities){
      activities.resize(molalities.size());
      for(unsigned i =0; i < molalities.size(); ++i)
          activities[i] =  IAC(molalities, i)*molalities[i];
    };
    virtual realT CalculateActivity(const rArray1d& molalities, int speciesIndex){
          return IAC(molalities,speciesIndex)*molalities[speciesIndex];
    };
    
    /// ion activity coefficient
    virtual realT IAC(const rArray1d& molalities, int i){ 
      return pow(10, Log10IAC(molalities, i) ); 
    };

    virtual realT Log10IAC(const rArray1d& molalities __attribute__((unused)), int speciesIndex __attribute__((unused))){return 0.0;};
};

class ExtendedDebyeHuckelActivityFunction:public ActivityFunction{
  public:

    ExtendedDebyeHuckelActivityFunction(const rArray1d& valences=rArray1d(),
                                        realT A=0.5, realT B=0.33,const rArray1d& ao=rArray1d(), realT Bdot=0.04);
    
    ExtendedDebyeHuckelActivityFunction(const sArray1d& species,realT Tc);

    ~ExtendedDebyeHuckelActivityFunction(){};
    virtual void CalculateActivities(rArray1d& activities,              
                                     const rArray1d& molalities){
      activities.resize(molalities.size());
      realT sqrtI = sqrt(GPChemistry::IonicStrength(molalities,  m_valences));
      for(unsigned i =0; i < molalities.size(); ++i){
        realT iac = pow(10,
                        sqrtI*(m_Bdot*sqrtI - m_A*m_valences[i]*m_valences[i]/(1.0+m_ao[i]*m_B*sqrtI) )
                        );
        activities[i] =  iac*molalities[i];
      }
    };

    virtual realT Log10IAC(const rArray1d& molalities, int i){ 
      realT sqrtI = sqrt(GPChemistry::IonicStrength(molalities,  m_valences));
      return sqrtI*(m_Bdot*sqrtI - m_A*m_valences[i]*m_valences[i]/(1.0+m_ao[i]*m_B*sqrtI) );
    };

    virtual realT LogIAC(const rArray1d& molalities, int i){ 
      realT sqrtI = sqrt(GPChemistry::IonicStrength(molalities,  m_valences));
      return ln10*sqrtI*(m_Bdot*sqrtI - m_A*m_valences[i]*m_valences[i]/(1.0+m_ao[i]*m_B*sqrtI) );
    };

    virtual realT LogIAC(const realT* molalities, int i, int nC __attribute__((unused)) ){
      realT sqrtI = sqrt(GPChemistry::IonicStrength(molalities,  m_valences));
      return ln10*sqrtI*(m_Bdot*sqrtI - m_A*m_valences[i]*m_valences[i]/(1.0+m_ao[i]*m_B*sqrtI) );

    }

    virtual realT dLogIACdmu(const rArray1d& molalities, int i){ 
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      realT b = 1.0+m_ao[i]*m_B*sqrtI;
      return ln10*(m_Bdot-m_A*m_valences[i]*m_valences[i]*(b-sqrtI)/(2*sqrtI*b*b));
    };

    virtual realT dLogIACdmu(const realT* molalities, int i, int nC __attribute__((unused)) ){
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      realT b = 1.0+m_ao[i]*m_B*sqrtI;
      return ln10*(m_Bdot-m_A*m_valences[i]*m_valences[i]*(b-sqrtI)/(2*sqrtI*b*b));
    }

  private:
    rArray1d m_valences;
    realT m_A;
    realT m_B;
    realT m_Bdot;
    rArray1d m_ao;
  

};

class DaviesActivityFunction:public ActivityFunction{
  public:

    DaviesActivityFunction(const rArray1d& valences=rArray1d(),realT A=0);
    DaviesActivityFunction(const sArray1d& species,realT Tc);
    DaviesActivityFunction(const DaviesActivityFunction& rhs);
    ~DaviesActivityFunction(){};

    virtual void CalculateActivities(rArray1d& activities,              
                                     const rArray1d& molalities){
      activities.resize(molalities.size());
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      realT c = -m_A*(sqrtI/(1.0+sqrtI) +0.3*I ) ;
      for(rArray1d::size_type i =0; i < molalities.size(); ++i){
        realT iac = pow(10.0, c*m_valences[i]*m_valences[i] ); 
        activities[i] =  iac*molalities[i];
      }
    };

    // 
    virtual realT Log10IAC(const rArray1d& molalities, int i){ 
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      return -m_A*m_valences[i]*m_valences[i]*(sqrtI/(1.0+sqrtI) -0.3*I  );
    };

    // natural log of the ion activity product
    virtual realT LogIAC(const rArray1d& molalities, int i){ 
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      return -m_A*ln10*m_valences[i]*m_valences[i]*(sqrtI/(1.0+sqrtI) -0.3*I );
    };

    virtual realT LogIAC(const realT* molalities, int i, int nC __attribute__((unused)) ){
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      return -m_A*ln10*m_valences[i]*m_valences[i]*(sqrtI/(1.0+sqrtI) -0.3*I );
    };

    virtual realT dLogIACdmu(const rArray1d& molalities, int i){ 
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      return -ln10*m_A*m_valences[i]*m_valences[i]*(1.0/(2*sqrtI*(1.0+sqrtI)*(1.0+sqrtI)) - 0.3); 
    };

    virtual realT dLogIACdmu(const realT* molalities, int i, int nC __attribute__((unused)) ){
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      return -ln10*m_A*m_valences[i]*m_valences[i]*(1.0/(2*sqrtI*(1.0+sqrtI)*(1.0+sqrtI)) - 0.3); 
    };
    
  private:
    rArray1d m_valences;
    realT m_A;

};

class SITActivityFunction:public ActivityFunction{
  public:

    SITActivityFunction(const sArray1d& species,realT Tc,
                        const std::string& background_cation,
                        const std::string& background_anion);

    SITActivityFunction(const SITActivityFunction& rhs);

    ~SITActivityFunction(){};
    virtual void CalculateActivities(rArray1d& activities,
                                   const rArray1d& molalities){
      activities.resize(molalities.size());
      for(unsigned i =0; i < molalities.size(); ++i){
        realT iac = pow(10,Log10IAC(molalities, i));
          activities[i] =  iac*molalities[i];
      }
    };

    virtual realT Log10IAC(const rArray1d& molalities, int i){
      realT sqrtI = sqrt(GPChemistry::IonicStrength(molalities,  m_valences));

//      unsigned indx = 0;
      realT rv = -m_A*m_valences[i]*m_valences[i]*sqrtI/(1.0+m_Ba_j*sqrtI);
      if(m_valences[i] > 0.0){
        rv += m_e[i]*molalities[m_anionIndx];
      } else if(m_valences[i] < 0.0) {
        rv += m_e[i]*molalities[m_cationIndx];
      }

      return rv;
    };

    virtual realT LogIAC(const rArray1d& molalities, int i){
      return ln10*Log10IAC(molalities, i);
    };

    // derivative of natural log of IAC for species i wrt the ionic strength
    virtual realT dLogIACdmu(const rArray1d& molalities, int i){
      realT I = GPChemistry::IonicStrength(molalities,  m_valences);
      realT sqrtI = sqrt(I);
      realT b = 1.0+m_Ba_j*sqrtI;
      return -ln10*m_A*m_valences[i]*m_valences[i]/(2*sqrtI*b*b);
    };

    // derivative of natural log of IAC for species i wrt the concentration of species j
    virtual realT dLogIACdm(const rArray1d& molalities, unsigned int i, unsigned int j){
      realT rv = 0.0;
      if(m_valences[i] > 0.0 && j == m_anionIndx){
        rv = ln10*m_e[i];
      } else if(m_valences[i] < 0.0 && j == m_cationIndx){
        rv = ln10*m_e[i];
      }
      return rv;
    };

  private:
    rArray1d m_valences;
    realT m_A;
    realT m_Ba_j;
    rArray1d m_e; // ionic interaction coefficients

    std::string m_background_cation;
    std::string m_background_anion;

    unsigned m_cationIndx; // index of background cation
    unsigned m_anionIndx;  // index of background anion


};

#endif /* CHEMISTRYMANAGER_H_ */
