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
 * @file ChemicalEquilibriumSolver.cpp
 * @author walsh24
 * @date July 25, 2011
 */
 
#include <cmath>
#include <vector>

#include "PhysicsSolvers/SolverFactory.h"
#include "FindEquilibriumNew.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"
#include "ObjectManagers/ChemistryManager.h"
#include "ObjectManagers/UnitManager.h"
#include "Utilities/ChemistryUtilities.h"

using namespace GPChemistry;

using namespace ODRF;


namespace{

  template <class T>
  void GetVectorAttribute(std::vector<T>& vals,TICPP::HierarchicalDataNode* hdn, std::string attribute,std::string sep = ",", bool doTrim = true){
    std::string vStr = hdn->GetAttributeString(attribute);
    if( !vStr.empty() ){
      sArray1d vStrs = Tokenize(vStr,sep);
      vals.resize(vStrs.size());
      for(unsigned i =0; i < vStrs.size(); ++i){
        if(doTrim) Trim(vStrs[i]," \t\n");
        vals[i] = fromString<T>(vStrs[i]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////////////

// Upate field with function

ChemicalEquilibriumSolver::ChemicalEquilibriumSolver( const std::string& name,
                                                      ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

ChemicalEquilibriumSolver::~ChemicalEquilibriumSolver()
{
}

/*
 * <ChemicalEquilibriumSolver name="CaUpdate"       * name of the solver
 *            objecttype="Face"                  * location of field (Default=Node)
 *            species="Scalar Scalar" />   * variableTypes (assumed scalar for all if omitted) 
 */
void ChemicalEquilibriumSolver::ReadXML(TICPP::HierarchicalDataNode* const hdn)
{  
  SolverBase::ReadXML(hdn);

  int rank(0);
  #if GPAC_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  std::string objectTypeStr = hdn->GetAttributeStringOrDefault("object", PhysicalDomainT::FiniteElementFaceManagerStr() );
  {
    std::string oldStr = hdn->GetAttributeStringOrDefault("objecttype", "" ); // throw error if using old syntax
    if(!oldStr.empty()) {
      throw GPException("UpdateFieldWithFunction: Attempting to set objecttype - use 'object' instead.");
    }
  }
  m_objectKey = PhysicalDomainT::GetObjectDataStructureConditionKey(objectTypeStr);

  m_regionName = hdn->GetAttributeStringOrDefault("regionname",""); // used for element regions only


  m_setNames = hdn->GetStringVector("setnames");
  if(m_setNames.empty()) m_setNames = hdn->GetStringVector("setname");
  
  // flags
  m_recordBrineSpecies = hdn->GetAttributeOrDefault("recordBrineSpecies", false);
  
  // Temperature
  ///////////////
  m_Tc = hdn->GetAttributeOrDefault("temperature", "25");
  
  // Species
  //////////////
  if (rank == 0)  std::cout << "Species" << "\n";
  HierarchicalDataNode* speciesNode = hdn->GetChild("Species");
  if(!speciesNode)  
      throw GPException("ChemicalEquilibriumSolver: Must have Species defined in the input file");
  {
    m_species.ReadXML(speciesNode);

  }

 
  // Elements
  /////////////////////
  if (rank == 0)  std::cout << "Elements" << std::endl;
  HierarchicalDataNode* elementNode = hdn->GetChild("Elements");
  if(!elementNode)  
      throw GPException("ChemicalEquilibriumSolver: Must have Elements defined in the input file");
  {
    m_elements.ReadXML(elementNode,m_species);
  }



  // Equilibrium data (solid and aqueous)
  ///////////////////////////////////////
  if (rank == 0)  std::cout << "Equilibrium Data" << std::endl;
  HierarchicalDataNode* eqNode = hdn->GetChild("EquilibriumData");
  if(!eqNode)  
      throw GPException("ChemicalEquilibriumSolver: Must have EquilibriumData defined in the input file");
  {
    m_equilibriumData.ReadXML(eqNode,m_species);

    // additional equilibrium constraints
	sArray1d eqCondStrs;
    GetVectorAttribute(eqCondStrs,eqNode, "equilibrium_constraints",",");


    // resize constraint equations
    localIndex nA = eqCondStrs.size()+1; // +1 = charge balance
    localIndex n = m_species.size;
    localIndex nE = m_elements.size;
    nA = nE + 1; // override user defined constraints - see below
    m_AA.resize2(nA,n);
    m_AA_e.resize2(nA,nE);
    m_bb.resize(nA);

/*
 * Here I wanted to have user defined constraint equations - looks like I didn't finish it
    for(unsigned ii =0; ii < nA-1; ++ii){
        sArray1d lhs_rhs = Tokenize(eqCondStrs[ii],"="); // lhs,rhs
        CHEM_UTILS::ParseSpeciesConstraintStr(lhs_rhs[0],m_elements ,m_species,  m_AA[ii]);
        CHEM_UTILS::ParseSpeciesConstraintStr(lhs_rhs[1], m_elements , m_species, m_AA[ii],false);
        CHEM_UTILS::ParseElementConstraintStr(lhs_rhs[0],m_elements, m_AA_e[ii]);
        CHEM_UTILS::ParseElementConstraintStr(lhs_rhs[1], m_elements, m_AA_e[ii],false);
    }

*/
    // have added this in in the mean time
    for(unsigned i =0; i < nA-1; ++i){
      for(unsigned j=0; j < n; ++j){
         m_AA(i,j) = m_elements.species[i][j];
      }
      for(unsigned j=0; j < nE; ++j){
          m_AA_e(i,j) = 0;
      }
      m_AA_e(i,i) = 1;

      m_elements.brineConcentrations[i] *= m_convertToMolPerL;
      m_bb[i] = m_elements.brineConcentrations[i];
    }

    // Charge balance
    for(unsigned j=0; j < n; ++j) m_AA(nA-1,j) = m_species.valences[j];
    m_bb[nA-1] = 0;
    realT maxValence = 1.0;  // rescale
    for(unsigned j=0; j < n; ++j){
  	  if( fabs( m_AA(nA-1,j) ) > maxValence) maxValence = fabs(m_AA(nA-1,j));
    }
    for(unsigned j=0; j < n; ++j) m_AA(nA-1,j) /= maxValence;
    m_bb[nA-1] /= maxValence;  // redundant but included in case alkalinity added in future.

  }


    
}

void ChemicalEquilibriumSolver::RegisterFields( PhysicalDomainT& domain )
{

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);

  // register elements
  for(sArray1d::size_type i =0; i < m_elements.size; ++i){
	  objectManager.AddKeylessDataField( FieldInfo::realField,  m_elements.names[i], true, true );
  }


  objectManager.AddKeylessDataField( FieldInfo::realField, "pH", true, true );

  // record data pointers
  for(sArray1d::size_type i =0; i < m_elements.size; ++i){
    Array1dT<realT>* ptr = &objectManager.GetFieldData<realT>(m_elements.names[i]);
    m_elementConcPtrs.push_back(ptr);
  }

  // brine species
  if(m_recordBrineSpecies){
    m_brineSpeciesConcPtrs.resize(m_species.size);
    for(sArray1d::size_type j =0; j < m_species.size; ++j){
      std::string sp = m_species.names[j];
      std::replace( sp.begin(), sp.end(), '+', 'p'); //+ -> p
      std::replace( sp.begin(), sp.end(), '-', 'm'); //- -> m
      std::replace( sp.begin(), sp.end(), '(', '_'); //( -> _
      std::replace( sp.begin(), sp.end(), ')', '_'); //) -> _

      std::string spFrtStr = sp+"_br";
      objectManager.AddKeylessDataField( FieldInfo::realField,  spFrtStr, true, true );
      m_brineSpeciesConcPtrs[j]
            = &objectManager.GetFieldData<realT>(spFrtStr);

    }
  }


  m_pHFieldPtr = &objectManager.GetFieldData<realT>("pH");

}

void ChemicalEquilibriumSolver::Initialize( PhysicalDomainT& domain , SpatialPartition& partition ){


  // Unit scale
  /////////////
  {
    using namespace GPUnits;
    UnitManager& um = UnitManager::Instance(); 
    m_convertToMolPerL = um.ConvertTo("mol/l",1);
    m_convertFromMolPerL=  1.0/m_convertToMolPerL;
  }


  // Reaction Fronts
  ///////////////////

  unsigned n = m_species.size;
  //unsigned ne = m_elements.size;
  //unsigned na = m_AA.Dimension(0); // aqueous species constraints including charge balance
  unsigned nk = m_species.eqConstants.size();
  
  // equilibrium data for fronts
  m_KK = rArray2d(nk,n);
  m_pp= rArray1d(nk);
  for(unsigned i =0; i < nk; ++i){
    std::string& kstr = m_species.eqConstants[i];
    for(unsigned j=0; j < n; ++j) m_KK(i,j) = m_equilibriumData.equations[kstr][j];  
    rArray1d& log10Kdata = m_equilibriumData.log10Kdata[kstr];
    m_pp[i] = log(10)* Log10K(log10Kdata, m_Tc); // convert log_10 data to log_e

    // rescale
    realT kMax = fabs(m_KK(i,0));
    for(unsigned j=0; j < n; ++j){
    	if(kMax < fabs(m_KK(i,j)) ) kMax = fabs(m_KK(i,j));
    }
    for(unsigned j=0; j < n; ++j) m_KK(i,j) /= kMax;
    m_pp[i] /= kMax;
  }


  
  // Initial guess for brine species
  //////////////////////////////////
  
  m_species.brineConcentrations = rArray1d (2*n+1,1e-7);  
  // ~ mean of elements
  for(unsigned i=0; i < m_elements.size; ++i){
    realT net = 0.0;
    for(unsigned ii=0; ii < n; ++ii) net += m_elements.species[i][ii];
    
    for(unsigned ii=0; ii < n; ++ii){
      if(m_elements.species[i][ii] > 0) m_species.brineConcentrations(ii) += m_elements.brineConcentrations(i)*m_elements.species[i][ii]/net;
    }
  }
  
  
  //m_species.brineConcentrations = rArray1d (2*n+1,1e-7);  
  for(unsigned i =0; i < n; ++i) m_species.brineConcentrations(i+n) = log(std::max(m_species.brineConcentrations(i),1e-64));

  // build rhs of brine system
  UpdateBrineConstraints();

  // Brine system
  m_brineSystem.Initialize(m_KK, m_pp, m_AA, m_bb, m_species.names,m_Tc);


  // pH
  m_H_index = m_species.indices["H+"];

}


inline
void ChemicalEquilibriumSolver::UpdateBrineConstraints(){

  unsigned ne = m_elements.size;
  unsigned na = m_AA.Dimension(0);

  for(unsigned i =0; i < na-1; ++i){ // last row in A x = b is reserved for charge balance
    realT val = 0.0;
	for(unsigned j=0; j < ne; ++j){
	  val += m_AA_e(i,j)*m_elements.brineConcentrations[j];
	}
	m_bb[i] = val;
  }
}



/**
 * 
 * 
**/

double ChemicalEquilibriumSolver::TimeStep( const realT& time ,
                                const realT& dt ,
                                const int cycleNumber,
                                PhysicalDomainT& domain ,
                                const sArray1d& namesOfSolverRegions ,
                                SpatialPartition& partition ,
                                FractunatorBase* const fractunator )
{  

  m_stabledt.m_maxdt = 0.9*std::numeric_limits<double>::max();
//  iArray1d& faceGhostRank       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();


  int rank(0);
	//realT t1;
  #if GPAC_MPI
	    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);
  if(m_setNames.empty()){
	for(localIndex i = 0; i < objectManager.DataLengths(); ++i){
	  UpdateEquilibrium(i);
	}
  } else {
	for( sArray1d::size_type i =0; i < m_setNames.size(); ++i){
	  lSet& subset = objectManager.GetSet(m_setNames[i]);
      for( lSet::const_iterator si=subset.begin() ; si!=subset.end() ; ++si )
		{
    	  UpdateEquilibrium(*si);
		}
	}
  }
	
  return dt;
}

void ChemicalEquilibriumSolver::UpdateEquilibrium(localIndex indx ){

	// Get elements species from object
	for(unsigned i =0; i < m_elements.size; ++i){
	     m_elements.brineConcentrations(i) = (*m_elementConcPtrs[i])[indx]*m_convertToMolPerL;
	}

	// Update brine equilibrium equations
	UpdateBrineConstraints();

    // Find Equilibrium
	const unsigned maxNumItrs = 1000;
	const realT tolx = 2*std::numeric_limits<realT>::epsilon();
	//std::cout << tolx << std::endl;
	//const int n = m_species.size;
	bool rootFound(false); unsigned numItrs(0); realT ff(0.0);

	// Brine system
	/////////////////
	//std::cout << "Brine System"  << std::endl;
	NewtonRaphson( m_brineSystem, m_species.brineConcentrations, maxNumItrs, tolx, rootFound, numItrs, ff);

	// copy solution back to object
	for(unsigned i =0; i < m_elements.size; ++i){
      realT val = 0.0;
      for(unsigned j =0; j < m_species.size;++j){
        val += m_elements.species[i][j]*m_species.brineConcentrations(j);
      }
	  (*m_elementConcPtrs[i])[indx] = val*m_convertFromMolPerL;
	}

	// pH
	(*m_pHFieldPtr)[indx] = -log10(m_species.brineConcentrations(m_H_index));
// (*m_pHFieldPtr)[indx] = -log10( exp(m_species.brineConcentrations(m_H_index+ m_species.size) ));
}


REGISTER_SOLVER( ChemicalEquilibriumSolver )
