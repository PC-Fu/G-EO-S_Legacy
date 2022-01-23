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

#include "QueryChemistryDatabase.h"
#include "ObjectManagers/UnitManager.h"
#include "Utilities/StringUtilities.h"
#include "PhysicsSolvers/SolverFactory.h"

//#include "Utilities/ReducedRowEchelonForm.h"

QueryChemistryDatabase::QueryChemistryDatabase(const std::string& name,
                                               ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

void QueryChemistryDatabase::ReadXML(TICPP::HierarchicalDataNode* const hdn)
{  
  SolverBase::ReadXML( hdn );

  int rank(0);
  #if GPAC_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  using namespace GPUnits;
//  UnitManager& theUnitManager = UnitManager::Instance();
  
  // get species
  //////////////
  if (rank == 0)  std::cout << "Species" << std::endl;
  HierarchicalDataNode* speciesNode = hdn->GetChild("Species");
  {
    m_species.names = speciesNode->GetStringVector("names", ",");
    m_species.size = m_species.names.size();
    if (rank == 0)  std::cout << "    ";
    for(unsigned i =0; i < m_species.size; ++i){
      if (rank == 0)  std::cout << m_species.names[i] << " ";
      //m_species.indices[m_species.names[i]] = i;
      Trim(m_species.names[i]);
    }
    if (rank == 0)  std::cout << std::endl;
    
  }

  // get solid phases
  ///////////////////
  if (rank == 0)  std::cout << "Solid Phases" << std::endl;
  HierarchicalDataNode* phasesNode = hdn->GetChild("SolidPhases");
  {
    m_solidphases.names = phasesNode->GetStringVector("names", ",");
    m_solidphases.size = m_solidphases.names.size();
    for(unsigned int i =0; i < m_solidphases.size; ++i){
      if (rank == 0)  std::cout << m_solidphases.names[i] << " ";
      Trim(m_solidphases.names[i]);
    }
  }

}


void QueryChemistryDatabase::Initialize( PhysicalDomainT& domain , SpatialPartition& partition  ){
  ChemistryManager& chemistryManager = ChemistryManager::Instance();
  std::vector<GPChemistry::ChemicalReaction>& solidPhaseReactionList = chemistryManager.m_phases;

  m_reactionSpecies.insert(m_species.names.begin(),m_species.names.end());
  
  for(unsigned i =0; i < m_solidphases.size; ++i){
  	std::vector<GPChemistry::ChemicalReaction>::iterator reactionItr = solidPhaseReactionList.begin();
  	bool notFound = true;
  	while( notFound && reactionItr != solidPhaseReactionList.end() ){
  	  notFound = !streq(m_solidphases.names[i],reactionItr->name);
  	  if(notFound) ++reactionItr;
  	}
  	
  	if(notFound){
  	  std::cout << "Warning solid phase " << m_solidphases.names[i] << " was not found in chemical database" << std::endl;
    }else{
  		m_reactionSpecies.insert(reactionItr->species.begin(),reactionItr->species.end()); 
  	}
  }
  
  chemistryManager.BuildReactionList(m_reactionSpecies, m_solutionReactions,m_solidPhaseReactions);

  
  std::cout << std::endl;
  std::cout << "Species" << std::endl;
  for(std::set<std::string>::iterator itr = m_reactionSpecies.begin(); 
      itr != m_reactionSpecies.end(); ++itr)
  {
    std::cout  << "    [" << *itr <<"]" << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Solution Reactions" << std::endl;
  for(std::vector<GPChemistry::ChemicalReaction>::size_type i =0; i < m_solutionReactions.size(); ++i)
  {
    std::cout  << "    " << m_solutionReactions[i].name << ": " << m_solutionReactions[i].formula << std::endl;
    std::cout  << "        " ;
    for(sArray1d::size_type j =0; j < m_solutionReactions[i].species.size(); ++j){
      std::cout  << m_solutionReactions[i].species[j] << ", ";	
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  
  
  std::cout << "Solid phase reactions" << std::endl;
  for(std::vector<GPChemistry::ChemicalReaction>::size_type i =0; i < m_solidPhaseReactions.size(); ++i)
  {
    std::cout  << "    " << m_solidPhaseReactions[i].name << ": " << m_solidPhaseReactions[i].formula << std::endl;
    std::cout  << "        " ;
    for(std::vector<GPChemistry::ChemicalReaction>::size_type j =0; j < m_solidPhaseReactions[i].species.size(); ++j){
      std::cout  << m_solidPhaseReactions[i].species[j] << ", ";	
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  exit(0);

}


REGISTER_SOLVER( QueryChemistryDatabase )
