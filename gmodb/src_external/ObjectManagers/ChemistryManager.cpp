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

#include "ChemistryManager.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <functional>

#include "DataStructures/Tables/Table.h"
#include "DataStructures/Tables/TableTypes.h"

using namespace GPChemistry;

namespace{
  bool XStartsWithY(const std::string& X, const std::string& Y){
    return (X.compare(0,Y.size(),Y) == 0);
  }
  bool XContainsY(const std::string& X, const std::string& Y){
    return (X.find(Y) != std::string::npos);
  }
  
  std::string FirstWord(const std::string& str){
    size_t indxA = str.find_first_not_of(" \t");
    size_t count = str.find_first_of(" \t",indxA) - indxA;
    return str.substr(indxA,count);
  }

  std::string FirstWordAfterEquals(const std::string& str){
    sArray1d sv = Tokenize(str,"=");
    if(sv.size() < 2) return "";
    Trim(sv[1]);
    size_t indx = sv[1].find_first_of(" \t");
    return sv[1].substr(0,indx);
  }
  
  bool isNumber(const std::string& str, realT& num){
  	std::istringstream iss( str );
  	iss >> num;
  	if(!iss) return false;
  	return(iss.rdbuf()->in_avail() == 0);
  }
  
  /// returns valence and changes element string to a standard form (Fe+++ -> Fe+3)
  realT getValenceAndAdjustStr(std::string& str){
  	Trim(str);
  	size_t indxA = str.find_first_of("+");
  	realT sgn = 1.0;
  	if(indxA ==  std::string::npos){
      indxA = str.find_first_of("-");
  	  sgn = -1.0;
  	  if(indxA ==  std::string::npos) return 0.0;
  	} 
  	realT val;
  	if(!isNumber ( str.substr(indxA), val ) ){
  	  int count = str.size() - indxA;
  	  val =  sgn*(count);
  	  str = str.substr(0,indxA+1); // includes sign
  	  if(int(fabs(val)) > 1.0) str += toString<int>(count);
  	}
  	
  	return val;
  }
  
  /// returns stoichiometry and removes coefficient from element string (if present)
  realT getStoichiometry(std::string& species){
    int indx = species.find_first_not_of("0123456789.-+");
    if(indx > 0) species.insert(indx," "); // needed to separate 4e- otherwise next step wont work

    std::istringstream iss( species );
    realT stoich = 1.0;
    iss >> stoich >> species; // will only change stoich and species if species starts with numeric expression

    return stoich;
  }
  
  void RemoveEmptyStrings(sArray1d& strVector){
  	std::vector<std::string>::iterator it = remove_if(strVector.begin(), strVector.end(),mem_fun_ref(&std::string::empty));
    // erase the removed elements
    strVector.erase(it, strVector.end());
  }
  
  template <class T>
  void GetVectorAttribute(std::vector<T>& vals,TICPP::HierarchicalDataNode* hdn, std::string attribute,std::string sep = ",", bool doTrim = true){
    std::string vStr = hdn->GetAttributeString(attribute);
    sArray1d vStrs = Tokenize(vStr,sep);
    vals.resize(vStrs.size());
    for(sArray1d::size_type i =0; i < vStrs.size(); ++i){
    	if(doTrim) Trim(vStrs[i]," \t\n");
    	vals[i] = fromString<T>(vStrs[i]);
    }	
  }

  
  enum READMODE
  {
    SOLUTION_MASTER_SPECIES,
    SOLUTION_SPECIES,
    PHASES,
    EXCHANGE_MASTER_SPECIES,
    EXCHANGE_SPECIES,
    SURFACE_MASTER_SPECIES,
    SURFACE_SPECIES,
    RATES,
    UNDEF
  };
  
  
    
  void ReadPhreeqcReactions(READMODE mode,std::vector<GPChemistry::ChemicalReaction>& reactions,
                            std::ifstream& database, 
                            std::map<READMODE,std::pair<std::streampos,std::streampos> >& modeRanges){
    
    std::streampos linePosition = modeRanges[mode].first;
//    std::streampos endpos = modeRanges[mode].second;
    std::string inputline; 
                            	
    database.clear();
    database.seekg(linePosition);
    ChemicalReaction reaction;
    std::string prevline = "";
    
    while ( linePosition != modeRanges[mode].second ){
      getline(database,inputline);
      RemoveComments(inputline, '#');
      
      std::string dummyStr;
      
      if( XContainsY(inputline,"=") ){
      	
      	// record old reaction
      	if(!reaction.name.empty()) reactions.push_back(reaction);
      	
      	// start new reaction record
      	reaction = ChemicalReaction();
      	
      	if(mode == PHASES){
          reaction.name = prevline; 
      	} else if(mode == SOLUTION_SPECIES) {
          reaction.name = FirstWordAfterEquals(inputline); 
      	} else if(mode == EXCHANGE_SPECIES || mode == SURFACE_SPECIES ){
          reaction.name = FirstWord(inputline); 
      	}
      	Trim(reaction.name," \t\r");
        reaction.formula = inputline; Trim(reaction.formula," \t\r");
      	
      } else {
      	std::istringstream linestream(inputline);
      	
      	Trim(inputline," \t");
      	
      	if(XStartsWithY(inputline,"-analytical_expression") ||
           XStartsWithY(inputline,"analytical_expression") || 
           XStartsWithY(inputline,"-analytic") || 
           XStartsWithY(inputline,"-ae") || 
           XStartsWithY(inputline,"ae") || 
           XStartsWithY(inputline,"-a_e") ||
           XStartsWithY(inputline,"a_e") ){

           linestream >> dummyStr;
           for(int i =0; i < 5; ++i) linestream >> reaction.log_K_coeffs[i];

        } else if(XStartsWithY(inputline,"-log_K")||
                  XStartsWithY(inputline,"log_K") ||
                  XStartsWithY(inputline,"log_k") ||
                  XStartsWithY(inputline,"-logk") ||
                  XStartsWithY(inputline,"logk")  ){ 

           linestream >> dummyStr >> reaction.log_K;

        } else if(XStartsWithY(inputline,"-delta_h")||
                  XStartsWithY(inputline,"delta_h") ||
                  XStartsWithY(inputline,"-deltah") ||
                  XStartsWithY(inputline,"deltah") ){

           linestream >> dummyStr >> reaction.delta_H;
        } else if(XStartsWithY(inputline,"-llnl_gamma") ){
           linestream >> dummyStr >> reaction.llnl_gamma;
        }
      }
      
      prevline =  inputline;
      linePosition = database.tellg();
    }
    
    // record final reaction
    if(!reaction.name.empty()) reactions.push_back(reaction);
    
  }
  
  void ReadPhreeqcMasterSpecies(READMODE mode, std::vector<GPChemistry::MasterSpecies>& species,
                                std::ifstream& database, 
                                std::map<READMODE,std::pair<std::streampos,std::streampos> >& modeRanges){
                                 	
    std::streampos linePosition = modeRanges[mode].first;
    std::streampos endpos = modeRanges[mode].second;
    std::string inputline;
    database.clear();
    database.seekg(linePosition,std::ios_base::beg);
    while ( linePosition != endpos && !database.eof() ){
    	
      getline(database,inputline);
      RemoveComments(inputline, '#');
      Trim(inputline," \t\n");
     
      // exit(0);
 
      if(inputline.size() > 1){
        // std::cout << species.size() << " " << inputline << std::endl;
        std::istringstream linestream(inputline);
        MasterSpecies masterSpecies;
        linestream >> masterSpecies.element;
        linestream >> masterSpecies.species;
         
        masterSpecies.alkalinity = 0.0;
        masterSpecies.gfw_formula = "";
        masterSpecies.gfw = 0.0;
         
        linestream >> masterSpecies.alkalinity;
        linestream >> masterSpecies.gfw_formula;
        linestream >> masterSpecies.gfw;
  
        species.push_back(masterSpecies);
        //std::cout << masterSpecies.element << " " << masterSpecies.species << std::endl;
      }
      linePosition = database.tellg();
    } 
  }
  
  sArray1d SplitAtCapitals(std::string& elementsStr){
    sArray1d rv;

    size_t pos = elementsStr.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ");
    size_t npos = elementsStr.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ",pos+1);
    while(pos < elementsStr.size() )
    {
      rv.push_back(elementsStr.substr(pos,npos-pos) );
    
      pos = npos;
      npos = elementsStr.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ",pos+1);
    }
    return rv;
  }

  // "Ar2" -> "Ar" and returns 2
  realT GetSubscriptFromElementStr(std::string& elementStr){
    size_t pos = elementStr.find_first_of("0123456789.");
    realT rv = 1.0;
    if(pos <  elementStr.size() ){
      rv = fromString<realT>(elementStr.substr(pos) );
      elementStr = elementStr.substr(0,pos);
    }
    return rv;
  }

}  // namespace


void ChemistryManager::ReadXML(TICPP::HierarchicalDataNode* hdn){
  
   
    std::string cdbFilename = hdn->GetAttributeString("chemistry_database");
    std::string sitFilename = hdn->GetAttributeString("sit_database");
//    bool m_isReducing = hdn->GetAttributeOrDefault<bool>("reducing_environment",false);

    m_rank = 0;
    #if GPAC_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    #endif

    if(!cdbFilename.empty()){
      if (m_rank  == 0)
        std::cout << "Reading Chemistry Database: " << cdbFilename << std::endl;

      ReadPhreeqcDatabase(cdbFilename);
    }


    if(!sitFilename.empty()){
      if (m_rank  == 0)
        std::cout << "Reading SIT Database: " << sitFilename << std::endl;
    
      ReadSITDatabase(sitFilename);
    }
}

void ChemistryManager::ReadPhreeqcDatabase(const std::string& chemistryDatabaseFilename)
{

  
  std::map<READMODE,std::pair<std::streampos,std::streampos> > modeRanges;

  std::string solmsKey("SOLUTION_MASTER_SPECIES");
  std::string solKey("SOLUTION_SPECIES");
  std::string phasesKey("PHASES");
  std::string exmsKey("EXCHANGE_MASTER_SPECIES");
  std::string exKey("EXCHANGE_SPECIES");
  std::string surfmsKey("SURFACE_MASTER_SPECIES");
  std::string surfKey("SURFACE_SPECIES");
  std::string ratesKey("RATES");

  std::string inputline,prevline;

  READMODE mode = UNDEF;

  std::ifstream database;

  database.open(chemistryDatabaseFilename.c_str());

  if(!database.is_open()){
  	throw GPException("Error ReadPhreeqcDatabase: Could not open file  " + chemistryDatabaseFilename);
  }

  // std::cout<<"Reading Chemistry Database...\n";

  std::streampos linePosition =  database.tellg();

  while ( !database.eof() )
  {
    getline(database,inputline);
    RemoveComments(inputline, '#');

    //// Switch modes
    if( XStartsWithY(inputline,solmsKey) )
    {
        modeRanges[mode].second = linePosition;
        mode = SOLUTION_MASTER_SPECIES;
        getline(database,inputline);
        modeRanges[mode].first = database.tellg();
    } 
    else if( XStartsWithY(inputline,solKey) )
    {
        modeRanges[mode].second = linePosition;
        mode = SOLUTION_SPECIES;
        getline(database,inputline);
        modeRanges[mode].first = database.tellg();
    } 
    else if( XStartsWithY(inputline,phasesKey) )
    {
        modeRanges[mode].second = linePosition;
        mode = PHASES;
        getline(database,inputline);
        modeRanges[mode].first = database.tellg();
    } 
    else if( XStartsWithY(inputline,exmsKey) )
    {
        modeRanges[mode].second = linePosition;
        mode = EXCHANGE_MASTER_SPECIES;
        getline(database,inputline);
        modeRanges[mode].first = database.tellg();
    } 
    else if( XStartsWithY(inputline,exKey) )
    {
        modeRanges[mode].second = linePosition;
        mode = EXCHANGE_SPECIES;
        getline(database,inputline);
        modeRanges[mode].first = database.tellg();
    } 
    else if( XStartsWithY(inputline,surfmsKey) )
    {
        modeRanges[mode].second = linePosition;
        mode = SURFACE_MASTER_SPECIES;
        getline(database,inputline);
        modeRanges[mode].first = database.tellg();
    }
    else if( XStartsWithY(inputline,surfKey) )
    {
        modeRanges[mode].second = linePosition;
        mode = SURFACE_SPECIES;
        getline(database,inputline);
        modeRanges[mode].first = database.tellg();
    }
    else if( XStartsWithY(inputline,ratesKey) )
    {
        modeRanges[mode].second = linePosition;
        mode = RATES;
        getline(database,inputline);
        modeRanges[mode].first = database.tellg();
    }
    linePosition =  database.tellg();
  }
 
  //// Master species
  ReadPhreeqcMasterSpecies(SOLUTION_MASTER_SPECIES, m_solutionMasterSpecies, database, modeRanges);
  ReadPhreeqcMasterSpecies(EXCHANGE_MASTER_SPECIES, m_exchangeMasterSpecies, database, modeRanges);
  ReadPhreeqcMasterSpecies(SURFACE_MASTER_SPECIES, m_surfaceMasterSpecies, database, modeRanges);
 
  std::cout << "    Finished reading master species" << std::endl;
    
    
  //// Solution reactions
  ReadPhreeqcReactions(SOLUTION_SPECIES, m_solutionReactions, database, modeRanges);
  ReadPhreeqcReactions(PHASES, m_phases, database, modeRanges);
  ReadPhreeqcReactions(EXCHANGE_SPECIES, m_exchangeReactions, database, modeRanges);
  ReadPhreeqcReactions(SURFACE_SPECIES, m_surfaceReactions, database, modeRanges);
  
  
  // RATES
  {
  	//TODO - requires BASIC language parser
  }
  
  
  // Record reaction species, valences and stoichiometry and add to species index map
  RecordReactionSpecies(m_solutionReactions);
  RecordReactionSpecies(m_phases);
  RecordReactionSpecies(m_exchangeReactions);
  RecordReactionSpecies(m_surfaceReactions);
  
  // number reaction species
  std::map<std::string,localIndex>::iterator itr;
  std::map<std::string,localIndex>::iterator iEnd = m_speciesIndexMap.end();
  localIndex count = 0;
  for(itr = m_speciesIndexMap.begin(); itr != iEnd; ++itr){
  	itr->second = count;
  	count++;
  }
  
  // record species indicies in reaction
  UpdateReactionSpeciesIndicies(m_solutionReactions);
  UpdateReactionSpeciesIndicies(m_phases);
  UpdateReactionSpeciesIndicies(m_exchangeReactions);
  UpdateReactionSpeciesIndicies(m_surfaceReactions);
  std::cout << "    Finished reading reactions" << std::endl;
  
  // Finish up
  database.close();

  std::cout << "Done Reading Chemistry Database.\n";

}

/**
 * Record ionic interaction coefficients for specific-ion interaction theory (SIT)
 *
 * SIT parameter file is a set of cation/anion pairs arranged as a comma separated list in
 * the following format:
 * Cation, Anion, ionic interaction coefficient
 */
void ChemistryManager::ReadSITDatabase(const std::string& chemistryDatabaseFilename)
{

  std::ifstream database;

  database.open(chemistryDatabaseFilename.c_str());

  if(!database.is_open()){
    throw GPException("Error ReadSITDatabase: Could not open file  " + chemistryDatabaseFilename);
  }

  std::string inputline;

  while ( getline(database,inputline) )
  {

    sArray1d CationAnionInteraction = Tokenize(inputline,",");
    if(CationAnionInteraction.size() < 3 && CationAnionInteraction.size() > 0){
      throw GPException("Error ReadSITDatabase: Failed to read line :\n" +  inputline + "\n");

    } else if(CationAnionInteraction.size() >= 3) {
      std::string cation = CationAnionInteraction[0];
      std::string anion = CationAnionInteraction[1];
      realT ionicInteractionCoefficient = fromString<realT>(CationAnionInteraction[2]);

      RemoveSpaces(cation);Trim(cation);
      RemoveSpaces(anion);Trim(anion);

      // convert to standard form and confirm that we have a cation/anion pair
      realT vC = getValenceAndAdjustStr(cation);
      realT vA = getValenceAndAdjustStr(anion);

      if(vC <= 0 && vA > 0){
          // oops - valences are reversed
          std::string temp = cation;
          cation = anion;
          anion = temp;

      } else if (vC < 0 || vA > 0){
        std::cout << "Warning ReadSITDatabase: Detected Ion pair of same sign in database: \n"  + inputline + "\n";
      }

      m_SITParameters[std::pair<std::string,std::string>(cation,anion)] = ionicInteractionCoefficient;

    }

  }

}

/// Record reaction species, valences, stoichiometry, and add species to reaction map
void ChemistryManager::RecordReactionSpecies(std::vector<GPChemistry::ChemicalReaction>& reactions){
  for(size_t i =0; i < reactions.size(); ++i){
  	std::string formula = reactions[i].formula;
    sArray1d reactProd = Tokenize(formula,"=");
    sArray1d react = TokenizeSeq(reactProd[0]," +");
    sArray1d prod = TokenizeSeq(reactProd[1]," +");
    
    Trim(react," \t\r");
    Trim(prod," \t\r");
    
    RemoveEmptyStrings(react);
    RemoveEmptyStrings(prod);
    
    int n = react.size() + prod.size(); 
    sArray1d speciesStrs(n);
    lArray1d speciesIndices(n,0); // added later 
    rArray1d valences(n);
  	rArray1d stoich_coeffs(n);
    
    for(size_t j =0; j < react.size(); ++j ){
  	  std::string& species = react[j];

  	  realT stoich = getStoichiometry(species);
  	  realT val = getValenceAndAdjustStr(species);
  	  
  	 // std::cout << react[j] << " " << stoich <<" " << species << " " << val <<std::endl;
  	  
  	  m_speciesIndexMap[species]; // add entry to map
  	  speciesStrs[j] = species;
  	  valences[j] = val;
  	  stoich_coeffs[j] = -stoich; // -ve stoich indicates reactants
    }
    
    for(size_t j =0; j < prod.size(); ++j ){
  	  std::string& species = prod[j];
  	  
  	  realT stoich = getStoichiometry(species);
  	  realT val = getValenceAndAdjustStr(species);
  	  
  	 // std::cout << prod[j] << " " << stoich << " " << species << " " << val <<std::endl;
  	  
  	  m_speciesIndexMap[species]; // add entry to map
  	  int jj = j + react.size();
  	  speciesStrs[jj] = species;
  	  valences[jj] = val;
  	  stoich_coeffs[jj] = stoich;
    }
    
    reactions[i].species = speciesStrs;
    reactions[i].speciesIndices.resize(n);
    reactions[i].valences = valences;
    reactions[i].stoich_coeffs = stoich_coeffs;
  }
  
}

/// Records the reaction species index for each species in each reaction after the 
/// species index map has been constructed. 
void ChemistryManager::UpdateReactionSpeciesIndicies(std::vector<GPChemistry::ChemicalReaction>& reactions){
  for(size_t i =0; i < reactions.size(); ++i){
  	for(size_t j =0; j < reactions[i].species.size(); ++j ){
  		reactions[i].speciesIndices[j] = m_speciesIndexMap[ reactions[i].species[j] ];
    }
  }
}

bool ChemistryManager::CheckIfReactionIsPossible(const GPChemistry::ChemicalReaction& reaction,
                               std::set<std::string>& species,
                               std::vector<GPChemistry::ChemicalReaction>& reactionList,
                               bool recordMasterSpeciesReactions){

  bool aReactionIsPossible = (recordMasterSpeciesReactions || reaction.species.size() > 2);
  if(!recordMasterSpeciesReactions && reaction.species.size() == 2){
  	aReactionIsPossible = !streq(reaction.species[0],reaction.species[1]);
  }
  	
  

  if(aReactionIsPossible){
  	
    bool forwardReactionIsPossible = true;
    bool backwardReactionIsPossible = true;
  
    for(unsigned i =0; i < reaction.species.size() && aReactionIsPossible; ++i)
    {
       const std::string& speciesStr = reaction.species[i];
       if(reaction.stoich_coeffs[i] < 0){
         if(forwardReactionIsPossible) forwardReactionIsPossible = isMember(speciesStr,species);
       } else {
         if(backwardReactionIsPossible) backwardReactionIsPossible = isMember(speciesStr,species);
       }
 
       aReactionIsPossible = forwardReactionIsPossible || backwardReactionIsPossible;
       
       // block 
       bool isReducingEnvironment = true;
       if(isReducingEnvironment && aReactionIsPossible){
         aReactionIsPossible = !( streq(reaction.species[i],"e-") || streq(reaction.species[i],"O2"));
       }
    }

    if(aReactionIsPossible){
      species.insert(reaction.species.begin(), reaction.species.end());
      reactionList.push_back(reaction);
    }
  }

  return aReactionIsPossible;
}

/// Creates a list of possible reactions from a set of species.
/// returns the full set of possible reaction species and the associated reactions 
void ChemistryManager::BuildReactionList(std::set<std::string>& reactionSpecies,
                                         std::vector<GPChemistry::ChemicalReaction>& activeSolutionReactions,
                                         std::vector<GPChemistry::ChemicalReaction>& activeSolidPhaseReactions)
{

  // 
  reactionSpecies.insert("H2O");
  reactionSpecies.insert("H+");

  activeSolutionReactions.resize(0);
  activeSolidPhaseReactions.resize(0);
  
  std::vector<GPChemistry::ChemicalReaction>& solutionReactions = m_solutionReactions;
  std::vector<GPChemistry::ChemicalReaction>& solidPhaseReactions = m_phases;

  unsigned numSolidPhaseReactions = solidPhaseReactions.size();
  unsigned numSolutionReactions = solutionReactions.size();
  std::vector<bool> solidPhaseReactionHasBeenRecorded(numSolidPhaseReactions,false );
  std::vector<bool> solutionReactionHasBeenRecorded(numSolutionReactions,false );

  // Check for valid reactions from the set of species until no new reactions are found
  bool noReactionsFound = false;
  while(!noReactionsFound)
  {
    noReactionsFound = true;
  
    // loop over solid phase reactions
    for(unsigned i =0 ; i < numSolidPhaseReactions; ++i)
    {
      if(!solidPhaseReactionHasBeenRecorded[i])
      {
        solidPhaseReactionHasBeenRecorded[i] = CheckIfReactionIsPossible(solidPhaseReactions[i],reactionSpecies,activeSolidPhaseReactions);
        if(solidPhaseReactionHasBeenRecorded[i]) noReactionsFound = false;
      }
    }

    // loop over aqueous reactions
    for(unsigned i =0 ; i < numSolutionReactions; ++i)
    {
      if(!solutionReactionHasBeenRecorded[i])
      {
        solutionReactionHasBeenRecorded[i] = CheckIfReactionIsPossible(solutionReactions[i],reactionSpecies,activeSolutionReactions);
        if(solutionReactionHasBeenRecorded[i]) noReactionsFound = false;
      }
    }
  }

}

void ChemistryManager::ResetLookupTables(){
#ifndef __xlC__

  // Element atomic weights
  m_atomicWeights = 
    CreateStlMap<std::string, realT >
      ("Ag",107.86820)
      ("Al",26.98154)
      ("Am",243.00000)
      ("Ar",39.94800)
      ("As",74.92160)
      ("Au",196.96655)
      ("B",10.81100)
      ("Ba",137.32700)
      ("Be",9.01218)
      ("Bi",208.98038)
      ("Br",79.90400)
      ("C",12.01070)
      ("Ca",40.07800)
      ("Cd",112.41100)
      ("Ce",140.11600)
      ("Cl",35.45270)
      ("Co",58.93320)
      ("Cr",51.99610)
      ("Cs",132.90545)
      ("Cu",63.54600)
      ("Dy",162.50000)
      ("Er",167.26000)
      ("Eu",151.96400)
      ("F",18.99840)
      ("Fe",55.84500)
      ("Fr",223.00000)
      ("Ga",69.72300)
      ("Gd",157.25000)
      ("H",1.00794)
      ("He",4.00260)
      ("Hf",178.49000)
      ("Hg",200.59000)
      ("Ho",164.93032)
      ("I",126.90447)
      ("In",114.81800)
      ("K",39.09830)
      ("Kr",83.80000)
      ("La",138.90550)
      ("Li",6.94100)
      ("Lu",174.96700)
      ("Mg",24.30500)
      ("Mn",54.93805)
      ("Mo",95.94000)
      ("N",14.00674)
      ("Na",22.98977)
      ("Nb",92.90638)
      ("Nd",144.24000)
      ("Ne",20.17970)
      ("Ni",58.69340)
      ("Np",237.00000)
      ("O",15.99940)
      ("P",30.97376)
      ("Pb",207.20000)
      ("Pd",106.42000)
      ("Pm",145.00000)
      ("Pr",140.90765)
      ("Pt",195.07800)
      ("Pu",244.00000)
      ("Ra",226.00000)
      ("Rb",85.46780)
      ("Re",186.20700)
      ("Rh",102.90550)
      ("Rn",222.00000)
      ("Ru",101.07000)
      ("S",32.06600)
      ("Sb",121.76000)
      ("Sc",44.95591)
      ("Se",78.96000)
      ("Si",28.08550)
      ("Sm",150.36000)
      ("Sn",118.71000)
      ("Sr",87.62000)
      ("Tb",158.92534)
      ("Tc",98.00000)
      ("Th",232.03810)
      ("Ti",47.86700)
      ("Tl",204.38330)
      ("Tm",168.93421)
      ("U",238.02890)
      ("V",50.94150)
      ("W",183.84000)
      ("Xe",131.29000)
      ("Y",88.90585)
      ("Yb",173.04000)
      ("Zn",65.39000)
      ("Zr",91.22400);

  // temperature (C) lookup table for debye_huckel_A parameter
  {
    rArray1d x(8), y(8);
    localIndex i = 0;
    x[i] = 0.01; y[i++] = 0.4939;
    x[i] = 25.;  y[i++] = 0.5144;
    x[i] = 60.;  y[i++] = 0.5465;
    x[i] = 100.; y[i++] = 0.5995;
    x[i] = 150.; y[i++] = 0.6855;
    x[i] = 200.; y[i++] = 0.7994;
    x[i] = 250.; y[i++] = 0.9593;
    x[i] = 300.; y[i++] = 1.218;
    m_DebyeHuckelA.SetGrid(&x);
    m_DebyeHuckelA.SetValues(y);
  }
  /// temperature (C) lookup table for debye_huckel_B parameter
  {
    rArray1d x(8), y(8);
    localIndex i = 0;
    x[i] = 0.01; y[i++] = 0.3253;
    x[i] = 25.;  y[i++] = 0.3288;
    x[i] = 60.;  y[i++] = 0.3346;
    x[i] = 100.; y[i++] = 0.3421;
    x[i] = 150.; y[i++] = 0.3525;
    x[i] = 200.; y[i++] = 0.3639;
    x[i] = 250.; y[i++] = 0.3766;
    x[i] = 300.; y[i++] = 0.3925;
    m_DebyeHuckelB.SetGrid(&x);
    m_DebyeHuckelB.SetValues(y);
  }
                                     
  /// temperature (C) lookup table for WATEQ extended Debye Huckel parameter (B dot)
  {
    rArray1d x(8), y(8);
    localIndex i = 0;
    x[i] = 0.01; y[i++] = 0.0394;
    x[i] = 25.;  y[i++] = 0.041;
    x[i] = 60.;  y[i++] = 0.0438;
    x[i] = 100.; y[i++] = 0.046;
    x[i] = 150.; y[i++] = 0.047;
    x[i] = 200.; y[i++] = 0.047;
    x[i] = 250.; y[i++] = 0.034;
    x[i] = 300.; y[i++] = 0.;
    m_WATEQBdot.SetGrid(&x);
    m_WATEQBdot.SetValues(y);
  }

  // ion size parameters for use in extended-debye-huckel equation
  m_ionSizeParameters =               
    CreateStlMap<std::string, realT >              
      ( "(UO2)3(CO3)6-6", 4.0000 )
      ( "Am(CO3)5-6", 4.0000 )              
      ( "Np(CO3)5-6", 4.0000 )              
      ( "Th(CO3)5-6", 4.0000 )              
      ( "U(CO3)5-6", 4.0000 )              
      ( "NpO2(CO3)3-5", 4.0000 )              
      ( "PuO2(CO3)3-5", 4.0000 )              
      ( "UO2(CO3)3-5", 4.0000 )              
      ( "Np(CO3)4-4", 4.0000 )              
      ( "P2O7-4", 4.0000 )              
      ( "Pd(SO4)3-4", 4.0000 )              
      ( "Pt(SO4)3-4", 4.0000 )              
      ( "PuO2(CO3)3-4", 4.0000 )              
      ( "Rh(SO4)3-4", 4.0000 )              
      ( "Ru(SO4)3-4", 4.0000 )              
      ( "Th(SO4)4-4", 4.0000 )              
      ( "U(CO3)4-4", 4.0000 )              
      ( "UO2(CO3)3-4", 4.0000 )              
      ( "Ag(CO3)2-3", 4.0000 )              
      ( "AgCl4-3", 4.0000 )              
      ( "Am(CO3)3-3", 4.0000 )              
      ( "AsO4-3", 4.0000 )              
      ( "Ce(PO4)2-3", 4.0000 )              
      ( "CrO4-3", 4.0000 )              
      ( "Dy(PO4)2-3", 4.0000 )              
      ( "Er(PO4)2-3", 4.0000 )              
      ( "Eu(PO4)2-3", 4.0000 )              
      ( "Gd(PO4)2-3", 4.0000 )              
      ( "HP2O7-3", 4.0000 )              
      ( "Ho(PO4)2-3", 4.0000 )              
      ( "KP2O7-3", 4.0000 )              
      ( "La(PO4)2-3", 4.0000 )              
      ( "Lu(PO4)2-3", 4.0000 )              
      ( "NaP2O7-3", 4.0000 )              
      ( "Nd(PO4)2-3", 4.0000 )              
      ( "NpO2(CO3)2-3", 4.0000 )              
      ( "PO4-3", 4.0000 )              
      ( "Pm(PO4)2-3", 4.0000 )              
      ( "Pr(PO4)2-3", 4.0000 )              
      ( "Rh(SO4)3-3", 4.0000 )              
      ( "Ru(SO4)3-3", 4.0000 )              
      ( "RuCl6-3", 4.0000 )              
      ( "Sm(PO4)2-3", 4.0000 )              
      ( "Tb(PO4)2-3", 4.0000 )              
      ( "Th(OH)4PO4-3", 4.0000 )              
      ( "Tm(PO4)2-3", 4.0000 )              
      ( "VO2(HPO4)2-3", 4.0000 )              
      ( "VO4-3", 4.0000 )              
      ( "Y(PO4)2-3", 4.0000 )              
      ( "Yb(PO4)2-3", 4.0000 )              
      ( "(UO2)11(CO3)6(OH)12-2", 4.0000 )              
      ( "AgCl3-2", 4.0000 )              
      ( "AsO3F-2", 4.0000 )              
      ( "AuCl3-2", 4.0000 )              
      ( "BeF4-2", 4.0000 )              
      ( "BeO2-2", 4.0000 )              
      ( "CO3-2", 4.5000 )              
      ( "CaP2O7-2", 4.0000 )              
      ( "Cd(CN)4-2", 4.0000 )              
      ( "Cd(CO3)2-2", 4.0000 )              
      ( "Cd(N3)4-2", 4.0000 )              
      ( "CdCl4-2", 4.0000 )              
      ( "CdI4-2", 4.0000 )              
      ( "CdO2-2", 4.0000 )              
      ( "CdP2O7-2", 4.0000 )              
      ( "CoO2-2", 4.0000 )              
      ( "Cr2O7-2", 4.0000 )              
      ( "CrO4-2", 4.0000 )              
      ( "CuCl3-2", 4.0000 )              
      ( "CuCl4-2", 4.0000 )              
      ( "CuO2-2", 4.0000 )              
      ( "Fe(CO3)2-2", 4.0000 )              
      ( "Fe(OH)4-2", 4.0000 )              
      ( "H2P2O7-2", 4.0000 )              
      ( "HAsO4-2", 4.0000 )              
      ( "HPO3-2", 4.0000 )              
      ( "HPO4-2", 4.0000 )              
      ( "HVO4-2", 4.0000 )              
      ( "HgCl4-2", 4.0000 )              
      ( "MgP2O7-2", 4.0000 )              
      ( "MnO2-2", 4.0000 )              
      ( "MnO4-2", 4.0000 )              
      ( "MoO4-2", 4.5000 )              
      ( "N2O2-2", 4.0000 )              
      ( "Na2P2O7-2", 4.0000 )              
      ( "NaHP2O7-2", 4.0000 )              
      ( "NiO2-2", 4.0000 )              
      ( "NiP2O7-2", 4.0000 )              
      ( "PO3F-2", 4.0000 )              
      ( "PbCl4-2", 4.0000 )              
      ( "PbI4-2", 4.0000 )              
      ( "PbP2O7-2", 4.0000 )              
      ( "Pd(SO4)2-2", 4.0000 )              
      ( "PdCl4-2", 4.0000 )              
      ( "Pt(SO4)2-2", 4.0000 )              
      ( "PtCl4-2", 4.0000 )              
      ( "PuO2(CO3)2-2", 4.0000 )              
      ( "Rh(SO4)2-2", 4.0000 )              
      ( "RhCl4-2", 4.0000 )              
      ( "Ru(OH)2Cl4-2", 4.0000 )              
      ( "Ru(SO4)2-2", 4.0000 )              
      ( "RuCl4-2", 4.0000 )              
      ( "RuCl5-2", 4.0000 )              
      ( "RuO4-2", 4.0000 )              
      ( "S-2", 5.0000 )              
      ( "S2-2", 4.0000 )              
      ( "S2O3-2", 4.0000 )              
      ( "S2O4-2", 5.0000 )              
      ( "S2O5-2", 4.0000 )              
      ( "S2O6-2", 4.0000 )              
      ( "S2O8-2", 4.0000 )              
      ( "S3-2", 4.0000 )              
      ( "S3O6-2", 4.0000 )              
      ( "S4-2", 4.0000 )              
      ( "S4O6-2", 4.0000 )              
      ( "S5-2", 4.0000 )              
      ( "S5O6-2", 4.0000 )              
      ( "SO3-2", 4.5000 )              
      ( "SO4-2", 4.0000 )              
      ( "SeO3-2", 4.0000 )              
      ( "SeO4-2", 4.0000 )              
      ( "SiF6-2", 4.0000 )              
      ( "SrP2O7-2", 4.0000 )              
      ( "TcO4-2", 4.0000 )              
      ( "Th(SO4)3-2", 4.0000 )              
      ( "UF6-2", 4.0000 )              
      ( "UO2(CO3)2-2", 4.0000 )              
      ( "UO2(N3)4-2", 4.0000 )              
      ( "UO2(SO4)2-2", 4.0000 )              
      ( "UO2F4-2", 4.0000 )              
      ( "UO4-2", 4.0000 )              
      ( "WO4-2", 5.0000 )              
      ( "Zn(CN)4-2", 4.0000 )              
      ( "Zn(SCN)4-2", 4.0000 )              
      ( "ZnI4-2", 4.0000 )              
      ( "ZnO2-2", 4.0000 )              
      ( "(UO2)2CO3(OH)3-", 4.0000 )              
      ( "(UO2)3(OH)7-", 4.0000 )              
      ( "Ag(HS)2-", 4.0000 )              
      ( "AgCO3-", 4.0000 )              
      ( "AgCl2-", 4.0000 )              
      ( "AgO-", 4.0000 )              
      ( "Al(SO4)2-", 4.0000 )              
      ( "AlF4-", 4.0000 )              
      ( "AlO2-", 4.0000 )              
      ( "Am(CO3)2-", 4.0000 )              
      ( "Am(SO4)2-", 4.0000 )              
      ( "Au(HS)2-", 4.0000 )              
      ( "AuCl2-", 4.0000 )              
      ( "AuCl4-", 4.0000 )              
      ( "B2O(OH)5-", 4.0000 )              
      ( "BF2(OH)2-", 4.0000 )              
      ( "BF3OH-", 4.0000 )              
      ( "BF4-", 4.0000 )              
      ( "BH4-", 4.0000 )              
      ( "BO2-", 4.0000 )              
      ( "BeF3-", 4.0000 )              
      ( "BiO2-", 4.0000 )              
      ( "Br-", 3.0000 )              
      ( "Br3-", 4.0000 )              
      ( "BrO-", 4.0000 )              
      ( "BrO3-", 3.5000 )              
      ( "BrO4-", 4.0000 )              
      ( "CN-", 3.0000 )              
      ( "CaPO4-", 4.0000 )              
      ( "Cd(CN)3-", 4.0000 )              
      ( "Cd(N3)3-", 4.0000 )              
      ( "Cd(SCN)3-", 4.0000 )              
      ( "CdBr3-", 4.0000 )              
      ( "CdCl3-", 4.0000 )              
      ( "CdI3-", 4.0000 )              
      ( "Ce(CO3)2-", 4.0000 )              
      ( "Ce(HPO4)2-", 4.0000 )              
      ( "Cl-", 3.0000 )              
      ( "ClO-", 4.0000 )              
      ( "ClO2-", 4.0000 )              
      ( "ClO3-", 3.5000 )              
      ( "ClO4-", 3.5000 )              
      ( "CrO2-", 4.0000 )              
      ( "CrO3Cl-", 4.0000 )              
      ( "CuCl2-", 4.0000 )              
      ( "CuCl3-", 4.0000 )              
      ( "Dy(CO3)2-", 4.0000 )              
      ( "Dy(HPO4)2-", 4.0000 )              
      ( "Dy(OH)4-", 4.0000 )              
      ( "Dy(SO4)2-", 4.0000 )              
      ( "Er(CO3)2-", 4.0000 )              
      ( "Er(HPO4)2-", 4.0000 )              
      ( "Er(OH)4-", 4.0000 )              
      ( "Er(SO4)2-", 4.0000 )              
      ( "Eu(CO3)2-", 4.0000 )              
      ( "Eu(HPO4)2-", 4.0000 )              
      ( "Eu(SO4)2-", 4.0000 )              
      ( "F-", 3.5000 )              
      ( "Fe(OH)3-", 4.0000 )              
      ( "Fe(SO4)2-", 4.0000 )              
      ( "FeO2-", 4.0000 )              
      ( "Formate", 3.5000 )              
      ( "GaO2-", 4.0000 )              
      ( "Gd(CO3)2-", 4.0000 )              
      ( "Gd(HPO4)2-", 4.0000 )              
      ( "Gd(OH)4-", 4.0000 )              
      ( "Gd(SO4)2-", 4.0000 )              
      ( "H2AsO3-", 4.0000 )              
      ( "H2AsO4-", 4.0000 )              
      ( "H2PO2-", 4.0000 )              
      ( "H2PO3-", 4.0000 )              
      ( "H2PO4-", 4.0000 )              
      ( "H2VO4-", 4.0000 )              
      ( "H3P2O7-", 4.0000 )              
      ( "HAsO3F-", 4.0000 )              
      ( "HBeO2-", 4.0000 )              
      ( "HCO3-", 4.0000 )              
      ( "HCdO2-", 4.0000 )              
      ( "HCoO2-", 4.0000 )              
      ( "HCrO4-", 4.0000 )              
      ( "HCuO2-", 4.0000 )              
      ( "HF2-", 4.0000 )              
      ( "HFeO2-", 4.0000 )              
      ( "HHfO3-", 4.0000 )              
      ( "HHgO2-", 4.0000 )              
      ( "HMnO2-", 4.0000 )              
      ( "HMoO4-", 4.0000 )              
      ( "HN2O2-", 4.0000 )              
      ( "HNiO2-", 4.0000 )              
      ( "HO2-", 4.0000 )              
      ( "HPO3F-", 4.0000 )              
      ( "HPbO2-", 4.0000 )              
      ( "HRuO5-", 4.0000 )              
      ( "HS-", 3.5000 )              
      ( "HS2O3-", 4.0000 )              
      ( "HS2O4-", 4.0000 )              
      ( "HSO3-", 4.0000 )              
      ( "HSO4-", 4.0000 )              
      ( "HSO5-", 4.0000 )              
      ( "HSe-", 4.0000 )              
      ( "HSeO3-", 4.0000 )              
      ( "HSeO4-", 4.0000 )              
      ( "HSiO3-", 4.0000 )              
      ( "HSnO2-", 4.0000 )              
      ( "HUO3-", 4.0000 )              
      ( "HUO4-", 4.0000 )              
      ( "HWO4-", 3.5000 )              
      ( "HZnO2-", 4.0000 )              
      ( "HZrO3-", 4.0000 )              
      ( "HgCl3-", 4.0000 )              
      ( "Ho(CO3)2-", 4.0000 )              
      ( "Ho(HPO4)2-", 4.0000 )              
      ( "Ho(SO4)2-", 4.0000 )              
      ( "I-", 3.0000 )              
      ( "I3-", 4.0000 )              
      ( "IO-", 4.0000 )              
      ( "IO3-", 4.0000 )              
      ( "IO4-", 3.5000 )              
      ( "InO2-", 4.0000 )              
      ( "KSO4-", 4.0000 )              
      ( "La(CO3)2-", 4.0000 )              
      ( "La(HPO4)2-", 4.0000 )              
      ( "La(SO4)2-", 4.0000 )              
      ( "LiSO4-", 4.0000 )              
      ( "Lu(CO3)2-", 4.0000 )              
      ( "Lu(HPO4)2-", 4.0000 )              
      ( "Lu(SO4)2-", 4.0000 )              
      ( "MnCl3-", 4.0000 )              
      ( "MnO4-", 3.5000 )              
      ( "N3-", 4.0000 )              
      ( "NH4SO4-", 4.0000 )              
      ( "NO2-", 3.0000 )              
      ( "NO3-", 3.0000 )              
      ( "NaCO3-", 4.0000 )              
      ( "NaSO4-", 4.0000 )              
      ( "NbO3-", 4.0000 )              
      ( "Nd(CO3)2-", 4.0000 )              
      ( "Nd(HPO4)2-", 4.0000 )              
      ( "Nd(OH)4-", 4.0000 )              
      ( "Nd(SO4)2-", 4.0000 )              
      ( "Ni(OH)3-", 4.0000 )              
      ( "NiHP2O7-", 4.0000 )              
      ( "NpO2(OH)2-", 4.0000 )              
      ( "NpO2CO3-", 4.0000 )              
      ( "OH-", 3.5000 )              
      ( "Pb(HS)3-", 4.0000 )              
      ( "Pb(OH)3-", 4.0000 )              
      ( "PbBr3-", 4.0000 )              
      ( "PbCl3-", 4.0000 )              
      ( "PbI3-", 4.0000 )              
      ( "PdCl3-", 4.0000 )              
      ( "Pm(CO3)2-", 4.0000 )              
      ( "Pm(HPO4)2-", 4.0000 )              
      ( "Pm(SO4)2-", 4.0000 )              
      ( "Pr(CO3)2-", 4.0000 )              
      ( "Pr(HPO4)2-", 4.0000 )              
      ( "Pr(SO4)2-", 4.0000 )              
      ( "PtCl3-", 4.0000 )              
      ( "PuO2CO3-", 4.0000 )              
      ( "ReO4-", 4.0000 )              
      ( "Rh(SO4)2-", 4.0000 )              
      ( "RhCl3-", 4.0000 )              
      ( "RhCl4-", 4.0000 )              
      ( "Ru(OH)2Cl3-", 4.0000 )              
      ( "Ru(SO4)2-", 4.0000 )              
      ( "RuCl3-", 4.0000 )              
      ( "RuCl4-", 4.0000 )              
      ( "RuO4-", 4.0000 )              
      ( "SCN-", 3.5000 )              
      ( "SbO2-", 4.0000 )              
      ( "ScO2-", 4.0000 )              
      ( "Sm(CO3)2-", 4.0000 )              
      ( "Sm(HPO4)2-", 4.0000 )              
      ( "Sm(OH)4-", 4.0000 )              
      ( "Sm(SO4)2-", 4.0000 )              
      ( "SnCl3-", 4.0000 )              
      ( "SnF3-", 4.0000 )              
      ( "Tb(CO3)2-", 4.0000 )              
      ( "Tb(HPO4)2-", 4.0000 )              
      ( "Tb(SO4)2-", 4.0000 )              
      ( "TcCO3(OH)3-", 4.0000 )              
      ( "TcO(OH)3-", 4.0000 )              
      ( "TcO4-", 4.0000 )              
      ( "Th(OH)3CO3-", 4.0000 )              
      ( "TlO2-", 4.0000 )              
      ( "Tm(CO3)2-", 4.0000 )              
      ( "Tm(HPO4)2-", 4.0000 )              
      ( "Tm(SO4)2-", 4.0000 )              
      ( "UF5-", 4.0000 )              
      ( "UO2(N3)3-", 4.0000 )              
      ( "UO2(SCN)3-", 4.0000 )              
      ( "UO2-", 4.0000 )              
      ( "UO2F3-", 4.0000 )              
      ( "UO2PO4-", 4.0000 )              
      ( "UO3-", 4.0000 )              
      ( "VO2F2-", 4.0000 )              
      ( "VO2HPO4-", 4.0000 )              
      ( "VO2SO4-", 4.0000 )              
      ( "Y(CO3)2-", 4.0000 )              
      ( "Y(HPO4)2-", 4.0000 )              
      ( "Y(OH)4-", 4.0000 )              
      ( "Y(SO4)2-", 4.0000 )              
      ( "Yb(CO3)2-", 4.0000 )              
      ( "Yb(HPO4)2-", 4.0000 )              
      ( "Yb(OH)4-", 4.0000 )              
      ( "Yb(SO4)2-", 4.0000 )              
      ( "ZnBr3-", 4.0000 )              
      ( "ZnCl3-", 4.0000 )              
      ( "ZnI3-", 4.0000 )              
      ( "(NH4)2Sb2S4", 3.0000 )              
      ( "AgCl", 3.0000 )              
      ( "AgF", 3.0000 )              
      ( "AgNO3", 3.0000 )              
      ( "AgOH", 3.0000 )              
      ( "AlF3", 3.0000 )              
      ( "Am(OH)3", 3.0000 )              
      ( "Ar", 3.0000 )              
      ( "AuCl", 3.0000 )              
      ( "B(OH)3", 3.0000 )              
      ( "BaCO3", 3.0000 )              
      ( "BeCl2", 3.0000 )              
      ( "BeF2", 3.0000 )              
      ( "BeO", 3.0000 )              
      ( "Br2", 3.0000 )              
      ( "CO2", 3.0000 )              
      ( "CaCO3", 3.0000 )              
      ( "CaCl2", 3.0000 )              
      ( "CaHPO4", 3.0000 )              
      ( "CaSO4", 3.0000 )              
      ( "Cd(CN)2", 3.0000 )              
      ( "Cd(N3)2", 3.0000 )              
      ( "Cd(OH)Cl", 3.0000 )              
      ( "Cd(SCN)2", 3.0000 )              
      ( "CdBr2", 3.0000 )              
      ( "CdCO3", 3.0000 )              
      ( "CdCl2", 3.0000 )              
      ( "CdF2", 3.0000 )              
      ( "CdI2", 3.0000 )              
      ( "CdO", 3.0000 )              
      ( "CdSO4", 3.0000 )              
      ( "CdSeO4", 3.0000 )              
      ( "Ce(OH)3", 3.0000 )              
      ( "CePO4", 3.0000 )              
      ( "CoBr2", 3.0000 )              
      ( "CoI2", 3.0000 )              
      ( "CoO", 3.0000 )              
      ( "CoSO4", 3.0000 )              
      ( "CoSeO4", 3.0000 )              
      ( "CsBr", 3.0000 )              
      ( "CsCl", 3.0000 )              
      ( "CsI", 3.0000 )              
      ( "CsOH", 3.0000 )              
      ( "CuCl", 3.0000 )              
      ( "CuCl2", 3.0000 )              
      ( "CuO", 3.0000 )              
      ( "CuSO4", 3.0000 )              
      ( "Dy(OH)3", 3.0000 )              
      ( "DyPO4", 3.0000 )              
      ( "Er(OH)3", 3.0000 )              
      ( "ErPO4", 3.0000 )              
      ( "Ethane", 3.0000 )              
      ( "Eu(OH)3", 3.0000 )              
      ( "EuF3", 3.0000 )              
      ( "EuPO4", 3.0000 )              
      ( "FeCO3", 3.0000 )              
      ( "FeCl2", 3.0000 )              
      ( "FeO", 3.0000 )              
      ( "FeSO4", 3.0000 )              
      ( "Formaldehyde", 3.0000 )              
      ( "Formic_acid", 3.0000 )              
      ( "Gd(OH)3", 3.0000 )              
      ( "GdPO4", 3.0000 )              
      ( "H2", 3.0000 )              
      ( "H2CrO4", 3.0000 )              
      ( "H2N2O2", 3.0000 )              
      ( "H2O", 3.0000 )              
      ( "H2O2", 3.0000 )              
      ( "H2PO3F", 3.0000 )              
      ( "H2S", 3.0000 )              
      ( "H2S2O3", 3.0000 )              
      ( "H2S2O4", 3.0000 )              
      ( "H2SO3", 3.0000 )              
      ( "H2SO4", 3.0000 )              
      ( "H2Se", 3.0000 )              
      ( "H2SeO3", 3.0000 )              
      ( "H3AsO4", 3.0000 )              
      ( "H3PO2", 3.0000 )              
      ( "H3PO3", 3.0000 )              
      ( "H3PO4", 3.0000 )              
      ( "H3VO4", 3.0000 )              
      ( "H4P2O7", 3.0000 )              
      ( "HAlO2", 3.0000 )              
      ( "HAsO2", 3.0000 )              
      ( "HAsS2", 3.0000 )              
      ( "HBiO2", 3.0000 )              
      ( "HBrO", 3.0000 )              
      ( "HClO", 3.0000 )              
      ( "HClO2", 3.0000 )              
      ( "HCrO2", 3.0000 )              
      ( "HF", 3.0000 )              
      ( "HFeO2", 3.0000 )              
      ( "HGaO2", 3.0000 )              
      ( "HIO", 3.0000 )              
      ( "HIO3", 3.0000 )              
      ( "HInO2", 3.0000 )              
      ( "HN3", 3.0000 )              
      ( "HNO2", 3.0000 )              
      ( "HNO3", 3.0000 )              
      ( "HNbO3", 3.0000 )              
      ( "HSbO2", 3.0000 )              
      ( "HScO2", 3.0000 )              
      ( "HTlO2", 3.0000 )              
      ( "HUO2", 3.0000 )              
      ( "He", 3.0000 )              
      ( "HfO2", 3.0000 )              
      ( "HgCl2", 3.0000 )              
      ( "HgO", 3.0000 )              
      ( "Ho(OH)3", 3.0000 )              
      ( "HoPO4", 3.0000 )              
      ( "KBr", 3.0000 )              
      ( "KCl", 3.0000 )              
      ( "KHSO4", 3.0000 )              
      ( "KI", 3.0000 )              
      ( "KOH", 3.0000 )              
      ( "Kr", 3.0000 )              
      ( "La(OH)3", 3.0000 )              
      ( "LaF3", 3.0000 )              
      ( "LaPO4", 3.0000 )              
      ( "LiCl", 3.0000 )              
      ( "LiOH", 3.0000 )              
      ( "Lu(OH)3", 3.0000 )              
      ( "LuPO4", 3.0000 )              
      ( "Methane", 3.0000 )              
      ( "Methanol", 3.0000 )              
      ( "MgCO3", 3.0000 )              
      ( "MgHPO4", 3.0000 )              
      ( "MgSO4", 3.0000 )              
      ( "Mn(NO3)2", 3.0000 )              
      ( "MnO", 3.0000 )              
      ( "MnSO4", 3.0000 )              
      ( "MnSeO4", 3.0000 )              
      ( "N2", 3.0000 )              
      ( "NH3", 3.0000 )              
      ( "NH4SbO2", 3.0000 )              
      ( "NaB(OH)4", 3.0000 )              
      ( "NaBr", 3.0000 )              
      ( "NaCl", 3.0000 )              
      ( "NaF", 3.0000 )              
      ( "NaHCO3", 3.0000 )              
      ( "NaHSiO3", 3.0000 )              
      ( "NaI", 3.0000 )              
      ( "NaOH", 3.0000 )              
      ( "Nd(OH)3", 3.0000 )              
      ( "NdPO4", 3.0000 )              
      ( "Ne", 3.0000 )              
      ( "Ni(NO3)2", 3.0000 )              
      ( "Ni(OH)2", 3.0000 )              
      ( "NiO", 3.0000 )              
      ( "NiSeO4", 3.0000 )              
      ( "Np(OH)4", 3.0000 )              
      ( "NpO2F", 3.0000 )              
      ( "NpO2OH", 3.0000 )              
      ( "O2", 3.0000 )              
      ( "O2(g)", 3.0000 )              
      ( "Pb(BrO3)2", 3.0000 )              
      ( "Pb(ClO3)2", 3.0000 )              
      ( "Pb(HS)2", 3.0000 )              
      ( "Pb(OH)2", 3.0000 )              
      ( "Pb(SCN)2", 3.0000 )              
      ( "PbBr2", 3.0000 )              
      ( "PbCl2", 3.0000 )              
      ( "PbF2", 3.0000 )              
      ( "PbHPO4", 3.0000 )              
      ( "PbI2", 3.0000 )              
      ( "PbO", 3.0000 )              
      ( "PdCl2", 3.0000 )              
      ( "PdO", 3.0000 )              
      ( "PdSO4", 3.0000 )              
      ( "Pm(OH)3", 3.0000 )              
      ( "PmPO4", 3.0000 )              
      ( "Pr(OH)3", 3.0000 )              
      ( "PrPO4", 3.0000 )              
      ( "PtCl2", 3.0000 )              
      ( "PtO", 3.0000 )              
      ( "PtSO4", 3.0000 )              
      ( "PuO2(OH)2", 3.0000 )              
      ( "PuO2CO3", 3.0000 )              
      ( "PuO2F2", 3.0000 )              
      ( "PuO2OH", 3.0000 )              
      ( "RbBr", 3.0000 )              
      ( "RbCl", 3.0000 )              
      ( "RbF", 3.0000 )              
      ( "RbI", 3.0000 )              
      ( "RbOH", 3.0000 )              
      ( "RhCl2", 3.0000 )              
      ( "RhCl3", 3.0000 )              
      ( "RhO", 3.0000 )              
      ( "RhSO4", 3.0000 )              
      ( "Rn", 3.0000 )              
      ( "Ru(Cl)3", 3.0000 )              
      ( "Ru(OH)2Cl2", 3.0000 )              
      ( "Ru(OH)2SO4", 3.0000 )              
      ( "Ru(OH)4", 3.0000 )              
      ( "RuCl2", 3.0000 )              
      ( "RuCl3", 3.0000 )              
      ( "RuO", 3.0000 )              
      ( "RuO4", 3.0000 )              
      ( "RuSO4", 3.0000 )              
      ( "SO2", 3.0000 )              
      ( "SiO2", 3.0000 )              
      ( "Sm(OH)3", 3.0000 )              
      ( "SmPO4", 3.0000 )              
      ( "Sn(OH)4", 3.0000 )              
      ( "Sn(SO4)2", 3.0000 )              
      ( "SnCl2", 3.0000 )              
      ( "SnF2", 3.0000 )              
      ( "SnO", 3.0000 )              
      ( "SrCO3", 3.0000 )              
      ( "SrHPO4", 3.0000 )              
      ( "SrSO4", 3.0000 )              
      ( "Tb(OH)3", 3.0000 )              
      ( "TbPO4", 3.0000 )              
      ( "TcCO3(OH)2", 3.0000 )              
      ( "TcO(OH)2", 3.0000 )              
      ( "Th(OH)4", 3.0000 )              
      ( "Th(SO4)2", 3.0000 )              
      ( "ThCl4", 3.0000 )              
      ( "ThF4", 3.0000 )              
      ( "Ti(OH)4", 3.0000 )              
      ( "TlCl", 3.0000 )              
      ( "TlF", 3.0000 )              
      ( "TlOH", 3.0000 )              
      ( "Tm(OH)3", 3.0000 )              
      ( "TmPO4", 3.0000 )              
      ( "U(SO4)2", 3.0000 )              
      ( "UF4", 3.0000 )              
      ( "UO2(H2PO4)2", 3.0000 )              
      ( "UO2(IO3)2", 3.0000 )              
      ( "UO2(N3)2", 3.0000 )              
      ( "UO2(SCN)2", 3.0000 )              
      ( "UO2", 3.0000 )              
      ( "UO2CO3", 3.0000 )              
      ( "UO2Cl2", 3.0000 )              
      ( "UO2F2", 3.0000 )              
      ( "UO2HPO4", 3.0000 )              
      ( "UO2OH", 3.0000 )              
      ( "UO2S2O3", 3.0000 )              
      ( "UO2SO3", 3.0000 )              
      ( "UO2SO4", 3.0000 )              
      ( "UO3", 3.0000 )              
      ( "VO(OH)3", 3.0000 )              
      ( "VO2F", 3.0000 )              
      ( "VO2H2PO4", 3.0000 )              
      ( "VOF2", 3.0000 )              
      ( "VOSO4", 3.0000 )              
      ( "Xe", 3.0000 )              
      ( "Y(OH)3", 3.0000 )              
      ( "YF3", 3.0000 )              
      ( "YPO4", 3.0000 )              
      ( "Yb(OH)3", 3.0000 )              
      ( "YbPO4", 3.0000 )              
      ( "Zn(N3)2", 3.0000 )              
      ( "Zn(OH)Cl", 3.0000 )              
      ( "Zn(SCN)2", 3.0000 )              
      ( "ZnBr2", 3.0000 )              
      ( "ZnCO3", 3.0000 )              
      ( "ZnCl2", 3.0000 )              
      ( "ZnHPO4", 3.0000 )              
      ( "ZnI2", 3.0000 )              
      ( "ZnO", 3.0000 )              
      ( "ZnSO4", 3.0000 )              
      ( "ZnSeO4", 3.0000 )              
      ( "ZrO2", 3.0000 )              
      ( "Acetic_acid", 3.0000 )              
      ( "(UO2)3(OH)5+", 4.0000 )              
      ( "(UO2)3O(OH)2(HCO3)+", 4.0000 )              
      ( "(UO2)4(OH)7+", 4.0000 )              
      ( "Ag+", 2.5000 )              
      ( "AlF2+", 4.0000 )              
      ( "AlHPO4+", 4.0000 )              
      ( "AlO+", 4.0000 )              
      ( "AlSO4+", 4.0000 )              
      ( "Am(OH)2+", 4.0000 )              
      ( "AmCO3+", 4.0000 )              
      ( "AmF2+", 4.0000 )              
      ( "AmO2+", 4.0000 )              
      ( "AmSO4+", 4.0000 )              
      ( "Au+", 4.0000 )              
      ( "BaB(OH)4+", 4.0000 )              
      ( "BaCl+", 4.0000 )              
      ( "BaF+", 4.0000 )              
      ( "BaHCO3+", 4.0000 )              
      ( "BaNO3+", 4.0000 )              
      ( "BaOH+", 4.0000 )              
      ( "BeCl+", 4.0000 )              
      ( "BeF+", 4.0000 )              
      ( "BeOH+", 4.0000 )              
      ( "BiO+", 4.0000 )              
      ( "CaB(OH)4+", 4.0000 )              
      ( "CaCl+", 4.0000 )              
      ( "CaF+", 4.0000 )              
      ( "CaHCO3+", 4.0000 )              
      ( "CaHSiO3+", 4.0000 )              
      ( "CaNO3+", 4.0000 )              
      ( "CaOH+", 4.0000 )              
      ( "CdBr+", 4.0000 )              
      ( "CdCN+", 4.0000 )              
      ( "CdCl+", 4.0000 )              
      ( "CdF+", 4.0000 )              
      ( "CdHCO3+", 4.0000 )              
      ( "CdI+", 4.0000 )              
      ( "CdN3+", 4.0000 )              
      ( "CdNO2+", 4.0000 )              
      ( "CdOH+", 4.0000 )              
      ( "CdSCN+", 4.0000 )              
      ( "Ce(OH)2+", 4.0000 )              
      ( "CeCO3+", 4.0000 )              
      ( "CeF2+", 4.0000 )              
      ( "CeHPO4+", 4.0000 )              
      ( "CeSO4+", 4.0000 )              
      ( "Co2(OH)3+", 4.0000 )              
      ( "CoCl+", 4.0000 )              
      ( "CoF+", 4.0000 )              
      ( "CoNO3+", 4.0000 )              
      ( "CoOH+", 4.0000 )              
      ( "CrCl2+", 4.0000 )              
      ( "CrO+", 4.0000 )              
      ( "Cs+", 2.5000 )              
      ( "Cu+", 4.0000 )              
      ( "CuCl+", 4.0000 )              
      ( "CuF+", 4.0000 )              
      ( "CuOH+", 4.0000 )              
      ( "Dy(OH)2+", 4.0000 )              
      ( "DyCO3+", 4.0000 )              
      ( "DyHPO4+", 4.0000 )              
      ( "DySO4+", 4.0000 )              
      ( "Er(OH)2+", 4.0000 )              
      ( "ErCO3+", 4.0000 )              
      ( "ErHPO4+", 4.0000 )              
      ( "ErSO4+", 4.0000 )              
      ( "Eu(OH)2+", 4.0000 )              
      ( "EuBr2+", 4.0000 )              
      ( "EuCO3+", 4.0000 )              
      ( "EuCl2+", 4.0000 )              
      ( "EuF2+", 4.0000 )              
      ( "EuHPO4+", 4.0000 )              
      ( "EuSO4+", 4.0000 )              
      ( "FeCl+", 4.0000 )              
      ( "FeF+", 4.0000 )              
      ( "FeF2+", 4.0000 )              
      ( "FeH2PO4+", 4.0000 )              
      ( "FeHCO3+", 4.0000 )              
      ( "FeO+", 4.0000 )              
      ( "FeOH+", 4.0000 )              
      ( "FeSO4+", 4.0000 )              
      ( "Fr+", 2.5000 )              
      ( "GaO+", 4.0000 )              
      ( "Gd(OH)2+", 4.0000 )              
      ( "GdCO3+", 4.0000 )              
      ( "GdF2+", 4.0000 )              
      ( "GdHPO4+", 4.0000 )              
      ( "GdSO4+", 4.0000 )              
      ( "H+", 9.0000 )              
      ( "HHfO2+", 4.0000 )              
      ( "HUO2+", 4.0000 )              
      ( "HZrO2+", 4.0000 )              
      ( "HgCl+", 4.0000 )              
      ( "HgF+", 4.0000 )              
      ( "HgOH+", 4.0000 )              
      ( "Ho(OH)2+", 4.0000 )              
      ( "HoCO3+", 4.0000 )              
      ( "HoHPO4+", 4.0000 )              
      ( "HoSO4+", 4.0000 )              
      ( "InO+", 4.0000 )              
      ( "K+", 3.0000 )              
      ( "La(OH)2+", 4.0000 )              
      ( "LaCO3+", 4.0000 )              
      ( "LaF2+", 4.0000 )              
      ( "LaHPO4+", 4.0000 )              
      ( "LaSO4+", 4.0000 )              
      ( "Li+", 6.0000 )              
      ( "Lu(OH)2+", 4.0000 )              
      ( "LuCO3+", 4.0000 )              
      ( "LuHPO4+", 4.0000 )              
      ( "LuSO4+", 4.0000 )              
      ( "MgB(OH)4+", 4.0000 )              
      ( "MgCl+", 4.0000 )              
      ( "MgF+", 4.0000 )              
      ( "MgHCO3+", 4.0000 )              
      ( "MgHSiO3+", 4.0000 )              
      ( "MgOH+", 4.0000 )              
      ( "Mn2(OH)3+", 4.0000 )              
      ( "MnCl+", 4.0000 )              
      ( "MnF+", 4.0000 )              
      ( "MnHCO3+", 4.0000 )              
      ( "MnNO3+", 4.0000 )              
      ( "MnOH+", 4.0000 )              
      ( "N2H5+", 4.0000 )              
      ( "NH4+", 2.5000 )              
      ( "Na+", 4.0000 )              
      ( "Nd(OH)2+", 4.0000 )              
      ( "NdCO3+", 4.0000 )              
      ( "NdHPO4+", 4.0000 )              
      ( "NdSO4+", 4.0000 )              
      ( "NiBr+", 4.0000 )              
      ( "NiCl+", 4.0000 )              
      ( "NiF+", 4.0000 )              
      ( "NiNO3+", 4.0000 )              
      ( "NiOH+", 4.0000 )              
      ( "NpO2+", 4.0000 )              
      ( "PH4+", 4.0000 )              
      ( "PbBr+", 4.0000 )              
      ( "PbBrO3+", 4.0000 )              
      ( "PbCl+", 4.0000 )              
      ( "PbClO3+", 4.0000 )              
      ( "PbF+", 4.0000 )              
      ( "PbI+", 4.0000 )              
      ( "PbNO3+", 4.0000 )              
      ( "PbOH+", 4.0000 )              
      ( "PbSCN+", 4.0000 )              
      ( "PdCl+", 4.0000 )              
      ( "PdOH+", 4.0000 )              
      ( "Pm(OH)2+", 4.0000 )              
      ( "PmCO3+", 4.0000 )              
      ( "PmHPO4+", 4.0000 )              
      ( "PmSO4+", 4.0000 )              
      ( "Pr(OH)2+", 4.0000 )              
      ( "PrCO3+", 4.0000 )              
      ( "PrHPO4+", 4.0000 )              
      ( "PrSO4+", 4.0000 )              
      ( "PtCl+", 4.0000 )              
      ( "PtOH+", 4.0000 )              
      ( "PuO2+", 4.0000 )              
      ( "PuO2F+", 4.0000 )              
      ( "PuO2OH+", 4.0000 )              
      ( "Rb+", 2.5000 )              
      ( "RhCl+", 4.0000 )              
      ( "RhCl2+", 4.0000 )              
      ( "RhO+", 4.0000 )              
      ( "RhOH+", 4.0000 )              
      ( "RhSO4+", 4.0000 )              
      ( "Ru(Cl)2+", 4.0000 )              
      ( "Ru(OH)2+", 4.0000 )              
      ( "Ru(OH)2Cl+", 4.0000 )              
      ( "RuCl+", 4.0000 )              
      ( "RuCl2+", 4.0000 )              
      ( "RuO+", 4.0000 )              
      ( "RuOH+", 4.0000 )              
      ( "RuSO4+", 4.0000 )              
      ( "ScO+", 4.0000 )              
      ( "Sm(OH)2+", 4.0000 )              
      ( "SmCO3+", 4.0000 )              
      ( "SmHPO4+", 4.0000 )              
      ( "SmSO4+", 4.0000 )              
      ( "Sn(OH)3+", 4.0000 )              
      ( "SnCl+", 4.0000 )              
      ( "SnF+", 4.0000 )              
      ( "SnOH+", 4.0000 )              
      ( "SrCl+", 4.0000 )              
      ( "SrF+", 4.0000 )              
      ( "SrHCO3+", 4.0000 )              
      ( "SrNO3+", 4.0000 )              
      ( "SrOH+", 4.0000 )              
      ( "Tb(OH)2+", 4.0000 )              
      ( "TbCO3+", 4.0000 )              
      ( "TbHPO4+", 4.0000 )              
      ( "TbSO4+", 4.0000 )              
      ( "TcOOH+", 4.0000 )              
      ( "Th(OH)3+", 4.0000 )              
      ( "ThCl3+", 4.0000 )              
      ( "ThF3+", 4.0000 )              
      ( "Tl+", 2.5000 )              
      ( "TlO+", 4.0000 )              
      ( "Tm(OH)2+", 4.0000 )              
      ( "TmCO3+", 4.0000 )              
      ( "TmHPO4+", 4.0000 )              
      ( "TmSO4+", 4.0000 )              
      ( "UF3+", 4.0000 )              
      ( "UO+", 4.0000 )              
      ( "UO2(H2PO4)(H3PO4)+", 4.0000 )              
      ( "UO2+", 4.0000 )              
      ( "UO2Br+", 4.0000 )              
      ( "UO2BrO3+", 4.0000 )              
      ( "UO2Cl+", 4.0000 )              
      ( "UO2ClO3+", 4.0000 )              
      ( "UO2F+", 4.0000 )              
      ( "UO2H2PO4+", 4.0000 )              
      ( "UO2IO3+", 4.0000 )              
      ( "UO2N3+", 4.0000 )              
      ( "UO2NO3+", 4.0000 )              
      ( "UO2OH+", 4.0000 )              
      ( "UO2OSi(OH)3+", 4.0000 )              
      ( "UO2SCN+", 4.0000 )              
      ( "V(OH)2+", 4.0000 )              
      ( "VO+", 4.0000 )              
      ( "VO2+", 4.0000 )              
      ( "VOF+", 4.0000 )              
      ( "VOH+", 4.0000 )              
      ( "VOOH+", 4.0000 )              
      ( "VSO4+", 4.0000 )              
      ( "Y(OH)2+", 4.0000 )              
      ( "YCO3+", 4.0000 )              
      ( "YF2+", 4.0000 )              
      ( "YHPO4+", 4.0000 )              
      ( "YSO4+", 4.0000 )              
      ( "Yb(OH)2+", 4.0000 )              
      ( "YbCO3+", 4.0000 )              
      ( "YbHPO4+", 4.0000 )              
      ( "YbSO4+", 4.0000 )              
      ( "ZnBr+", 4.0000 )              
      ( "ZnCl+", 4.0000 )              
      ( "ZnClO4+", 4.0000 )              
      ( "ZnF+", 4.0000 )              
      ( "ZnH2PO4+", 4.0000 )              
      ( "ZnHCO3+", 4.0000 )              
      ( "ZnI+", 4.0000 )              
      ( "ZnN3+", 4.0000 )              
      ( "ZnOH+", 4.0000 )              
      ( "(PuO2)2(OH)2+2", 4.5000 )              
      ( "(UO2)2(OH)2+2", 4.5000 )              
      ( "(UO2)3(OH)4+2", 4.5000 )              
      ( "(VO)2(OH)2+2", 4.5000 )              
      ( "Ag+2", 4.5000 )              
      ( "AlF+2", 4.5000 )              
      ( "AlH2PO4+2", 4.5000 )              
      ( "AlOH+2", 4.5000 )              
      ( "Am+2", 4.5000 )              
      ( "AmCl+2", 4.5000 )              
      ( "AmF+2", 4.5000 )              
      ( "AmH2PO4+2", 4.5000 )              
      ( "AmN3+2", 4.5000 )              
      ( "AmNO2+2", 4.5000 )              
      ( "AmNO3+2", 4.5000 )              
      ( "AmO2+2", 4.5000 )              
      ( "AmOH+2", 4.5000 )              
      ( "Ba+2", 5.0000 )              
      ( "Be+2", 8.0000 )              
      ( "BiOH+2", 4.5000 )              
      ( "Ca+2", 6.0000 )              
      ( "Cd(NH3)+2", 4.5000 )              
      ( "Cd(NH3)2+2", 4.5000 )              
      ( "Cd(NH3)4+2", 4.5000 )              
      ( "Cd+2", 5.0000 )              
      ( "Ce(OH)2+2", 4.5000 )              
      ( "Ce+2", 4.5000 )              
      ( "CeCl+2", 4.5000 )              
      ( "CeF+2", 4.5000 )              
      ( "CeH2PO4+2", 4.5000 )              
      ( "CeHCO3+2", 4.5000 )              
      ( "CeNO3+2", 4.5000 )              
      ( "CeOH+2", 4.5000 )              
      ( "Co+2", 6.0000 )              
      ( "CoOH+2", 4.5000 )              
      ( "Cr+2", 6.0000 )              
      ( "CrBr+2", 4.5000 )              
      ( "CrCl+2", 4.5000 )              
      ( "CrOH+2", 4.5000 )              
      ( "Cu(NH3)2+2", 4.5000 )              
      ( "Cu(NH3)3+2", 4.5000 )              
      ( "Cu+2", 6.0000 )              
      ( "CuNH3+2", 4.5000 )              
      ( "Dy+2", 4.5000 )              
      ( "DyCl+2", 4.5000 )              
      ( "DyF+2", 4.5000 )              
      ( "DyH2PO4+2", 4.5000 )              
      ( "DyHCO3+2", 4.5000 )              
      ( "DyNO3+2", 4.5000 )              
      ( "DyOH+2", 4.5000 )              
      ( "Er+2", 4.5000 )              
      ( "ErCl+2", 4.5000 )              
      ( "ErF+2", 4.5000 )              
      ( "ErH2PO4+2", 4.5000 )              
      ( "ErHCO3+2", 4.5000 )              
      ( "ErNO3+2", 4.5000 )              
      ( "ErOH+2", 4.5000 )              
      ( "Eu+2", 4.5000 )              
      ( "EuBr+2", 4.5000 )              
      ( "EuBrO3+2", 4.5000 )              
      ( "EuCl+2", 4.5000 )              
      ( "EuF+2", 4.5000 )              
      ( "EuH2PO4+2", 4.5000 )              
      ( "EuHCO3+2", 4.5000 )              
      ( "EuIO3+2", 4.5000 )              
      ( "EuNO3+2", 4.5000 )              
      ( "EuOH+2", 4.5000 )              
      ( "Fe+2", 6.0000 )              
      ( "FeCl+2", 4.5000 )              
      ( "FeF+2", 4.5000 )              
      ( "FeH2PO4+2", 4.5000 )              
      ( "FeNO2+2", 4.5000 )              
      ( "FeNO3+2", 4.5000 )              
      ( "FeOH+2", 4.5000 )              
      ( "GaOH+2", 4.5000 )              
      ( "Gd+2", 4.5000 )              
      ( "GdCl+2", 4.5000 )              
      ( "GdF+2", 4.5000 )              
      ( "GdH2PO4+2", 4.5000 )              
      ( "GdHCO3+2", 4.5000 )              
      ( "GdNO3+2", 4.5000 )              
      ( "GdOH+2", 4.5000 )              
      ( "HfO+2", 4.5000 )              
      ( "Hg+2", 5.0000 )              
      ( "Hg2+2", 4.0000 )              
      ( "Ho+2", 4.5000 )              
      ( "HoF+2", 4.5000 )              
      ( "HoH2PO4+2", 4.5000 )              
      ( "HoHCO3+2", 4.5000 )              
      ( "HoNO3+2", 4.5000 )              
      ( "HoOH+2", 4.5000 )              
      ( "InCl+2", 4.5000 )              
      ( "InF+2", 4.5000 )              
      ( "InOH+2", 4.5000 )              
      ( "La+2", 4.5000 )              
      ( "LaCl+2", 4.5000 )              
      ( "LaF+2", 4.5000 )              
      ( "LaH2PO4+2", 4.5000 )              
      ( "LaHCO3+2", 4.5000 )              
      ( "LaNO3+2", 4.5000 )              
      ( "LaOH+2", 4.5000 )              
      ( "LuCl+2", 4.5000 )              
      ( "LuF+2", 4.5000 )              
      ( "LuH2PO4+2", 4.5000 )              
      ( "LuHCO3+2", 4.5000 )              
      ( "LuNO3+2", 4.5000 )              
      ( "LuOH+2", 4.5000 )              
      ( "Mg+2", 8.0000 )              
      ( "Mn+2", 6.0000 )              
      ( "N2H6+2", 4.5000 )              
      ( "Nd+2", 4.5000 )              
      ( "NdCl+2", 4.5000 )              
      ( "NdF+2", 4.5000 )              
      ( "NdH2PO4+2", 4.5000 )              
      ( "NdHCO3+2", 4.5000 )              
      ( "NdNO3+2", 4.5000 )              
      ( "NdOH+2", 4.5000 )              
      ( "Ni(NH3)2+2", 4.5000 )              
      ( "Ni(NH3)6+2", 4.5000 )              
      ( "Ni+2", 6.0000 )              
      ( "NpF2+2", 4.5000 )              
      ( "NpO2+2", 4.5000 )              
      ( "Pb+2", 4.5000 )              
      ( "Pb3(OH)4+2", 4.5000 )              
      ( "Pd+2", 4.5000 )              
      ( "Pm+2", 4.5000 )              
      ( "PmCl+2", 4.5000 )              
      ( "PmF+2", 4.5000 )              
      ( "PmH2PO4+2", 4.5000 )              
      ( "PmHCO3+2", 4.5000 )              
      ( "PmNO3+2", 4.5000 )              
      ( "PmOH+2", 4.5000 )              
      ( "Pr+2", 4.5000 )              
      ( "PrCl+2", 4.5000 )              
      ( "PrF+2", 4.5000 )              
      ( "PrH2PO4+2", 4.5000 )              
      ( "PrHCO3+2", 4.5000 )              
      ( "PrNO3+2", 4.5000 )              
      ( "PrOH+2", 4.5000 )              
      ( "Pt+2", 4.5000 )              
      ( "PuF2+2", 4.5000 )              
      ( "PuO2+2", 4.5000 )              
      ( "PuOH+2", 4.5000 )              
      ( "Ra+2", 5.0000 )              
      ( "Rh+2", 4.5000 )              
      ( "RhCl+2", 4.5000 )              
      ( "RhOH+2", 4.5000 )              
      ( "Ru(OH)2+2", 4.5000 )              
      ( "Ru+2", 4.5000 )              
      ( "RuCl+2", 4.5000 )              
      ( "RuOH+2", 4.5000 )              
      ( "ScOH+2", 4.5000 )              
      ( "Sm+2", 4.5000 )              
      ( "SmCl+2", 4.5000 )              
      ( "SmF+2", 4.5000 )              
      ( "SmH2PO4+2", 4.5000 )              
      ( "SmHCO3+2", 4.5000 )              
      ( "SmNO3+2", 4.5000 )              
      ( "SmOH+2", 4.5000 )              
      ( "Sn(OH)2+2", 4.5000 )              
      ( "Sn+2", 6.0000 )              
      ( "SnSO4+2", 4.5000 )              
      ( "Sr+2", 5.0000 )              
      ( "Tb+2", 4.5000 )              
      ( "TbCl+2", 4.5000 )              
      ( "TbF+2", 4.5000 )              
      ( "TbH2PO4+2", 4.5000 )              
      ( "TbHCO3+2", 4.5000 )              
      ( "TbNO3+2", 4.5000 )              
      ( "TbOH+2", 4.5000 )              
      ( "TcO+2", 4.5000 )              
      ( "ThCl2+2", 4.5000 )              
      ( "ThF2+2", 4.5000 )              
      ( "ThSO4+2", 4.5000 )              
      ( "TlCl+2", 4.5000 )              
      ( "TlOH+2", 4.5000 )              
      ( "Tm+2", 4.5000 )              
      ( "TmCl+2", 4.5000 )              
      ( "TmF+2", 4.5000 )              
      ( "TmH2PO4+2", 4.5000 )              
      ( "TmHCO3+2", 4.5000 )              
      ( "TmNO3+2", 4.5000 )              
      ( "TmOH+2", 4.5000 )              
      ( "U(NO3)2+2", 4.5000 )              
      ( "U(SCN)2+2", 4.5000 )              
      ( "UF2+2", 4.5000 )              
      ( "UO+2", 4.5000 )              
      ( "UO2+2", 4.5000 )              
      ( "UO2H3PO4+2", 4.5000 )              
      ( "UOH+2", 4.5000 )              
      ( "USO4+2", 4.5000 )              
      ( "V+2", 6.0000 )              
      ( "VO+2", 4.5000 )              
      ( "VOH+2", 4.5000 )              
      ( "YCl+2", 4.5000 )              
      ( "YF+2", 4.5000 )              
      ( "YH2PO4+2", 4.5000 )              
      ( "YHCO3+2", 4.5000 )              
      ( "YNO3+2", 4.5000 )              
      ( "YOH+2", 4.5000 )              
      ( "Yb+2", 4.5000 )              
      ( "YbCl+2", 4.5000 )              
      ( "YbF+2", 4.5000 )              
      ( "YbH2PO4+2", 4.5000 )              
      ( "YbHCO3+2", 4.5000 )              
      ( "YbNO3+2", 4.5000 )              
      ( "YbOH+2", 4.5000 )              
      ( "Zn(NH3)+2", 4.5000 )              
      ( "Zn(NH3)2+2", 4.5000 )              
      ( "Zn(NH3)3+2", 4.5000 )              
      ( "Zn(NH3)4+2", 4.5000 )              
      ( "Zn+2", 6.0000 )              
      ( "ZrO+2", 4.5000 )              
      ( "(UO2)2OH+3", 5.0000 )              
      ( "Al+3", 9.0000 )              
      ( "Am+3", 9.0000 )              
      ( "Au+3", 5.0000 )              
      ( "Bi+3", 9.0000 )              
      ( "Cd2OH+3", 5.0000 )              
      ( "Ce+3", 9.0000 )              
      ( "CeOH+3", 5.0000 )              
      ( "Co+3", 9.0000 )              
      ( "Cr+3", 9.0000 )              
      ( "Dy+3", 9.0000 )              
      ( "Er+3", 9.0000 )              
      ( "Eu+3", 9.0000 )              
      ( "Fe+3", 9.0000 )              
      ( "Ga+3", 9.0000 )              
      ( "Gd+3", 9.0000 )              
      ( "HfOH+3", 5.0000 )              
      ( "Ho+3", 9.0000 )              
      ( "In+3", 9.0000 )              
      ( "La+3", 9.0000 )              
      ( "Lu+3", 9.0000 )              
      ( "Mn+3", 9.0000 )              
      ( "Mn2OH+3", 5.0000 )              
      ( "Nd+3", 9.0000 )              
      ( "Ni2OH+3", 5.0000 )              
      ( "NpF+3", 5.0000 )              
      ( "NpOH+3", 5.0000 )              
      ( "Pb2OH+3", 5.0000 )              
      ( "Pm+3", 9.0000 )              
      ( "Pr+3", 9.0000 )              
      ( "Pu+3", 9.0000 )              
      ( "PuF+3", 5.0000 )              
      ( "PuOH+3", 5.0000 )              
      ( "Rh+3", 9.0000 )              
      ( "Ru+3", 9.0000 )              
      ( "Sc+3", 9.0000 )              
      ( "Sm+3", 9.0000 )              
      ( "SnOH+3", 5.0000 )              
      ( "Tb+3", 9.0000 )              
      ( "ThCl+3", 5.0000 )              
      ( "ThF+3", 5.0000 )              
      ( "ThOH+3", 5.0000 )              
      ( "Tl+3", 9.0000 )              
      ( "Tm+3", 9.0000 )              
      ( "U+3", 9.0000 )              
      ( "UBr+3", 5.0000 )              
      ( "UCl+3", 5.0000 )              
      ( "UF+3", 5.0000 )              
      ( "UI+3", 5.0000 )              
      ( "UNO3+3", 5.0000 )              
      ( "UOH+3", 5.0000 )              
      ( "USCN+3", 5.0000 )              
      ( "V+3", 9.0000 )              
      ( "Y+3", 9.0000 )              
      ( "Yb+3", 9.0000 )              
      ( "ZrOH+3", 9.0000 )              
      ( "Al2(OH)2+4", 5.5000 )              
      ( "Am+4", 11.0000 )              
      ( "Cd4(OH)4+4", 5.5000 )              
      ( "Ce+4", 11.0000 )              
      ( "Ce3(OH)5+4", 5.5000 )              
      ( "Co4(OH)4+4", 5.5000 )              
      ( "Cr2(OH)2+4", 5.5000 )              
      ( "Fe2(OH)2+4", 5.5000 )              
      ( "Hf+4", 11.0000 )              
      ( "La2(OH)2+4", 5.5000 )              
      ( "Mg4(OH)4+4", 5.5000 )              
      ( "Nd2(OH)2+4", 5.5000 )              
      ( "Ni4(OH)4+4", 5.5000 )              
      ( "Np+4", 11.0000 )              
      ( "Pb4(OH)4+4", 5.5000 )              
      ( "Pb6(OH)8+4", 5.5000 )              
      ( "Pu+4", 11.0000 )              
      ( "Ru4(OH)12+4", 5.5000 )              
      ( "Sn+4", 11.0000 )              
      ( "Th+4", 11.0000 )              
      ( "U+4", 11.0000 )              
      ( "V2(OH)2+4", 5.5000 )              
      ( "Y2(OH)2+4", 5.5000 )              
      ( "Zr+4", 11.0000 )              
      ( "Al3(OH)4+5", 6.0000 )              
      ( "Cr3(OH)4+5", 6.0000 )              
      ( "Fe3(OH)4+5", 6.0000 )              
      ( "Th2(OH)3+5", 6.0000 )              
      ( "Ce2(OH)2+6", 6.0000 )              
      ( "La5(OH)9+6", 6.0000 )              
      ( "Th2(OH)2+6", 6.0000 )
      ( "Al13O4(OH)24+7", 6.0000 );
#endif
}


/////////////////////////////////////////////////////////////////////////////////////////////////////

namespace GPChemistry
{
  void GetElementCountFromFormula(std::string formulaOrig, std::map<std::string, realT>& elementCount)
  {
  
    std::string formula = formulaOrig; // keep copy for error reporting
    sArray1d elVect;
    Trim(formula);
    
    // electron 
    ///////////
    if( streq(formula,"e-") || streq(formula,"e") ){
      elementCount["e"] += 1;
      return;	
    }
    
    // remove valence
    //////////////////
    size_t valenceStart = formula.find_first_of("-+");
    formula = formula.substr(0,valenceStart);
    
    // remove brackets
    //////////////////
    size_t nextClosingBracket = formula.find_first_of(")");
    while(nextClosingBracket < formula.size()){

      // find matching opening bracket
      size_t openingBracket = formula.find_last_of("(", nextClosingBracket);
      if(openingBracket == std::string::npos){
      	throw GPException("Error: GetElementCountFromFormula - mismatched brackets in formula: "+ formulaOrig);
      }

      // find multiple
      size_t afterNumber = formula.find_first_not_of("0123456789.", nextClosingBracket+1);
      size_t l = afterNumber-nextClosingBracket-1;
      realT multiple = 1.0;
      if(l > 0) multiple = fromString<realT>( formula.substr(nextClosingBracket+1,l) );
      
      // loop over elements in bracket string, identify coefficients and apply multiple
      std::string bracketStr = formula.substr(openingBracket+1, 
                                              nextClosingBracket-openingBracket-1);
   
      elVect = SplitAtCapitals(bracketStr);
      bracketStr = "";
      for(unsigned i =0; i < elVect.size();++i)
      {
         std::string& elStr = elVect[i];
         realT subscript = GetSubscriptFromElementStr(elStr);
         bracketStr += elStr+toString(subscript*multiple);
      }
   
      // replace substring 
      formula.replace(openingBracket, afterNumber-openingBracket, bracketStr);

      // search for back bracket
      nextClosingBracket = formula.find_first_of(")");
    }

    // Count elements
    //////////////////

    elVect = SplitAtCapitals(formula);

    for(unsigned i =0; i < elVect.size();++i){
      std::string& element = elVect[i];
 
      realT subscript = GetSubscriptFromElementStr(element);

      elementCount[element] += subscript;
    }
  } 
  
  void ReadXMLEquilibriumData(TICPP::HierarchicalDataNode* hdn,
                            std::map< std::string, unsigned>& speciesIndicies, 
                            std::map< std::string, rArray1d >& equations,
                            std::map< std::string, rArray1d >& log10Kdata){
                            	
    int rank(0);
    #if GPAC_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif
  

    if (rank == 0)  std::cout << "Equilibrium Data" << std::endl;
    TICPP::HierarchicalDataNode* eqNode = hdn->GetChild("EquilibriumData");
    if(!eqNode)	
        throw GPException("ReactionFrontSolver: Must have EquilibriumData defined in the input file");
    {
      sArray1d analyticStrs, eqStrs;
      GetVectorAttribute(eqStrs,eqNode, "equilibrium_equations",";");
      GetVectorAttribute(analyticStrs,eqNode, "logKdata",";");
  	 
  	 
      for(sArray1d::size_type j = 0; j < eqStrs.size(); ++j){
    	sArray1d subStrs = Tokenize(eqStrs[j],"=");
    	Trim(subStrs[0]);
  	   	std::string kstr =  subStrs[0];
  	 	sArray1d numDenom = Tokenize(subStrs[1],"/");
  	 	
  	    rArray1d eqVector(speciesIndicies.size(),0.0);
  	    
  	      // numerator
  	    sArray1d nums = Tokenize(numDenom[0],"[");
  	    for(sArray1d::size_type i =0; i < nums.size(); ++i){
  	      realT p = 1.0;
  	      if(nums[i].size() > 0){
  	        sArray1d sp_pow = Tokenize(nums[i],"]"); // check for powers
  	        
  	        Trim(sp_pow[0]);
  	        if(sp_pow[0].size() > 0  && !streq(sp_pow[0],"H2O") ){
  	          
  	          if(!isMember(sp_pow[0],speciesIndicies))
  	   	        throw GPException("Equilibrium equation component " + sp_pow[0] + " not found in species list");
  	   	        
  	          int indx = speciesIndicies[sp_pow[0]];
  	        
  	          // get power
  	          if(sp_pow.size() > 1){
  	            Trim(sp_pow[1]);
  	            if (sp_pow[1].size() > 0){
  	        	  std::replace(sp_pow[1].begin(),sp_pow[1].end(),'^',' '); // remove "^"
  	         	  p = fromString<realT>(sp_pow[1]); 
  	            }
  	          }
  	          eqVector[indx] += p;
  	        }
  	      }
  	    }
  	    
  	    // denominator
  	 	if(numDenom.size() > 1){
  	      sArray1d denoms = Tokenize(numDenom[1],"[");
  	      for(sArray1d::size_type i =0; i < denoms.size(); ++i){
  	        realT p = 1.0;
  	        if(denoms[i].size() > 0){
  	          sArray1d sp_pow = Tokenize(denoms[i],"]"); // check for powers
  	        
  	          Trim(sp_pow[0]);
  	          if(sp_pow[0].size() > 0 && !streq(sp_pow[0],"H2O") ){
  	           
  	            if(!isMember(sp_pow[0],speciesIndicies) )
  	   	          throw GPException("Equilibrium equation component " + sp_pow[0] + " not found in species list");
  	   	        
  	            int indx = speciesIndicies[sp_pow[0]];
  	        
  	            // get power
  	            if(sp_pow.size() > 1){
  	              Trim(sp_pow[1]);
  	              if (sp_pow[1].size() > 0){
  	        	    std::replace(sp_pow[1].begin(),sp_pow[1].end(),'^',' '); // remove "^"
  	        	    p = fromString<realT>(sp_pow[1]); 
  	              }
  	            }
  	            eqVector[indx] -= p;
  	          }
  	        }
  	  	  }
  	    }
  	 	
  	    equations[kstr] = eqVector;
  	    
  	    // print data
  	    if (rank == 0){
  	      std::cout << "    " << kstr << ":\t";
  	      for(sArray1d::size_type i =0; i < eqVector.size(); ++i) std::cout << " "<< eqVector[i];
  	      std::cout << std::endl;
  	    }
  	  }
  	 
  	  if (rank == 0)  std::cout << "Log K data" << std::endl;
  	  for(sArray1d::size_type j = 0; j < analyticStrs.size(); ++j){
  	 	sArray1d aSubStr = Tokenize(analyticStrs[j],",");
  	 	std::string kstr =  aSubStr[0];
  	 	Trim(kstr);
  	    rArray1d analyticVector(5,0.0);
  	    int iend = (6<=aSubStr.size())? 6: aSubStr.size();
  	    for(int i = 1; i < iend ; ++i){
  	    	analyticVector[i-1] = fromString<realT>(aSubStr[i]);
  	    }
  	    log10Kdata[kstr] = analyticVector;
  	    
  	    // print data
  	    if (rank == 0){
  	      std::cout << "    " << kstr << ":\t";
  	      for(rArray1d::size_type i =0; i < analyticVector.size(); ++i) std::cout << " "<< analyticVector[i];
  	      std::cout << std::endl;
  	    }
  	  }
    }
  }
  
  /// returns valence of element in formula, string is unaffected. 
  realT GetValence(std::string formula){
  	return getValenceAndAdjustStr(formula);
  }
  
}



ExtendedDebyeHuckelActivityFunction::
ExtendedDebyeHuckelActivityFunction(const rArray1d& valences,realT A, realT B, const rArray1d& ao, realT Bdot):
m_valences(valences),
m_A(A),
m_B(B),
m_Bdot(Bdot),
m_ao(ao)
{
/* empty */
}

ExtendedDebyeHuckelActivityFunction::
ExtendedDebyeHuckelActivityFunction(const sArray1d& species,realT Tc)
{
  int rank(0);
  #if GPAC_MPI
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  unsigned n = species.size();
  m_valences = GetValences(species);
  ChemistryManager& chemistryManager = ChemistryManager::Instance();
  m_A = chemistryManager.m_DebyeHuckelA.Lookup(&Tc);
  m_B = chemistryManager.m_DebyeHuckelB.Lookup(&Tc);
  m_Bdot = chemistryManager.m_WATEQBdot.Lookup(&Tc);
  /*
  if(rank == 0){
	  std::cout << "A " << m_A << std::endl;
	  std::cout << "B " << m_B << std::endl;
	  std::cout << "Bdot " << m_Bdot << std::endl;
  }
  */
  if(rank == 0) std::cout << "Extended Debye Huckel Activity Function: \n";
  m_ao.resize(n);
  for(unsigned i =0 ; i < n; ++i){
    if(isMember(species[i],chemistryManager.m_ionSizeParameters)){
      m_ao[i] = chemistryManager.m_ionSizeParameters[species[i]];
      if(rank == 0) std::cout << species[i] << " Ion size parameter: " << m_ao[i] << "\n";
    } else {
      m_ao[i] = 4.0;
      if(rank == 0) std::cout << "ExtendedDebyeHuckelActivityFunction: Warning - did not find species " << species[i] << " in ion size parameter lookup table. \n"
                << "Setting ion size parameter to " << m_ao[i] << std::endl;
    }
  }
}

DaviesActivityFunction::
DaviesActivityFunction(const rArray1d& valences,realT A):
m_valences(valences),
m_A(A)
{
/* empty */
}

DaviesActivityFunction::
DaviesActivityFunction(const sArray1d& species,realT Tc)
{
  m_valences = GetValences(species);
  ChemistryManager& chemistryManager = ChemistryManager::Instance();
  m_A = chemistryManager.m_DebyeHuckelA.Lookup(&Tc);
}

DaviesActivityFunction::
DaviesActivityFunction(const DaviesActivityFunction& rhs)
{
  m_valences = rhs.m_valences;
  m_A = rhs.m_A;
}



SITActivityFunction::
SITActivityFunction(const sArray1d& species,realT Tc,
                    const std::string& background_cation,
                    const std::string& background_anion):
m_Ba_j(1.5),
m_background_cation(background_cation),
m_background_anion(background_anion),
m_cationIndx(999999),
m_anionIndx(999999)
{
  m_valences = GetValences(species);
  ChemistryManager& chemistryManager = ChemistryManager::Instance();
  m_A = chemistryManager.m_DebyeHuckelA.Lookup(&Tc);


  unsigned numSpecies = species.size();
  m_e = std::vector<realT> (numSpecies,0.0);

  std::pair< std::string,std::string> cationAnionPair;

  sArray1d anList;
  anList += "Cl-","ClO4-","NO3-";

  sArray1d catList;
  catList += "Na+","K+","Li+";

  for(unsigned i =0; i < numSpecies; ++i){
    std::string speciesA = species[i];
    if(m_valences[i] != 0.0){

      bool isCation = (m_valences[i] > 0.0);
      if(isCation){
        cationAnionPair.first  = species[i];
        cationAnionPair.second = background_anion;
        if(background_cation == species[i]) m_cationIndx = i;
      }else {
        cationAnionPair.first  = background_cation;
        cationAnionPair.second = species[i];
        if(background_anion == species[i]) m_anionIndx = i;
      }

      if( !isMember(cationAnionPair,chemistryManager.m_SITParameters) ){

        bool matchFound = false;
        if(isCation){
          for(unsigned j =0 ; j < anList.size()&& !matchFound; ++j){
            cationAnionPair.second = anList[j];
            matchFound = isMember(cationAnionPair,chemistryManager.m_SITParameters);
          }
        } else {
          for(unsigned j =0 ; j < catList.size()&& !matchFound; ++j){
            cationAnionPair.first = catList[j];
            matchFound = isMember(cationAnionPair,chemistryManager.m_SITParameters) ;
          }
        }

        if(!matchFound)
          throw GPException("Error SITActivityFunction: Ion pair not found for " + species[i] + ".");


      }

      m_e[i] = chemistryManager.m_SITParameters[cationAnionPair];
    }
  }

  if(m_cationIndx == 999999 )
    throw GPException("Error SITActivityFunction: Background cation " + m_background_cation + " not found in species list.");

  if(m_anionIndx  == 999999)
    throw GPException("Error SITActivityFunction: Background anion " + m_background_anion + " not found in species list.");
}

SITActivityFunction::
SITActivityFunction(const SITActivityFunction& rhs)
{
  m_valences = rhs.m_valences;
  m_A = rhs.m_A;
  m_Ba_j = rhs.m_Ba_j;
  m_background_cation = rhs.m_background_cation;
  m_background_anion = rhs.m_background_anion;
}


