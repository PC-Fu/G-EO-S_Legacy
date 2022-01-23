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
 * @file ChemistryUtilities.cpp
 * @author walsh24
 * @date July 25, 2011
 */


#include "ChemistryUtilities.h"
#include "Utilities/StringUtilities.h"
#include "ObjectManagers/ChemistryManager.h"
#include "ObjectManagers/FunctionManager.h"

#include "Common/intrinsic_typedefs.h"
#include "Common/Common.h"

using namespace CHEM_UTILS;

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

void CHEM_UTILS::SpeciesData::ReadXML( TICPP::HierarchicalDataNode* speciesNode ){


	int rank(0);

	#if GPAC_MPI
	    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif


    GetVectorAttribute(this->names,speciesNode, "names");
    this->size = this->names.size();
    if (rank == 0)  std::cout << "    ";
    for(unsigned i =0; i < this->size; ++i){
      //if (rank == 0)  std::cout << this->names[i] << " ";
      this->indices[this->names[i]] = i;
    }
    if (rank == 0)  std::cout << "\n";
    
    // valences    
    this->valences.resize(this->size);
    this->ionicStrengthCoefficients.resize(this->size);
    for(unsigned spIndx=0; spIndx < this->size; ++spIndx){
      this->valences[spIndx] = GPChemistry::GetValence(this->names[spIndx]);
      if (rank == 0)   std::cout<< "    " << this->names[spIndx] << ": " << this->valences[spIndx] << std::endl ;
      this->ionicStrengthCoefficients[spIndx] = 0.5*this->valences[spIndx]*this->valences[spIndx];
    }
    
    GetVectorAttribute(this->eqConstants,speciesNode, "equilibrium_constants");
    
    this->brineConcentrations.resize(this->size);

}


void CHEM_UTILS::ElementsData::ReadXML( TICPP::HierarchicalDataNode* elementNode, const CHEM_UTILS::SpeciesData& speciesData ){

	int rank(0);

	#if GPAC_MPI
	    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

     GetVectorAttribute(this->names,elementNode, "names");
     this->size = this->names.size();

     for(unsigned i = 0; i < this->size; ++i)
       this->indices[this->names[i]] = i;
     
  

     // get element count in derived species
     this->species.resize(this->size,rArray1d(speciesData.size,0.0));
     
     for(unsigned i =0; i < speciesData.size; ++i){
       
       std::string formula = speciesData.names[i];
       std::map<std::string, realT> elMap;
       GPChemistry::GetElementCountFromFormula(formula, elMap);

       std::map<std::string, realT>::iterator
         itr = elMap.begin(),
         iend = elMap.end();
       for(;itr != iend; ++itr){
         const std::string& el = itr->first;
         realT& count = itr->second;
         if(isMember(el,this->indices) ){  // ignores H,O
           unsigned elIndex = this->indices[el];
           this->species[elIndex][i] = count;
         }
       } 
     }

      // print data
     for(unsigned i =0; i < this->size; ++i){
       if (rank == 0){
         std::cout << "    " << this->names[i] << ":\t";
         for(unsigned j = 0; j < speciesData.size; ++j)
             std::cout << " " << this->species[i][j];
         std::cout << std::endl;
       }
     }


     GetVectorAttribute(this->brineConcentrations, elementNode,"brine_concentrations");
     if(this->brineConcentrations.size()  == 0)
         this->brineConcentrations= rArray1d(this->size,0.0);

     if(this->brineConcentrations.size() != this->size )
         throw GPException("ElementData: Number of brine concentrations must match number of elements. ");
}


void CHEM_UTILS::EquilibriumData::ReadXML( TICPP::HierarchicalDataNode* eqNode, const CHEM_UTILS::SpeciesData& speciesData ){

	int rank(0);

	#if GPAC_MPI
	    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif

     sArray1d analyticStrs, eqStrs;
     GetVectorAttribute(eqStrs,eqNode, "equilibrium_equations",";");
     GetVectorAttribute(analyticStrs,eqNode, "logKdata",";");
     
     
     for(unsigned j = 0; j < eqStrs.size(); ++j){
       sArray1d subStrs = Tokenize(eqStrs[j],"=");
       Trim(subStrs[0]);
       std::string kstr =  subStrs[0];
       sArray1d numDenom = Tokenize(subStrs[1],"/");
       
        rArray1d eqVector(speciesData.size,0.0);
        
        // numerator
        sArray1d nums = Tokenize(numDenom[0],"[");
        for(unsigned i =0; i < nums.size(); ++i){
          realT p = 1.0;
          if(nums[i].size() > 0){
            sArray1d sp_pow = Tokenize(nums[i],"]"); // check for powers
            
            Trim(sp_pow[0]);
            if(sp_pow[0].size() > 0  && !streq(sp_pow[0],"H2O") ){
              
              if(!isMember(sp_pow[0],speciesData.indices))
                 throw GPException("Equilibrium equation component " + sp_pow[0] + " not found in species list");
                 
              int indx = speciesData.indices.at(sp_pow[0]);
            
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
          for(unsigned i =0; i < denoms.size(); ++i){
            realT p = 1.0;
            if(denoms[i].size() > 0){
              sArray1d sp_pow = Tokenize(denoms[i],"]"); // check for powers
            
              Trim(sp_pow[0]);
              if(sp_pow[0].size() > 0 && !streq(sp_pow[0],"H2O") ){
               
                if(!isMember(sp_pow[0],speciesData.indices) )
                   throw GPException("Equilibrium equation component " + sp_pow[0] + " not found in species list");
                 
                int indx = speciesData.indices.at(sp_pow[0]);
            
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
       
        this->equations[kstr] = eqVector;
        
       // print data
        if (rank == 0){
          std::cout << "    " << kstr << ":\t";
          for(unsigned i =0; i < eqVector.size(); ++i) std::cout << " "<< eqVector[i];
          std::cout << std::endl;
        }
     }
     
     if (rank == 0)  std::cout << "Log K data" << std::endl;
     for(unsigned j = 0; j < analyticStrs.size(); ++j){
       sArray1d aSubStr = Tokenize(analyticStrs[j],",");
       std::string kstr =  aSubStr[0];
       Trim(kstr);
        rArray1d analyticVector(5,0.0);
        unsigned iend = (6<=aSubStr.size())? 6: aSubStr.size();
        for(unsigned i = 1; i < iend ; ++i){
          analyticVector[i-1] = fromString<realT>(aSubStr[i]);
        }
        this->log10Kdata[kstr] = analyticVector;
        
       // print data
      /** 
        if (rank == 0){
          realT logK =  Log10K(analyticVector, m_Tc); // calculate value at default temperature
          
          std::cout << "    " << kstr << ":\t" << logK << "\t(";
          for(unsigned i =0; i < analyticVector.size(); ++i) std::cout << " "<< analyticVector[i];
          std::cout << ")" << std::endl;
        }
      **/
     }
}

/**

Builds a set of species constraint equations from strings in the form:
q_Ca+2 + q_CaHCO3+ + q_CaCO3 = 3*q_CaHCO3+ + 3*q_CaCO3 + 3*q_CO3-2 + 3*q_HCO3- + 3*q_H2CO3 + 3*q_CO2
q_Ca+2 + q_CaHCO3+ + q_CaCO3 = 3*Q_C

Lowercase prefix indicates species constraints
Uppercase prefix indicates element constraints

Currently no distinction for different prefixes. In future may expand constrain equations to allow input of different types based on prefix, i.e. q_Ca+2 (flux of Ca+2), A_C (amount of C).

**/
void CHEM_UTILS::ParseSpeciesConstraintStr(const std::string& theStr,const CHEM_UTILS::ElementsData& elementsData, const CHEM_UTILS::SpeciesData& speciesData , realT* eqVector,bool isLHS){

  int rank(0);

  #if GPAC_MPI
	    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  sArray1d lhs = TokenizeSeq(theStr," +"); // note there must be a space before the plus
  unsigned n = speciesData.size;
  for(unsigned iii = 0; iii < lhs.size(); ++iii){
    realT coef = 1.0;
    std::string spStr;
    char preChar = ' ';
    sArray1d coef_sp = Tokenize(lhs[iii],"*");
    if(coef_sp.size() > 1){
      size_t indx = coef_sp[1].find("_");
      spStr = coef_sp[1].substr(indx+1);
      Trim(coef_sp[0]);
      coef =  fromString<realT>(coef_sp[0]);
      if (rank == 0)  std::cout << coef <<std::endl;
      if(indx > 0) preChar = coef_sp[1][indx-1];
    } else {
      size_t indx = coef_sp[0].find("_");
      spStr = coef_sp[0].substr(indx+1);
      if(indx > 0) preChar = coef_sp[0][indx-1];
    }
    Trim(spStr);
    if(!streq(spStr,"0")){
      if( islower(preChar) ){ // lowercase indicates individual species constraints
        if( isMember(spStr,speciesData.indices) ){
          int indx = speciesData.indices.at(spStr);
          eqVector[indx] += (isLHS)? coef:-coef;
        } else {
          throw GPException("Species constraint component " + spStr + " not found in species list. Did you mean to use an uppercase prefix instead (element constraint)?");
        }
      } else {  // uppercase indicates element constraints
        if( isMember(spStr,elementsData.indices) ){
          int e_indx = elementsData.indices.at(spStr);
          realT coefB = (isLHS)? coef:-coef;
          for(unsigned ii=0; ii < n; ++ii){
            eqVector[ii] +=coefB*elementsData.species[e_indx][ii];
          }
        } else {
          throw GPException("Species constraint component " + spStr + " not found in element list. Did you mean to use a lowercase prefix instead (species constraint)?");
        }
      }
    }
  }
}

/**

Builds a set of element constraint equations from strings in the form:
Q_Ca = 3*Q_C

Lowercase prefix indicates species constraints
Uppercase prefix indicates element constraints

Currently no distinction for different prefixes. In future may expand constrain equations to allow input of different types based on prefix, i.e. Q_Ca (flux of Ca), A_C (amount of C).

**/
void CHEM_UTILS::ParseElementConstraintStr(const std::string& theStr,const CHEM_UTILS::ElementsData& elementsData, realT* eqVector,bool isLHS){

  int rank(0);

  #if GPAC_MPI
		    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  sArray1d lhs = TokenizeSeq(theStr," +"); // note there must be a space before the plus

  for(unsigned iii = 0; iii < lhs.size(); ++iii){
    realT coef = 1.0;
    std::string spStr;
    char preChar = ' ';
    sArray1d coef_sp = Tokenize(lhs[iii],"*");
    if(coef_sp.size() > 1){
      size_t indx = coef_sp[1].find("_");
      spStr = coef_sp[1].substr(indx+1);
      Trim(coef_sp[0]);
      coef =  fromString<realT>(coef_sp[0]);
      //if (rank == 0)  std::cout << coef <<std::endl;
      if(indx > 0) preChar = coef_sp[1][indx-1];
    } else {
      size_t indx = coef_sp[0].find("_");
      spStr = coef_sp[0].substr(indx+1);
      if(indx > 0) preChar = coef_sp[0][indx-1];
    }
    Trim(spStr);
    if(!streq(spStr,"0")){
      if( islower(preChar) ){ // lowercase indicates individual species constraints
         throw GPException("ParseElementConstraintStr: Element constraint parsing does not support species constraints. Did you mean to use an uppercase prefix instead (element constraint)?");

      } else {  // uppercase indicates element constraints
        if( isMember(spStr,elementsData.indices) ){
          int e_indx = elementsData.indices.at(spStr);
          eqVector[e_indx] += (isLHS)? coef:-coef;
        } else {
          throw GPException("ParseElementConstraintStr: Front condition component " + spStr + " not found in element list.");
        }
      }
    }
  }
}






