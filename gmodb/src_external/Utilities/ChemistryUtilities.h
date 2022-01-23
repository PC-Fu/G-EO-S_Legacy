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
 * @file ChemistryUtilities.h
 * @author walsh24
 * @date July 25, 2011
 */

#ifndef CHEMISTRY_UTILITIES_H_
#define CHEMISTRY_UTILITIES_H_

#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

#include "Common/Common.h"

namespace CHEM_UTILS{

   struct SpeciesData{
     sArray1d names;
     unsigned size;
     std::map<std::string,int> indices;
     iArray1d valences;
     rArray1d ionicStrengthCoefficients; // 0.5*valence^2
     sArray1d eqConstants;
     rArray1d brineConcentrations;  // Concentrations followed by their natural log
     rArray1d frontConcentrations;  // Concentrations followed by their natural log

     public:
       void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;

   };

   struct ElementsData{
     sArray1d names;
     unsigned size;
     std::map<std::string,int> indices;
     std::vector<rArray1d > species;
     rArray1d brineConcentrations;  // Net element concentrations

     public:
       void ReadXML( TICPP::HierarchicalDataNode* hdn, const SpeciesData& speciesData ) ;
   };

   struct EquilibriumData{
     std::map< std::string, rArray1d > equations;
     std::map< std::string, rArray1d > log10Kdata;

     public:
       void ReadXML( TICPP::HierarchicalDataNode* hdn, const SpeciesData& speciesData  ) ;
   };

   void ParseSpeciesConstraintStr(const std::string& theStr,const CHEM_UTILS::ElementsData& elementsData,
		                   const CHEM_UTILS::SpeciesData&  speciesData, realT* eqVector,bool isLHS=true);

   void ParseElementConstraintStr(const std::string& theStr, const CHEM_UTILS::ElementsData& elementsData , realT* eqVector,bool isLHS=true);

}


#endif /* CHEMISTRY_UTILITIES_H_ */
