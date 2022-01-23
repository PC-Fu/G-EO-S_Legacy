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
#ifndef geochemicalreactionsolver_H
#define geochemicalreactionsolver_H

/**
 * @file GeochemicalReactionSolver.h
 * @author walsh24
 */

#include "Common/Common.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/FunctionManager.h"
#include "ObjectManagers/FaceManagerT.h"
#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"
#include "PhysicsSolvers/SolverBase.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "PhysicsSolvers/PhysicsSolverStrings.h"
#include "OneDReactionFront.h"


#include "ObjectManagers/ChemistryManager.h"


#include "Utilities/Functions.h"
#include "Utilities/Utilities.h"


class GeochemicalReactionSolver : public SolverBase
{
public:
  GeochemicalReactionSolver( const std::string& name,
                             ProblemManagerT* const pm );
  virtual ~GeochemicalReactionSolver(){};
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;
  void RegisterFields( PhysicalDomainT& domain ){};
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);
  void InitializeCommunications( PartitionBase& partition ){};

  virtual void TimeStep(const realT& time,
                        const realT& dt,
                        const int cycleNumber,
                        PhysicalDomainT& domain,
                        const sArray1d& namesOfSolverRegions,
                        SpatialPartition& partition,
                        FractunatorBase* const fractunator);

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "GeochemicalReactionSolver";};

private:

   void BuildEquilibriumEquations();
   
   // Data
   struct{
     sArray1d names;
     unsigned size;
     rArray1d brineConcentrations;
   } m_species;
   
   // fixme need to remove H,e,O from elements
   struct{
     sArray1d names;
     std::map<std::string, unsigned > indexMap;
     rArray2d concentrationEquations; 
     unsigned size;
   } m_elements;

   // fixme need to make solid phases only those phases listed in input file
   struct{
     sArray1d names;
     unsigned size;
     sArray1d reactionRateFunctionNames;
     sArray1d reactionRateFunctionVariables;
     std::vector< Function* > reactionRateFunctionPtrs;

     std::vector<GPChemistry::ChemicalReaction> reactions;
     rArray2d equilibriumEquations; 
     rArray1d logK;
     rArray2d reactionElements; // returns elements produced when dissolved 
     rArray1d volumeFractions;
   } m_solidphases;

   
   std::vector<GPChemistry::ChemicalReaction>  m_solutionReactions;
   std::vector<GPChemistry::ChemicalReaction>  m_solidPhaseReactions;
   std::set<std::string> m_reactionSpecies;
   
   std::set<std::string> m_solutionSpecies;
   std::set<std::string> m_solutionMasterSpecies;
   
   std::map<std::string, unsigned > m_speciesIndexMap;
   sArray1d m_indexSpeciesMap; // change to m_species.names when working
   std::vector< std::map<std::string, realT> > m_speciesElements;
   rArray1d m_valences;
   
   realT m_temperature;

   // aqueous species equilibrium, logK
   rArray2d m_equilibriumEquations; 
   rArray1d m_LogK;


   // equilibrium system      
   ODRF::EquilibriumSystem m_brineSystem;  // Brine equilibrium

   // face set   
   const lSet* m_faceSet;
   std::string m_faceSetName;


   // conversion factors
   realT m_convertToMolPerL;
   realT m_convertFromMolPerL;
   realT m_convertFromMolPerLSMsqrd;
   
   // parallel communication
   sArray1d synced_face_fields; 

   realT m_traceConcentration;
 
};

#endif
