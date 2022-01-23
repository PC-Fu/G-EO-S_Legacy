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
 * @file FindEquilibrium.h
 * @author walsh24
 * @date July 25, 2011
 */

#ifndef FINDEQUILIBRIUM_H_
#define FINDEQUILIBRIUM_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

#include "Common/Common.h"

#include "PhysicsSolvers/Chemistry/OneDReactionFront.h"

#include "Utilities/ChemistryUtilities.h"


/// 
class ChemicalEquilibriumSolver : public SolverBase
{
public:
  ChemicalEquilibriumSolver( const std::string& name,
                             ProblemManagerT* const pm );
  virtual ~ChemicalEquilibriumSolver();
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition  ) {}


  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber, PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "ChemicalEquilibriumSolver";};

private:

   sArray1d m_setNames;
   
   PhysicalDomainT::ObjectDataStructureKeys m_objectKey;
   std::string m_regionName; // only used if setting field in an element region
   
   ODRF::EquilibriumSystem m_brineSystem;  // Brine equilibrium

   // Data
   CHEM_UTILS::SpeciesData m_species;

   CHEM_UTILS::ElementsData m_elements;

   CHEM_UTILS::EquilibriumData m_equilibriumData;

   // Temperature   
   realT m_Tc; // Temperature in Celcius
 
   // equilibrium equations
   rArray2d m_KK; // equilibrium data for aqueous phases
   rArray1d m_pp;

   rArray2d m_AA; // constraint equations for aqueous species
   rArray1d m_bb;
   rArray2d m_AA_e; // constraint equations for aqueous elements (used to form m_bb)

   // conversion factors
   realT m_convertFromMolPerL;
   realT m_convertToMolPerL;

   // data pointers
   std::vector<Array1dT<realT>* > m_elementConcPtrs; // element concentrations

   //pH
   Array1dT<realT>* m_pHFieldPtr; // pH
   localIndex m_H_index;

   // flags
   bool m_recordBrineSpecies;
   std::vector<Array1dT<realT>* > m_brineSpeciesConcPtrs; // only used if brine species are recorded

   void UpdateBrineConstraints(void);
   void UpdateEquilibrium(localIndex indx );
};

#endif /* UPDATEFIELDWFUNCTION_H_ */
