
#ifndef QCD_H_
#define QCD_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "ObjectManagers/ChemistryManager.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"

#include "Common/Common.h"
#include "Utilities/Functions.h"
#include "Utilities/Utilities.h"

/// 
class QueryChemistryDatabase : public SolverBase
{
public:
  QueryChemistryDatabase(const std::string& name,
                         ProblemManagerT* const pm );
  virtual ~QueryChemistryDatabase(){};
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn) ;
  void RegisterFields( PhysicalDomainT& domain  ){};
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);
  void InitializeCommunications( PartitionBase& partition  ){};


  double TimeStep( const realT& time,
                 const realT& dt ,
                 const int cycleNumber,
                 PhysicalDomainT& domain ,
                 const sArray1d& namesOfSolverRegions ,
                 SpatialPartition& partition ,
                 FractunatorBase* const fractunator ){ return dt; };

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "QueryChemistryDatabase";};

private:
   
   // Data
   struct{
     sArray1d names;
     unsigned size;
   } m_species;


   struct{
     sArray1d names;
     unsigned size;
   } m_solidphases;

   
   std::vector<GPChemistry::ChemicalReaction>  m_solutionReactions;
   std::vector<GPChemistry::ChemicalReaction>  m_solidPhaseReactions;
   std::set<std::string> m_reactionSpecies;
 
};


#endif
