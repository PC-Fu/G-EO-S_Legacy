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
 * @file SteadyStateParallelPlateFlowSolver_TwoD.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef REACTIONFRONTMATERIALUPDATE_H_
#define REACTIONFRONTMATERIALUPDATE_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"
#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include <set>

#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif


class ReactionFrontMaterialUpdate : public SolverBase
{
public:
  ReactionFrontMaterialUpdate( const std::string& name,
                               ProblemManagerT* const pm );
  virtual ~ReactionFrontMaterialUpdate();
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);
  void InitializeCommunications( PartitionBase & partition  ){};

  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 
  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "ReactionFrontMaterialUpdate";};

protected:
  void UpdateElementModuli(ElementRegionT& elementRegion, localIndex element, realT distance);

private:
  sArray1d m_frontNames;
  localIndex m_numFronts;
  rArray1d m_bulkModuli;
  rArray1d m_shearModuli;
  const lSet* m_faceSet;
  std::string m_faceSetName;

  rArray1d m_frontDistances;
  
};


#endif /* REACTIONFRONTMATERIALUPDATE_H_ */
