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
 * @file GroupLabeler.h
 * @author walsh24
 * @date Nov 21, 2013
 */

#ifndef GroupLabeler_H_
#define GroupLabeler_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"
#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include "src_external/Utilities/DisjointSet.h"


#include <set>

#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

/*
 * Label connected groups (currently faces, based on edge connectivity) that pass a particular test (currently memberField > 0)
 * In the future hope to extend to arbitrary objects and connecting objects with arbitrary membership tests
 */

class GroupLabeler : public SolverBase
{
public:
  GroupLabeler( const std::string& name,
                               ProblemManagerT* const pm );
  virtual ~GroupLabeler();
  
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
  static const char* SolverName(){return "GroupLabeler";};


private:

  std::string m_isMemberFieldName; // > 0 if member  <=0 otherwise (ultimately will replace this with user defined member test function)
  std::string m_groupIdFieldName; // name of field to assign group id to
};


#endif /* GroupLabeler_H_ */
