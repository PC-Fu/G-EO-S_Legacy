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
 * @file UpdateFieldWithFunction.h
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

//namespace{

  /// Base class for non-linear systems of equations
  class NonLinearSystemBase{

    public:
      NonLinearSystemBase(){};
      virtual ~NonLinearSystemBase(){};
      virtual void GetFunction(const rArray1d& X,rArray1d& F) = 0;
      virtual void GetJacobian(const rArray1d& X,rArray2d& gradF) = 0;

  };
	
  /// 
  class EquilibriumSystem: public NonLinearSystemBase{
    public:
      EquilibriumSystem():
        m_K(),m_p(),m_A(),m_b(){};
      EquilibriumSystem(const rArray2d& K, const rArray1d& p, 
                        const rArray2d& A, const rArray1d& b):
        m_K(K),m_p(p),m_A(A),m_b(b){};
      ~EquilibriumSystem(){};
      void GetFunction(const rArray1d& X,rArray1d& F);
      void GetJacobian(const rArray1d& X,rArray2d& gradF);
    private:
      rArray2d m_K;
      rArray1d m_p;
      rArray2d m_A;
      rArray1d m_b;
  };
	
//}

/// 
class FindEquilibrium : public SolverBase
{
public:
  FindEquilibrium( const std::string& name,
                   ProblemManagerT* const pm );
  virtual ~FindEquilibrium();
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition  ) {}


  double TimeStep( const realT& time,
                 const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions,
                 SpatialPartition& partition,
                 FractunatorBase* const fractunator );

  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "FindEquilibrium";};

private:

   sArray1d m_species;
   Array1dT<FieldType> m_variable_types;
   
   PhysicalDomainT::ObjectDataStructureKeys m_objectKey;
   std::string m_regionName; // only used if setting field in an element region
   
   EquilibriumSystem m_equilibriumSystem;
  
};

#endif /* UPDATEFIELDWFUNCTION_H_ */
