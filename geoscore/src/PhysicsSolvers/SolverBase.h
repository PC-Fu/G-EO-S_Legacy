//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2014, Lawrence Livermore National Security, LLC.
//
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//
//  Randolph Settgast		Stuart Walsh
//  Scott Johnson		Pengcheng Fu
//  Joshua White
//
//  LLNL-CODE-656616
//  GEOS-CORE, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GEOS-CORE. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
//
//  Please also read "Additional BSD Notice" below.
//
//  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//  * Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the 
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Additional BSD Notice
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file SolverBase.h
 * @author settgast1
 * @date Feb 10, 2011
 */

#ifndef SOLVERBASE_H_
#define SOLVERBASE_H_


#include "Common/typedefs.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "MPI_Communications/SpatialPartition.h"
#include "ObjectManagers/PhysicalDomainT.h"
//#include "ObjectManagers/ProblemManagerT.h"

#if GPAC_MPI
  class Epetra_MpiComm;
#else
class Epetra_SerialComm;
#endif

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_SolverMap_CrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "AztecOO.h"

#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Thyra_OperatorVectorClientSupport.hpp"
#include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
#include "Thyra_AztecOOLinearOpWithSolve.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraLinearOpBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultInverseLinearOp.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_MLPreconditionerFactory.hpp"


#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace ML_Epetra
{ class MultiLevelPreconditioner; }

class AztecOO;

class PhysicalDomainT;
class ProblemManagerT;
class FractunatorBase;

namespace EpetraBlock
{
  enum ID
  {
    displacement,
    fluidMass,
    numBlockDof
  };
}

/**
 * @author settgast
 * @note class to hold the epetra system matrices and vectors.
 */
class Epetra_System
{

public:

  Epetra_System():
    m_nblocks(1),
    m_blockID(),
    m_solverNames(),
    m_solverNameMap()
  {
    for( int i=0 ; i<EpetraBlock::numBlockDof ; ++i )
    {
      m_blockID[i] = -1;
    }


    m_rowMap.resize(m_nblocks);
    m_solution.resize(m_nblocks);
    m_rhs.resize(m_nblocks);
    m_sparsity.resize2(m_nblocks,m_nblocks);
    m_matrix.resize2(m_nblocks,m_nblocks);
//    m_scratch.resize2(m_nblocks,m_nblocks);

#if USECPP11!=1
    for( int i=0 ; i<m_nblocks ; ++i )
    {
      m_rowMap[i] = NULL;
      m_solution[i] = NULL;
      m_rhs[i] = NULL;
      for( int j=0 ; j<m_nblocks ; ++j )
      {
        m_sparsity[i][j] = NULL;
        m_matrix[i][j] = NULL;
      }
    }
#endif

  }


  ~Epetra_System()
  {
#if USECPP11!=1
    for( int i=0 ; i<m_nblocks ; ++i )
    {
      delete m_rowMap[i];
      delete m_solution[i];
      delete m_rhs[i];
      for( int j=0 ; j<m_nblocks ; ++j )
      {
        delete m_sparsity[i][j];
        delete m_matrix[i][j];
      }
    }
#endif

  }

  void SetNumBlocks( const int numBlocks )
  {
    m_nblocks = numBlocks;
    m_rowMap.resize(m_nblocks);
    m_solution.resize(m_nblocks);
    m_rhs.resize(m_nblocks);
    m_sparsity.resize2(m_nblocks,m_nblocks);
    m_matrix.resize2(m_nblocks,m_nblocks);
//    m_scratch.resize2(m_nblocks,m_nblocks);
  }

  void SetBlockID( const int block, const EpetraBlock::ID id, const std::string& name )
  {
    if( block>=m_nblocks )
    {
      throw GPException("SolverBase.h:Epetra_System::SetBlockID(): block index is larger than number of blocks\n");
    }
    else if( block<0 )
    {
      throw GPException("SolverBase.h:Epetra_System::SetBlockID(): block index is less than 0\n");
    }

    if( m_blockID[id]==-1 )
    {
      m_blockID[id] = block;
    }
    else
    {
      throw GPException("error in Epetra_System::SetBlockID(). EpetraBlock::ID ("+toString(id)+") has already been used to register block "+ toString(block)+ ". Ya betta check yoself befo ya wreck yoself.\n" );
    }

    if( m_solverNames[id].empty() )
    {
      m_solverNames[id] = name;
    }
    else
    {
      throw GPException("error in Epetra_System::SetBlockID(). EpetraBlock::ID ("+toString(id)+") has already been used to register solvername "+ m_solverNames[id]+ ". Ya betta check yoself befo ya wreck yoself.\n" );
    }

    std::map<std::string, int>::iterator iterSolverNameMap = m_solverNameMap.find(name);
    if( iterSolverNameMap==m_solverNameMap.end() )
    {
      m_solverNameMap[name] = block;
    }
    else
    {
      throw GPException("error in Epetra_System::SetBlockID(). Solver Name ("+name+") has already been used to register block "+ toString(block)+ ". Ya betta check yoself befo ya wreck yoself.\n" );
    }
  }

  int GetBlockID( EpetraBlock::ID id ) const
  {
    return m_blockID[id];
  }

  std::string GetSolverName( EpetraBlock::ID id ) const
  {
    return m_solverNames[id];
  }


#if USECPP11==1
  std::shared_ptr<Epetra_Map>&
#else
  Epetra_Map*&
#endif
  GetRowMap( const EpetraBlock::ID dofID )
  {
    if( m_blockID[dofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetRowMap():m_blockID isn't set \n");
    }
    return m_rowMap[ m_blockID[dofID] ];
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FEVector>&
#else
  Epetra_FEVector*&
#endif
  GetSolutionVector( const EpetraBlock::ID dofID )
  {
    if( m_blockID[dofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetSolutionVector():m_blockID isn't set \n");
    }
    return m_solution[ m_blockID[dofID] ];
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FEVector>&
#else
  Epetra_FEVector*&
#endif
  GetResidualVector( const EpetraBlock::ID dofID )
  {
    if( m_blockID[dofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetResidualVector():m_blockID isn't set \n");
    }
    return m_rhs[ m_blockID[dofID] ];
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FECrsGraph>&
#else
  Epetra_FECrsGraph*&
#endif
  GetSparsity( const EpetraBlock::ID rowDofID,
               const EpetraBlock::ID colDofID )
  {
    if( m_blockID[rowDofID]==-1 && m_blockID[colDofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetSparsity():m_blockID isn't set \n");
    }
    return m_sparsity[ m_blockID[rowDofID] ][ m_blockID[colDofID] ] ;
  }

#if USECPP11==1
  std::shared_ptr<Epetra_FECrsMatrix>&
#else
  Epetra_FECrsMatrix*&
#endif
  GetMatrix( const EpetraBlock::ID rowDofID,
             const EpetraBlock::ID colDofID )
  {
    if( m_blockID[rowDofID]==-1 && m_blockID[colDofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetMatrix():m_blockID isn't set \n");
    }
    return m_matrix[ m_blockID[rowDofID] ][ m_blockID[colDofID] ] ;
  }
/*
#if USECPP11==1
  std::shared_ptr<Epetra_FECrsMatrix>
#else
  Epetra_FECrsMatrix*&
#endif
  GetScratch( const EpetraBlock::ID rowDofID, const EpetraBlock::ID colDofID )
  {
    if( m_blockID[rowDofID]==-1 && m_blockID[colDofID]==-1 )
    {
      throw GPException("SolverBase.h:Epetra_System::GetScratch():m_blockID isn't set \n");
    }
    return m_scratch[ m_blockID[rowDofID] ][ m_blockID[colDofID] ] ;
  }
*/
#if USECPP11==1
  Array1dT<std::shared_ptr<Epetra_Map> >         m_rowMap;
  Array1dT<std::shared_ptr<Epetra_FEVector> >    m_solution;
  Array1dT<std::shared_ptr<Epetra_FEVector> >    m_rhs;
  Array2dT<std::shared_ptr<Epetra_FECrsGraph> >  m_sparsity;
  Array2dT<std::shared_ptr<Epetra_FECrsMatrix> > m_matrix;
//  Array2dT<std::shared_ptr<Epetra_FECrsMatrix> > m_scratch;
#else
  Array1dT<Epetra_Map*>         m_rowMap;
  Array1dT<Epetra_FEVector*>    m_solution;
  Array1dT<Epetra_FEVector*>    m_rhs;
  Array2dT<Epetra_FECrsGraph*>  m_sparsity;
  Array2dT<Epetra_FECrsMatrix*> m_matrix;
#endif



private:
  int m_nblocks;
  int m_blockID[EpetraBlock::numBlockDof];
  std::string m_solverNames[EpetraBlock::numBlockDof];
  std::map<std::string, int> m_solverNameMap;

  Epetra_System( Epetra_System& );
  Epetra_System& operator=(Epetra_System&);

};

class SolverBase
{
public:
  SolverBase(const std::string& name,
             ProblemManagerT* const pm );
  virtual ~SolverBase();

  virtual double TimeStep( const realT& time,
                         const realT& dt,
                         const int cycleNumber,
                         PhysicalDomainT& domain,
                         const sArray1d& namesOfSolverRegions,
                         SpatialPartition& partition,
                         FractunatorBase* const fractunator ) = 0;

  virtual realT UpdateTimeStepMid(realT dt) { return dt; }
  
  virtual void PostProcess( PhysicalDomainT& domain,
                            SpatialPartition& partition,
                            const sArray1d& namesOfSolverRegions);

  virtual void SetMaxStableTimeStep( PhysicalDomainT& domain,
                                     const sArray1d& namesOfSolverRegions,
                                     SpatialPartition& partition);

  virtual void Initialize( PhysicalDomainT& domain, SpatialPartition& partition ) = 0;

  virtual void InitializeCommunications( PartitionBase& partition ) = 0;

  /// returns name of specific solver instances
  std::string Name(){return m_name;};

  virtual void RegisterFields( PhysicalDomainT& domain ) = 0;

  virtual void ReadXML( TICPP::HierarchicalDataNode* const hdn ) ;

  virtual void WriteSilo( SiloFile& siloFile ) const;

  virtual void ReadSilo( const SiloFile& siloFile );

  StableTimeStep m_stabledt;
  realT m_courant;



  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> m_syncedFields;

  std::string m_name;


#if USECPP11==1
  void SetRowMapPtr( std::shared_ptr<Epetra_Map>& rowMap )            { m_rowMap=rowMap; }
  void SetSparsityPtr( std::shared_ptr<Epetra_FECrsGraph>& sparsity ) { m_sparsity = sparsity; }
  void SetMatrixPtr( std::shared_ptr<Epetra_FECrsMatrix>& matrix )    { m_matrix = matrix; }
  void SetSolutionPtr( std::shared_ptr<Epetra_FEVector>& solution )   { m_solution = solution; }
  void SetRhsPtr( std::shared_ptr<Epetra_FEVector>& rhs )             { m_rhs = rhs; }
  void SetEpetraSystemPtr( std::shared_ptr<Epetra_System>& epetraSystem ) { m_epetraSystem = epetraSystem; }
#else
  void SetRowMapPtr( Epetra_Map* const rowMap )            { m_rowMap=rowMap; }
  void SetSparsityPtr( Epetra_FECrsGraph* const sparsity ) { m_sparsity = sparsity; }
  void SetMatrixPtr( Epetra_FECrsMatrix* const matrix )    { m_matrix = matrix; }
  void SetSolutionPtr( Epetra_FEVector* const solution )   { m_solution = solution; }
  void SetRhsPtr( Epetra_FEVector* const rhs )             { m_rhs = rhs; }
#endif


  void Solve ( PhysicalDomainT&  domain,
               SpatialPartition& partition,
               const realT time,
               const realT dt );


  void SolveBlock( PhysicalDomainT&  domain,
                   SpatialPartition& partition,
                   Epetra_System& epetraSystem,
                   const realT time,
                   const realT dt );

  // TODO these should be pure virtuals;
  virtual void SetInitialGuess( const PhysicalDomainT& domain,
                                realT* const local_solution ) {}

  virtual void SetupMLPreconditioner( const PhysicalDomainT& domain, ML_Epetra::MultiLevelPreconditioner* MLPrec ){}

  virtual void SetLinearSolverParameters( AztecOO& solver );


  virtual realT CheckSolution( const realT* const local_solution,
                               const PhysicalDomainT& domain,
                               const localIndex dofOffset )
  { return 1.0; }

  virtual realT CheckSolutionBlock( const PhysicalDomainT& domain)
  { return 1.0; }


  virtual void PropagateSolution( const realT* const local_solution,
                                  const realT scalingFactor,
                                  PhysicalDomainT& domain,
                                  const localIndex dofOffset ) {}

  virtual void PropagateSolutionBlock( const realT scalingFactor,
                                       PhysicalDomainT& domain) {}


  virtual void PostSyncConsistency( PhysicalDomainT& domain, SpatialPartition& partition ) {}

  const std::string& TrilinosIndexString() const
  {
    return m_TrilinosIndexStr;
  }

  Epetra_System* m_system;

protected:

  std::string m_TrilinosIndexStr;


  struct
  {
    realT krylov_tol;       // Solver convergence criteria
    int   m_maxIters;       // Maximum number of solver iterations
    int   m_kspace;         // Number of krylov vectors before GMRES restart
    realT m_ilut_fill;      // Fill factor for ILUT preconditioner
    realT m_ilut_drop;      // Drop tolerance for ILUT preconditioner
    bool  m_useMLPrecond;   // Use ML preconditioner
    bool  m_useRowScaling;  // Use row scaling
    bool  m_useBicgstab;    // Use bicgstab instead of gmres
    bool  m_verbose;        // print extra info
    bool  m_useDirectSolver; // Use Direct solver
    bool  m_useNewtonSolve;  // Use Newton-Raphson iterations
    int   m_maxIterNewton;   // Maximum number of Newton-Raphson iterations
    realT m_tolNewton;       // tolerance for Newton convergence
  }
  m_numerics;

#if GPAC_MPI
  const Epetra_MpiComm* epetra_comm;
#else
  const Epetra_SerialComm* epetra_comm;
#endif


#if USECPP11==1

  std::shared_ptr<Epetra_System> m_epetraSystem;

  std::shared_ptr<Epetra_Map>         m_rowMap;
  std::shared_ptr<Epetra_FECrsGraph>  m_sparsity;
  std::shared_ptr<Epetra_FECrsMatrix> m_matrix;
  std::shared_ptr<Epetra_FEVector>    m_solution;
  std::shared_ptr<Epetra_FEVector>    m_rhs;

#else

  Epetra_Map*         m_rowMap;
  Epetra_FECrsGraph*  m_sparsity;
  Epetra_FECrsMatrix* m_matrix;
  Epetra_FEVector*    m_solution;
  Epetra_FEVector*    m_rhs;

#endif


private:

  std::string m_lsp;


  int m_timeIntegrationFlag;

  virtual void SetTimeIntegrationFlag( )  {}

  virtual void TimeStepDerived( const realT& time,
                                 const realT& dt,
                                 PhysicalDomainT& domain,
                                 const sArray1d& namesOfSolverRegions,
                                 SpatialPartition& partition ) {};


  virtual void WriteSiloDerived( SiloFile& siloFile ) const {}
  virtual void ReadSiloDerived( const SiloFile& siloFile ) {}

};


#endif /* SOLVERBASE_H_ */
