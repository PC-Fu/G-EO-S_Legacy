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
 * @file SolverBase.cpp
 * @author settgast1
 * @date Feb 10, 2011
 */

#include "SolverBase.h"
#include "ObjectManagers/ProblemManagerT.h"
#include "Amesos.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Utilities/TrilinosUtilities.h"

SolverBase::SolverBase( const std::string& name,
                        ProblemManagerT* const pm ):
m_stabledt(),
m_courant(1.0),
m_syncedFields(),
m_name(name),
epetra_comm(&(pm->m_epetraComm))
#if USECPP11!=1
,
m_system(NULL),
m_rowMap(NULL),
m_sparsity(NULL),
m_matrix(NULL),
m_solution(NULL),
m_rhs(NULL),
m_TrilinosIndexStr()
#endif

{
  m_stabledt.m_maxdt = std::numeric_limits<double>::max();
}

SolverBase::~SolverBase()
{
#if USECPP11!=1
//  delete m_rowMap;
//  delete m_sparsity;
//  delete m_matrix;
//  delete m_solution;
//  delete m_rhs;
#endif
}


void SolverBase::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  m_numerics.krylov_tol        = hdn->GetAttributeOrDefault("tol",1.0e-6);
  m_numerics.m_kspace          = hdn->GetAttributeOrDefault("kspace",300);
  m_numerics.m_ilut_fill       = hdn->GetAttributeOrDefault("ilut_fill",3.0);
  m_numerics.m_ilut_drop       = hdn->GetAttributeOrDefault("ilut_drop",0.0);
  m_numerics.m_maxIters        = hdn->GetAttributeOrDefault<int>("maxSolverIterations",10000);
  m_numerics.m_useMLPrecond    = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);
  m_numerics.m_useRowScaling   = hdn->GetAttributeOrDefault<bool>("useRowScaling",false);
  m_numerics.m_useBicgstab     = hdn->GetAttributeOrDefault<bool>("useBicgstab",false);
  m_numerics.m_verbose         = hdn->GetAttributeOrDefault<bool>("verbose",false);
  m_numerics.m_useDirectSolver = hdn->GetAttributeOrDefault<bool>("useDirectSolver",false);
  m_numerics.m_useNewtonSolve  = hdn->GetAttributeOrDefault<bool>("useNewtonSolve",false);
  m_numerics.m_maxIterNewton   = hdn->GetAttributeOrDefault<int>("maxIterNewton",30);
  m_numerics.m_tolNewton       = hdn->GetAttributeOrDefault<realT>("tolNewton",1e-06);

  m_courant = hdn->GetAttributeOrDefault<realT>("courant",0.5);
  m_lsp = hdn->GetAttributeStringOrDefault("lsp","Klu");
}

void SolverBase::PostProcess( PhysicalDomainT& domain,
                              SpatialPartition& partition,
                              const sArray1d& namesOfSolverRegions)
{
  //Do nothing.  To be implemented by each solver as needed;
  //For operations that we don't need at each step but need for plotting.
}


void SolverBase::SetMaxStableTimeStep( PhysicalDomainT& domain,
                                       const sArray1d& namesOfSolverRegions,
                                       SpatialPartition& partition)
{
  TimeStep( 0.0, 0.0, 0, domain, namesOfSolverRegions, partition, NULL );
}


void SolverBase::WriteSilo( SiloFile& siloFile ) const
{
  siloFile.DBWriteWrapper( "m_stabledt__m_maxdt", m_stabledt.m_maxdt );

  std::map<std::string, sArray1d> syncedFields;
  for( std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d>::const_iterator i=m_syncedFields.begin() ;
      i!= m_syncedFields.end() ; ++i )
  {
    if( !(i->second.empty()) )
    {
      syncedFields[ PhysicalDomainT::GetObjectDataStructureName( i->first ) ] = i->second;
    }
  }

  siloFile.DBWriteWrapper( "syncedFields", syncedFields );

  siloFile.DBWriteWrapper( "m_name", m_name );

}

void SolverBase::ReadSilo( const SiloFile& siloFile )
{
  siloFile.DBReadWrapper( "m_stabledt__m_maxdt", m_stabledt.m_maxdt );

  std::map<std::string, sArray1d> syncedFields;
  siloFile.DBReadWrapper( "syncedFields", syncedFields );


  for( std::map<std::string, sArray1d>::const_iterator i=syncedFields.begin() ;
      i!= syncedFields.end() ; ++i )
  {
    m_syncedFields[ PhysicalDomainT::GetObjectDataStructureKey(i->first)] = i->second;
  }


  siloFile.DBReadWrapper( "m_name", m_name );

}



void SolverBase::Solve ( PhysicalDomainT&  domain,
                         SpatialPartition& partition,
                         const realT time,
                         const realT dt )
{

  // initial guess for solver
  int dummy;
  double* local_solution = NULL;

  m_solution->ExtractView(&local_solution,&dummy);

  SetInitialGuess( domain, local_solution  );

  if(m_numerics.m_useRowScaling)
  {
      printf("Scaling matrix/rhs in place\n");
      Epetra_Vector scaling(m_matrix->RowMap());
      m_matrix->InvRowSums(scaling); 
      m_matrix->LeftScale(scaling);

      Epetra_MultiVector tmp (*m_rhs);
      m_rhs->Multiply(1.0,scaling,tmp,0.0);
  }

#if USECPP11==1
  Epetra_LinearProblem problem( m_matrix.get(),
                                m_solution.get(),
                                m_rhs.get() );
#else
  Epetra_LinearProblem problem( m_matrix,
                                m_solution,
                                m_rhs );

#endif

  // @annavarapusr1: Needed to use direct solver without changing it for everyone else
  if(m_numerics.m_useDirectSolver) // If Chandra's test problems, use direct solver
  {
    Amesos_BaseSolver* Solver;
    Amesos Factory;

    std::string SolverType = "Klu";

    if( !(m_lsp.compare("Slu")) )
    {
      SolverType = "Amesos_Superludist";
    }
    Solver = Factory.Create(SolverType, problem);

    int ierr = Solver->SymbolicFactorization();
    if (ierr > 0)
      std::cerr << "ERROR!" << std::endl;

    ierr = Solver->NumericFactorization();
    if (ierr > 0)
      std::cerr << "ERROR!" << std::endl;

    Solver->Solve();

    if( Solver!=NULL )
      delete Solver;
  }
  else
  {
    AztecOO solver(problem);

    SetLinearSolverParameters( solver );

    ML_Epetra::MultiLevelPreconditioner* MLPrec = NULL;
    //SetupMLPreconditioner( domain, MLPrec );

    if( m_numerics.m_useMLPrecond )
    {
      Teuchos::ParameterList MLList;
      ML_Epetra::SetDefaults("SA",MLList);
      //MLList.set("aggregation: type", "Uncoupled");
      //MLList.set("aggregation: type", "MIS");
      MLList.set("prec type", "MGW");
      MLList.set("smoother: type","ILU");
      MLList.set("ML output",0);
      MLList.set("PDE equations",3);

#if USECPP11==1
      MLPrec = new ML_Epetra::MultiLevelPreconditioner(*m_matrix.get(), MLList);
#else
      MLPrec = new ML_Epetra::MultiLevelPreconditioner(*m_matrix, MLList);
#endif
      solver.SetPrecOperator(MLPrec);
    }
    else // use ILUT preconditioner with domain decomp
    {
      solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
      solver.SetAztecOption(AZ_subdomain_solve,AZ_ilut);
      solver.SetAztecParam(AZ_ilut_fill,m_numerics.m_ilut_fill);
      solver.SetAztecParam(AZ_drop,m_numerics.m_ilut_drop);
    }

    solver.Iterate(m_numerics.m_maxIters,
                   m_numerics.krylov_tol);

    if(m_numerics.m_useMLPrecond)
    {
      delete MLPrec; MLPrec=NULL;
    }
  }

  // copy vector solution into geos data structures

  realT scalingFactor = CheckSolution( local_solution, domain, 0 );
  PropagateSolution( local_solution, scalingFactor, domain, 0 );

  // re-sync ghost nodes

  partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);

  // copy vector solution into geos data structures

  PostSyncConsistency( domain, partition );
}


void SolverBase::SetLinearSolverParameters( AztecOO& solver )
{
    if(m_numerics.m_useBicgstab)
    {
      solver.SetAztecOption(AZ_solver,AZ_bicgstab);
    }
    else
    {
      solver.SetAztecOption(AZ_solver,AZ_gmres);
    }
    solver.SetAztecOption(AZ_conv,AZ_r0);
    solver.SetAztecOption(AZ_kspace,m_numerics.m_kspace);
    solver.SetAztecOption(AZ_output,AZ_all);
}


void SolverBase::SolveBlock ( PhysicalDomainT&  domain,
                              SpatialPartition& partition,
                              Epetra_System& epetraSystem,
                              const realT time,
                              const realT dt )
{

		// code fragment to write out linear system to matlab
		// for debugging purposes

  EpetraExt::RowMatrixToMatlabFile("umatrix00.dat",*epetraSystem.m_matrix[0][0]);
  EpetraExt::RowMatrixToMatlabFile("umatrix01.dat",*epetraSystem.m_matrix[0][1]);
  EpetraExt::RowMatrixToMatlabFile("umatrix10.dat",*epetraSystem.m_matrix[1][0]);
  EpetraExt::RowMatrixToMatlabFile("umatrix11.dat",*epetraSystem.m_matrix[1][1]);

  EpetraExt::MultiVectorToMatlabFile("urhs0.dat",*epetraSystem.m_rhs[0]);
  EpetraExt::MultiVectorToMatlabFile("urhs1.dat",*epetraSystem.m_rhs[1]);

		// there are several flags to control solver behavior.
		// these should be compared in a scaling study.
		//
		// 1. whether to use inner solvers or just the 
		//    sub-block preconditioners directly. false
		//    is probably better.
		// 2. whether to use a block diagonal or a full
		//    block triangular preconditioner.  false is
		//    probably better.
		// 3. whether to perform an explicit row scaling
		//    of the linear system before solving.  note
		//    that the matrix and rhs are modified in place
		//    by this operation.  true is probably better.
		// 4. whether to use BiCGstab or GMRES for the 
		//    krylov solver.  GMRES is generally more robust,
		//    BiCGstab sometimes shows better parallel performance.
		//    false is probably better.

  const bool use_inner_solver  = false;
  const bool use_diagonal_prec = false;
  const bool use_row_scaling   = true;
  const bool use_bicgstab      = false;

		// perform an explicit row scaling of the linear system,
		// R*A*x = R*b, where R is a diagonal scaling matrix.
		// we will use inverse row sums for the scaling.

		// todo: recall that the matrix row stretches across two
		// blocks.  to compute the true row sum we would have
		// to add across both blocks.  for now we just use the row
		// sum from the diagonal block for the scaling. 

  if(use_row_scaling)
  {

    for(unsigned b=0; b<2; ++b)
    {
      Epetra_Vector scale_one(epetraSystem.m_matrix[b][b]->RowMap());
      Epetra_Vector scale_two(epetraSystem.m_matrix[b][b]->RowMap());

      epetraSystem.m_matrix[b][0]->InvRowSums(scale_one); 
      epetraSystem.m_matrix[b][1]->InvRowSums(scale_two);

      scale_one.Reciprocal(scale_one);
      scale_two.Reciprocal(scale_two);  // not ideal, could choke if 1/0 or 1/NaN appears
      scale_one.Update(1.0,scale_two,1.0);
      scale_one.Reciprocal(scale_one);

      for(unsigned c=0; c<2; ++c)
        epetraSystem.m_matrix[b][c]->LeftScale(scale_one);

      Epetra_MultiVector tmp (*epetraSystem.m_rhs[b]);
      epetraSystem.m_rhs[b]->Multiply(1.0,scale_one,tmp,0.0);
    }

		// legacy version, where we only did the row sum
		// for the diagonal blocks.  keeping this here
		// because it may be more robust against 1/0 errors.

    /*
    for( int a=0 ; a<2 ; ++a )
    {
      Epetra_Vector scaling(epetraSystem.m_matrix[a][a]->RowMap());
      scaling.Scale(0.0);

#if 0
        TrilinosUtilities::RowSum(*(epetraSystem.m_matrix[a][a]),scaling);
#else
        TrilinosUtilities::RowSum(*(epetraSystem.m_matrix[a][0]),*(epetraSystem.m_matrix[a][1]),scaling);
#endif
      scaling.Reciprocal(scaling);
      for(unsigned b=0; b<2; ++b)
      {
        epetraSystem.m_matrix[a][b]->LeftScale(scaling);
      }
      Epetra_MultiVector tmp (*epetraSystem.m_rhs[a]);
      epetraSystem.m_rhs[a]->Multiply(1.0,scaling,tmp,0.0);

    }
    */

  }


		// try adding a small diagonal perturbation
		// to the system to see if that improves
		// convergence behavior

  /*
  {
    for(unsigned b=0; b<1; ++b) // just A block for now
    {
      Epetra_Vector perturb(epetraSystem.m_matrix[b][b]->RowMap());
      Epetra_Vector ones(epetraSystem.m_matrix[b][b]->RowMap());
                    ones.PutScalar(1.0);

      epetraSystem.m_matrix[b][b]->ExtractDiagonalCopy(perturb); 
      perturb.Update(-1.0e-3,ones,1.0);
      epetraSystem.m_matrix[b][b]->ReplaceDiagonalValues(perturb); 
    }
  }
  */

		// set initial guess to zero.  this is not strictly
		// necessary but is good for comparing solver performance.

  //epetraSystem.m_solution[0]->PutScalar(0.0);
  //epetraSystem.m_solution[1]->PutScalar(0.0);

		// we want to use thyra to wrap epetra operators and vectors
		// for individual blocks.  this is an ugly conversion, but 
		// it is basically just window dressing.
		//
		// note the use of Teuchos::RCP reference counted pointers.
		// The general syntax is usually one of:
		//
		//   RCP<T> Tptr = rcp(new T) 
		//   RCP<T> Tptr = nonMemberConstructor();
		//   RCP<T> Tptr (t_ptr,false)
		//
		// where "false" implies the RCP does not own the object and
		// should not attempt to delete it when finished.


  RCP<const Thyra::LinearOpBase<double> >  matrix_block[2][2];
  RCP<Thyra::MultiVectorBase<double> >     lhs_block[2];
  RCP<Thyra::MultiVectorBase<double> >     rhs_block[2];


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  std::cout<<"breakpoint 1"<<std::endl;
  //MPI_Barrier(MPI_COMM_WORLD);
  for(unsigned i=0; i<2; ++i)
  for(unsigned j=0; j<2; ++j)
  {
    std::cout<<"breakpoint 2, rank = "<<rank<<", epetraSystem.m_matrix["<<i<<"]["<<j<<"]->NumMyRows() = "<<epetraSystem.m_matrix[i][j]->NumMyRows()<<std::endl;
    RCP<Epetra_Operator> mmm (&*epetraSystem.m_matrix[i][j],false);
    std::cout<<"breakpoint 3, rank = "<<rank<<std::endl;
    matrix_block[i][j] = Thyra::epetraLinearOp(mmm);
    std::cout<<"breakpoint 4, rank = "<<rank<<std::endl;
  }
  

  for(unsigned i=0; i<2; ++i)
  { 
    RCP<Epetra_MultiVector> lll (&*epetraSystem.m_solution[i],false);
    RCP<Epetra_MultiVector> rrr (&*epetraSystem.m_rhs[i],false);

    lhs_block[i] = Thyra::create_MultiVector(lll,matrix_block[i][i]->domain());
    rhs_block[i] = Thyra::create_MultiVector(rrr,matrix_block[i][i]->range());
  }

		// now use thyra to create an operator representing
		// the full block 2x2 system

  RCP<const Thyra::LinearOpBase<double> > matrix = Thyra::block2x2(matrix_block[0][0],
                                                                   matrix_block[0][1],
                                                                   matrix_block[1][0],
                                                                   matrix_block[1][1]);

		// creating a representation of the blocked
		// rhs is a little uglier. (todo: check if there is
		// a cleaner way to do this.)

  RCP<Thyra::ProductMultiVectorBase<double> > rhs;

//  std::cout<<"breakpoint 6"<<std::endl;

  {
    Teuchos::Array<RCP<Thyra::MultiVectorBase<double> > > mva;
    Teuchos::Array<RCP<const Thyra::VectorSpaceBase<double> > > mvs;

    for(unsigned i=0; i<2; ++i)
    {
      mva.push_back(rhs_block[i]);
      mvs.push_back(rhs_block[i]->range());
    }

    RCP<const Thyra::DefaultProductVectorSpace<double> > vs = Thyra::productVectorSpace<double>(mvs);

    rhs = Thyra::defaultProductMultiVector<double>(vs,mva);
  }

		// do the identical operation for the lhs

  RCP<Thyra::ProductMultiVectorBase<double> > lhs;

  {
    Teuchos::Array<RCP<Thyra::MultiVectorBase<double> > > mva;
    Teuchos::Array<RCP<const Thyra::VectorSpaceBase<double> > > mvs;

    for(unsigned i=0; i<2; ++i)
    {
      mva.push_back(lhs_block[i]);
      mvs.push_back(lhs_block[i]->range());
    }

    RCP<const Thyra::DefaultProductVectorSpace<double> > vs = Thyra::productVectorSpace<double>(mvs);

    lhs = Thyra::defaultProductMultiVector<double>(vs,mva);
  }


		// from randy:
		// create diagonal schur complement approximation
		// S = C - B2 diag(A)inv B1

#if 0
  RCP<const Thyra::LinearOpBase<double> > SchurEstimate;
  RCP<const Thyra::LinearOpBase<double> >  invDiag00LinOp;
  Epetra_CrsMatrix SchurEstimate00( Copy, *(epetraSystem.m_rowMap[1]),0 );
  if(1)
  {
    /*
    Epetra_FECrsGraph sparsity_00_diag( Copy, *(epetraSystem.m_rowMap[0]), 1 );

    iArray1d nonEmptyDOF;
    for( int a=0; a<epetraSystem.m_matrix[0][0]->NumMyRows(); ++a )
    {
      nonEmptyDOF.push_back(a);
    }
    sparsity_00_diag.InsertIndices( nonEmptyDOF.size(),
                                    nonEmptyDOF.data(),
                                    nonEmptyDOF.size(),
                                    nonEmptyDOF.data() );
    */

    Epetra_CrsMatrix invDiag00( Copy, *(epetraSystem.m_rowMap[0]), 0 );



    Epetra_Vector diagVec( *(epetraSystem.m_rowMap[0]) );
    Epetra_Vector invDiagVec( *(epetraSystem.m_rowMap[0]) );
    epetraSystem.m_matrix[0][0]->ExtractDiagonalCopy( diagVec );
    invDiagVec.Reciprocal( diagVec );

//    std::cout<<diagVec<<std::endl;
//    std::cout<<invDiagVec<<std::endl;

    invDiag00.ReplaceDiagonalValues( invDiagVec );


    Epetra_CrsMatrix scratch( Copy, *(epetraSystem.m_rowMap[0]),0 );

    EpetraExt::MatrixMatrix::Multiply( invDiag00, false,
                                       *(epetraSystem.m_matrix[1][0]), false,
                                       scratch );

    EpetraExt::MatrixMatrix::Multiply( *(epetraSystem.m_matrix[0][1]), false,
                                       scratch, false,
                                       SchurEstimate00 );

    EpetraExt::MatrixMatrix::Add( *(epetraSystem.m_matrix[1][1]), false, 1.0,
                                  SchurEstimate00, -1.0 );

    SchurEstimate00.FillComplete();
    
    // JAW write for checking
    EpetraExt::RowMatrixToMatlabFile("schur_est.dat",SchurEstimate00);

//    invDiag00LinOp = Thyra::epetraLinearOp( RCP<Epetra_Operator>( &invDiag00, false ) );
    RCP<Epetra_Operator> mmm (&SchurEstimate00,false);
    SchurEstimate = Thyra::epetraLinearOp(mmm);

//    std::cout<<*(epetraSystem.m_matrix[1][1])<<std::endl;
//    std::cout<<SchurEstimate00<<std::endl;
  }
//  SchurEstimate = matrix_block[1][1];

//  SchurEstimate = Thyra::add(matrix_block[1][1],Thyra::scale(-1.0,Thyra::multiply(matrix_block[0][1],invDiag00LinOp,matrix_block[1][0])));
#endif


		// for the preconditioner, we need two approximate inverses,
		// one for the (0,0) block and one for the approximate
		// schur complement.  for now, we will use the (1,1) block
		// as our schur complement approximation, though we should
		// explore better approaches later. 

		// we store both "sub operators" in a 1x2 array:

  RCP<const Thyra::LinearOpBase<double> > sub_op[2];

		// each implicit "inverse" is based on an inner krylov solver,
		// with their own sub-preconditioners.  this leads to a very
		// accurate approximation of the inverse operator, but can be
		// overly expensive.  the other option is to ditch the inner
		// krylov solver, and just use the sub-preconditioners directly.

		// the implicit inverse for each diagonal block is built in 
		// three steps:
		//   1.  define solver parameters
		//   2.  build a solver factory
		//   3.  build the inner solver operator
	

  for(unsigned i=0; i<2; ++i) // loop over diagonal blocks
  {
    RCP<Teuchos::ParameterList> list = rcp(new Teuchos::ParameterList("solver_list"),true);

      list->set("Linear Solver Type","AztecOO");
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Max Iterations",m_numerics.m_maxIters);
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Tolerance",1e-3*m_numerics.krylov_tol);
      if(use_bicgstab)
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","BiCGStab");
      else 
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","GMRES");
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",0);//int(m_numerics.m_verbose));
      //list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Convergence Test","rhs");

      //if(i==0)
      //{
        if(m_numerics.m_useMLPrecond)
        {
          list->set("Preconditioner Type","ML");
          list->sublist("Preconditioner Types").sublist("ML").set("Base Method Defaults","SA");
          list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("aggregation: type","Uncoupled");
          list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("smoother: type","ILU");
          list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("ML output",0);
          list->sublist("Preconditioner Types").sublist("ML").sublist("ML Settings").set("PDE equations",(i==0?3:1));
        }
        else
        {
          list->set("Preconditioner Type","Ifpack");
          list->sublist("Preconditioner Types").sublist("Ifpack").set("Prec Type","ILU");
        }
      //}
      //else
      //{
      //  list->set("Preconditioner Type","None");
      //}

    Stratimikos::DefaultLinearSolverBuilder builder;

      builder.setParameterList(list);

    if(use_inner_solver)
    {
      RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > strategy = createLinearSolveStrategy(builder);

      //if(i==0)
        sub_op[i] = Thyra::inverse(*strategy,matrix_block[i][i]);
      //else
      //{
      //  RCP<const Thyra::LinearOpBase<double> > BAinvBt = Thyra::multiply(matrix_block[0][1],sub_op[0],matrix_block[1][0]);
      //  RCP<const Thyra::LinearOpBase<double> > schur = Thyra::add(matrix_block[1][1],Thyra::scale(-1.0,BAinvBt));
      //  sub_op[i] = Thyra::inverse(*strategy,schur);
      //}
    }
    else
    {
      RCP<const Thyra::PreconditionerFactoryBase<double> > strategy = createPreconditioningStrategy(builder);
      RCP<Thyra::PreconditionerBase<double> > tmp;

      //if(i==0) 
        tmp = prec(*strategy,matrix_block[i][i]);
      //else     
      //  tmp = prec(*strategy,SchurEstimate);
 
     sub_op[i] = tmp->getUnspecifiedPrecOp();
    }
  }
 

		// create zero operators for off diagonal blocks

  RCP<const Thyra::LinearOpBase<double> > zero_01
    = rcp(new Thyra::DefaultZeroLinearOp<double>(matrix_block[0][0]->range(),
                                                 matrix_block[1][1]->domain())); 

  RCP<const Thyra::LinearOpBase<double> > zero_10
    = rcp(new Thyra::DefaultZeroLinearOp<double>(matrix_block[1][1]->range(),
                                                 matrix_block[0][0]->domain())); 

		// now build the block preconditioner

  RCP<const Thyra::LinearOpBase<double> > preconditioner;

  if(use_diagonal_prec)
  {
    preconditioner = Thyra::block2x2(sub_op[0],zero_01,zero_10,sub_op[1]);
  }
  else
  {
    RCP<const Thyra::LinearOpBase<double> > eye_00
      = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(matrix_block[0][0]->range())); 
 
    RCP<const Thyra::LinearOpBase<double> > eye_11
      = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(matrix_block[1][1]->range())); 

    /*
    RCP<const Thyra::LinearOpBase<double> > P_A,P_B1,P_B2,P_S;

    P_A  = Thyra::block2x2(sub_op[0],zero_01,zero_10,eye_11);
    P_B1 = Thyra::block2x2(eye_00,Thyra::scale(-1.0,matrix_block[0][1]),zero_10,eye_11);
    P_B2 = Thyra::block2x2(eye_00,zero_01,Thyra::scale(-1.0,matrix_block[1][0]),eye_11);
    P_S  = Thyra::block2x2(eye_00,zero_01,zero_10,sub_op[1]);
    
    preconditioner = Thyra::multiply(P_A,P_B1,P_S);
    //preconditioner = Thyra::multiply(P_S,P_B2,P_A);
    */

    RCP<const Thyra::LinearOpBase<double> > mAinvB1, mB2Ainv; 
    
    mAinvB1 = Thyra::scale(-1.0, Thyra::multiply(sub_op[0],matrix_block[0][1]) );
    mB2Ainv = Thyra::scale(-1.0, Thyra::multiply(matrix_block[1][0],sub_op[0]) );

    RCP<const Thyra::LinearOpBase<double> > Linv,Dinv,Uinv,Eye;
    
    //Eye = Thyra::block2x2(eye_00,zero_01,zero_10,eye_11);
    Dinv = Thyra::block2x2(sub_op[0],zero_01,zero_10,sub_op[1]);
    Linv = Thyra::block2x2(eye_00,zero_01,mB2Ainv,eye_11);
    Uinv = Thyra::block2x2(eye_00,mAinvB1,zero_10,eye_11);
    
    //preconditioner = Eye;
    //preconditioner = Dinv;
    //preconditioner = Thyra::multiply(Uinv,Dinv);
    //preconditioner = Thyra::multiply(Dinv,Linv);
    preconditioner = Thyra::multiply(Uinv,Dinv,Linv);
  }


		// define solver strategy for blocked system. this is
		// similar but slightly different from the sub operator
		// construction, since now we have a user defined preconditioner

  {
    RCP<Teuchos::ParameterList> list = rcp(new Teuchos::ParameterList("list"),true);

      list->set("Linear Solver Type","AztecOO");
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Max Iterations",m_numerics.m_maxIters);
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").set("Tolerance",m_numerics.krylov_tol);
      if(use_bicgstab)
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","BiCGStab");
      else 
        list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver","GMRES");
      list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency",1);
      //list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Convergence Test","rhs");
      //list->sublist("Linear Solver Types").sublist("AztecOO").sublist("Forward Solve").sublist("AztecOO Settings").set("Orthogonalization","Modified");

      list->set("Preconditioner Type","None"); // will use user-defined P

    Stratimikos::DefaultLinearSolverBuilder builder;

      builder.setParameterList(list);

    RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > strategy = createLinearSolveStrategy(builder);

    RCP<Thyra::LinearOpWithSolveBase<double> > solver = strategy->createOp();

    Thyra::initializePreconditionedOp<double>(*strategy,
                                               matrix,
                                               Thyra::rightPrec<double>(preconditioner),
                                               solver.ptr());


    double true_residual_initial;
    double true_residual_final;
    int    iteration_count;

    // JAW check true residual
    RCP<Thyra::VectorBase<double> > Ax = Thyra::createMember(matrix->range());
    RCP<Thyra::VectorBase<double> > r  = Thyra::createMember(matrix->range());
    {      
      Thyra::apply(*matrix, Thyra::NOTRANS,*lhs,Ax.ptr());
      Thyra::V_VmV<double>(r.ptr(),*rhs,*Ax);
      true_residual_initial = Thyra::norm(*r);
    }


    Thyra::SolveStatus<double> status = 
      solver->solve(Thyra::NOTRANS,*rhs,lhs.ptr());  // actual solve!
      iteration_count = status.extraParameters->get<int>("Iteration Count");


    // JAW true check residual
    {      
      Thyra::apply(*matrix, Thyra::NOTRANS,*lhs,Ax.ptr());
      Thyra::V_VmV<double>(r.ptr(),*rhs,*Ax);
      true_residual_final = Thyra::norm(*r);
    }

    if(partition.m_rank == 0)
    {
      FILE* fp = fopen("solver_profile.txt","a");
      fprintf(fp,"%d %.9e %.9e\n", iteration_count, true_residual_initial, true_residual_final);
      fclose(fp);
    }


		// write system out for checking

    EpetraExt::RowMatrixToMatlabFile("matrix00.dat",*epetraSystem.m_matrix[0][0]);
    EpetraExt::RowMatrixToMatlabFile("matrix01.dat",*epetraSystem.m_matrix[0][1]);
    EpetraExt::RowMatrixToMatlabFile("matrix10.dat",*epetraSystem.m_matrix[1][0]);
    EpetraExt::RowMatrixToMatlabFile("matrix11.dat",*epetraSystem.m_matrix[1][1]);

    EpetraExt::MultiVectorToMatlabFile("rhs0.dat",*epetraSystem.m_rhs[0]);
    EpetraExt::MultiVectorToMatlabFile("rhs1.dat",*epetraSystem.m_rhs[1]);

    EpetraExt::MultiVectorToMatlabFile("sol0.dat",*epetraSystem.m_solution[0]);
    EpetraExt::MultiVectorToMatlabFile("sol1.dat",*epetraSystem.m_solution[1]);
  }

  // to do: copy vector solution into geos data structures
  realT scalingFactor = CheckSolutionBlock( domain );
  PropagateSolutionBlock( scalingFactor, domain );

  // re-sync ghost nodes
  partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);

  // copy vector solution into geos data structures
  PostSyncConsistency( domain, partition );

}





/* ................ SCRATCH CODE .....................................

		// make identity and zero operators

  RCP<const Thyra::LinearOpBase<double> > Eye_11
    = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(A_11->range())); 

  RCP<const Thyra::LinearOpBase<double> > Zero_01
    = Teuchos::rcp(new Thyra::DefaultZeroLinearOp<double>(A_00->range(),A_11->domain())); 

  RCP<const Thyra::LinearOpBase<double> > Zero_10
    = Teuchos::rcp(new Thyra::DefaultZeroLinearOp<double>(A_11->range(),A_00->domain())); 

*/

/*
  EpetraExt::RowMatrixToMatlabFile("matrix00.dat",*epetraSystem.m_matrix[0][0]);
  EpetraExt::RowMatrixToMatlabFile("matrix11.dat",*epetraSystem.m_matrix[1][1]);


		// test diagonal blocks

  for(unsigned b=0; b<1; ++b)
  {
    printf("Checking singularity of (%d,%d) block ...\n",b,b);

    Epetra_LinearProblem problem( &*epetraSystem.m_matrix[b][b],
                                  &*epetraSystem.m_solution[b],
                                  &*epetraSystem.m_rhs[b] );


    epetraSystem.m_rhs[b]->PutScalar(1.0);
    epetraSystem.m_solution[b]->PutScalar(0.0);

    AztecOO solver(problem);

      solver.SetAztecOption(AZ_solver,AZ_gmres);
      solver.SetAztecOption(AZ_conv,AZ_r0);
      solver.SetAztecOption(AZ_kspace,m_numerics.m_kspace);
      solver.SetAztecOption(AZ_output,AZ_all);

    SetLinearSolverParameters( solver );
*/
/*
    ML_Epetra::MultiLevelPreconditioner* MLPrec = NULL;

    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type","ILU");
    MLList.set("ML output",0);
    MLList.set("PDE equations",1);//(b==1?1:3));

    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*epetraSystem.m_matrix[b][b], MLList);
    solver.SetPrecOperator(MLPrec);
*/
/*
    solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
    solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
    solver.SetAztecParam(AZ_ilut_fill,m_numerics.m_ilut_fill);

    solver.Iterate(100,1e-6);

    std::ofstream sol_output;
    char solNamei[100] = "solution.vec";
    sol_output.open( solNamei );
    epetraSystem.m_solution[0]->Print( sol_output );
    sol_output.close();

//    delete MLPrec;
  }
  return;
*/


/*
		// comms for ML

  ML_Comm * ml_comm;
  ML_Comm_Create(&ml_comm);

		// create ML block operator

  ML_Operator * matrix = ML_Operator_Create(ml_comm);

  ML_Operator_BlkMatInit(matrix, ml_comm, 2, 2, ML_DESTROY_EVERYTHING);

  EpetraExt::CrsMatrix_SolverMap map;

  for(unsigned i=0; i<2; ++i)
  for(unsigned j=0; j<2; ++j)
  {
    Epetra_RowMatrix * epetra_block = ML_Epetra::ModifyEpetraMatrixColMap(*epetraSystem.m_matrix[i][j],map);
    ML_Operator * ml_block = ML_Operator_Create(ml_comm);
    ML_Operator_WrapEpetraMatrix(epetra_block,ml_block);
    ML_Operator_BlkMatInsert(matrix,ml_block,i,j);
  }

  ML_Operator_BlkMatFinalize(matrix);

		// extract diagonal blocks
		// approximate S with A_11 for now

  ML_Operator * matrix_00 = ML_Operator_BlkMatExtract(matrix,0,0);
  ML_Operator * matrix_11 = ML_Operator_BlkMatExtract(matrix,1,1);

  Epetra_RowMatrix* A[2];
  A[0] = dynamic_cast<Epetra_RowMatrix*>((Epetra_CrsMatrix *) matrix_00->data);
  A[1] = dynamic_cast<Epetra_RowMatrix*>((Epetra_CrsMatrix *) matrix_11->data);

		// set up ML parameter lists

  Teuchos::ParameterList TopList;
  Teuchos::ParameterList SubList[2];

  ML_Epetra::SetDefaults("SA",TopList);
  ML_Epetra::SetDefaults("SA",SubList[0]);
  ML_Epetra::SetDefaults("SA",SubList[1]);

  TopList.set("aggregation: type", "Uncoupled");
  TopList.set("smoother: type","ILU");
  
  SubList[0].set("PDE equations",3);
  SubList[0].set("aggregation: type", "Uncoupled");
  SubList[0].set("smoother: type","ILU");

  SubList[1].set("PDE equations",1);
  SubList[1].set("aggregation: type", "Uncoupled");
  SubList[1].set("smoother: type","ILU");

		// create composite AMG preconditioner

  ML_Epetra::MultiLevelPreconditioner * ml_pre = 
      new ML_Epetra::MultiLevelPreconditioner(matrix, TopList);
      //new ML_Epetra::MultiLevelPreconditioner(matrix, TopList, A, SubList, 2);

		// convert ML matrix to Epetra Block Matrix

  ML_Epetra::RowMatrix e_matrix(matrix,epetra_comm);


  RCP<Epetra_Vector> RHS = Teuchos::rcp(new Epetra_Vector(e_matrix.OperatorRangeMap()));
  RHS->PutScalar(7.0);
  RCP<Epetra_Vector> srcRHS = Teuchos::rcp(new Epetra_Vector(*RHS));
  e_matrix.Apply(*srcRHS,*RHS);

  // set initial guess 
  Epetra_Vector LHS(e_matrix.OperatorDomainMap()); 
  LHS.PutScalar(0.0);

  Epetra_LinearProblem Problem(&e_matrix, &LHS, &*RHS);
  AztecOO solver(Problem);
  
  solver.SetPrecOperator(ml_pre);
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 1);
  solver.Iterate(100, 1e-8);

  delete ml_pre;
  ML_Operator_Destroy(&matrix);
  ML_Comm_Destroy(&ml_comm);
*/

//const Epetra_RowMatrix & rm = ml_precOp->RowMatrix();


		// ML preconditioners for A00 and A11
/*
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("SA",MLList);
    MLList.set("aggregation: type", "Uncoupled");
    MLList.set("smoother: type","ILU");
    MLList.set("ML output",0);
    MLList.set("PDE equations",3);

    const RCP<Epetra_Operator> tmp 
      = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*epetraSystem.m_matrix[0][0], MLList));

    const RCP<const Thyra::PreconditionerBase<double> > ml_A00
      = Thyra::epetraLinearOp(tmp);
*/

/*
  const RCP<Thyra::MLPreconditionerFactory> ml_A00_factory
    = Teuchos::rcp(new Thyra::MLPreconditionerFactory());

  const RCP<Thyra::MLPreconditionerFactory> ml_A11_factory
    = Teuchos::rcp(new Thyra::MLPreconditionerFactory());


  const RCP<Teuchos::ParameterList> ml_A00_list 
    = Teuchos::rcp(new Teuchos::ParameterList("ml_A00_list"),true);

  const RCP<Teuchos::ParameterList> ml_A11_list 
    = Teuchos::rcp(new Teuchos::ParameterList("ml_A11_list"),true);


  ML_Epetra::SetDefaults("SA",*ml_A00_list);
    ml_A00_list->set("aggregation: type", "Uncoupled");
    ml_A00_list->set("smoother: type","ILU");
    ml_A00_list->set("ML output",0);
    ml_A00_list->set("PDE equations",3);

  ML_Epetra::SetDefaults("SA",*ml_A11_list);
    ml_A11_list->set("aggregation: type", "Uncoupled");
    ml_A11_list->set("smoother: type","ILU");
    ml_A11_list->set("ML output",0);
    ml_A11_list->set("PDE equations",1);

 
  ml_A00_factory->setParameterList(ml_A00_list);
  ml_A11_factory->setParameterList(ml_A11_list);



  RCP<Thyra::PreconditionerBase<double> > ml_A00 
    = ml_A00_factory->createPrec();

  Thyra::initializePrec(A_00,&*ml_A00,Thyra::ESupportSolveUse()); 
 */

/*
		// old way

  const RCP<Teuchos::ParameterList> inv_A_00_list 
    = Teuchos::rcp(new Teuchos::ParameterList("inv_A_00_list"),true);


    inv_A_00_list->sublist("Forward Solve").set("Max Iterations",100);
    inv_A_00_list->sublist("Forward Solve").set("Tolerance",1e-10);
    inv_A_00_list->sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Solver", "GMRES");
    inv_A_00_list->sublist("Forward Solve").sublist("AztecOO Settings").set("Size of Krylov Subspace",100);
    inv_A_00_list->sublist("Forward Solve").sublist("AztecOO Settings").set("Output Frequency", 0);
    inv_A_00_list->sublist("Forward Solve").sublist("AztecOO Settings").set("Aztec Preconditioner", "none");

  const RCP<Thyra::LinearOpWithSolveFactoryBase<double> > inv_A_00_factory
    = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());

    inv_A_00_factory->setParameterList(inv_A_00_list);

  const RCP<const Thyra::LinearOpBase<double> > inv_A_00
    = Thyra::inverse(*inv_A_00_factory,A_00);

		// do same for A_11 block

  const RCP<const Thyra::LinearOpBase<double> > inv_A_11
    = Thyra::inverse(*inv_A_00_factory,A_11);
*/
/*
#if 0

  Amesos_BaseSolver* Solver;
  Amesos Factory;

  std::string SolverType = "Klu";

//  if( !(m_lsp.compare("Slu")) )
//  {
//    SolverType = "Amesos_Superludist";
//  }

  Solver = Factory.Create(SolverType, problem);
//  if (Solver == 0) {
//      cerr << "Specified solver is not available" << endl;
//  }

//  Amesos_Superludist* Solver = new Amesos_Superludist(problem);

  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  Teuchos::ParameterList List;
  List.set("MaxProcs", Comm.NumProc());
  Solver->SetParameters(List);

  int ierr = Solver->SymbolicFactorization();
  if (ierr > 0)
    cerr << "ERROR!" << endl;

  ierr = Solver->NumericFactorization();
  if (ierr > 0)
    cerr << "ERROR!" << endl;

  Solver->Solve();


  if( Solver!=NULL )
    delete Solver;
#else
*/

/*
#if 0

  //Amesos_BaseSolver* Solver;
  Amesos Factory;

  std::string SolverType = "Klu";

  if( !(m_lsp.compare("Slu")) )
  {
    SolverType = "Amesos_Superludist";
  }

//  Solver = Factory.Create(SolverType, problem);
//  if (Solver == 0) {
//      cerr << "Specified solver is not available" << endl;
//  }

  Amesos_Superludist* Solver = new Amesos_Superludist(problem);

  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  Teuchos::ParameterList List;
  List.set("MaxProcs", Comm.NumProc());
  Solver->SetParameters(List);

  int ierr = Solver->SymbolicFactorization();
  if (ierr > 0)
    cerr << "ERROR!" << endl;

  ierr = Solver->NumericFactorization();
  if (ierr > 0)
    cerr << "ERROR!" << endl;

  Solver->Solve();


  if( Solver!=NULL )
    delete Solver;
#else
*/
