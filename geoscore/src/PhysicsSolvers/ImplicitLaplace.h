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
#ifndef IMPLICIT_LAPLACE_H
#define IMPLICIT_LAPLACE_H

/**
 * @file ImplicitLaplace.h
 * @author white230
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
#include "PhysicsSolverStrings.h"

#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_LinearProblem.h"

#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_RCP.hpp"
#include "AztecOO.h"
#include "ml_MultiLevelPreconditioner.h"

/**
 * Physics solver for Laplace's equation 
 */

template <int dim>
class ImplicitLaplaceSolver : public SolverBase
{
  public:
    ImplicitLaplaceSolver( const std::string& name,
                           ProblemManagerT* const pm );

    ~ImplicitLaplaceSolver();

    static const char* SolverName();

    virtual void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

    virtual void RegisterFields(PhysicalDomainT &domain);

    virtual void InitializeCommunications(PartitionBase& partition);

    virtual double TimeStep(const realT& time,
                          const realT& dt,
                          const int cycleNumber,
                          PhysicalDomainT& domain,
                          const sArray1d& namesOfSolverRegions,
                          SpatialPartition& partition,
                          FractunatorBase* const fractunator);

  private:

    void ReadXML(TICPP::HierarchicalDataNode* const hdn);

    void SetupSystem (PhysicalDomainT& domain, SpatialPartition& partition);
    void Assemble    (PhysicalDomainT& domain, SpatialPartition& partition);
    void Solve       (PhysicalDomainT& domain, SpatialPartition& partition);

    #if GPAC_MPI
      const Epetra_MpiComm & epetra_comm;
    #else
      const Epetra_SerialComm & epetra_comm;
    #endif

    const int this_mpi_process;
    const int n_mpi_processes;
    const bool     verbose;

    struct
    {
      double diffusion;
      double source;
    }
    equation_data;

    struct
    {
      double krylov_tol;
    }
    numerics;

    typedef std::map<ElementManagerT::RegKeyType, ElementRegionT > RegionMap;

    Epetra_Map*         row_map;
    Epetra_FECrsGraph*  sparsity;
    Epetra_FECrsMatrix* matrix;
    Epetra_FEVector*    solution;
    Epetra_FEVector*    rhs;

    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;
};
 

#endif

