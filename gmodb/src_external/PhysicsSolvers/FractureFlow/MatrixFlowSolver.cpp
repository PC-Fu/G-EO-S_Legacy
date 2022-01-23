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
//  LLNL-CODE-656690
//  GMOD-B, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GMOD-B. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
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
 * @file MatrixFlowSolver.cpp
 * @author hao1
 * @date May 10, 2014
 */

#include "MatrixFlowSolver.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"
#include "ElementLibrary/FiniteElementBase.h"
#include "ElementLibrary/FiniteElementUtilities.h"

// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

using namespace BoundaryConditionFunctions;
using namespace PS_STR;

MatrixFlowSolver::MatrixFlowSolver(const std::string &name, ProblemManagerT* const pm):
      SolverBase(name, pm),
      m_useMLPrecond(false),
      this_mpi_process(pm->m_epetraComm.MyPID()),
      n_mpi_processes(pm->m_epetraComm.NumProc()),
      m_epetra_comm(pm->m_epetraComm),
      row_map(),
      sparsity(),
      matrix(),
      solution(),
      rhs(),
      syncedFields(),
      m_TrilinosIndexStr(),
      m_dt(0.0)
{
  //  ++m_instances; 
  //  m_TrilinosIndexStr = "FFS_SS_" +  toString<int>(m_instances) + "_GlobalDof";
  m_TrilinosIndexStr = "MFS_SS_GlobalDof";
}

MatrixFlowSolver::~MatrixFlowSolver() {}


void MatrixFlowSolver::ReadXML( TICPP::HierarchicalDataNode* hdn )
{ 

  m_matrix_regionName = hdn->GetStringVector("matrixRegionName");

  m_matrix_regionNumber = m_matrix_regionName.size();

  if(m_matrix_regionName.empty())
    throw GPException("ERROR: 'matrixRegionName' is not set");


  std::string tempStr = hdn->GetAttributeString("timeIntegration");

  if(tempStr.empty())
    throw GPException("ERROR: 'timeIntegration' is not set");

  if( tempStr == "Implicit"){

    m_timeIntegrationOption =  Implicit;

  } else if(tempStr == "Explicit") {

    m_timeIntegrationOption = Explicit;

  } else {

    throw GPException("Invalid value input into MatrixFlowSolver::IntToTimeIntegration()");

  }

  if(m_timeIntegrationOption == Implicit) {

    // Linear solver parameters
    m_numerics.m_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-5);
    m_numerics.m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);

    m_useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);

    // Nonlinear solver parameters
    m_NRNumerics.m_tol = hdn->GetAttributeOrDefault<realT>("nrTol", 1e-5);
    m_NRNumerics.m_maxIters = hdn->GetAttributeOrDefault<realT>("nrMaxIterations", 10);

  }

  // fluid and rock properties

  m_fluidRockProperty.m_compress = hdn->GetAttributeOrDefault("fluidCompress", "5.0e-10 1/Pa");

  m_fluidRockProperty.m_mu = hdn->GetAttributeOrDefault("fluidViscosity", "0.001 Pa.s");

  m_fluidRockProperty.m_rho_o = hdn->GetAttributeOrDefault("rhoRef","1000.0 kg/m^3");

  m_fluidRockProperty.m_press_o = hdn->GetAttributeOrDefault("pressRef","1.0e5 Pa");

  m_fluidRockProperty.m_rockCompress = hdn->GetAttributeOrDefault("rockCompress", "3.5e-10 1/Pa");

  m_fieldName = hdn->GetStringVector("fieldName");

  if(m_fieldName.empty())
    throw GPException("ERROR: fieldName is not set!");

}

void MatrixFlowSolver::RegisterFields( PhysicalDomainT& domain )
{

  ElementManagerT& elementManager = domain.m_feElementManager;
  NodeManagerT& nodeManager = domain.m_feNodeManager;

  nodeManager.AddKeyedDataField<FieldInfo::pressure>();
  domain.m_feFaceManager.AddKeylessDataField<int>( "flowFaceType", true, true );
  //domain.m_feFaceManager.AddKeylessDataField<realT>("totalLeakedVolume", true, true); This should have been done in the coupled HF solver
  domain.m_feFaceManager.AddKeylessDataField<realT>("faceToMatrixLeakOffRate", true, true);
  nodeManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);
  rArray1d &totalLeakedVolume = domain.m_feFaceManager.GetFieldData<realT>("totalLeakedVolume");

  if (true) // These are needed for coupling with the porous media flow model
  {
    domain.m_feNodeManager.AddKeylessDataField<realT>("fractureNodalPressure",false,true); //We dont need it for restart because the parallel flow solver will calculate it before invoking the porous media solver.
    domain.m_feNodeManager.AddKeylessDataField<int>("isFractureNode",false,false);
  }

  totalLeakedVolume = 0.0;

  for(int i =0; i < m_matrix_regionNumber; ++i){

    std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.find(m_matrix_regionName[i]);

    if(it == elementManager.m_ElementRegions.end())
      throw GPException("Fracture region " + m_matrix_regionName[i] + " is not found!");

    ElementRegionT& elemRegion = it->second;
    m_elemRegion.push_back(&elemRegion);

    elemRegion.AddKeylessDataField<realT>(PermeabilityStr,true,true);
    elemRegion.AddKeylessDataField<realT>(PorosityStr,true,true);

  }

}

void MatrixFlowSolver::InitializeCommunications( PartitionBase& partition )
{

  syncedFields.clear();

  syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::pressure>::Name());

}


void MatrixFlowSolver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{

  //common initialization

  NodeManagerT& nodeManager = domain.m_feNodeManager;

  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  int num, index, id;
  int m, n, p, q, num1, num2;

  int node_num = nodeManager.m_numNodes;

  num = m_elemRegion.size();

  m_elem_is_ghost.resize(num);

  m_femCoefs.resize(num);

  m_volume.resize(num);

  m_density.resize(node_num);

  m_flowFaceNodeFlag.resize(node_num); 

  for(int i = 0; i < num; i++) {

    m_elem_is_ghost[i] = &(m_elemRegion[i]->GetFieldData<FieldInfo::ghostRank>());
    m_femCoefs[i].resize(m_elemRegion[i]->m_numElems);

    m = m_elemRegion[i]->m_numNodesPerElem * (m_elemRegion[i]->m_numNodesPerElem - 1) / 2;

    for(unsigned int j = 0; j < m_elemRegion[i]->m_numElems; j++) {
      m_femCoefs[i][j].resize(m);
    }

    m_volume[i].resize(m_elemRegion[i]->m_numElems);

  }


  iArray1d &trilinos_index = nodeManager.GetFieldData<int>(m_TrilinosIndexStr);

  trilinos_index = 0;

  // initialize Dirichlet boundary conditions

  SetDirichletBCs(domain);

  // calculate some geometric parameters 
  // (e.g. element volume, connection area)

  if(m_timeIntegrationOption == Explicit) {

    m_volume2.resize(node_num);

  }

  MakeParameters(domain);

  if(m_timeIntegrationOption == Implicit) {

    m_oldDensity.resize(node_num);

    CalculateInterpCoefs(domain);

    const iArray1d& node_is_ghost  = nodeManager.GetFieldData<FieldInfo::ghostRank>();

    // define trilinos index

    int n_local_rows = 0;

    for(int i = 0; i < node_num; i++) {

      if(i == nodeManager.GetParentIndex(i)) {

        if(node_is_ghost[i] < 0 && trilinos_index[i] != -1) {

          n_local_rows++;

        }
      }

    }

    std::vector<int> gather(n_mpi_processes);

    m_epetra_comm.GatherAll(&n_local_rows,
                            &gather.front(),
                            1);

    int first_local_row = 0;
    int n_global_rows = 0;

    for(int proc=0; proc<n_mpi_processes; ++proc)
    {
      n_global_rows += gather[proc];
      if(proc<this_mpi_process)
        first_local_row += gather[proc];
    }

    unsigned local_count = 0;

    for(int i = 0; i < node_num; i++) {

      if(node_is_ghost[i] < 0 && i == nodeManager.GetParentIndex(i)) {

        if(trilinos_index[i] != -1) {
          trilinos_index[i] = first_local_row + local_count;
          local_count++;

        }

      }
      else {

        trilinos_index[i] = -INT_MAX;

      }

    }

    assert(static_cast<int>(local_count) == n_local_rows);

    std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedIndexFields;
    syncedIndexFields[PhysicalDomainT::FiniteElementNodeManager].push_back(m_TrilinosIndexStr);

    partition.SetBufferSizes(syncedIndexFields, CommRegistry::matrixFlowSolver);
    partition.SynchronizeFields(syncedIndexFields, CommRegistry::matrixFlowSolver);

    // create epetra map

    row_map = Teuchos::rcp(new Epetra_Map(n_global_rows, n_local_rows,0,m_epetra_comm));

    // set up sparsity graph

    sparsity = Teuchos::rcp(new Epetra_FECrsGraph(Copy,*row_map,0));

    // loop over all non-fracture elements, make node pairs

    std::vector<int> dofIndex(2);

    Array1dT<lArray1d> node_pairs(node_num);

    localIndex lcnode0, lcnode1, tmp;
    bool not_found;

    for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

      const ElementRegionT* elemRegion = m_elemRegion[i];

      // loop over all elements in the region
      for( localIndex k=0 ; k<elemRegion->m_numElems; k++) {

        if((*m_elem_is_ghost[i])[k] < 0)
        {
          // get the elementToNodeMap for element k
          const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(k);

          // loop over all nodes in elementToNodeMap
          for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem-1 ; a++) {
            // get local index of the node from elementToNodeMap
            lcnode0 = nodeManager.GetParentIndex(elementToNodeMap[a]);

            for( unsigned int b=a+1 ; b<elemRegion->m_numNodesPerElem ; b++) {
              // get local index of the node from elementToNodeMap
              lcnode1 = nodeManager.GetParentIndex(elementToNodeMap[b]);

              if(lcnode1 < lcnode0) {
                tmp = lcnode1;
                lcnode1 = lcnode0;
                lcnode0 = tmp;
              }

              not_found = 1;

              for(unsigned int j = 0; j < node_pairs[lcnode0].size(); j++) {

                if(node_pairs[lcnode0][j] == lcnode1) {
                  not_found = 0;
                  break;
                }

              }

              if(not_found) {

                node_pairs[lcnode0].push_back(lcnode1);

                if(trilinos_index[lcnode0] != -1 && trilinos_index[lcnode1] != -1) {

                  dofIndex[0] = trilinos_index[lcnode0];
                  dofIndex[1] = trilinos_index[lcnode1];


                  sparsity->InsertGlobalIndices(dofIndex.size(),
                                                &dofIndex.front(),
                                                dofIndex.size(),
                                                &dofIndex.front());

                }
                else if(trilinos_index[lcnode0] != -1) {

                  dofIndex[0] = trilinos_index[lcnode0];

                  sparsity->InsertGlobalIndices(1,
                                                &dofIndex.front(),
                                                1,
                                                &dofIndex.front());

                }
                else if(trilinos_index[lcnode1] != -1) {

                  dofIndex[0] = trilinos_index[lcnode1];

                  sparsity->InsertGlobalIndices(1,
                                                &dofIndex.front(),
                                                1,
                                                &dofIndex.front());

                }

              }

            }

          }

        }

      }

    }

    sparsity->GlobalAssemble();
    sparsity->OptimizeStorage();

    matrix   = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*sparsity));
    solution = Teuchos::rcp(new Epetra_FEVector(*row_map));
    rhs      = Teuchos::rcp(new Epetra_FEVector(*row_map));

  }



  partition.SetBufferSizes(syncedFields, CommRegistry::matrixFlowSolver);

}

void MatrixFlowSolver::DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set, realT time){

  MultiVarDirichletBoundaryCondition* bc = dynamic_cast<MultiVarDirichletBoundaryCondition*> (bcBase);

  if(bc) {

    if(!bc->isClamped()) {

      const rArray1d &value = bc->GetValues(time);

      if(object.GetObjectType() == ObjectDataStructureBaseT::NodeManager && bc->m_objectKey == PhysicalDomainT::FiniteElementNodeManager) {

        rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

        lSet::const_iterator kf=set.begin();
        localIndex id;

        for(localIndex i =0; i < set.size() ; ++i, ++kf) {

          id = domain.m_feNodeManager.GetParentIndex(*kf);
          pressure[id] = value[0];

        }

      }

    }

  }

}


void MatrixFlowSolver::Assemble(PhysicalDomainT&  domain,
                                SpatialPartition& partition,
                                const realT& time)
{

  matrix->PutScalar(0.0);
  rhs->PutScalar(0.0);
  solution->PutScalar(0.0);

  int nvar = 1;

  Epetra_IntSerialDenseVector  elem_index (nvar);
  Epetra_SerialDenseVector     elem_rhs     (nvar);
  Epetra_SerialDenseMatrix     elem_matrix  (nvar,nvar);

  int num = 2;

  Epetra_IntSerialDenseVector  face_index(num);
  Epetra_SerialDenseVector     face_rhs(num);
  Epetra_SerialDenseMatrix     face_matrix(num, num);

  int index;
  realT pv, rho, rhoav;
  realT term0, term1, term2, term3, term4, term5;

  realT m_compress = m_fluidRockProperty.m_compress;
  realT m_mu = m_fluidRockProperty.m_mu;

  const iArray1d &trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);
  const rArray1d &pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  // loop over elements (accumulation term)

  localIndex lcnode0, lcnode1;

  for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

    const ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {

      if((*m_elem_is_ghost[i])[j] < 0) {

        const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(j);
        // loop over all nodes in elementToNodeMap

        for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem; a++) {

          // get local index of the node from elementToNodeMap
          lcnode0 = domain.m_feNodeManager.GetParentIndex(elementToNodeMap[a]);

          if(trilinos_index[lcnode0] != -1 && m_flowFaceNodeFlag[lcnode0] == 0) {

            elem_matrix(0,0) = m_compress * m_density[lcnode0] * m_volume[i][j];
            elem_index(0) = trilinos_index[lcnode0];

            elem_rhs(0) =  (m_density[lcnode0] - m_oldDensity[lcnode0]) * m_volume[i][j];

            matrix->SumIntoGlobalValues(elem_index, elem_matrix);
            rhs->SumIntoGlobalValues(elem_index, elem_rhs);

          }

        }

      }

    }

  }

  for(localIndex i = 0; i < domain.m_feNodeManager.m_numNodes; i++) {

    if(i != domain.m_feNodeManager.GetParentIndex(i))
      continue;

    if(trilinos_index[i] != -1 && m_flowFaceNodeFlag[i] == 1) {

      elem_matrix(0,0) = 1.0;
      elem_index(0) = trilinos_index[i];

      elem_rhs(0) =  0.0;

      matrix->SumIntoGlobalValues(elem_index, elem_matrix);
      rhs->SumIntoGlobalValues(elem_index, elem_rhs);

    }
  }

  //connection

  for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

    const ElementRegionT* elemRegion = m_elemRegion[i];

    // loop over all elements in the region
    for( localIndex k=0 ; k<elemRegion->m_numElems; k++) {

      if((*m_elem_is_ghost[i])[k] < 0) {

        // get the elementToNodeMap for element k
        const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(k);

        num = 0;

        // loop over all nodes in elementToNodeMap
        for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem-1 ; a++) {

          // get local index of the node from elementToNodeMap
          lcnode0 = domain.m_feNodeManager.GetParentIndex(elementToNodeMap[a]);

          for( unsigned int b=a+1 ; b<elemRegion->m_numNodesPerElem ; b++) {

            // get local index of the node from elementToNodeMap
            lcnode1 = domain.m_feNodeManager.GetParentIndex(elementToNodeMap[b]);

            rhoav = (m_density[lcnode0] + m_density[lcnode1]) * 0.5;
            term0 = pressure[lcnode0] - pressure[lcnode1];

            term1 = rhoav * (m_femCoefs[i])[k][num] / m_mu * m_dt;

            term2 = term0 * 0.5 * m_density[lcnode0] * m_compress * (m_femCoefs[i])[k][num] / m_mu * m_dt;

            term3 = term0 * 0.5 * m_density[lcnode1] * m_compress * (m_femCoefs[i])[k][num] / m_mu * m_dt;


            if((trilinos_index[lcnode0] != -1 && m_flowFaceNodeFlag[lcnode0] ==0) && (trilinos_index[lcnode1] != -1 && m_flowFaceNodeFlag[lcnode1] == 0)) {

              face_matrix(0,0) = term1 + term2;
              face_matrix(1,1) = term1 - term3;

              face_matrix(0,1) =  -term1 + term3;
              face_matrix(1,0) =  -term1 - term2;

              face_index[0] = trilinos_index[lcnode0] * nvar;
              face_index[1] = trilinos_index[lcnode1] * nvar;

              matrix->SumIntoGlobalValues(face_index, face_matrix);

              elem_index(0) = trilinos_index[lcnode0];
              elem_rhs(0) =  term1 * term0;

              rhs->SumIntoGlobalValues(elem_index, elem_rhs);

              elem_index(0) = trilinos_index[lcnode1];
              elem_rhs(0) = -term1 * term0;

              rhs->SumIntoGlobalValues(elem_index, elem_rhs);


            }
            else if(trilinos_index[lcnode0] != -1 && m_flowFaceNodeFlag[lcnode0] == 0) {

              elem_matrix(0,0) = term1 + term2;
              elem_index[0] = trilinos_index[lcnode0];

              matrix->SumIntoGlobalValues(elem_index, elem_matrix);

              elem_rhs(0) =  term1 * term0;
              rhs->SumIntoGlobalValues(elem_index, elem_rhs);

            }
            else if(trilinos_index[lcnode1] != -1 && m_flowFaceNodeFlag[lcnode1] == 0) {

              elem_matrix(0,0) = term1 - term3;
              elem_index[0] = trilinos_index[lcnode1];

              matrix->SumIntoGlobalValues(elem_index, elem_matrix);

              elem_rhs(0) =  -term1 * term0;
              rhs->SumIntoGlobalValues(elem_index, elem_rhs);

            }

            num++;

          }
        }

      }

    }

  }


  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

}


realT MatrixFlowSolver::ExplicitSolve(const realT& dt,
                                      PhysicalDomainT&  domain,
                                      SpatialPartition& partition)
{

  realT rho, rhoav;
  realT term0, term1;

  realT m_compress = m_fluidRockProperty.m_compress;
  realT m_mu = m_fluidRockProperty.m_mu;

  NodeManagerT& nodeManager = domain.m_feNodeManager;

  const iArray1d &trilinos_index = nodeManager.GetFieldData<int>(m_TrilinosIndexStr);
  const iArray1d& node_is_ghost  = nodeManager.GetFieldData<FieldInfo::ghostRank>();
  rArray1d &pressure  = nodeManager.GetFieldData<FieldInfo::pressure>();

  int node_num = nodeManager.m_numNodes;

  rArray1d temp_solution(node_num);

  temp_solution = 0.0;

  localIndex lcnode0, lcnode1;

  int num;

  //connection

  for(localIndex i = 0; i < m_elemRegion.size(); i++) {

    const ElementRegionT* elemRegion = m_elemRegion[i];

    // loop over all elements in the region
    for( localIndex k=0 ; k<elemRegion->m_numElems; k++) {

      if((*m_elem_is_ghost[i])[k] < 0) {

        // get the elementToNodeMap for element k
        const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(k);

        num = 0;

        // loop over all nodes in elementToNodeMap
        for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem-1 ; a++) {

          // get local index of the node from elementToNodeMap
          lcnode0 = nodeManager.GetParentIndex(elementToNodeMap[a]);

          for( unsigned int b=a+1 ; b<elemRegion->m_numNodesPerElem ; b++) {

            // get local index of the node from elementToNodeMap
            lcnode1 = nodeManager.GetParentIndex(elementToNodeMap[b]);

            rhoav = (m_density[lcnode0] + m_density[lcnode1]) * 0.5;

            term0 = pressure[lcnode0] - pressure[lcnode1];

            term1 = rhoav * (m_femCoefs[i])[k][num] / m_mu;

            if(trilinos_index[lcnode0] != -1 && m_flowFaceNodeFlag[lcnode0] == 0 && trilinos_index[lcnode1] != -1 &&  m_flowFaceNodeFlag[lcnode1] == 0) {
              temp_solution[lcnode0] +=  - term1 * term0;
              temp_solution[lcnode1] +=  term1 * term0;

            }
            else if(trilinos_index[lcnode0] != -1 && m_flowFaceNodeFlag[lcnode0] == 0) {

              temp_solution[lcnode0] +=  - term1 * term0;

            }
            else if(trilinos_index[lcnode1] != -1 && m_flowFaceNodeFlag[lcnode1] == 1) {

              temp_solution[lcnode1] +=  term1 * term0;

            }

            num++;

          }
        }

      }

    }

  }

  realT out = -1e20;
  localIndex id;

  for(int i = 0; i < node_num; i++) {

    if(i != nodeManager.GetParentIndex(i))
      continue;

    if(node_is_ghost[i] < 0 && trilinos_index[i] >= 0 && m_flowFaceNodeFlag[i] == 0) {

      rho = temp_solution[i] / m_volume2[i];

      if(fabs(rho) > out)
        out = fabs(rho);

      rho = m_density[i] + rho * dt;

      pressure[i] = log(rho / m_fluidRockProperty.m_rho_o) / m_compress + m_fluidRockProperty.m_press_o;

      //Fu note:
      // Many sources define the fluid compressibility as   dP = dV/V / Cf;
      // However, the equation above yields practically the same result as long as density is smaller than 1.1.

    }

  }

  // re-sync ghost nodes
  partition.SynchronizeFields(syncedFields, CommRegistry::fractureFlowSolver);

  return out;

}


realT MatrixFlowSolver::ImplicitSolve(PhysicalDomainT&  domain,
                                      SpatialPartition& partition)
{

  realT out, tmp;

  int nvar = 1;

  int index;
  double* local_solution = NULL;

  rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();
  const iArray1d& node_is_ghost  = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();
  const iArray1d &trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);

  solution->ExtractView(&local_solution,&index);


  Epetra_LinearProblem problem(&(*matrix),
                               &(*solution),
                               &(*rhs));

  Teuchos::ParameterList MLList;

  ML_Epetra::SetDefaults("SA",MLList);

  ML_Epetra::MultiLevelPreconditioner* MLPrec = NULL;
  if(m_useMLPrecond){
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*matrix, MLList);
  }


  // Later we need to test performance of different solvers/preconditioners 
  // for fracture flow problems 

  AztecOO solver(problem);

  if(m_useMLPrecond) solver.SetPrecOperator(MLPrec);

  solver.SetAztecOption(AZ_solver,AZ_bicgstab);
  solver.SetAztecOption(AZ_precond,AZ_dom_decomp);
  solver.SetAztecOption(AZ_subdomain_solve,AZ_ilu);
  solver.SetAztecOption(AZ_conv,AZ_rhs);

  solver.Iterate(m_numerics.m_maxIters,m_numerics.m_tol);

  if(m_useMLPrecond){
    delete MLPrec;
  }

  solution->MinValue(&out);
  if(out < 0.0) out = -out;

  solution->MaxValue(&tmp);
  if(tmp < 0.0) tmp = -tmp;

  if(out < tmp)
    out = tmp;

  long long int rowTmp;
  int lid;

  for(unsigned int i = 0; i < domain.m_feNodeManager.m_numNodes; i++) {

    if(i != domain.m_feNodeManager.GetParentIndex(i))
      continue;

    if(node_is_ghost[i] < 0 && trilinos_index[i] >= 0) {

      rowTmp = static_cast<long long int>(trilinos_index[i]);
      lid = row_map->LID(rowTmp);
      pressure[i] -= local_solution[lid];

    }

  }

  // re-sync ghost nodes
  partition.SynchronizeFields(syncedFields, CommRegistry::matrixFlowSolver);

  return out;

}

double MatrixFlowSolver::TimeStep(const realT& time,
                                const realT& dt,
                                const int cycleNumber,
                                PhysicalDomainT& domain,
                                const sArray1d& namesOfSolverRegions,
                                SpatialPartition& partition,
                                FractunatorBase* const fractunator)
{

  m_dt = dt;

  UpdateDirichletBCs(domain, time);

  UpdateFlowFaceBCs(domain, time);

  UpdateFluidRockProps(domain);

  if(m_timeIntegrationOption == Implicit) {

    MakeOldAcc(domain);

    bool converge = 0;
    realT out;

    for(int i = 0; i < m_NRNumerics.m_maxIters; i++) {

      Assemble(domain,partition, time);

      out = ImplicitSolve(domain,partition);

      //for debugging purpose
      printf("Current time=%g NR iter=%d Max residual=%g\n",time, i, out);

      if(out < m_NRNumerics.m_tol) {
        converge = 1;
        break;
      }

      UpdateFluidRockProps(domain);

    }

    if(!converge) {

      printf("Max. Residual = %e\n", out);
      throw GPException("Error: hit max. NR iterations -- need timestep cutback (not implemented: ");
    }

  }
  else if(m_timeIntegrationOption == Explicit) {

    realT resid = ExplicitSolve(dt, domain, partition);

    // for debugging  purpose
    printf("time=%g cycle=%d dt=%g resid=%g\n",time, cycleNumber, dt, resid);

  } 

  return dt;
}

void MatrixFlowSolver::MakeParameters(PhysicalDomainT& domain) {

  R1Tensor la, lb;
  localIndex index;
  realT w;
  int kf, eg;
  localIndex num;

  NodeManagerT& nodeManager = domain.m_feNodeManager;

  //loop over all local pairs of all the elements

  for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];

    const rArray1d &permeability = elemRegion->GetFieldData<realT>(PermeabilityStr);
    const rArray1d &porosity = elemRegion->GetFieldData<realT>(PorosityStr);

    for(localIndex j=0; j < elemRegion->m_numElems; j++) {

      if((*m_elem_is_ghost[i])[j] < 0) {

        m_volume[i][j] = 0.0;

        for(unsigned int a=0; a < elemRegion->m_numIntegrationPointsPerElem; a++) {

          m_volume[i][j] += elemRegion->m_detJ(j,a);

        }

        m_volume[i][j] /= realT(elemRegion->m_numNodesPerElem);

        num = 0;

        for(unsigned int m=0; m < elemRegion->m_numNodesPerElem-1; m++) {

          for(unsigned int n=m+1; n < elemRegion->m_numNodesPerElem; n++) {

            (m_femCoefs[i])[j][num] = 0.0;

            for(unsigned int a=0; a < elemRegion->m_numIntegrationPointsPerElem; a++) {

              w = 0.0;

              for(unsigned int dir = 0; dir < elemRegion->m_ElementDimension; dir++) {

                w += elemRegion->m_dNdX(j)(a,m)[dir] * elemRegion->m_dNdX(j)(a,n)[dir];
              }

              (m_femCoefs[i])[j][num] -= w * elemRegion->m_detJ(j,a);

            }

            num++;
          }

        }


        m_volume[i][j] *= porosity[j];

        for(unsigned int k = 0; k < (m_femCoefs[i])[j].size(); k++) {

          (m_femCoefs[i])[j][k] *= permeability[j];

        }

      }

    }

  }

  if(m_timeIntegrationOption == Explicit) {

    m_volume2 = 0.0;

    for(localIndex i = 0; i < m_elemRegion.size(); i++) {

      const ElementRegionT* elemRegion = m_elemRegion[i];

      for(localIndex j=0; j< elemRegion->m_numElems; j++) {

        const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(j);
        // loop over all nodes in elementToNodeMap

        for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem; a++) {

          // get local index of the node from elementToNodeMap
          index = domain.m_feNodeManager.GetParentIndex(elementToNodeMap[a]);

          m_volume2[index] += m_volume[i][j];

        }

      }

    }

    //deallocate memory
    m_volume.clear();

  }

}

void MatrixFlowSolver::CalculateInterpCoefs(PhysicalDomainT& domain)
{

  // Are always parent faces put before child faces?

  FaceManagerT& faceManager = domain.m_feFaceManager;
  Array1dT<R1Tensor>& faceCenter = faceManager.GetFieldData<R1Tensor>("FaceCenter");
  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  R1Tensor la, elemCenter, faceNormal;
  localIndex index, childIndex;

  m_interpCoefs.resize(faceManager.m_numFaces);

  rArrayPtrs permeabilities(m_elemRegion.size());

  for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

    permeabilities[m_elemRegion[i]->m_regionNumber] = &(m_elemRegion[i]->GetFieldData<realT>(PermeabilityStr));

  }

  realT m_mu = m_fluidRockProperty.m_mu;

  for(localIndex i = 0; i < faceManager.m_numFaces; i++) 
  {

    index = faceManager.GetParentIndex(i);

    if(i != index) continue;

    for(localIndex j = 0; j < faceManager.m_toElementsRelation[i].size(); j++) 
    {

      const std::pair< ElementRegionT*, localIndex >& elem = faceManager.m_toElementsRelation[i][j];

      if (std::find(m_matrix_regionName.begin(), m_matrix_regionName.end(),  elem.first -> m_regionName) != m_matrix_regionName.end())
        //Fu: We have to consider the situation where the element attached to this face is not ina porous medium flow region.
      {
        m_interpCoefs[i].push_back(InterpCoefs(elem.first, elem.second));
      }

    }

    if(flowFaceType[i] == 0) {

      if(faceManager.m_childIndices[i].size() > 0) {

        for(localIndex k = 0; k < faceManager.m_childIndices[i].size(); k++) {

          // This can handle both 2D and 3D parent/child faces
          if(faceManager.m_childIndices[i][k] != i) {

            childIndex = faceManager.m_childIndices[i][k];
            break;

          }

        }

        const std::pair< ElementRegionT*, localIndex >& elem = faceManager.m_toElementsRelation[childIndex][0];
        if (std::find(m_matrix_regionName.begin(), m_matrix_regionName.end(),  elem.first -> m_regionName) != m_matrix_regionName.end())
        {
          m_interpCoefs[i].push_back(InterpCoefs(elem.first, elem.second));
        }
      }

    }

  }

  for(localIndex i = 0; i < faceManager.m_numFaces; i++) 
  {

    if(m_interpCoefs[i].size() > 0) 
    {

      faceNormal = faceManager.FaceNormal(domain.m_feNodeManager, i);

      for(localIndex j = 0; j < m_interpCoefs[i].size(); j++)
      {

        const ElementRegionT *elemRegion = m_interpCoefs[i][j].m_elemRegion;
        const localIndex &elemIndex = m_interpCoefs[i][j].m_elemIndex;

        la = elemRegion->GetElementCenter(elemIndex, domain.m_feNodeManager) - faceCenter[i];

        if(Dot(faceNormal, la) < 0)
          faceNormal *= -1.0;

        unsigned int ndim = elemRegion->m_finiteElement->Dim();
        unsigned int ndofs = elemRegion->m_numNodesPerElem;

        Basis *basis = elemRegion->m_elementBasis;

        Array1dT<R1Tensor> nodeCoords(ndofs);

        const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(elemIndex);

        // loop over all nodes in elementToNodeMap
        for(unsigned int a = 0 ; a < elemRegion->m_numNodesPerElem; a++) {

          localIndex nodeIndex = domain.m_feNodeManager.GetParentIndex(elementToNodeMap[a]);

          nodeCoords[a] = (*domain.m_feNodeManager.m_refposition)[nodeIndex];

        }

        Array1dT<R1Tensor> dNdX;

        if(!(elemRegion->m_elementType).compare("linear")) {

          dNdX.resize(elemRegion->m_numNodesPerElem);

          for(unsigned int a = 0 ; a < elemRegion->m_numNodesPerElem; a++) {

            for(unsigned int dir = 0; dir < ndim; dir++) {

              dNdX[a][dir] = elemRegion->m_dNdX(elemIndex)(0,a)[dir];

            }

          }
        }
        else {

          FiniteElementUtilities::InterpdNdX(faceCenter[i], nodeCoords, basis, ndofs, ndim, dNdX);

        }

        for(unsigned int a = 0 ; a < elemRegion->m_numNodesPerElem; a++) {

          (m_interpCoefs[i][j].m_coefs)[a] = 0.0;

          for(unsigned int dir = 0; dir < ndim; dir++) {

            (m_interpCoefs[i][j].m_coefs)[a] += dNdX[a][dir] * faceNormal[dir];

          }

          (m_interpCoefs[i][j].m_coefs)[a] *= -(*permeabilities[elemRegion->m_regionNumber])[elemIndex] / m_mu;

        }

      }

    }

  }

}

void MatrixFlowSolver::UpdateFluidRockProps(PhysicalDomainT& domain) {

  const rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  for(unsigned int i=0; i<domain.m_feNodeManager.m_numNodes; i++) {
    if(i != domain.m_feNodeManager.GetParentIndex(i))
      continue;

    m_density[i] = m_fluidRockProperty.m_rho_o * exp(m_fluidRockProperty.m_compress * (pressure[i] - m_fluidRockProperty.m_press_o));


  }

}

void MatrixFlowSolver::SetDirichletBCs(PhysicalDomainT& domain)
{

  ApplyBoundaryCondition<realT>(this, &MatrixFlowSolver::SetDirichletNodeBoundaryCondition, domain, domain.m_feNodeManager, "MultiVars", 0.0);

}

void MatrixFlowSolver::SetDirichletNodeBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time) {

  if(object.GetObjectType() == ObjectDataStructureBaseT::NodeManager && bcBase->m_objectKey == PhysicalDomainT::FiniteElementNodeManager) {

    MultiVarDirichletBoundaryCondition* bc = dynamic_cast<MultiVarDirichletBoundaryCondition*> (bcBase);

    if(bc) {

      bc->CheckVars(m_fieldName);

      iArray1d &trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);

      lSet::const_iterator kf=set.begin();
      localIndex id;

      for(localIndex i =0; i < set.size() ; ++i, ++kf) {

        id = domain.m_feNodeManager.GetParentIndex(*kf);
        trilinos_index[id] = -1;

      }

    }

  }

}

void MatrixFlowSolver::UpdateDirichletBCs(PhysicalDomainT& domain, const realT &time)
{

  ApplyBoundaryCondition<realT>(this, &MatrixFlowSolver::DirichletBoundaryCondition, domain, domain.m_feNodeManager, "MultiVars", 0.0);

}

void MatrixFlowSolver::UpdateFlowFaceBCs(PhysicalDomainT& domain, const realT &time)
{

  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  const rArray1d* fractureNodalPressure = domain.m_feNodeManager.GetFieldDataPointer<realT>("fractureNodalPressure");

  rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  m_flowFaceNodeFlag = 0; 

  for(localIndex i = 0; i < flowFaceType.size(); i++) 
  {
    if(flowFaceType[i] == 0) {

      for(localIndex j = 0; j < domain.m_feFaceManager.m_toNodesRelation[i].size(); j++)
      {
        localIndex id = domain.m_feNodeManager.GetParentIndex(domain.m_feFaceManager.m_toNodesRelation[i][j]);

        m_flowFaceNodeFlag[id] = 1;


      }

    }

  }

  if(fractureNodalPressure) {

    for(localIndex i = 0; i < flowFaceType.size(); i++)
    {
      if(flowFaceType[i] == 0) {

        for(localIndex j = 0; j < domain.m_feFaceManager.m_toNodesRelation[i].size(); j++) {

          localIndex id = domain.m_feNodeManager.GetParentIndex(domain.m_feFaceManager.m_toNodesRelation[i][j]);

          pressure[id] = (*fractureNodalPressure)[id];

        }

      }

    }

  }

}

void MatrixFlowSolver::MakeOldAcc(PhysicalDomainT& domain) {

  for(unsigned int i=0; i< domain.m_feNodeManager.m_numNodes; i++) {

    if(i != domain.m_feNodeManager.GetParentIndex(i))
      continue;

    m_oldDensity[i] = m_density[i];

  }

}

void MatrixFlowSolver::CalculateLeakoff(const realT &dt, PhysicalDomainT& domain)
{

  rArray1d &totalLeakedVolume = domain.m_feFaceManager.GetFieldData<realT>("totalLeakedVolume");

  const iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  const rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  const iArray1d& face_is_ghost  = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();

  rArray1d& faceToMatrixLeakOffRate = domain.m_feFaceManager.GetFieldData<realT>("faceToMatrixLeakOffRate");


  for( localIndex kf=0 ; kf<domain.m_feFaceManager.DataLengths() ; ++kf)
  {

    if(face_is_ghost[kf] < 0) {

      if (flowFaceType[kf] == 0)
      {
        localIndex face0, face1;
        if (domain.m_feFaceManager.m_toEdgesRelation[kf].size() == 2 ) //2D
        {
          face0 = domain.m_feFaceManager.m_childIndices[kf][0];
          face1 = domain.m_feFaceManager.m_childIndices[kf][1];
        }
        else
        {
          face0 = kf;
          face1 = domain.m_feFaceManager.m_childIndices[kf][0];
        }

        localIndex pf = domain.m_feFaceManager.GetParentIndex(kf);

        faceToMatrixLeakOffRate[kf] = 0.0;
        for(localIndex i = 0; i < m_interpCoefs[pf].size(); i++) {

          realT leakOffFlux = 0.0;
          const localIndex elemIndex = m_interpCoefs[pf][i].m_elemIndex;
          const ElementRegionT *elemRegion = m_interpCoefs[pf][i].m_elemRegion;

          const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(elemIndex);

          // loop over all nodes in elementToNodeMap
          for(unsigned int a = 0 ; a < elemRegion->m_numNodesPerElem; a++) {

            localIndex nodeIndex = domain.m_feNodeManager.GetParentIndex(elementToNodeMap[a]);

            leakOffFlux += pressure[nodeIndex] * (m_interpCoefs[pf][i].m_coefs)[a];
          }
          if (i==0)
          {
            totalLeakedVolume[face0] += leakOffFlux * dt;
          }
          else
          {
            totalLeakedVolume[face1] += leakOffFlux * dt;
          }
          // dt here is actually the matrix flow update interval

          faceToMatrixLeakOffRate[kf] += leakOffFlux;
        }
      }
    }
  }

}

void MatrixFlowSolver::PostProcess(PhysicalDomainT & domain,
                                   SpatialPartition& partition,
                                   const sArray1d& namesOfSolverRegions)
{
  rArray1d& faceFluidPressure = domain.m_feFaceManager.GetFieldData<FieldInfo::pressure>();
  rArray1d& nodalPressure = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();
  for (localIndex i=0; i<domain.m_feNodeManager.DataLengths(); ++i)
  {
    nodalPressure[i] = nodalPressure[domain.m_feNodeManager.GetParentIndex(i)];
  }

}

/// Register solver in the solver factory
REGISTER_SOLVER( MatrixFlowSolver )

