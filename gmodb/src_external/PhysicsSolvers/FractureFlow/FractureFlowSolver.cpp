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
 * @file FractureFlowSolver.cpp
 * @author hao1
 * @date Oct. 21, 2013
 */

#include "FractureFlowSolver.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"
#include "ElementLibrary/FiniteElementBase.h"
#include "ElementLibrary/FiniteElementUtilities.h"
#include "ObjectManagers/UnitManager.h"

// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

using namespace BoundaryConditionFunctions;
using namespace PS_STR;

FractureFlowSolver::FractureFlowSolver(const std::string &name, ProblemManagerT* const pm):
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

  ++m_instances; 
  m_TrilinosIndexStr = "FFS_SS_" +  toString<int>(m_instances) + "_GlobalDof";
}

FractureFlowSolver::~FractureFlowSolver() {}

/*

      <FractureFlowSolver name="fracMatrixFlowSolver"
                          coupleMatrixFractureFlow="1"
                          fractureRegionName="EB200" 
                          updateAperture="0" 
                          fieldName="Pressure"
                          rockCompress="3.5e-10 1/Pa"
                          downVector="0.0 0.0 0.0"
      />


*/

void FractureFlowSolver::ReadXML( TICPP::HierarchicalDataNode* hdn )
{ 

  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("updateAperture",false);
  m_doRockPropUpdate = hdn->GetAttributeOrDefault<bool>("updateRockProps",false);
  m_doMFCoupling = hdn->GetAttributeOrDefault<bool>("coupleMatrixFractureFlow",false);

  // Linear solver parameters
  m_numerics.m_tol = hdn->GetAttributeOrDefault<realT>("tol",1e-5);
  m_numerics.m_maxIters = hdn->GetAttributeOrDefault<int>("maxSolverIterations",1000);
  
  m_useMLPrecond = hdn->GetAttributeOrDefault<bool>("useMLPreconditioner",false);

  // Nonlinear solver parameters
  m_NRNumerics.m_tol = hdn->GetAttributeOrDefault<realT>("nrTol", 1e-5);
  m_NRNumerics.m_maxIters = hdn->GetAttributeOrDefault<realT>("nrMaxIterations", 10);

  // equations
  m_eqt.ncomp = 1;
  m_eqt.thermal = hdn->GetAttributeOrDefault<bool>("thermal",false);
  m_eqt.nvar = m_eqt.thermal ? m_eqt.ncomp + 1 : m_eqt.ncomp;


  // fluid and rock properties

  m_fluidRockProperty.m_compress = hdn->GetAttributeOrDefault("fluidCompress", "5.0e-10 1/Pa");

  m_fluidRockProperty.m_mu = hdn->GetAttributeOrDefault("fluidViscosity", "0.001 Pa.s");

  m_fluidRockProperty.m_rho_o = hdn->GetAttributeOrDefault("rhoRef","1000.0 kg/m^3");

  m_fluidRockProperty.m_press_o = hdn->GetAttributeOrDefault("pressRef","1.0e5 Pa");

  //  m_fluidRockProperty.m_rockCompress = hdn->GetAttributeOrDefault("rockCompress", "3.5e-10 1/Pa");

  m_fluidRockProperty.m_minAperture = hdn->GetAttributeOrDefault(MinimumApertureStr,"1 um");

  m_fluidRockProperty.m_maxAperture = hdn->GetAttributeOrDefault(MaximumApertureStr,"4 mm");

  m_fluidRockProperty.m_maxPorosity = hdn->GetAttributeOrDefault("maxrockPorosity","0.999");

  m_fluidRockProperty.m_minPorosity = hdn->GetAttributeOrDefault("minRockPorosity","0.001");

  m_gravFactor = hdn->GetAttributeOrDefault("gravFactor", "9.81 m^2/s");

  m_downVector = hdn->GetAttributeTensorOrDefault("downVector", R1Tensor(0.0, 0.0, -1.0));

  m_fracture_regionName = hdn->GetStringVector("fractureRegionName");
  m_fracture_regionNumber = m_fracture_regionName.size();

  if(m_fracture_regionName.empty())
      throw GPException("ERROR: 'fractureRegionName' is not set");

  m_fieldName = hdn->GetStringVector("fieldName");

  if(m_fieldName.empty())
      throw GPException("ERROR: fieldName is not set!");

  m_matrix_regionName = hdn->GetStringVector("matrixRegionName");
  m_matrix_regionNumber = m_matrix_regionName.size();

  if(m_matrix_regionName.empty())
      throw GPException("ERROR: 'matrixRegionName' is not set");

  if(m_doRockPropUpdate) {

    m_fluidRockProperty.m_rockCompress = hdn->GetAttributeVector<realT>("rockCompress");

    if((m_fluidRockProperty.m_rockCompress).size()!= m_matrix_regionNumber)
      throw GPException("ERROR: rockCompress properties are not set correctly!");

    m_fluidRockProperty.m_rockPermCoef = hdn->GetAttributeVector<realT>("rockPermCoef");

    if((m_fluidRockProperty.m_rockPermCoef).size()!= m_matrix_regionNumber)
      throw GPException("ERROR: rockPermCoef properties are not set correctly!");

  }

  if(m_eqt.thermal) {

    //    m_fluidRockProperty.m_fluidCond = hdn->GetAttributeOrDefault("fluidConductivity", "0.5 W/m/K");
    m_fluidRockProperty.m_fluidHeatCap = hdn->GetAttributeOrDefault("fluidHeatCap", "1000.0 J/kg/K");

    //m_fluidRockProperty.m_rockDensity = hdn->GetAttributeOrDefault("rockDensity", "2500.0 kg/m^3");
   	//m_fluidRockProperty.m_rockCond = hdn->GetAttributeOrDefault("rockConductivity", "2.0 W/m/K");
    //m_fluidRockProperty.m_rockHeatCap = hdn->GetAttributeOrDefault("rockHeatCap", "1500.0 J/kg/K");
    UnitManager& theUnitManager = UnitManager::Instance();
    realT defaultDensity = theUnitManager.Convert("2500.0 kg/m^3");
    realT defaultConductivity = theUnitManager.Convert("2.0 W/m/K");
    realT defaultHeatCap = theUnitManager.Convert("1500.0 J/kg/K");

    m_fluidRockProperty.m_rockDensity = hdn->GetAttributeVectorOrDefault<realT>("rockDensity", std::vector<realT>(1,defaultDensity));

    if((m_fluidRockProperty.m_rockDensity).size()!= m_matrix_regionNumber)
      throw GPException("ERROR: rockDensity properties are not set correctly!");

    m_fluidRockProperty.m_rockCond = hdn->GetAttributeVectorOrDefault<realT>("rockConductivity", std::vector<realT>(1,defaultConductivity));

    if((m_fluidRockProperty.m_rockCond).size()!= m_matrix_regionNumber)
      throw GPException("ERROR: rockCond properties are not set correctly!");

    m_fluidRockProperty.m_rockHeatCap = hdn->GetAttributeVectorOrDefault<realT>("rockHeatCap", std::vector<realT>(1,defaultHeatCap));

    if((m_fluidRockProperty.m_rockHeatCap).size()!= m_matrix_regionNumber)
      throw GPException("ERROR: rockHeatCap properties are not set correctly!");

    m_thermalProductionBCSetName = hdn->GetStringVector("thermalProductionBCSetName");

  }

  m_dtMin = hdn->GetAttributeOrDefault<realT>("dtMin",10.0);
  m_dtMax = hdn->GetAttributeOrDefault<realT>("dtMax",5.0e6);


}

void FractureFlowSolver::RegisterFields( PhysicalDomainT& domain )
{

  ElementManagerT& elementManager = domain.m_feElementManager;
  NodeManagerT& nodeManager = domain.m_feNodeManager;

  nodeManager.AddKeyedDataField<FieldInfo::pressure>();
  domain.m_feFaceManager.AddKeylessDataField<int>( "flowFaceType", true, true );
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  if(m_eqt.thermal)
    nodeManager.AddKeylessDataField<realT>("Temperature",true,true);

  nodeManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);

  localIndex num;

  num = (elementManager.m_ElementRegions).size();

  m_elemRegionIndex.resize(num);
  m_elemRegionIndex = -1;

  for(localIndex i =0; i < m_fracture_regionNumber; ++i){

    std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.find(m_fracture_regionName[i]);

    if(it == elementManager.m_ElementRegions.end())
      throw GPException("Fracture region " + m_fracture_regionName[i] + " is not found!");

    ElementRegionT& elemRegion = it->second;
    m_elemRegionIndex[elemRegion.m_regionNumber] = i;

    m_elemRegion.push_back(&elemRegion);

    elemRegion.AddKeylessDataField<realT>(ApertureStr,true,true);

  }

  if(m_doMFCoupling) {

	  localIndex index;

    for(localIndex i =0; i < m_matrix_regionNumber; ++i) {

      std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.find(m_matrix_regionName[i]);

      if(it == elementManager.m_ElementRegions.end())
	throw GPException("Matrix region " + m_matrix_regionName[i] + " is not found!");

      ElementRegionT& elemRegion = it->second;

      if(m_elemRegionIndex[elemRegion.m_regionNumber] >=0)
	throw GPException("Matrix region " + m_matrix_regionName[i] + " is already defined!");

	index = i + m_fracture_regionNumber;
	m_elemRegionIndex[elemRegion.m_regionNumber] = index;
	m_elemRegion.push_back(&elemRegion);

	elemRegion.AddKeylessDataField<realT>(PermeabilityStr,true,true);
	elemRegion.AddKeylessDataField<realT>(PorosityStr,true,true);

	elemRegion.AddKeylessDataField<realT>("InitPermeability",true,true);
	elemRegion.AddKeylessDataField<realT>("InitPorosity",true,true);
	elemRegion.AddKeylessDataField<realT>("InitPressure",true,true);

	elemRegion.AddKeylessDataField<realT>("SolidVolume",true,true);

    }
    
  }

}

void FractureFlowSolver::InitializeCommunications( PartitionBase& partition )
{

    syncedFields.clear();

    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::pressure>::Name());

    if(m_eqt.thermal)
      syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back("Temperature");

}


void FractureFlowSolver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{

  NodeManagerT& nodeManager = domain.m_feNodeManager;
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");


  int num, index;
  int m;

  int nvar = m_eqt.nvar;
		
  int node_num = nodeManager.m_numNodes;
	    
  num = m_elemRegion.size();

  m_elem_is_ghost.resize(num);
  m_accumIndex.resize(num);
  m_volume.resize(num);


  m_femCoefs.resize(num);

  m_density.resize(node_num);
  m_oldDensity.resize(node_num);
  m_elev.resize(node_num);

  index = 0;

  for(int i = 0; i < num; i++) {

    m_elem_is_ghost[i] = &(m_elemRegion[i]->GetFieldData<FieldInfo::ghostRank>());
    m_volume[i].resize(m_elemRegion[i]->m_numElems);

    m_femCoefs[i].resize(m_elemRegion[i]->m_numElems);

    m = m_elemRegion[i]->m_numNodesPerElem * (m_elemRegion[i]->m_numNodesPerElem - 1) / 2;

    for(unsigned int j = 0; j < m_elemRegion[i]->m_numElems; j++) {
      m_femCoefs[i][j].resize(m);
    }

    m_accumIndex[i] = index;
    index += m_elemRegion[i]->m_numElems;

    for(localIndex j = 0; j < m_elemRegion[i]->m_numElems; j++) {
      m_elemIndex.push_back(std::make_pair(i, j));
    }
  }

  if(m_eqt.thermal) {

    m_femCoefs2.resize(num);
    m_oldTemperature.resize(node_num);
    m_rockHeatCap.resize(num);

    for(int i = 0; i < num; i++) {

      m_rockHeatCap[i].resize(m_elemRegion[i]->m_numElems);
      m_femCoefs2[i].resize(m_elemRegion[i]->m_numElems);
      for(unsigned int j = 0; j < m_elemRegion[i]->m_numElems; j++) {
	m_femCoefs2[i][j].resize(m);
      }

    }

  }

  iArray1d &trilinos_index = nodeManager.GetFieldData<int>(m_TrilinosIndexStr);
  const iArray1d& node_is_ghost  = nodeManager.GetFieldData<FieldInfo::ghostRank>();

  trilinos_index = 0;

  // initialize Dirichlet boundary conditions

  SetDirichletBCs(domain);

  // calculate some geometric parameters 
  // (e.g. element volume, connection area)

  MakeFractureParameters(domain);

  MakeMatrixParameters(domain);

  // define trilinos index

  // we need to examine the effect of element order on trilinos linear solver performance later !!

  // local rows

  int n_local_rows = 0;

  for(int i = 0; i < node_num; i++) {

    if(node_is_ghost[i] < 0 && trilinos_index[i] != -1) {

      n_local_rows++;

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

    if(node_is_ghost[i] < 0) {

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

  partition.SetBufferSizes(syncedIndexFields, CommRegistry::fractureFlowSolver);
  partition.SynchronizeFields(syncedIndexFields, CommRegistry::fractureFlowSolver);

  SetSrcFluxBCs(domain);

  // create epetra map

  row_map = Teuchos::rcp(new Epetra_Map(n_global_rows * nvar, n_local_rows * nvar,0,m_epetra_comm));

  // set up sparsity graph

  sparsity = Teuchos::rcp(new Epetra_FECrsGraph(Copy,*row_map,0));
	
  // loop over all non-fracture elements, make node pairs

  std::vector<int> dofIndex(2*nvar);

  std::vector<int> dofIndex2(nvar);

  Array1dT<lArray1d> node_pairs(node_num);

  localIndex lcnode0, lcnode1, tmp;
  bool not_found;

  for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

    const ElementRegionT* elemRegion = m_elemRegion[i];

    // loop over all elements in the region
    for( localIndex k=0 ; k<elemRegion->m_numElems; k++) {

      /*
      bool isActiveElement = false;
      if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
      {
        isActiveElement = true;
      }
      else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
      {
        if (flowFaceType[elemRegion->m_toFacesRelation[k][0]] == 0) isActiveElement = true;
      }

      if((*m_elem_is_ghost[i])[k] < 0 && isActiveElement )

      */

      if(i < m_fracture_regionNumber && flowFaceType[elemRegion->m_toFacesRelation[k][0]] != 0)
	continue;

      if((*m_elem_is_ghost[i])[k] < 0)
      {

	// get the elementToNodeMap for element k
	const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(k);

	// loop over all nodes in elementToNodeMap
	for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem-1 ; a++) {
	  // get local index of the node from elementToNodeMap
	  lcnode0 = elementToNodeMap[a];

	  for( unsigned int b=a+1 ; b<elemRegion->m_numNodesPerElem ; b++) {
	    // get local index of the node from elementToNodeMap
	    lcnode1 = elementToNodeMap[b];

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

		for(int kk = 0; kk < m_eqt.nvar; kk++) {

		  dofIndex[kk] = trilinos_index[lcnode0] * m_eqt.nvar + kk;
		  dofIndex[kk + m_eqt.nvar] = trilinos_index[lcnode1] * m_eqt.nvar + kk;

		}
		sparsity->InsertGlobalIndices(dofIndex.size(),
					      &dofIndex.front(),
					      dofIndex.size(),
					      &dofIndex.front());

	      }
	      else if(trilinos_index[lcnode0] != -1) {

		for(int kk = 0; kk < m_eqt.nvar; kk++) {

		  dofIndex[kk] = trilinos_index[lcnode0] * m_eqt.nvar + kk;

		}

		sparsity->InsertGlobalIndices(m_eqt.nvar,
					      &dofIndex.front(),
					      m_eqt.nvar,
					      &dofIndex.front());

	      }
	      else if(trilinos_index[lcnode1] != -1) {

		for(int kk = 0; kk < m_eqt.nvar; kk++) {

		  dofIndex[kk] = trilinos_index[lcnode1] * m_eqt.nvar + kk;

		}

		sparsity->InsertGlobalIndices(m_eqt.nvar,
					      &dofIndex.front(),
					      m_eqt.nvar,
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

  partition.SetBufferSizes(syncedFields, CommRegistry::fractureFlowSolver);

 if(m_doRockPropUpdate) {

   CalculateInterpCoefs(domain);

 }

}

void FractureFlowSolver::DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set, realT time){

  MultiVarDirichletBoundaryCondition* bc = dynamic_cast<MultiVarDirichletBoundaryCondition*> (bcBase);

  if(bc) {

    if(!bc->isClamped()) {

      const rArray1d &value = bc->GetValues(time);

      if(object.GetObjectType() == ObjectDataStructureBaseT::NodeManager && bc->m_objectKey == PhysicalDomainT::FiniteElementNodeManager) {

	rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

	lSet::const_iterator kf=set.begin();
    
	for(localIndex i =0; i < set.size() ; ++i, ++kf) {

	  pressure[*kf] = value[0];

	}

	if(m_eqt.thermal) {

	  kf=set.begin();

	  rArray1d &temperature = domain.m_feNodeManager.GetFieldData<realT>("Temperature");

	  for(localIndex i =0; i < set.size() ; ++i, ++kf) {

	    temperature[*kf] = value[1];

	  }

	}

      }

    }

  }

}

void FractureFlowSolver::SrcFluxBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time)
{

  realT q;
  localIndex index;

  int nvar = m_eqt.nvar;

  Epetra_IntSerialDenseVector  elem_index(nvar);
  Epetra_SerialDenseVector     elem_rhs(nvar);

  MultiVarSrcFluxBoundaryCondition* bc = dynamic_cast<MultiVarSrcFluxBoundaryCondition*> (bcBase);

  if(bc) {

    const rArray1d &fluxValue = bc->GetValues(time);

    if(object.GetObjectType() == ObjectDataStructureBaseT::NodeManager && bc->m_objectKey == PhysicalDomainT::FiniteElementNodeManager) {

      const iArray1d& node_is_ghost = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();

      const iArray1d &trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);

      lSet::const_iterator kf=set.begin();

      index = 0;

      for( localIndex i =0; i < set.size() ; ++i, ++kf ){

	if(node_is_ghost[*kf] < 0 && trilinos_index[*kf] >=0) {

	    elem_index[0] = trilinos_index[*kf] * nvar;
	    q = fluxValue[0];

	    if(bc->isAllocedByWeight())
	      q = fluxValue[0] * bc->allocFac()[index++];

	    elem_rhs[0] = -m_dt * q;

	    if(m_eqt.thermal) {

	      elem_index[1] = trilinos_index[*kf] * nvar + 1;
	      elem_rhs[1] = -m_dt * q * fluxValue[1];

	    }
	  
	    rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	}

      }

    }

  }

}

void FractureFlowSolver::Assemble(PhysicalDomainT&  domain,
				  SpatialPartition& partition,
				  const realT& time)
{
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  matrix->PutScalar(0.0);
  rhs->PutScalar(0.0);
  solution->PutScalar(0.0);

  int nvar = m_eqt.nvar;

  Epetra_IntSerialDenseVector  elem_index (nvar);
  Epetra_SerialDenseVector     elem_rhs     (nvar);
  Epetra_SerialDenseMatrix     elem_matrix  (nvar,nvar);

  int num = 2*nvar;

  Epetra_IntSerialDenseVector  face_index(num);
  Epetra_SerialDenseVector     face_rhs(num);
  Epetra_SerialDenseMatrix     face_matrix(num, num);

  realT rhoav;
  realT term0, term1, term2, term3, term4, term5;

  realT m_compress = m_fluidRockProperty.m_compress;
  realT m_mu = m_fluidRockProperty.m_mu;

  const iArray1d &trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);
  const rArray1d &pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  // loop over elements (accumulation term)

  localIndex lcnode0, lcnode1;

  if(m_eqt.thermal) {

    realT term6, term7;

    realT m_fluidHeatCap = m_fluidRockProperty.m_fluidHeatCap;
    //realT m_rockDensity = m_fluidRockProperty.m_rockDensity;
     //realT rockDensity = m_fluidRockProperty.m_rockDensity[0];

    const rArray1d &temperature = domain.m_feNodeManager.GetFieldData<realT>("Temperature");

    // fracture and matrix elements (ignore rock compresibility

    for(unsigned int i = 0; i < m_elemRegion.size(); i++) {


      const ElementRegionT* elemRegion = m_elemRegion[i];

      for(localIndex j=0; j< elemRegion->m_numElems; j++) {

	/*
        bool isActiveElement = false;
        if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
        {
          isActiveElement = true;
        }
        else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
        {
          if (flowFaceType[elemRegion->m_toFacesRelation[j][0]] == 0) isActiveElement = true;
        }

	if((*m_elem_is_ghost[i])[j] < 0 && isActiveElement) {
	*/
     if(i < m_fracture_regionNumber && flowFaceType[elemRegion->m_toFacesRelation[j][0]] != 0)
    		continue;  


        const rArray1d* solidVolumePtr;
        
	if(i >= m_fracture_regionNumber) {
          solidVolumePtr = &(elemRegion->GetFieldData<realT>("SolidVolume"));
        }

	if((*m_elem_is_ghost[i])[j] < 0) {

	  const realT rockDensity = (i >= m_fracture_regionNumber)? m_fluidRockProperty.m_rockDensity[i-m_fracture_regionNumber] : 0.0;// set to 0 for fluid regions


	  const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(j);
	  // loop over all nodes in elementToNodeMap

	  for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem; a++) {

	    // get local index of the node from elementToNodeMap
	    lcnode0 = elementToNodeMap[a];

	    if(trilinos_index[lcnode0] != -1) {

	      elem_matrix(0,0) = m_compress * m_density[lcnode0] * m_volume[i][j];
	      elem_index(0) = trilinos_index[lcnode0] * nvar;

	      elem_rhs(0) =  (m_density[lcnode0] - m_oldDensity[lcnode0]) * m_volume[i][j]; 

	      elem_matrix(1,1) = m_density[lcnode0] * m_fluidHeatCap * m_volume[i][j];

	      elem_matrix(0,1) = 0.0;

	      elem_matrix(1,0) = m_compress * m_density[lcnode0] * temperature[lcnode0] * m_fluidHeatCap * m_volume[i][j];

	      elem_index(1) = trilinos_index[lcnode0] * nvar + 1;

	      elem_rhs(1) =  (m_density[lcnode0] * temperature[lcnode0] - m_oldDensity[lcnode0] * m_oldTemperature[lcnode0]) * m_fluidHeatCap * m_volume[i][j];

	      if(i >= m_fracture_regionNumber) {

		elem_rhs(1) +=  (temperature[lcnode0] - m_oldTemperature[lcnode0]) * rockDensity * m_rockHeatCap[i][j] * (*solidVolumePtr)[j];

		elem_matrix(1,1) += rockDensity * m_rockHeatCap[i][j] * (*solidVolumePtr)[j];

	      }

	      matrix->SumIntoGlobalValues(elem_index, elem_matrix);
	      rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	    }

	  }

	}

      }

    }
  
    //connection

    for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

      const ElementRegionT* elemRegion = m_elemRegion[i];

      // loop over all elements in the region
      for( localIndex k=0 ; k<elemRegion->m_numElems; k++) {

	/*
        bool isActiveElement = false;
        if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
        {
          isActiveElement = true;
        }
        else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
        {
          if (flowFaceType[elemRegion->m_toFacesRelation[k][0]] == 0) isActiveElement = true;
        }

	if((*m_elem_is_ghost[i])[k] < 0 && isActiveElement) {
	*/

	if(i < m_fracture_regionNumber && flowFaceType[elemRegion->m_toFacesRelation[k][0]] != 0)
	  continue;

	if((*m_elem_is_ghost[i])[k] < 0) {

	  // get the elementToNodeMap for element k
	  const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(k);

	  num = 0;

	  // loop over all nodes in elementToNodeMap
	  for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem-1 ; a++) {

	    // get local index of the node from elementToNodeMap
	    lcnode0 = elementToNodeMap[a];

	    for( unsigned int b=a+1 ; b<elemRegion->m_numNodesPerElem ; b++) {

	      // get local index of the node from elementToNodeMap
	      lcnode1 = elementToNodeMap[b];

	      rhoav = (m_density[lcnode0] + m_density[lcnode1]) * 0.5;
	      term0 = pressure[lcnode0] - pressure[lcnode1] - rhoav * (m_elev[lcnode0] - m_elev[lcnode1]);

	      term1 = rhoav * (m_femCoefs[i])[k][num] / m_mu * m_dt;
 
	      term2 = term0 * 0.5 * m_density[lcnode0] * m_compress * (m_femCoefs[i])[k][num] / m_mu * m_dt;

	      term3 = term0 * 0.5 * m_density[lcnode1] * m_compress * (m_femCoefs[i])[k][num] / m_mu * m_dt;

	      term4 = 0.5 * m_density[lcnode0] * m_compress * (m_elev[lcnode0] - m_elev[lcnode1]);

	      term5 = 0.5 * m_density[lcnode1] * m_compress * (m_elev[lcnode0] - m_elev[lcnode1]);

	      if(trilinos_index[lcnode0] != -1 && trilinos_index[lcnode1] != -1) {

		face_matrix(0,0) = (1.0 - term4) * term1 + term2;
		face_matrix(nvar,nvar) = (1.0 + term5) * term1 - term3;
  	    	 
		face_matrix(0,nvar) = (-1.0 - term5) * term1 + term3;
		face_matrix(nvar,0) = (-1.0 + term4) * term1 - term2;
        
		face_index[0] = trilinos_index[lcnode0] * nvar;
		face_index[1] = trilinos_index[lcnode0] * nvar + 1;
		face_index[nvar] = trilinos_index[lcnode1] * nvar;
		face_index[nvar+1] = trilinos_index[lcnode1] * nvar + 1;  	 

		face_matrix(0,1) = 0.0;
		face_matrix(0,nvar+1) = 0.0;

		face_matrix(nvar,1) = 0.0;
		face_matrix(nvar,nvar+1) = 0.0;

		if(term1 * term0 >= 0.0) {

		  term6 = term1 * term0 * m_fluidHeatCap * temperature[lcnode0];

		  face_matrix(1,0) = ((1.0 - term4) * term1 + term2) * m_fluidHeatCap * temperature[lcnode0];
		  face_matrix(1,1) = term1 * term0 * m_fluidHeatCap;
		  face_matrix(1,nvar) = ((-1.0 - term5) * term1 + term3) * m_fluidHeatCap * temperature[lcnode0];
		  face_matrix(1,nvar+1) = 0.0;

		  face_matrix(nvar+1,0) = ((-1.0 + term4) * term1 - term2) * m_fluidHeatCap * temperature[lcnode0];
		  face_matrix(nvar+1,1) = -term1 * term0 * m_fluidHeatCap;
		  face_matrix(nvar+1,nvar) = ((1.0 + term5) * term1 - term3) * m_fluidHeatCap * temperature[lcnode0];
		  face_matrix(nvar+1,nvar+1) = 0.0;

		}
		else {

		  term6 = term1 * term0 * m_fluidHeatCap * temperature[lcnode1];

		  face_matrix(1,0) = ((1.0 - term4) * term1 + term2) * m_fluidHeatCap * temperature[lcnode1];
		  face_matrix(1,1) = 0.0;
		  face_matrix(1,nvar+1) = term1 * term0 * m_fluidHeatCap;
		  face_matrix(1,nvar) = ((-1.0 - term5) * term1 + term3) * m_fluidHeatCap * temperature[lcnode1];

		  face_matrix(nvar+1,0) = ((-1.0 + term4) * term1 - term2) * m_fluidHeatCap * temperature[lcnode1];
		  face_matrix(nvar+1,1) = 0.0;
		  face_matrix(nvar+1,nvar+1) = -term1 * term0 * m_fluidHeatCap;
		  face_matrix(nvar+1,nvar) = ((1.0 + term5) * term1 - term3) * m_fluidHeatCap * temperature[lcnode1];

		}


		term7 =  (m_femCoefs2[i])[k][num] * (temperature[lcnode0] - temperature[lcnode1]) * m_dt;

		face_matrix(1,1) += (m_femCoefs2[i])[k][num] * m_dt;
		face_matrix(1,nvar+1) += -(m_femCoefs2[i])[k][num] * m_dt;

		face_matrix(nvar+1,1) += -(m_femCoefs2[i])[k][num] * m_dt;
		face_matrix(nvar+1,nvar+1) += (m_femCoefs2[i])[k][num] * m_dt;


		matrix->SumIntoGlobalValues(face_index, face_matrix);

		elem_index(0) = trilinos_index[lcnode0] * nvar;
		elem_rhs(0) =  term1 * term0;
		elem_index(1) = trilinos_index[lcnode0] * nvar + 1;
		elem_rhs(1) =  term6 + term7;

		rhs->SumIntoGlobalValues(elem_index, elem_rhs);

		elem_index(0) = trilinos_index[lcnode1] * nvar;
		elem_rhs(0) = -term1 * term0;
		elem_index(1) = trilinos_index[lcnode1] * nvar + 1;
		elem_rhs(1) = -term6 - term7;

		rhs->SumIntoGlobalValues(elem_index, elem_rhs);


	      }
	      else if(trilinos_index[lcnode0] != -1) {

		elem_index[0] = trilinos_index[lcnode0] * nvar;
		elem_index[1] = trilinos_index[lcnode0] * nvar + 1;

		elem_matrix(0,0) = (1.0 - term4) * term1 + term2;
		elem_matrix(0,1) = 0.0;

		if(term1 * term0 >= 0.0) {

		  term6 = term1 * term0 * m_fluidHeatCap * temperature[lcnode0];
		  elem_matrix(1,0) = ((1.0 - term4) * term1 + term2) * m_fluidHeatCap * temperature[lcnode0];
		  elem_matrix(1,1) = term1 * term0 * m_fluidHeatCap;

		}
		else {

		  term6 = term1 * term0 * m_fluidHeatCap * temperature[lcnode1];
		  elem_matrix(1,0) = ((1.0 - term4) * term1 + term2) * m_fluidHeatCap * temperature[lcnode1];
		  elem_matrix(1,1) = 0.0;

		}

		term7 =  (m_femCoefs2[i])[k][num] * (temperature[lcnode0] - temperature[lcnode1]) * m_dt;
		elem_matrix(1,1) += (m_femCoefs2[i])[k][num] * m_dt;

		matrix->SumIntoGlobalValues(elem_index, elem_matrix);

		elem_rhs(0) =  term1 * term0;
		elem_rhs(1) =  term6 + term7;

		rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	      }
	      else if(trilinos_index[lcnode1] != -1) {


		elem_index[0] = trilinos_index[lcnode1] * nvar;
		elem_index[1] = trilinos_index[lcnode1] * nvar + 1;

		elem_matrix(0,0) = (1.0 + term5) * term1 - term3;
		elem_matrix(0,1) = 0.0;

		if(term1 * term0 >= 0.0) {

		  term6 = term1 * term0 * m_fluidHeatCap * temperature[lcnode0];
		  elem_matrix(1,0) = ((1.0 + term5) * term1 - term3) * m_fluidHeatCap * temperature[lcnode0];
		  elem_matrix(1,1) = 0.0;

		}
		else {

		  term6 = term1 * term0 * m_fluidHeatCap * temperature[lcnode1];
		  elem_matrix(1,0) = ((1.0 + term5) * term1 - term3) * m_fluidHeatCap * temperature[lcnode1];
		  elem_matrix(1,1) = -term1 * term0 * m_fluidHeatCap;

		}

		term7 =  (m_femCoefs2[i])[k][num] * (temperature[lcnode0] - temperature[lcnode1]) * m_dt;
		elem_matrix(1,1) += (m_femCoefs2[i])[k][num] * m_dt;

		matrix->SumIntoGlobalValues(elem_index, elem_matrix);

		elem_rhs(0) =  -term1 * term0;
		elem_rhs(1) =  -term6 - term7;

		rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	      }

	      num++;

	    }
	  }

	}

      }

    }

  }
  else {

    // fracture and matrix elements (ignore rock compresibility

    for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

      const ElementRegionT* elemRegion = m_elemRegion[i];

      for(localIndex j=0; j< elemRegion->m_numElems; j++) {

	/*
        bool isActiveElement = false;
        if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
        {
          isActiveElement = true;
        }
        else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
        {
          if (flowFaceType[elemRegion->m_toFacesRelation[j][0]] == 0) isActiveElement = true;
        }

	if((*m_elem_is_ghost[i])[j] < 0 && isActiveElement) {
	*/

	if(i < m_fracture_regionNumber && flowFaceType[elemRegion->m_toFacesRelation[j][0]] != 0)
	  continue;

	if((*m_elem_is_ghost[i])[j] < 0){

	  const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(j);
	  // loop over all nodes in elementToNodeMap

	  for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem; a++) {

	    // get local index of the node from elementToNodeMap
	    lcnode0 = elementToNodeMap[a];

	    if(trilinos_index[lcnode0] != -1) {

	      elem_matrix(0,0) = m_compress * m_density[lcnode0] * m_volume[i][j];
	      elem_index(0) = trilinos_index[lcnode0] * nvar;

	      elem_rhs(0) =  (m_density[lcnode0] - m_oldDensity[lcnode0]) * m_volume[i][j]; 

	      matrix->SumIntoGlobalValues(elem_index, elem_matrix);
	      rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	    }

	  }

	}

      }

    }
  
    //connection

    for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

      const ElementRegionT* elemRegion = m_elemRegion[i];

      // loop over all elements in the region
      for( localIndex k=0 ; k<elemRegion->m_numElems; k++) {

	/*
        bool isActiveElement = false;
        if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
        {
          isActiveElement = true;
        }
        else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
        {
          if (flowFaceType[elemRegion->m_toFacesRelation[k][0]] == 0) isActiveElement = true;
        }

	if((*m_elem_is_ghost[i])[k] < 0 && isActiveElement) {

	*/


	if(i < m_fracture_regionNumber && flowFaceType[elemRegion->m_toFacesRelation[k][0]] != 0)
	  continue;

	if((*m_elem_is_ghost[i])[k] < 0) {

	  // get the elementToNodeMap for element k
	  const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(k);

	  num = 0;

	  // loop over all nodes in elementToNodeMap
	  for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem-1 ; a++) {

	    // get local index of the node from elementToNodeMap
	    lcnode0 = elementToNodeMap[a];

	    for( unsigned int b=a+1 ; b<elemRegion->m_numNodesPerElem ; b++) {

	      // get local index of the node from elementToNodeMap
	      lcnode1 = elementToNodeMap[b];

	      rhoav = (m_density[lcnode0] + m_density[lcnode1]) * 0.5;
	      term0 = pressure[lcnode0] - pressure[lcnode1] - rhoav * (m_elev[lcnode0] - m_elev[lcnode1]);

	      term1 = rhoav * (m_femCoefs[i])[k][num] / m_mu * m_dt;
 
	      term2 = term0 * 0.5 * m_density[lcnode0] * m_compress * (m_femCoefs[i])[k][num] / m_mu * m_dt;

	      term3 = term0 * 0.5 * m_density[lcnode1] * m_compress * (m_femCoefs[i])[k][num] / m_mu * m_dt;

	      term4 = 0.5 * m_density[lcnode0] * m_compress * (m_elev[lcnode0] - m_elev[lcnode1]);

	      term5 = 0.5 * m_density[lcnode1] * m_compress * (m_elev[lcnode0] - m_elev[lcnode1]);

	      if(trilinos_index[lcnode0] != -1 && trilinos_index[lcnode1] != -1) {

		face_matrix(0,0) = (1.0 - term4) * term1 + term2;
		face_matrix(1,1) = (1.0 + term5) * term1 - term3;
  	    	 
		face_matrix(0,1) = (-1.0 - term5) * term1 + term3;
		face_matrix(1,0) = (-1.0 + term4) * term1 - term2;
        
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
	      else if(trilinos_index[lcnode0] != -1) {

		elem_matrix(0,0) = (1.0 - term4) * term1 + term2;
		elem_index[0] = trilinos_index[lcnode0];

		matrix->SumIntoGlobalValues(elem_index, elem_matrix);

		elem_rhs(0) =  term1 * term0;
		rhs->SumIntoGlobalValues(elem_index, elem_rhs);

	      }
	      else if(trilinos_index[lcnode1] != -1) {

		elem_matrix(0,0) = (1.0 + term5) * term1 - term3;
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

  }

  //update srcflux boundary conditions
  UpdateSrcFluxBCs(domain, time);

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

}

realT FractureFlowSolver::Solve(PhysicalDomainT&  domain,
				SpatialPartition& partition)
{

  realT out, tmp;

  int nvar = m_eqt.nvar;

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

    if(node_is_ghost[i] < 0 && trilinos_index[i] >= 0) {

      rowTmp = static_cast<long long int>(trilinos_index[i]*nvar);
      lid = row_map->LID(rowTmp);
      pressure[i] -= local_solution[lid];

    }

  }

  if(m_eqt.thermal) {

    rArray1d &temperature = domain.m_feNodeManager.GetFieldData<realT>("Temperature");

    for(unsigned int i = 0; i < domain.m_feNodeManager.m_numNodes; i++) {

      if(node_is_ghost[i] < 0 && trilinos_index[i] >= 0) {

	rowTmp = static_cast<long long int>(trilinos_index[i]*nvar+1);
	lid = row_map->LID(rowTmp);
	temperature[i] -= local_solution[lid];

      }

    }

  }

  // re-sync ghost nodes
  partition.SynchronizeFields(syncedFields, CommRegistry::fractureFlowSolver);

  return out;

}

double FractureFlowSolver::TimeStep(const realT& time,
				  const realT& dt,
				  const int cycleNumber,
				  PhysicalDomainT& domain,
				  const sArray1d& namesOfSolverRegions,
				  SpatialPartition& partition,
				  FractunatorBase* const fractunator)
{

  m_dt = dt;

  UpdateDirichletBCs(domain, time);

  UpdateFluidProps(domain);

  MakeOldAcc(domain);

  bool converge = 0;
  realT out;

  for(int i = 0; i < m_NRNumerics.m_maxIters; i++) {

    Assemble(domain,partition, time);
    out = Solve(domain,partition);

    //for debugging purpose
    printf("Current time=%g NR iter=%d Max residual=%g\n",time, i, out); 

    if(out < m_NRNumerics.m_tol) {
      if (i > m_NRNumerics.m_maxIters * 0.5)
      {
        m_stabledt.m_maxdt = std::max( dt * 0.5, m_dtMin);
      }
      else if (i < std::max( int(m_NRNumerics.m_maxIters * 0.03), 4 ) )
      {
        m_stabledt.m_maxdt = std::min( dt * 2.0, m_dtMax);
      }
      else if (i < std::max( int(m_NRNumerics.m_maxIters * 0.1), 8 ) )
      {
        m_stabledt.m_maxdt = std::min( dt * 1.5, m_dtMax);
      }
      else if (i < std::max( int(m_NRNumerics.m_maxIters * 0.2), 8 ) )
      {
        m_stabledt.m_maxdt = std::min( dt * 1.25, m_dtMax);
      }
      else if (i < std::max( int(m_NRNumerics.m_maxIters * 0.35), 16 ) )
      {
        m_stabledt.m_maxdt = std::min( dt * 1.1, m_dtMax);
      }

      converge = 1;
      break;
    }

    UpdateFluidProps(domain);

  }

  if (dt < std::numeric_limits<realT>::min()) m_stabledt.m_maxdt = m_dtMin;

  if(!converge) {

    printf("Max. Residual = %e\n", out);
    throw GPException("Error: hit max. NR iterations -- need timestep cutback (not implemented: ");

  }

  if(m_doApertureUpdate)
    MakeFractureParameters(domain);

  if(m_doRockPropUpdate) {

    UpdateRockProps(domain);
    MakeMatrixParameters(domain);

  }

  return dt;
}

void FractureFlowSolver::MakeFractureParameters(PhysicalDomainT& domain) {

  R1Tensor la, lb;
  localIndex index;
  realT w;
  localIndex num;

  NodeManagerT& nodeManager = domain.m_feNodeManager;
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  for(unsigned int i = 0; i < m_fracture_regionNumber; i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j=0; j < elemRegion->m_numElems; j++) {

      /*
      bool isActiveElement = false;
      if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
      {
        isActiveElement = true;
      }
      else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
      {
        if (flowFaceType[elemRegion->m_toFacesRelation[j][0]] == 0) isActiveElement = true;
      }

     if((*m_elem_is_ghost[i])[j] < 0 && isActiveElement) {

      */

      if(flowFaceType[elemRegion->m_toFacesRelation[j][0]] != 0)
	continue;

      //      if((*m_elem_is_ghost[i])[j] < 0) {

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

	//      }

    }

  }


  Array1dT<rArray1d *> aperture(m_fracture_regionNumber);

  for(unsigned int i = 0; i < m_fracture_regionNumber; i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];
    aperture[i] = &(elemRegion->GetFieldData<realT>(ApertureStr));

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {

      /*
      bool isActiveElement = false;
      if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
      {
        isActiveElement = true;
      }
      else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
      {
        if (flowFaceType[elemRegion->m_toFacesRelation[j][0]] == 0) isActiveElement = true;
      }
      */

      if(flowFaceType[elemRegion->m_toFacesRelation[j][0]] != 0)
	continue;

      //      if (isActiveElement) {

       if((*aperture[i])[j] > m_fluidRockProperty.m_maxAperture)
   (*aperture[i])[j] = m_fluidRockProperty.m_maxAperture;

       if((*aperture[i])[j] < m_fluidRockProperty.m_minAperture)
   (*aperture[i])[j] = m_fluidRockProperty.m_minAperture;

       //      }

    }

  }

  for(unsigned int i = 0; i < m_fracture_regionNumber; i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {
      /*
      bool isActiveElement = false;
      if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
      {
        isActiveElement = true;
      }
      else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
      {
        if (flowFaceType[elemRegion->m_toFacesRelation[j][0]] == 0) isActiveElement = true;
      }

      if((*m_elem_is_ghost[i])[j] < 0 && isActiveElement) {

      */

      if(flowFaceType[elemRegion->m_toFacesRelation[j][0]] != 0)
	continue;

      //          if(isActiveElement) {

	m_volume[i][j] *= (*aperture[i])[j];

	for(unsigned int k = 0; k < (m_femCoefs[i])[j].size(); k++) {

       	  (m_femCoefs[i])[j][k] *= (*aperture[i])[j] * (*aperture[i])[j] * (*aperture[i])[j] / 12.0;


	}

	//	  }

    }
  }

}

void FractureFlowSolver::MakeMatrixParameters(PhysicalDomainT& domain) {

  R1Tensor la, lb;
  localIndex index;
  realT w;
  localIndex num;

  NodeManagerT& nodeManager = domain.m_feNodeManager;
  iArray1d& flowFaceType = domain.m_feFaceManager.GetFieldData<int>("flowFaceType");

  for(unsigned int i = 0; i < m_matrix_regionNumber; i++) {

    index = i + m_fracture_regionNumber;

    ElementRegionT* elemRegion = m_elemRegion[index];

    for(localIndex j=0; j < elemRegion->m_numElems; j++) {
      /*
      bool isActiveElement = false;
      if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
      {
        isActiveElement = true;
      }
      else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
      {
        if (flowFaceType[elemRegion->m_toFacesRelation[j][0]] == 0) isActiveElement = true;
      }

      */

      //      if((*m_elem_is_ghost[index])[j] < 0 && isActiveElement) {


	m_volume[index][j] = 0.0;

	for(unsigned int a=0; a < elemRegion->m_numIntegrationPointsPerElem; a++) {
 
	  m_volume[index][j] += elemRegion->m_detJ(j,a);

	}

	m_volume[index][j] /= realT(elemRegion->m_numNodesPerElem);

	num = 0;

	for(unsigned int m=0; m < elemRegion->m_numNodesPerElem-1; m++) {

	  for(unsigned int n=m+1; n < elemRegion->m_numNodesPerElem; n++) {

	    (m_femCoefs[index])[j][num] = 0.0;

	    for(unsigned int a=0; a < elemRegion->m_numIntegrationPointsPerElem; a++) {

	      w = 0.0;

	      for(unsigned int dir = 0; dir < elemRegion->m_ElementDimension; dir++) {

		w += elemRegion->m_dNdX(j)(a,m)[dir] * elemRegion->m_dNdX(j)(a,n)[dir];
	      } 

	      (m_femCoefs[index])[j][num] -= w * elemRegion->m_detJ(j,a);

	    }

	    num++;
	  }

	}

	//      }

    }

  }


  for(localIndex i = 0; i < m_matrix_regionNumber; i++) {

    index = i + m_fracture_regionNumber;
    ElementRegionT* elemRegion = m_elemRegion[index];

    const rArray1d &permeability = elemRegion->GetFieldData<realT>(PermeabilityStr);
    const rArray1d &porosity = elemRegion->GetFieldData<realT>(PorosityStr);
    rArray1d &solidVolume = elemRegion->GetFieldData<realT>("SolidVolume");

    const realT rockCond = (m_eqt.thermal)? m_fluidRockProperty.m_rockCond[i]: 0.0;

    for(localIndex j=0; j < elemRegion->m_numElems; j++) {
      /*
      bool isActiveElement = false;
      if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
      {
        isActiveElement = true;
      }
      else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGeometryID == "TRSH")
      {
        if (flowFaceType[elemRegion->m_toFacesRelation[j][0]] == 0) isActiveElement = true;
      }
      */

      //      if((*m_elem_is_ghost[index])[j] < 0 && isActiveElement) {

        solidVolume[j] = m_volume[index][j]*(1-porosity[j]);
	m_volume[index][j] *= porosity[j];

	for(unsigned int k = 0; k < (m_femCoefs[index])[j].size(); k++) {

          if(m_eqt.thermal){
	    (m_femCoefs2[index])[j][k] = (m_femCoefs[index])[j][k] * rockCond;  // moved here to avoid /0 when permeability =0
          }

	  (m_femCoefs[index])[j][k] *= permeability[j];

	}

	//      }

    }
  }


  if(m_eqt.thermal) {

	// This part was for using uniform thermal properties in the rock/matrix
    //realT m_rockCond = m_fluidRockProperty.m_rockCond;
    //realT m_rockHeatCapCoef = m_fluidRockProperty.m_rockHeatCap;
	  //realT rockCond = m_fluidRockProperty.m_rockCond[0];
	  //realT rockHeatCap = m_fluidRockProperty.m_rockHeatCap[0];


    for(localIndex i = 0; i < m_matrix_regionNumber; i++) {

      //const realT& rockCond = m_fluidRockProperty.m_rockCond[i];
      const realT& rockHeatCap = m_fluidRockProperty.m_rockHeatCap[i];

      index = i + m_fracture_regionNumber;
      ElementRegionT* elemRegion = m_elemRegion[index];

      const rArray1d &porosity = elemRegion->GetFieldData<realT>(PorosityStr);

      const rArray1d &permeability = elemRegion->GetFieldData<realT>(PermeabilityStr);
      for(localIndex j=0; j < elemRegion->m_numElems; j++) {
	/*
        bool isActiveElement = false;
        if (elemRegion->m_elementGeometryID == "C3D4" || elemRegion->m_elementGeometryID == "C3D8" || elemRegion->m_elementGeometryID == "STRI" || elemRegion->m_elementGeometryID == "CPE4")
        {
          isActiveElement = true;
        }
        else if (elemRegion->m_elementGeometryID == "S4R" || elemRegion->m_elementGvoeometryID == "TRSH")
        {
          if (flowFaceType[elemRegion->m_toFacesRelation[j][0]] == 0) isActiveElement = true;
        }
	*/

	//	if((*m_elem_is_ghost[index])[j] < 0 && isActiveElement) {

	  m_rockHeatCap[index][j] = rockHeatCap * (1.0 - porosity[j]);

	/*  for(unsigned int k = 0; k < (m_femCoefs[index])[j].size(); k++) {

	    (m_femCoefs2[index])[j][k] = (m_femCoefs[index])[j][k] * rockCond / permeability[j];
            // SW moved above before *permeability
	  }*/

	  //	}

      }
    }
  

    iArray1d flags(nodeManager.m_numNodes);
    flags = 0;

    for(localIndex i = 0; i < m_thermalProductionBCSetName.size(); i++) {

      lSet& set = nodeManager.GetSet(m_thermalProductionBCSetName[i]);

      lSet::const_iterator kn =set.begin();
    
      for(localIndex j =0; j < set.size() ; ++j, ++kn) {

	  flags[*kn] = 1;

      }

    }
  
    for(unsigned int i = 0; i < m_elemRegion.size(); i++) {

      ElementRegionT* elemRegion = m_elemRegion[i];

      for(localIndex j=0; j < elemRegion->m_numElems; j++) {

	const localIndex* const elemToNodeMap = elemRegion->m_toNodesRelation[j];

	num = 0;

	for(unsigned int m=0; m < elemRegion->m_numNodesPerElem-1; m++) {

	  localIndex p = elemToNodeMap[m];

	  for(unsigned int n=m+1; n < elemRegion->m_numNodesPerElem; n++) {

	    localIndex q = elemToNodeMap[n];

	    if(flags[p] == 1 || flags[q] == 1) {
	      (m_femCoefs2[i])[j][num] = 0.0;
	    }

	    num++;

	  }

	}

      }

    }

  }

}

void FractureFlowSolver::CalculateInterpCoefs(PhysicalDomainT& domain) {

  int index;
  NodeManagerT& nodeManager = domain.m_feNodeManager;

  m_interpCoefs.resize(m_matrix_regionNumber);

  for(localIndex i = 0; i < m_matrix_regionNumber; i++) {

    index = i + m_fracture_regionNumber;
    ElementRegionT* elemRegion = m_elemRegion[index];
    m_interpCoefs[i].resize2(elemRegion->m_numElems, elemRegion->m_numNodesPerElem);

    for(localIndex j=0; j < elemRegion->m_numElems; j++) {

      //consider using more accurate interpolation later
      for(unsigned int a = 0 ; a < elemRegion->m_numNodesPerElem; a++) {
	m_interpCoefs[i](j, a) = 1.0 / realT(elemRegion->m_numNodesPerElem);
      }
    }

  }

  const rArray1d& pressure  = nodeManager.GetFieldData<FieldInfo::pressure>();

  //define initial porosity, permeability and pressure

  for(localIndex i = 0; i < m_matrix_regionNumber; i++) {

      index = i + m_fracture_regionNumber;
      ElementRegionT* elemRegion = m_elemRegion[index];

      const rArray1d &porosity = elemRegion->GetFieldData<realT>(PorosityStr);
      const rArray1d &permeability = elemRegion->GetFieldData<realT>(PermeabilityStr);

      rArray1d &initPorosity = elemRegion->GetFieldData<realT>("InitPorosity");
      rArray1d &initPermeability = elemRegion->GetFieldData<realT>("InitPermeability");
     rArray1d &initPressure = elemRegion->GetFieldData<realT>("InitPressure");

      for(localIndex j=0; j < elemRegion->m_numElems; j++) {

	const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(j);
	initPressure[j] = 0.0;

	// loop over all nodes in elementToNodeMap
	for(unsigned int a = 0 ; a < elemRegion->m_numNodesPerElem; a++) {

	  localIndex nodeIndex = nodeManager.GetParentIndex(elementToNodeMap[a]);

	  initPressure[j] += pressure[nodeIndex] * m_interpCoefs[i](j ,a);

	}
         
	//update porosity and permeability

	initPorosity[j] = porosity[j];
	initPermeability[j] = permeability[j];
      
      }
 
  }

}


void FractureFlowSolver::UpdateFluidProps(PhysicalDomainT& domain) {

  const rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  for(unsigned int i=0; i<domain.m_feNodeManager.m_numNodes; i++) {

    m_density[i] = m_fluidRockProperty.m_rho_o * exp(m_fluidRockProperty.m_compress * (pressure[i] - m_fluidRockProperty.m_press_o));

  }

}


void FractureFlowSolver::UpdateRockProps(PhysicalDomainT& domain) {

  const rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  realT value;
  localIndex index;

  for(localIndex i = 0; i < m_matrix_regionNumber; i++) {

      index = i + m_fracture_regionNumber;
      ElementRegionT* elemRegion = m_elemRegion[index];

      rArray1d &porosity = elemRegion->GetFieldData<realT>(PorosityStr);
      rArray1d &permeability = elemRegion->GetFieldData<realT>(PermeabilityStr);

      const rArray1d &initPorosity = elemRegion->GetFieldData<realT>("InitPorosity");
      const rArray1d &initPermeability = elemRegion->GetFieldData<realT>("InitPermeability");
      const rArray1d &initPressure = elemRegion->GetFieldData<realT>("InitPressure");

      for(localIndex j=0; j < elemRegion->m_numElems; j++) {

	const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(j);
	value = 0.0;

	// loop over all nodes in elementToNodeMap
	for(unsigned int a = 0 ; a < elemRegion->m_numNodesPerElem; a++) {

	  localIndex nodeIndex = domain.m_feNodeManager.GetParentIndex(elementToNodeMap[a]);

	  value += pressure[nodeIndex] * m_interpCoefs[i](j ,a);

	}
         
	//update porosity and permeability

	porosity[j] = initPorosity[j] * exp((m_fluidRockProperty.m_rockCompress)[i] * (value - initPressure[j]));

	permeability[j] = initPermeability[j] * exp((m_fluidRockProperty.m_rockPermCoef)[i] * (value - initPressure[j]));
      
      }
 
  }

}

void FractureFlowSolver::SetDirichletBCs(PhysicalDomainT& domain)
{

    ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SetDirichletNodeBoundaryCondition, domain, domain.m_feNodeManager, "MultiVars", 0.0);

}

void FractureFlowSolver::SetDirichletNodeBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time) {

  if(object.GetObjectType() == ObjectDataStructureBaseT::NodeManager && bcBase->m_objectKey == PhysicalDomainT::FiniteElementNodeManager) {

    MultiVarDirichletBoundaryCondition* bc = dynamic_cast<MultiVarDirichletBoundaryCondition*> (bcBase);

    if(bc) {

      bc->CheckVars(m_fieldName);

      iArray1d &trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);

      lSet::const_iterator kf=set.begin();
    
      for(localIndex i =0; i < set.size() ; ++i, ++kf) {

	trilinos_index[*kf] = -1;

      }

    }

  }

}

void FractureFlowSolver::UpdateDirichletBCs(PhysicalDomainT& domain, const realT &time)
{

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::DirichletBoundaryCondition, domain, domain.m_feNodeManager, "MultiVars", 0.0);

}


void FractureFlowSolver::SetSrcFluxBCs(PhysicalDomainT& domain)
{

  // flux into node-based control volume

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SetSrcFluxBoundaryCondition, domain, domain.m_feNodeManager, "MultiVarSrcFlux", 0.0);                                       

}


void FractureFlowSolver::SetSrcFluxBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time)
{

  realT value, totalValue = 0.0;
  realT sum;
  localIndex m, n;

  MultiVarSrcFluxBoundaryCondition* bc = dynamic_cast<MultiVarSrcFluxBoundaryCondition*> (bcBase);

  if(bc) {

      bc->CheckVars(m_fieldName);

      if(bc->isAllocedByWeight()) {

	if(bc->m_objectKey == PhysicalDomainT::FiniteElementNodeManager) {

	  const iArray1d &trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);

	  const iArray1d& node_is_ghost = domain.m_feNodeManager.GetFieldData<FieldInfo::ghostRank>();

	  lSet::const_iterator kf=set.begin();

	  for( localIndex i =0; i < set.size() ; ++i, ++kf ){

	    if(node_is_ghost[*kf] < 0 && trilinos_index[*kf] >=0) {

	      value = 0.0;

	      for(std::set<std::pair< ElementRegionT*, localIndex > >::const_iterator elem = domain.m_feNodeManager.m_toElementsRelation[*kf].begin(); elem != domain.m_feNodeManager.m_toElementsRelation[*kf].end() ; ++elem) {

		m = m_elemRegionIndex[(elem->first)->m_regionNumber];
		n = elem->second;

		value += m_volume[m][n];

	      }

	      totalValue += value;
	      bc->allocFac().push_back(value);

	    }   

	  }

	  MPI_Allreduce(&totalValue, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	  for(localIndex i = 0; i < bc->allocFac().size(); ++i) {

	    bc->allocFac()[i] /= sum;

	  }

	}

      }

  }

}


void FractureFlowSolver::UpdateSrcFluxBCs(PhysicalDomainT& domain, const realT &time)
{

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolver::SrcFluxBoundaryCondition, domain, domain.m_feNodeManager, "MultiVarSrcFlux", 0.0);                                       
}


void FractureFlowSolver::MakeOldAcc(PhysicalDomainT& domain) {

  for(unsigned int i=0; i<domain.m_feNodeManager.m_numNodes; i++) {

      m_oldDensity[i] = m_density[i];

  }

  if(m_eqt.thermal) {

    const rArray1d &temperature = domain.m_feNodeManager.GetFieldData<realT>("Temperature");

    for(unsigned int i=0; i<domain.m_feNodeManager.m_numNodes; i++) {

      m_oldTemperature[i] = temperature[i];

    }

  }

}

/// Register solver in the solver factory
REGISTER_SOLVER( FractureFlowSolver )



/* Note for Yue
The face correponding to an flow element: faceID = elemRegion->m_toFacesRelation[elemID][0]
For 3D problems, if flowFaceType of a face is zero, this face should have one and only one child
childFaceID = m_feFaceManager.m_childIndices[faceID][0]

There are two solid elements connected to a flow face element (unless the flow face element is at the domain boundary, i.e. it is not a fracture)
These two elements can be accessed via:
The first one: faceManager.m_toElementsRelation[faceID][0]
The the second one: faceManager.m_toElementsRelation[childFaceID][0]
*/
