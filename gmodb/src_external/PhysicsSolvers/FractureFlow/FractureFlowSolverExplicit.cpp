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
 * @file FractureFlowSolverExplicit.cpp
 * @author hao1
 * @date Oct. 21, 2013
 */

#include "FractureFlowSolverExplicit.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Common/Common.h"
#include "Utilities/Utilities.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"


// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

using namespace BoundaryConditionFunctions;
using namespace PS_STR;

FractureFlowSolverExplicit::FractureFlowSolverExplicit(const std::string &name, ProblemManagerT* const pm):
SolverBase(name, pm),
syncedFields(),
m_TrilinosIndexStr(),
m_dt(0.0)
{

  ++m_instances; 
  m_TrilinosIndexStr = "FFS_SS_" +  toString<int>(m_instances) + "_GlobalDof";
}

FractureFlowSolverExplicit::~FractureFlowSolverExplicit() {}

/*

      <FractureFlowSolverExplicit name="fracMatrixFlowSolverExplicit"
                          coupleMatrixFractureFlow="1"
                          fractureRegionName="EB200" 
                          updateAperture="0" 
                          fieldName="Pressure"
                          rockCompress="3.5e-10 1/Pa"
                          downVector="0.0 0.0 0.0"
      />


*/

void FractureFlowSolverExplicit::ReadXML( TICPP::HierarchicalDataNode* hdn )
{ 

  m_doApertureUpdate = hdn->GetAttributeOrDefault<bool>("updateAperture",false);
  m_doMFCoupling = hdn->GetAttributeOrDefault<bool>("coupleMatrixFractureFlow",false);

  // equations (flow only)
  m_eqt.ncomp = 1;
  m_eqt.nvar = 1;

  // fluid and rock properties

  m_fluidRockProperty.m_compress = hdn->GetAttributeOrDefault("fluidCompress", "5.0e-10 1/Pa");

  m_fluidRockProperty.m_mu = hdn->GetAttributeOrDefault("fluidViscosity", "0.001 Pa.s");

  m_fluidRockProperty.m_rho_o = hdn->GetAttributeOrDefault("rhoRef","1000.0 kg/m^3");

  m_fluidRockProperty.m_press_o = hdn->GetAttributeOrDefault("pressRef","1.0e5 Pa");

  m_fluidRockProperty.m_rockCompress = hdn->GetAttributeOrDefault("rockCompress", "3.5e-10 1/Pa");

  m_fluidRockProperty.m_minAperture = hdn->GetAttributeOrDefault(MinimumApertureStr,"1 um");

  m_fluidRockProperty.m_maxAperture = hdn->GetAttributeOrDefault(MaximumApertureStr,"4 mm");

  m_fluidRockProperty.m_maxPorosity = hdn->GetAttributeOrDefault("maxrockPorosity","0.999");

  m_fluidRockProperty.m_minPorosity = hdn->GetAttributeOrDefault("minRockPorosity","0.001");

  m_fracture_regionName = hdn->GetStringVector("fractureRegionName");

  m_fracture_regionNumber = m_fracture_regionName.size();

  if(m_fracture_regionName.empty())
      throw GPException("ERROR: 'fractureRegionName' is not set");

  m_fieldName = hdn->GetStringVector("fieldName");

  if(m_fieldName.empty())
      throw GPException("ERROR: fieldName is not set!");

}

void FractureFlowSolverExplicit::RegisterFields( PhysicalDomainT& domain )
{

  ElementManagerT& elementManager = domain.m_feElementManager;
  NodeManagerT& nodeManager = domain.m_feNodeManager;

  nodeManager.AddKeyedDataField<FieldInfo::pressure>();    

  nodeManager.AddKeylessDataField<int>(m_TrilinosIndexStr,true,true);

  int num;

  num = (elementManager.m_ElementRegions).size();

  m_elemRegionIndex.resize(num);
  m_elemRegionIndex = -1;

  for(int i =0; i < m_fracture_regionNumber; ++i){

    std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.find(m_fracture_regionName[i]);

    if(it == elementManager.m_ElementRegions.end())
      throw GPException("Fracture region " + m_fracture_regionName[i] + " is not found!");

    ElementRegionT& elemRegion = it->second;
    m_elemRegionIndex[elemRegion.m_regionNumber] = i;

    m_elemRegion.push_back(&elemRegion);

    elemRegion.AddKeylessDataField<realT>(ApertureStr,true,true);

  }

  if(m_doMFCoupling) {

    int index;

    m_matrix_regionNumber = 0;

    std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.begin();

    for(localIndex i =0; i < elementManager.m_ElementRegions.size(); i++, it++){
      ElementRegionT& elemRegion = it->second;

      if(m_elemRegionIndex[elemRegion.m_regionNumber] == -1) {

	index = m_matrix_regionNumber + m_fracture_regionNumber;
	m_elemRegionIndex[elemRegion.m_regionNumber] = index;

	m_elemRegion.push_back(&elemRegion);

	elemRegion.AddKeylessDataField<realT>(PermeabilityStr,true,true);
	elemRegion.AddKeylessDataField<realT>(PorosityStr,true,true);

	m_matrix_regionNumber++;

      }
    }
  }

}

void FractureFlowSolverExplicit::InitializeCommunications( PartitionBase& partition )
{

    syncedFields.clear();

    syncedFields[PhysicalDomainT::FiniteElementNodeManager].push_back(Field<FieldInfo::pressure>::Name());

}


void FractureFlowSolverExplicit::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{

  NodeManagerT& nodeManager = domain.m_feNodeManager;

  int m, num;

  int node_num = nodeManager.m_numNodes;
	    
  num = m_elemRegion.size();

  m_elem_is_ghost.resize(num);

  m_femCoefs.resize(num);

  m_density.resize(node_num);
  m_volume.resize(node_num);

  for(int i = 0; i < num; i++) {

    m_elem_is_ghost[i] = &(m_elemRegion[i]->GetFieldData<FieldInfo::ghostRank>());
    m_femCoefs[i].resize(m_elemRegion[i]->m_numElems);

    m = m_elemRegion[i]->m_numNodesPerElem * (m_elemRegion[i]->m_numNodesPerElem - 1) / 2;

    for(localIndex j = 0; j < m_elemRegion[i]->m_numElems; j++) {
      m_femCoefs[i][j].resize(m);
    }

  }

  iArray1d &trilinos_index = nodeManager.GetFieldData<int>(m_TrilinosIndexStr);

  trilinos_index = 0;

  // initialize Dirichlet boundary conditions

  SetDirichletBCs(domain);

  // calculate some geometric parameters 
  // (e.g. element volume, connection area)

  MakeParameters(domain);

  partition.SetBufferSizes(syncedFields, CommRegistry::fractureFlowSolver);

}

void FractureFlowSolverExplicit::DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set, realT time){
 
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

      }

    }

  }

}



realT FractureFlowSolverExplicit::Solve(const realT& dt,
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

  rArray1d solution(node_num);

  solution = 0.0;

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
	  lcnode0 = elementToNodeMap[a];

	  for( unsigned int b=a+1 ; b<elemRegion->m_numNodesPerElem ; b++) {

	    // get local index of the node from elementToNodeMap
	    lcnode1 = elementToNodeMap[b];

	    rhoav = (m_density[lcnode0] + m_density[lcnode1]) * 0.5;

	    term0 = pressure[lcnode0] - pressure[lcnode1];

	    term1 = rhoav * (m_femCoefs[i])[k][num] / m_mu;
 
	    if(trilinos_index[lcnode0] != -1 && trilinos_index[lcnode1] != -1) {
	      solution[lcnode0] +=  - term1 * term0;
	      solution[lcnode1] +=  term1 * term0;

	    }
	    else if(trilinos_index[lcnode0] != -1) {

	      solution[lcnode0] +=  - term1 * term0;

	    }
	    else if(trilinos_index[lcnode1] != -1) {

	      solution[lcnode1] +=  term1 * term0;

	    }

	    num++;

	  }
	}

      }

    }

  }

  realT out = -1e20;

  for(int i = 0; i < node_num; i++) {

    if(node_is_ghost[i] < 0 && trilinos_index[i] >= 0) {

      rho = solution[i] / m_volume[i];

      if(fabs(rho) > out)
	out = fabs(rho);

      rho = m_density[i] + rho * dt;

      pressure[i] = log(rho / m_fluidRockProperty.m_rho_o) / m_compress + m_fluidRockProperty.m_press_o;

    }

  }

  // re-sync ghost nodes
  partition.SynchronizeFields(syncedFields, CommRegistry::fractureFlowSolver);

  return out;

}

double FractureFlowSolverExplicit::TimeStep(const realT& time,
					  const realT& dt,
					  const int cycleNumber,
					  PhysicalDomainT& domain,
					  const sArray1d& namesOfSolverRegions,
					  SpatialPartition& partition,
					  FractunatorBase* const fractunator)
{


  m_dt = dt;

  if(m_doApertureUpdate)
    MakeParameters(domain);

  UpdateDirichletBCs(domain, time);

  UpdateFluidRockProps(domain);

  realT resid = Solve(dt, domain, partition);

  // for debugging  purpose
  printf("time=%g cycle=%d dt=%g resid=%g\n",time, cycleNumber, dt, resid); 


  return dt;
}

void FractureFlowSolverExplicit::MakeParameters(PhysicalDomainT& domain) {

  R1Tensor la, lb;
  localIndex index;
  realT w;
  int kf, eg;
  localIndex num;

  NodeManagerT& nodeManager = domain.m_feNodeManager;

  //loop over all local pairs of all the elements

  rSArray2d volume; 
  volume.resize(m_elemRegion.size());

  for(localIndex i = 0; i < m_elemRegion.size(); i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];

    volume[i].resize(m_elemRegion[i]->m_numElems);

    for(localIndex j=0; j < elemRegion->m_numElems; j++) {

	volume[i][j] = 0.0;

	for(localIndex a=0; a < elemRegion->m_numIntegrationPointsPerElem; a++) {
 
	  volume[i][j] += elemRegion->m_detJ(j,a);

	}

	volume[i][j] /= realT(elemRegion->m_numNodesPerElem);

	num = 0;

	for(localIndex m=0; m < elemRegion->m_numNodesPerElem-1; m++) {

	  for(localIndex n=m+1; n < elemRegion->m_numNodesPerElem; n++) {

	    (m_femCoefs[i])[j][num] = 0.0;

	    for(localIndex a=0; a < elemRegion->m_numIntegrationPointsPerElem; a++) {

	      w = 0.0;

	      for(localIndex dir = 0; dir < elemRegion->m_ElementDimension; dir++) {

		w += elemRegion->m_dNdX(j)(a,m)[dir] * elemRegion->m_dNdX(j)(a,n)[dir];
	      } 

	      (m_femCoefs[i])[j][num] -= w * elemRegion->m_detJ(j,a);

	    }

	    num++;
	  }

	}

    }

  }


  Array1dT<rArray1d *> aperture(m_fracture_regionNumber);

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];
    aperture[i] = &(elemRegion->GetFieldData<realT>(ApertureStr));

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {

      if((*aperture[i])[j] > m_fluidRockProperty.m_maxAperture)  
	(*aperture[i])[j] = m_fluidRockProperty.m_maxAperture;

      if((*aperture[i])[j] < m_fluidRockProperty.m_minAperture)  
	(*aperture[i])[j] = m_fluidRockProperty.m_minAperture;

    }

  }

  // fracture elements

  for(int i = 0; i < m_fracture_regionNumber; i++) {

    ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {

      if((*m_elem_is_ghost[i])[j] < 0) {

	volume[i][j] *= (*aperture[i])[j];

	for(localIndex k = 0; k < (m_femCoefs[i])[j].size(); k++) {

	  (m_femCoefs[i])[j][k] *= (*aperture[i])[j] * (*aperture[i])[j] * (*aperture[i])[j] / 12.0;

	}

      }

    }
  }

  for(int i = 0; i < m_matrix_regionNumber; i++) {

    index = i + m_fracture_regionNumber;
    ElementRegionT* elemRegion = m_elemRegion[index];

    const rArray1d &permeability = elemRegion->GetFieldData<realT>(PermeabilityStr);
    const rArray1d &porosity = elemRegion->GetFieldData<realT>(PorosityStr);

    for(localIndex j=0; j < elemRegion->m_numElems; j++) {

      if((*m_elem_is_ghost[index])[j] < 0) {

	volume[index][j] *= porosity[j];
       
	for(localIndex k = 0; k < (m_femCoefs[index])[j].size(); k++) {

	  (m_femCoefs[index])[j][k] *= permeability[j];

	}

      }

    }

  }

  const iArray1d &trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_TrilinosIndexStr);

  localIndex lcnode;
  m_volume = 0.0;

  for(localIndex i = 0; i < m_elemRegion.size(); i++) {

    const ElementRegionT* elemRegion = m_elemRegion[i];

    for(localIndex j=0; j< elemRegion->m_numElems; j++) {

      const localIndex* const elementToNodeMap = elemRegion->ElementToNodeMap(j);
      // loop over all nodes in elementToNodeMap

      for( unsigned int a=0 ; a<elemRegion->m_numNodesPerElem; a++) {

	// get local index of the node from elementToNodeMap
	lcnode = elementToNodeMap[a];

	if(trilinos_index[lcnode] != -1) {

	  m_volume[lcnode] += volume[i][j];

	}

      }

    }

  }

}

void FractureFlowSolverExplicit::UpdateFluidRockProps(PhysicalDomainT& domain) {

  const rArray1d& pressure  = domain.m_feNodeManager.GetFieldData<FieldInfo::pressure>();

  for(localIndex i=0; i<domain.m_feNodeManager.m_numNodes; i++) {

    m_density[i] = m_fluidRockProperty.m_rho_o * exp(m_fluidRockProperty.m_compress * (pressure[i] - m_fluidRockProperty.m_press_o));

  }

}

void FractureFlowSolverExplicit::SetDirichletBCs(PhysicalDomainT& domain)
{

    ApplyBoundaryCondition<realT>(this, &FractureFlowSolverExplicit::SetDirichletNodeBoundaryCondition, domain, domain.m_feNodeManager, "MultiVars", 0.0);

}

void FractureFlowSolverExplicit::SetDirichletNodeBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time) {

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

void FractureFlowSolverExplicit::UpdateDirichletBCs(PhysicalDomainT& domain, const realT &time)
{

  ApplyBoundaryCondition<realT>(this, &FractureFlowSolverExplicit::DirichletBoundaryCondition, domain, domain.m_feNodeManager, "MultiVars", 0.0);

}


/// Register solver in the solver factory
REGISTER_SOLVER( FractureFlowSolverExplicit )
