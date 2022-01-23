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
 * @file ImplicitMechanicsSolver.cpp
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#include "LagrangeSmallStrainLinearElastic.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "PhysicsSolvers/PhysicsSolverStrings.h"


#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/LagrangeBasis.h"
#include "ElementLibrary/GaussQuadrature.h"


#include "Constitutive/CohesiveZone/InitiallyRigidCohesiveZone.h"


#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"



LagrangeSmallStrainLinearElastic::LagrangeSmallStrainLinearElastic(  const std::string& name,
                                                                          ProblemManagerT* const pm ):
LagrangeSmallStrain(name,pm)
{}


LagrangeSmallStrainLinearElastic::~LagrangeSmallStrainLinearElastic()
{
  // TODO Auto-generated destructor stub
}





realT LagrangeSmallStrainLinearElastic::CalculateElementResidualAndDerivative( const MaterialBaseParameterData& matParams,
                                                                                   const FiniteElementBase& fe,
                                                                                   const Array2dT<R1Tensor>& dNdX,
                                                                                   const realT* const detJ,
                                                                                   const Epetra_SerialDenseVector& dof_np1,
                                                                                   Epetra_SerialDenseMatrix& dRdU,
                                                                                   Epetra_SerialDenseVector& R )
{
  realT maxForce = 0;
  const realT G = matParams.init_shearModulus;
  realT lambda = matParams.Lame;

  if( LagrangeSolverBase::m_2dOption==LagrangeSolverBase::PlaneStress )
  {
    lambda = 2*lambda*G / ( lambda + 2*G );
  }

  dRdU.Scale(0);
  R.Scale(0);


  R1Tensor dNdXa;
  R1Tensor dNdXb;

  for( unsigned int q=0 ; q<fe.n_quadrature_points() ; ++q )
  {
    for( unsigned int a=0 ; a<fe.dofs_per_element() ; ++a )
    {
      dNdXa = dNdX(q,a);

      for( unsigned int b=0 ; b<fe.dofs_per_element() ; ++b )
      {
        dNdXb = dNdX(q,b);

        if( dim==3 )
        {
          dRdU(a*dim+0,b*dim+0) -= ( (dNdXa[1]*dNdXb[1]+ dNdXa[2]*dNdXb[2])*G + dNdXa[0]*dNdXb[0]*(2*G + lambda) ) * detJ[q];
          dRdU(a*dim+0,b*dim+1) -= ( dNdXa[1]*dNdXb[0]*G + dNdXa[0]*dNdXb[1]*lambda ) * detJ[q];
          dRdU(a*dim+0,b*dim+2) -= ( dNdXa[2]*dNdXb[0]*G + dNdXa[0]*dNdXb[2]*lambda ) * detJ[q];

          dRdU(a*dim+1,b*dim+0) -= ( dNdXa[0]*dNdXb[1]*G + dNdXa[1]*dNdXb[0]*lambda ) * detJ[q];
          dRdU(a*dim+1,b*dim+1) -= ( (dNdXa[0]*dNdXb[0] + dNdXa[2]*dNdXb[2])*G + dNdXa[1]*dNdXb[1]*(2*G + lambda) ) * detJ[q];
          dRdU(a*dim+1,b*dim+2) -= ( dNdXa[2]*dNdXb[1]*G + dNdXa[1]*dNdXb[2]*lambda ) * detJ[q];

          dRdU(a*dim+2,b*dim+0) -= ( dNdXa[0]*dNdXb[2]*G + dNdXa[2]*dNdXb[0]*lambda ) * detJ[q];
          dRdU(a*dim+2,b*dim+1) -= ( dNdXa[1]*dNdXb[2]*G + dNdXa[2]*dNdXb[1]*lambda ) * detJ[q];
          dRdU(a*dim+2,b*dim+2) -= ( (dNdXa[0]*dNdXb[0] + dNdXa[1]*dNdXb[1])*G + dNdXa[2]*dNdXb[2]*(2*G + lambda) ) * detJ[q];
        }
        else if( dim==2 )
        {
          dRdU(a*dim+0,b*dim+0) -= ( dNdXa[1]*dNdXb[1]*G + dNdXa[0]*dNdXb[0]*(2*G + lambda) ) * detJ[q];
          dRdU(a*dim+0,b*dim+1) -= ( dNdXa[1]*dNdXb[0]*G + dNdXa[0]*dNdXb[1]*lambda ) * detJ[q];

          dRdU(a*dim+1,b*dim+0) -= ( dNdXa[0]*dNdXb[1]*G + dNdXa[1]*dNdXb[0]*lambda ) * detJ[q];
          dRdU(a*dim+1,b*dim+1) -= ( dNdXa[0]*dNdXb[0]*G + dNdXa[1]*dNdXb[1]*(2*G + lambda) ) * detJ[q];
        }
      }
    }
  }


  for( unsigned int a=0 ; a<fe.dofs_per_element()*dim ; ++a )
  {
    realT nodeForce = 0;
    for( unsigned int b=0 ; b<fe.dofs_per_element()*dim ; ++b )
    {
      nodeForce += dRdU(a,b) * dof_np1(b);
      R(a) += dRdU(a,b) * dof_np1(b);
    }
    if( fabs(nodeForce) > maxForce )
    {
      maxForce = fabs(nodeForce);
    }
  }

  return maxForce;
}


void LagrangeSmallStrainLinearElastic::ApplyThermalStress( ElementRegionT& elemRegion,
                         const NodeManagerT& nodeManager,
                         const localIndex& elementID,
                         Epetra_SerialDenseVector& rhs,
                         const bool updateTemperatureOnly)
{

  rArray1d& temperature = elemRegion.GetFieldData<realT>("temperature");
  rArray1d& CTE = elemRegion.GetFieldData<realT>("linearCTE");

  if (m_useNodalTemperature > 0)
  {
    const rArray1d& nodalTemp = nodeManager.GetFieldData<realT>("Temperature");
    temperature[elementID] = 0.0;
    for (localIndex i = 0; i<elemRegion.m_toNodesRelation.Dimension(1); ++i)
    {
      temperature[elementID] += nodalTemp[elemRegion.m_toNodesRelation[elementID][i]];
    }
    temperature[elementID] /= elemRegion.m_toNodesRelation.Dimension(1);
  }

  if (!updateTemperatureOnly)
  {
    const localIndex paramIndex = elemRegion.m_mat->NumParameterIndex0() > 1 ? elementID : 0 ;

    const MaterialBaseParameterData& matParams = *(elemRegion.m_mat->ParameterData(paramIndex));

    const realT K =  matParams.Lame + matParams.init_shearModulus * 2.0 / 3.0;

    R2SymTensor thermalStress;
    thermalStress *= 0.0;

    for (int i = 0; i < dim; ++i) thermalStress(i,i) = (temperature[elementID] - m_refTemperature) * CTE [elementID] * K * 3;  //The factor 3 is because CTE is the linear CTE


    Array1dT<R1Tensor> fNode(elemRegion.m_numNodesPerElem);
    elemRegion.CalculateNodalForceFromStress(elementID, nodeManager, thermalStress, fNode);

    for (localIndex i = 0; i < elemRegion.m_numNodesPerElem; ++i)
    {
      for (int j = 0; j < dim; ++j)
      {
        rhs(i*dim+j) = fNode[i][j];
      }
    }
  }

}



void LagrangeSmallStrainLinearElastic::PostSyncConsistency( PhysicalDomainT& domain, SpatialPartition& partition )
{

  CalculateElementStresses( domain.m_feNodeManager, domain.m_feElementManager );



  if(true)
    if( domain.m_feNodeManager.HasField<R1Tensor>("cohesiveForce") &&
        domain.m_feFaceManager.HasField<R1Tensor>("cohesiveTraction") &&
        domain.m_feFaceManager.HasField<int>("ruptureState") &&
        domain.m_feFaceManager.m_cohesiveZone )
    {


      InitiallyRigidCohesiveZone czone;

      Epetra_IntSerialDenseVector  faceLocalDofIndex;
      Epetra_SerialDenseVector     face_rhs;
      Epetra_SerialDenseMatrix     face_matrix;

      Epetra_SerialDenseVector     face_dof_np1;

      // determine ghost elements

      //const iArray1d& isGhost = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();


      // apply cohesive forces
      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
      Array1dT<R1Tensor>& cohesiveForce = domain.m_feNodeManager.GetFieldData<R1Tensor>("cohesiveForce");

      const iArray1d& ruptureState = domain.m_feFaceManager.GetFieldData<int>("ruptureState");

      cohesiveForce = 0.0;



      // begin element loop, skipping ghost elements

      for(localIndex kf = 0 ; kf < domain.m_feFaceManager.m_numFaces ; ++kf)
      {
        if( ruptureState[kf] == 2 && !(childFaceIndex[kf].empty()) )
        {


          const localIndex faceIndex[2] = { kf, childFaceIndex[kf][0] };

          const R1Tensor N[2] = { domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[0] ),
                                  domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, faceIndex[1] )};

          R1Tensor Nbar = N[0];
          Nbar -= N[1];
          Nbar.Normalize();


//          if(isGhost[kf] < 0)
          {

            R1Tensor gap = domain.m_feFaceManager.CalculateGapVector( domain.m_feNodeManager, kf );
            R1Tensor cohesiveTraction;
            R2Tensor cStiffness;

            domain.m_feFaceManager.m_cohesiveZone->UpdateCohesiveZone( kf, gap, Nbar,
                                                                       domain.m_feFaceManager.m_toElementsRelation[faceIndex[0]][0],
                                                                       domain.m_feFaceManager.m_toElementsRelation[faceIndex[1]][0],
                                                                       cohesiveTraction, cStiffness );

            const localIndex faceID = faceIndex[0];
            const realT area = domain.m_feFaceManager.SurfaceArea( domain.m_feNodeManager, faceID );

            R1Tensor cohesiveNodalForce;
            cohesiveNodalForce.cA( area/domain.m_feFaceManager.m_toNodesRelation[kf].size(), cohesiveTraction );

            face_rhs.Scale(0.0);

            for( unsigned int a=0 ; a<domain.m_feFaceManager.m_toNodesRelation[kf].size() ; ++a )
            {
              const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[0]][a];
              const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[faceIndex[1]][a];
              for( int i=1 ; i<2 ; ++i )
              {
                cohesiveForce[localNodeIndex1][i] = cohesiveNodalForce[i];
                cohesiveForce[localNodeIndex2][i] = -cohesiveNodalForce[i];
              }
            }
          }
        }
      }
    }  // end element
}



void LagrangeSmallStrainLinearElastic::CalculateElementStresses( const NodeManagerT& nodeManager,
                                                                      ElementManagerT& elementManager )
{
  const Array1dT<R1Tensor>& disp = nodeManager.GetFieldData<FieldInfo::displacement>();

  for( std::map< ElementManagerT::RegKeyType, ElementRegionT >::iterator iReg = elementManager.m_ElementRegions.begin();
      iReg != elementManager.m_ElementRegions.end(); ++iReg)
  {
    ElementRegionT& elemRegion = iReg->second;

    const FiniteElementBase& fe = *(elemRegion.m_finiteElement);

    realT lambda = 0.0, G = 0.0, TwoG = 0.0;
//    std:cout<<elemRegion.m_mat->NumParameterIndex0();

    for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
    {

      const localIndex paramIndex = elemRegion.m_mat->NumParameterIndex0() > 1 ? element : 0 ;
      const MaterialBaseParameterData& param = *(elemRegion.m_mat->ParameterData(paramIndex));
      G = param.init_shearModulus;
      TwoG = 2*G;
      lambda = param.Lame;
      if( LagrangeSolverBase::m_2dOption==LagrangeSolverBase::PlaneStress )
      {
        lambda = 2*lambda*G / ( lambda + 2*G );
      }

//      if(elem_is_ghost[element] < 0)
      {
        const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[element];

        for( unsigned int q=0 ; q<fe.n_quadrature_points() ; ++q )
        {
          R2SymTensor stress ;
          for( unsigned int a=0 ; a<fe.dofs_per_element() ; ++a )
          {
            const localIndex b = localNodeIndices[a];
            const R1Tensor& u = disp[b];
            const R1Tensor& dNdXb = elemRegion.m_dNdX[element](q,a);

            const realT u0_x_dNdXb0 = u[0]*dNdXb[0];
            const realT u1_x_dNdXb1 = u[1]*dNdXb[1];

            if( dim==3 )
            {
              const realT u2_x_dNdXb2 = u[2]*dNdXb[2];

              stress(0,0) += ( u1_x_dNdXb1 + u2_x_dNdXb2 )*lambda + u0_x_dNdXb0*( TwoG + lambda );
              stress(1,1) += ( u0_x_dNdXb0 + u2_x_dNdXb2 )*lambda + u1_x_dNdXb1*( TwoG + lambda );
              stress(2,2) += ( u0_x_dNdXb0 + u1_x_dNdXb1 )*lambda + u2_x_dNdXb2*( TwoG + lambda );
              stress(1,2) += ( u[2]*dNdXb[1] + u[1]*dNdXb[2] )*G;
              stress(0,2) += ( u[2]*dNdXb[0] + u[0]*dNdXb[2] )*G;
              stress(0,1) += ( u[1]*dNdXb[0] + u[0]*dNdXb[1] )*G;
            }
            else if( dim==2 )
            {
              stress(0,0) += u1_x_dNdXb1*lambda + u0_x_dNdXb0*( TwoG + lambda );
              stress(1,1) += u0_x_dNdXb0*lambda + u1_x_dNdXb1*( TwoG + lambda );
              stress(0,1) += ( u[1]*dNdXb[0] + u[0]*dNdXb[1] )*G;
              if( LagrangeSolverBase::m_2dOption==LagrangeSolverBase::PlaneStress )
              {
                stress(2,2) = 0.0;
              }

            }

          }
          const realT pressure = stress.Trace() / 3.0;
          stress.PlusIdentity( -pressure );
          elemRegion.m_mat->StateData(element,q)->devStress = stress;
          elemRegion.m_mat->StateData(element,q)->pressure = pressure;
        }
      }
    }
  }

  //Take care of thermal stress
  for(localIndex i = 0; i < m_thermalRegionNames.size(); ++i)
  {
    std::map<std::string, ElementRegionT>::iterator it = elementManager.m_ElementRegions.find(m_thermalRegionNames[i]);

    ElementRegionT& elemRegion = it -> second;
    const FiniteElementBase& fe = *(elemRegion.m_finiteElement);
    rArray1d& temperature = elemRegion.GetFieldData<realT>("temperature");
    rArray1d& CTE = elemRegion.GetFieldData<realT>("linearCTE");

    for(localIndex element = 0; element < elemRegion.m_numElems; ++element)
    {
      const localIndex paramIndex = elemRegion.m_mat->NumParameterIndex0() > 1 ? element : 0 ;
      const MaterialBaseParameterData& matParams = *(elemRegion.m_mat->ParameterData(paramIndex));
      const realT K =  matParams.Lame + matParams.init_shearModulus * 2.0 / 3.0;

      const realT tStress = (temperature[element] - m_refTemperature) * CTE [element] * K * 3;

      const localIndex* const localNodeIndices = elemRegion.m_toNodesRelation[element];

      for( unsigned int q=0 ; q<fe.n_quadrature_points() ; ++q )
      {
        elemRegion.m_mat->StateData(element,q)->pressure -= tStress;
      }
    }
  }
}


void LagrangeSmallStrainLinearElastic::UpdateContactDataStructures( PhysicalDomainT& domain, const realT time)
{
  if(dim==3)
  {
    std::cout<<"Calculating neighbor lists"<<std::endl;
    const bool resort = domain.m_externalFaces.RecalculateNeighborList(
        domain.m_feNodeManager, domain.m_discreteElementSurfaceNodes,
        domain.m_discreteElementManager, 0.0, false, true);
    if (resort)
    {
      std::cout<<"Updating contact Manager"<<std::endl;
      domain.m_contactManager.Update(domain.m_externalFaces.m_neighborList);
    }
    {
      Array1dT<Array1dT<R1Tensor> > xs;
      xs.resize(domain.m_externalFaces.DataLengths());
      std::cout<<"Updating contact data structures"<<std::endl;
      domain.m_externalFaces.UpdateGeometricContactProperties(0.0, domain, xs, false);
    }
    if(m_numerics.m_useNewtonSolve == false || time==0)
    {
      const iArray1d& activeC = domain.m_contactManager.GetFieldData<int>("active");
      iArray1d& activeInit = domain.m_contactManager.GetFieldData<int>("activeInit");

      activeInit = activeC;
    }
  }

//  std::cout<<"Calculating neighbor lists"<<std::endl;
//  const bool resort = domain.m_externalFaces.RecalculateNeighborList(
//      domain.m_feNodeManager, domain.m_discreteElementSurfaceNodes,
//      domain.m_discreteElementManager, 0.0, false);
//  if (resort)
//  {
//    std::cout<<"Updating contact Manager"<<std::endl;
//    domain.m_contactManager.Update(domain.m_externalFaces.m_neighborList);
//  }
//  {
//    Array1dT<Array1dT<R1Tensor> > xs;
//    xs.resize(domain.m_externalFaces.DataLengths());
//    std::cout<<"Updating contact data structures"<<std::endl;
//    domain.m_externalFaces.UpdateGeometricContactProperties(0.0, domain, xs, false);
//  }

}

//@annavarapusr1: Interface stiffness specific functions

void LagrangeSmallStrainLinearElastic::InsertGlobalIndices( PhysicalDomainT& domain)
{
  // @annavarapusr1: Specify sparsity indices for the global matrix. In addition to the entries specified above, the contact face-pair
  // introduces non-zero entries for the Face-Sibling cross indices and these need to be specified here.
  // If this is not done, the local matrix entries do not get assembled in the global matrix.

  unsigned int totContFaces;

//  if(dim==2 or dim==3)
//  {
//    totContFaces = domain.m_feFaceManager.m_numFaces;
//  }

  if(dim==2)
  {
    totContFaces = domain.m_feFaceManager.m_numFaces;
  }
  else
  {
    totContFaces = domain.m_contactManager.DataLengths();
  }

//  unsigned int totContFaces = domain.m_contactManager.DataLengths();

  for (localIndex iContFace =  0; iContFace < totContFaces; iContFace++)
  {
    bool contActiv = false;

//    if(dim==2 or dim==3)
//    {
//      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");
//
//      if (!(childFaceIndex[iContFace].empty()))
//      {
//        contActiv = true;
//      }
//    }

    if(dim==2)
    {
      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");

      if (!(childFaceIndex[iContFace].empty()))
      {
        contActiv = true;
      }
    }
    else
    {
      const iArray1d& active = domain.m_contactManager.GetFieldData<int>("activeInit");

      if(active[iContFace]!=0)
      {
        contActiv = true;
      }
    }

    if(contActiv)
    {
      iArray1d rowDofIndex;
      iArray1d colDofIndex;

      const bool nitsche_active = domain.m_contactManager.m_nitsche_active;
      const bool nitsche_symmetry_active = domain.m_contactManager.m_nitsche_symmetry_active;

      GetRowDofIndexColDofIndex(iContFace, domain, nitsche_active, rowDofIndex, colDofIndex);

      m_sparsity->InsertGlobalIndices(rowDofIndex.size(), &rowDofIndex.front(),
                                      colDofIndex.size(), &colDofIndex.front());
      if(nitsche_symmetry_active)
      {
        m_sparsity->InsertGlobalIndices(colDofIndex.size(), &colDofIndex.front(),
                                        rowDofIndex.size(), &rowDofIndex.front());
      }
    }
  }
}


void LagrangeSmallStrainLinearElastic::GetContactStiffnessContribution( PhysicalDomainT& domain)
{
  rArray1d gauss(2);
  gauss[0] = -1 / sqrt(3);
  gauss[1] = 1 / sqrt(3);

  const bool nitsche_active = domain.m_contactManager.m_nitsche_active;
  const bool nitsche_symmetry_active = domain.m_contactManager.m_nitsche_symmetry_active;

  Epetra_IntSerialDenseVector rowDofIndexPenalty;
  Epetra_IntSerialDenseVector colDofIndexPenalty;

  Epetra_IntSerialDenseVector rowDofIndexNitsche;
  Epetra_IntSerialDenseVector colDofIndexNitsche;

  Epetra_SerialDenseMatrix faceMatrixPenalty;

  Epetra_SerialDenseMatrix faceMatrixNitsche;
  Epetra_SerialDenseMatrix faceMatrixNitscheSym;

  Epetra_SerialDenseVector faceRhs;

  unsigned int totContFaces;

//  if(dim==2 or dim==3)
//  {
//    totContFaces = domain.m_feFaceManager.m_numFaces;
//  }

  if(dim==2)
  {
    totContFaces = domain.m_feFaceManager.m_numFaces;
  }
  else
  {
    totContFaces = domain.m_contactManager.DataLengths();
  }

//    unsigned int totContFaces = domain.m_contactManager.DataLengths();

  for (localIndex iContFace =  0; iContFace < totContFaces; iContFace++)
  {
    bool contActiv = false;
    localIndex kf1, kf2;

//
//    if(dim==2 or dim==3)
//    {
//      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");
//
//      if (!(childFaceIndex[iContFace].empty()))
//      {
//        contActiv = true;
////        kf1 = childFaceIndex[iContFace][0];
////        kf2 = childFaceIndex[iContFace][1];
//        if(dim==2)
//        {
//          kf1 = childFaceIndex[iContFace][0];
//          kf2 = childFaceIndex[iContFace][1];
//        }
//        if(dim==3)
//        {
//          kf1 = iContFace;
//          kf2 = childFaceIndex[iContFace][0];
//        }
//      }
//    }

    if(dim==2)
    {
      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");

      if (!(childFaceIndex[iContFace].empty()))
      {
        contActiv = true;
        kf1 = childFaceIndex[iContFace][0];
        kf2 = childFaceIndex[iContFace][1];
      }
    }
    else
    {
      const iArray1d& active = domain.m_contactManager.GetFieldData<int>("activeInit");
      const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
      const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");

      if(active[iContFace]!=0)
      {
        contActiv = true;
        bool fe;
        kf1 = domain.m_externalFaces.FaceIndex(f1[iContFace], fe);
        kf2 = domain.m_externalFaces.FaceIndex(f2[iContFace], fe);
      }
    }

//    const iArray1d& active = domain.m_contactManager.GetFieldData<int>("active");
//    const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
//    const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
//
//    if(active[iContFace]!=0)
//    {
//      contActiv = true;
//      bool fe;
//      kf1 = domain.m_externalFaces.FaceIndex(f1[iContFace], fe);
//      kf2 = domain.m_externalFaces.FaceIndex(f2[iContFace], fe);
//    }

    if(contActiv)
    {
      const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();
      const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
      const ElementRegionT* elemRegionPar = ftoe[kf1][0].first;
      const ElementRegionT* elemRegionSib = ftoe[kf2][0].first;

      const localIndex numNodesParentEle = elemRegionPar->m_numNodesPerElem;
      const localIndex numNodesSiblingEle = elemRegionSib->m_numNodesPerElem;
      const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

      rArray1d xe, psi, eta;
      lArray1d localParentFaceNodes;
      rArray2d P_par(dim,dim), P_sib(dim,dim);
      rArray2d alphaPar(dim,dim), alphaSib(dim,dim);

      std::string s_previous_base = "uJumpPl_Previous_gp", s_current_base  = "uJumpPl_Current_gp";
      iArray1d stickGP_par(numNodes), stickGP_sib(numNodes);
      iArray1d openingGP_par(numNodes), openingGP_sib(numNodes);
      realT tol = domain.m_contactManager.m_traction_n_tol;

      // Declaration of Nitsche specific variables
      iArray1d FaceConnToElemConnPar(numNodes), FaceConnToElemConnSib(numNodes);
      rArray2d normalVoigtPar(dim,0.5*dim*(dim+1)), normalVoigtSib(dim,0.5*dim*(dim+1));
      rArray2d DPar(0.5*dim*(dim+1),0.5*dim*(dim+1)), DSib(0.5*dim*(dim+1),0.5*dim*(dim+1));
      rArray2d BPar(0.5*dim*(dim+1),dim*numNodesParentEle), BSib(0.5*dim*(dim+1),dim*numNodesSiblingEle);
      rArray2d nParDotDParBPar(dim,dim*numNodesParentEle), nSibDotDSibBSib(dim,dim*numNodesSiblingEle);
      rArray2d nParDotDParBParPermute(dim,dim*numNodesParentEle), nSibDotDSibBSibPermute(dim,dim*numNodesSiblingEle);
      realT gamPar, gamSib;

      faceMatrixPenalty.Reshape(2 * dim * numNodes, 2 * dim * numNodes); faceMatrixPenalty.Scale(0.0);
      faceRhs.Resize(2 * dim * numNodes); faceRhs.Scale(0.0);

      // Remark: Seems like this logic assumes faces are aligned. Should probably be changed for non-conforming meshes
      GetParentFaceNodesAndCoordsAndInterfaceGaussPoints(iContFace, domain, gauss, localParentFaceNodes, xe, psi, eta);

      // Returns global indices for the Face-pair that is under contact and face IDs for split faces
      GetRowDofIndexColDofIndex(iContFace, domain, false, rowDofIndexPenalty, colDofIndexPenalty);

      GetTransformationTensor(domain, kf1, P_par);
      GetTransformationTensor(domain, kf2, P_sib);

      GetInitialAlphaTensor(domain, kf1, P_par, alphaPar);
      GetInitialAlphaTensor(domain, kf2, P_sib, alphaSib);

      if(nitsche_active)
      {
        faceMatrixNitsche.Reshape(2 * dim * numNodes, dim *(numNodesParentEle + numNodesSiblingEle)); faceMatrixNitsche.Scale(0.0);

        faceMatrixNitscheSym.Reshape(dim *(numNodesParentEle + numNodesSiblingEle), 2 * dim * numNodes); faceMatrixNitscheSym.Scale(0.0);

        GetRowDofIndexColDofIndex(iContFace, domain, true, rowDofIndexNitsche, colDofIndexNitsche);

        GetFaceConnToElemConnMap(domain, kf1, false, FaceConnToElemConnPar);
        GetFaceConnToElemConnMap(domain, kf2, true, FaceConnToElemConnSib);

        GetNormalVoigt(domain, kf1, normalVoigtPar);
        GetNormalVoigt(domain, kf2, normalVoigtSib);

        GetElasticityTensorVoigt(domain, kf1, DPar);
        GetElasticityTensorVoigt(domain, kf2, DSib);

        const rArray1d& nitscheGamma = domain.m_externalFaces.GetFieldData<realT>("nitscheGamma");

        gamPar = nitscheGamma(faceToExternalFaceMap(kf1));
        gamSib = nitscheGamma(faceToExternalFaceMap(kf2));
        if(fabs(((gamPar+gamSib)-1))>1e-12)
        {
          std::cout<<"WARNING:: Nitsche weights don't sum to unity!"<<std::endl;
        }

        if(elemRegionPar->m_elementGeometryID.compare(0,4,"STRI") ==0 or elemRegionPar->m_elementGeometryID.compare(0,4,"C3D4") == 0)
        {
          GetShapeFunctionDerivativeMatrixConstStr(domain, FaceConnToElemConnPar, kf1, BPar);
          GetShapeFunctionDerivativeMatrixConstStr(domain, FaceConnToElemConnSib, kf2, BSib);

          // GetDBDotN includes the Nitsche weights so should be renamed to GetGammaDBDotN?
          GetDBDotN(DPar, BPar, normalVoigtPar, numNodesParentEle, gamPar, domain, P_par, nParDotDParBPar);
          GetDBDotN(DSib, BSib, normalVoigtSib, numNodesSiblingEle, gamSib, domain, P_sib, nSibDotDSibBSib);

          GetPermutedDBDotN(domain, nParDotDParBPar, FaceConnToElemConnPar, kf1, nParDotDParBParPermute);
          GetPermutedDBDotN(domain, nSibDotDSibBSib, FaceConnToElemConnSib, kf2, nSibDotDSibBSibPermute);
        }
      }

      //RHS contribution from penalty and Nitsche terms
      for (localIndex a = 0; a < numNodes; ++a)
      {
        localIndex edgeLocalIndexRow, sibLocalIndexRow;

        GetLocalIndexOnInterface(iContFace, a, domain, localParentFaceNodes, edgeLocalIndexRow, sibLocalIndexRow);

        const int aDof = 2 * a * dim;

        for (int i = 0; i < dim; ++i)
        {
          for (localIndex iGp = 0; iGp < numNodes; ++iGp)
          {
            realT jcob_sub = 0;
            rArray2d N_sub(dim,dim*numNodes);
            rArray1d tracPar(dim), tracSib(dim);

            stickGP_par(iGp) = 1; stickGP_sib(iGp) = 1;
            openingGP_par(iGp) = 0; openingGP_sib(iGp) = 0;

            GetJacobianAndShapeFunctionsOnInterface(numNodes, xe, psi(iGp), eta(iGp), jcob_sub, N_sub);

            if(nitsche_active)
            {
              if(elemRegionPar->m_elementGeometryID.compare(0,4,"CPE4") ==0 or elemRegionPar->m_elementGeometryID.compare(0,4,"C3D8") ==0 )
              {
                GetShapeFunctionDerivativeMatrixQuadHex(domain, FaceConnToElemConnPar, kf1, psi(iGp), eta(iGp), BPar);
                GetShapeFunctionDerivativeMatrixQuadHex(domain, FaceConnToElemConnSib, kf2, psi(iGp), eta(iGp), BSib);

                // GetDBDotN includes the Nitsche weights so should be renamed to GetGammaDBDotN?
                GetDBDotN(DPar, BPar, normalVoigtPar, numNodesParentEle, gamPar, domain, P_par,  nParDotDParBPar);
                GetDBDotN(DSib, BSib, normalVoigtSib, numNodesSiblingEle, gamSib, domain, P_sib, nSibDotDSibBSib);

                GetPermutedDBDotN(domain, nParDotDParBPar, FaceConnToElemConnPar, kf1, nParDotDParBParPermute);
                GetPermutedDBDotN(domain, nSibDotDSibBSib, FaceConnToElemConnSib, kf2, nSibDotDSibBSibPermute);

              }
            }

            if(domain.m_contactManager.m_sliding_law == 0)
            {
              GetTractionInterface(domain, N_sub, alphaPar, alphaSib, nParDotDParBParPermute, nSibDotDSibBSibPermute, localParentFaceNodes, kf1, kf2, iContFace, tracPar, tracSib);

              GetUpdatedTractionsOpeningMode(P_par, iGp, tol, openingGP_par, tracPar);
              GetUpdatedTractionsOpeningMode(P_sib, iGp, tol, openingGP_sib, tracSib);
            }
            else
            {
              GetTractionInterfaceSliding(N_sub, alphaPar, alphaSib, P_par, P_sib, nParDotDParBParPermute, nSibDotDSibBSibPermute, localParentFaceNodes, kf1, kf2, iContFace,
                                          iGp, s_previous_base, s_current_base, domain, stickGP_par, stickGP_sib, openingGP_par, openingGP_sib, tracPar, tracSib);
            }

            faceRhs(aDof+i)     -= jcob_sub*N_sub(i,dim*edgeLocalIndexRow+i)*tracPar(i);
            faceRhs(aDof+dim+i) -= jcob_sub*N_sub(i,dim*sibLocalIndexRow+i)*tracSib(i);
          }
        }
      }

      // Penalty/Nitsche stiffness for the face-pair
      for (localIndex a = 0; a < numNodes; ++a)
      {
        localIndex edgeLocalIndexRow, sibLocalIndexRow;

        GetLocalIndexOnInterface(iContFace, a, domain, localParentFaceNodes, edgeLocalIndexRow, sibLocalIndexRow);

        const int aDof = 2 * a * dim;

        for (localIndex b = 0; b < numNodesParentEle; ++b)
        {
          localIndex edgeLocalIndexCol, sibLocalIndexCol;
          if(b<numNodes)
          {
            GetLocalIndexOnInterface(iContFace, b, domain, localParentFaceNodes, edgeLocalIndexCol, sibLocalIndexCol);
          }
          const int bDof = 2 * b * dim;

          for (int i = 0; i < dim; ++i)
          {
            for (int j = 0; j < dim; ++j)
            {
              for (localIndex iGp = 0; iGp < numNodes; ++iGp)
              {
                realT jcob_sub = 0;
                rArray2d N_sub(dim,dim*numNodes), alphaParTimesNSub(dim,dim*numNodes), alphaSibTimesNSub(dim,dim*numNodes);

                GetJacobianAndShapeFunctionsOnInterface(numNodes, xe, psi(iGp), eta(iGp), jcob_sub, N_sub);

                alphaPar = 0; alphaSib = 0; nParDotDParBPar = 0; nSibDotDSibBSib = 0;

                GetTrialTangentModulus(domain, kf1, P_par, FaceConnToElemConnPar, psi(iGp), eta(iGp), DPar, normalVoigtPar, gamPar, alphaPar, nParDotDParBPar);
                GetTrialTangentModulus(domain, kf2, P_sib, FaceConnToElemConnSib, psi(iGp), eta(iGp), DSib, normalVoigtSib, gamSib, alphaSib, nSibDotDSibBSib);

                nParDotDParBParPermute = 0; nSibDotDSibBSibPermute = 0;
                GetPermutedDBDotN(domain, nParDotDParBPar, FaceConnToElemConnPar, kf1, nParDotDParBParPermute);
                GetPermutedDBDotN(domain, nSibDotDSibBSib, FaceConnToElemConnSib, kf2, nSibDotDSibBSibPermute);

                // tracPar and tracSib are only needed here to calculate slip directions in 3D.
                rArray1d tracPar(dim), tracSib(dim);
                if(domain.m_contactManager.m_sliding_law==1 or domain.m_contactManager.m_sliding_law==2)
                {
                  std::stringstream ss; ss << iGp;
                  std::string FieldNamePrevious = s_previous_base + ss.str(), FieldNameCurrent  = s_current_base + ss.str();
                  GetTractionInterfaceSlidingTrial(domain, N_sub, alphaPar, alphaSib, P_par, P_sib, nParDotDParBParPermute, nSibDotDSibBSibPermute, localParentFaceNodes,
                                                   kf1, kf2, iContFace, FieldNamePrevious, tracPar, tracSib);
                }

                UpdateTangentModulusForOpeningAndSliding(domain, P_par, FaceConnToElemConnPar, kf1, tracPar, openingGP_par(iGp), stickGP_par(iGp), true,
                                                         nParDotDParBPar, alphaPar, nParDotDParBParPermute);
                UpdateTangentModulusForOpeningAndSliding(domain, P_sib, FaceConnToElemConnSib, kf2, tracSib, openingGP_sib(iGp), stickGP_sib(iGp), false,
                                                         nSibDotDSibBSib, alphaSib, nSibDotDSibBSibPermute);

                if(b<numNodes)
                {
                  // Pre-Multiply interface tensor with shape-function matrix: \alpha_XYZ*N

                  MultiplyArray(alphaPar,N_sub,alphaParTimesNSub);
                  MultiplyArray(alphaSib,N_sub,alphaSibTimesNSub);

                  // Get face_matrix = N'*alpha_XYZ*N
                  for (int iCol = 0; iCol < dim; ++iCol)
                  {
                    faceMatrixPenalty(aDof + i, bDof + j)             -= N_sub(i, dim * edgeLocalIndexRow + iCol) * jcob_sub * alphaParTimesNSub(iCol, dim * edgeLocalIndexCol + j);
                    faceMatrixPenalty(aDof + i, bDof + dim + j)       += N_sub(i, dim * edgeLocalIndexRow + iCol) * jcob_sub * alphaSibTimesNSub(iCol, dim * sibLocalIndexCol + j);
                    faceMatrixPenalty(aDof + dim + i, bDof + j)       += N_sub(i, dim * sibLocalIndexRow + iCol) * jcob_sub * alphaParTimesNSub(iCol, dim * edgeLocalIndexCol + j);
                    faceMatrixPenalty(aDof + dim + i, bDof + dim + j) -= N_sub(i, dim * sibLocalIndexRow + iCol) * jcob_sub * alphaSibTimesNSub(iCol, dim * sibLocalIndexCol + j);
                  }
                }

                if(nitsche_active)
                {
                  for (int iCol = 0; iCol<dim; ++iCol)
                  {
                    faceMatrixNitsche(aDof+i, bDof+j)          += N_sub(i,dim*edgeLocalIndexRow+iCol)*jcob_sub*nParDotDParBParPermute(iCol,dim*b+j);
                    faceMatrixNitsche(aDof+i, bDof+dim+j)      -= N_sub(i,dim*edgeLocalIndexRow+iCol)*jcob_sub*nSibDotDSibBSibPermute(iCol,dim*b+j);
                    faceMatrixNitsche(aDof+dim+i, bDof+j)      -= N_sub(i,dim*sibLocalIndexRow+iCol)*jcob_sub*nParDotDParBParPermute(iCol,dim*b+j);
                    faceMatrixNitsche(aDof+dim+i, bDof+dim+j)  += N_sub(i,dim*sibLocalIndexRow+iCol)*jcob_sub*nSibDotDSibBSibPermute(iCol,dim*b+j);
                  }

                  if(nitsche_symmetry_active)
                  {
                    for (int iCol = 0; iCol<dim; ++iCol)
                    {
                      faceMatrixNitscheSym(bDof+j, aDof+i)          += N_sub(i,dim*edgeLocalIndexRow+iCol)*jcob_sub*nParDotDParBParPermute(iCol,dim*b+j);
                      faceMatrixNitscheSym(bDof+dim+j, aDof+i)      -= N_sub(i,dim*edgeLocalIndexRow+iCol)*jcob_sub*nSibDotDSibBSibPermute(iCol,dim*b+j);
                      faceMatrixNitscheSym(bDof+j, aDof+dim+i)      -= N_sub(i,dim*sibLocalIndexRow+iCol)*jcob_sub*nParDotDParBParPermute(iCol,dim*b+j);
                      faceMatrixNitscheSym(bDof+dim+j, aDof+dim+i)  += N_sub(i,dim*sibLocalIndexRow+iCol)*jcob_sub*nSibDotDSibBSibPermute(iCol,dim*b+j);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // Assemble into Global matrices
      m_matrix->SumIntoGlobalValues(rowDofIndexPenalty, colDofIndexPenalty, faceMatrixPenalty);
      if(nitsche_active)
      {
        m_matrix->SumIntoGlobalValues(rowDofIndexNitsche, colDofIndexNitsche, faceMatrixNitsche);

        if(nitsche_symmetry_active)
        {
          m_matrix->SumIntoGlobalValues(colDofIndexNitsche, rowDofIndexNitsche, faceMatrixNitscheSym);
        }
      }
      if(nitsche_active)
      {
        m_rhs->SumIntoGlobalValues(rowDofIndexNitsche, faceRhs);
      }
      else
      {
        m_rhs->SumIntoGlobalValues(rowDofIndexPenalty, faceRhs);
      }
    }
  }
}

// @annavarapusr1: Function definitions for interface contact specific quantities

void LagrangeSmallStrainLinearElastic::GetRowDofIndexColDofIndex( const localIndex i,
                                                                  const PhysicalDomainT& domain,
                                                                  const bool nitsche_active,
                                                                  iArray1d& rowDofIndex,
                                                                  iArray1d& colDofIndex)
{
  localIndex kf1, kf2;

//  if(dim==2 or dim==3)
//  {
//    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");
//
//    if (!(childFaceIndex[i].empty()))
//    {
////      kf1 = childFaceIndex[i][0];
////      kf2 = childFaceIndex[i][1];
//      if(dim==2)
//      {
//        kf1 = childFaceIndex[i][0];
//        kf2 = childFaceIndex[i][1];
//      }
//      if(dim==3)
//      {
//        kf1 = i;
//        kf2 = childFaceIndex[i][0];
//      }
//    }
//  }

  if(dim==2)
  {
    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");

    if (!(childFaceIndex[i].empty()))
    {
      kf1 = childFaceIndex[i][0];
      kf2 = childFaceIndex[i][1];
    }
  }
  else
  {
    const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
    const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
    bool fe;
    kf1 = domain.m_externalFaces.FaceIndex(f1[i], fe);
    kf2 = domain.m_externalFaces.FaceIndex(f2[i], fe);
  }

//  const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
//  const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
//  bool fe;
//  const localIndex kf1 = domain.m_externalFaces.FaceIndex(f1[i], fe);
//  const localIndex kf2 = domain.m_externalFaces.FaceIndex(f2[i], fe);

  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();
  const iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);

  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegionPar = ftoe[kf1][0].first;
  const ElementRegionT* elemRegionSib = ftoe[kf2][0].first;

  const localIndex numNodesParentEle = elemRegionPar->m_numNodesPerElem;
  const localIndex numNodesSiblingEle = elemRegionSib->m_numNodesPerElem;

  rowDofIndex.resize(2*dim*numNodes);
  iArray1d ParNodIDNotOnFace(numNodesParentEle-numNodes), SibNodIDNotOnFace(numNodesSiblingEle-numNodes);

  if(nitsche_active)
  {
    colDofIndex.resize(dim*(numNodesParentEle+numNodesSiblingEle));

    GetNodIndicesNotOnFace(domain, kf1, true, ParNodIDNotOnFace);
    GetNodIndicesNotOnFace(domain, kf2, true, SibNodIDNotOnFace);
  }
  else
  {
    colDofIndex.resize(2*dim*numNodes);
  }

  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    const int aDof = 2*a*dim;
    const localIndex aa = a == 0 ? a : numNodes - a;
    const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[kf1][a];
    const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[kf2][aa];

    for( int iDim=0 ; iDim<dim ; ++iDim )
    {
      rowDofIndex[aDof+iDim]       = dim*trilinos_index[localNodeIndex1]+iDim;
      rowDofIndex[aDof+dim+iDim]   = dim*trilinos_index[localNodeIndex2]+iDim;

      colDofIndex[aDof+iDim]       = dim*trilinos_index[localNodeIndex1]+iDim;
      colDofIndex[aDof+dim+iDim]   = dim*trilinos_index[localNodeIndex2]+iDim;
    }
  }

  if(nitsche_active)
  {
    for( localIndex a=numNodes ; a<numNodesParentEle ; ++a )
    {
      const int aDof = 2*a*dim;
      int iNodNotOnFace = a-numNodes;
      for( int jDim=0 ; jDim<dim ; ++jDim )
      {
        colDofIndex[aDof+jDim] = dim*trilinos_index[ParNodIDNotOnFace(iNodNotOnFace)]+jDim;
      }
    }
    for( localIndex a=numNodes ; a<numNodesSiblingEle ; ++a )
    {
      const int aDof = 2*a*dim;
      int iNodNotOnFace = a-numNodes;
      for( int jDim=0 ; jDim<dim ; ++jDim )
      {
        colDofIndex[aDof+dim+jDim] = dim*trilinos_index[SibNodIDNotOnFace(iNodNotOnFace)]+jDim;
      }
    }
  }
}

void LagrangeSmallStrainLinearElastic::GetRowDofIndexColDofIndex( const localIndex i,
                                                                  const PhysicalDomainT& domain,
                                                                  const bool nitsche_active,
                                                                  Epetra_IntSerialDenseVector& rowDofIndex,
                                                                  Epetra_IntSerialDenseVector& colDofIndex)
{
  localIndex kf1, kf2;
//
//  if(dim==2 or dim==3)
//  {
//    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");
//
//    if (!(childFaceIndex[i].empty()))
//    {
////      kf1 = childFaceIndex[i][0];
////      kf2 = childFaceIndex[i][1];
//      if(dim==2)
//      {
//        kf1 = childFaceIndex[i][0];
//        kf2 = childFaceIndex[i][1];
//      }
//      if(dim==3)
//      {
//        kf1 = i;
//        kf2 = childFaceIndex[i][0];
//      }
//    }
//  }

  if(dim==2)
  {
    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");

    if (!(childFaceIndex[i].empty()))
    {
      kf1 = childFaceIndex[i][0];
      kf2 = childFaceIndex[i][1];
    }
  }
  else
  {
    const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
    const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
    bool fe;
    kf1 = domain.m_externalFaces.FaceIndex(f1[i], fe);
    kf2 = domain.m_externalFaces.FaceIndex(f2[i], fe);
  }

//  const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
//  const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
//  bool fe;
//  const localIndex kf1 = domain.m_externalFaces.FaceIndex(f1[i], fe);
//  const localIndex kf2 = domain.m_externalFaces.FaceIndex(f2[i], fe);

  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();
  const iArray1d& trilinos_index = domain.m_feNodeManager.GetFieldData<int>(m_trilinosIndexStr);

  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegionPar = ftoe[kf1][0].first;
  const ElementRegionT* elemRegionSib = ftoe[kf2][0].first;

  const localIndex numNodesParentEle = elemRegionPar->m_numNodesPerElem;
  const localIndex numNodesSiblingEle = elemRegionSib->m_numNodesPerElem;

  rowDofIndex.Resize(2*dim*numNodes);
  iArray1d ParNodIDNotOnFace(numNodesParentEle-numNodes), SibNodIDNotOnFace(numNodesSiblingEle-numNodes);

  if(nitsche_active)
  {
    colDofIndex.Resize(dim*(numNodesParentEle+numNodesSiblingEle));

    GetNodIndicesNotOnFace(domain, kf1, true, ParNodIDNotOnFace);
    GetNodIndicesNotOnFace(domain, kf2, true, SibNodIDNotOnFace);
  }
  else
  {
    colDofIndex.Resize(2*dim*numNodes);
  }

  for( localIndex a=0 ; a<numNodes ; ++a )
  {
    const int aDof = 2*a*dim;
    const localIndex aa = a == 0 ? a : numNodes - a;
    const localIndex localNodeIndex1 = domain.m_feFaceManager.m_toNodesRelation[kf1][a];
    const localIndex localNodeIndex2 = domain.m_feFaceManager.m_toNodesRelation[kf2][aa];

    for( int iDim=0 ; iDim<dim ; ++iDim )
    {
      rowDofIndex[aDof+iDim]       = dim*trilinos_index[localNodeIndex1]+iDim;
      rowDofIndex[aDof+dim+iDim]   = dim*trilinos_index[localNodeIndex2]+iDim;

      colDofIndex[aDof+iDim]       = dim*trilinos_index[localNodeIndex1]+iDim;
      colDofIndex[aDof+dim+iDim]   = dim*trilinos_index[localNodeIndex2]+iDim;
    }
  }

  if(nitsche_active)
  {
    for( localIndex a=numNodes ; a<numNodesParentEle ; ++a )
    {
      const int aDof = 2*a*dim;
      int iNodNotOnFace = a-numNodes;
      for( int jDim=0 ; jDim<dim ; ++jDim )
      {
        colDofIndex[aDof+jDim] = dim*trilinos_index[ParNodIDNotOnFace(iNodNotOnFace)]+jDim;
      }
    }
    for( localIndex a=numNodes ; a<numNodesSiblingEle ; ++a )
    {
      const int aDof = 2*a*dim;
      int iNodNotOnFace = a-numNodes;
      for( int jDim=0 ; jDim<dim ; ++jDim )
      {
        colDofIndex[aDof+dim+jDim] = dim*trilinos_index[SibNodIDNotOnFace(iNodNotOnFace)]+jDim;
      }
    }
  }
}

void LagrangeSmallStrainLinearElastic::GetNodIndicesNotOnFace( const PhysicalDomainT& domain,
                                                               const localIndex kf,
                                                               const bool globFlag,
                                                               iArray1d& NodIDNotOnFace)
{
  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf].size();
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[kf][0].first;
  const localIndex numNodesParentEle = elemRegion->m_numNodesPerElem;
  const localIndex EleID = domain.m_feFaceManager.m_toElementsRelation[kf][0].second;

  lArray1d locConn(numNodesParentEle);
  for (localIndex iNod=0; iNod<numNodesParentEle; ++iNod)
  {
    locConn(iNod) = elemRegion->m_toNodesRelation[EleID][iNod];
  }

  int count = 0;
  for(localIndex iNod=0; iNod<numNodesParentEle; ++iNod)
  {
    bool NodIDFoundFlag = 0;
    for(localIndex jNod=0; jNod<numNodes; ++jNod)
    {
      if(locConn(iNod)==domain.m_feFaceManager.m_toNodesRelation[kf][jNod])
      {
        NodIDFoundFlag = 1;
      }
    }
    if(NodIDFoundFlag==0)
    {
      if(globFlag)
        NodIDNotOnFace(count) = locConn(iNod);
      else
        NodIDNotOnFace(count) = iNod;
      count = count + 1;
    }
  }
}
// GetLocalIndexOnInterface: Gets corresponding local ID for node "a" in childFaceIndex["kf"][0] and childFaceIndex[kf][1]

void LagrangeSmallStrainLinearElastic::GetLocalIndexOnInterface( const localIndex iC,
                                                                 const localIndex a,
                                                                 const PhysicalDomainT& domain,
                                                                 const lArray1d& localParentFaceNodes,
                                                                 localIndex& edgeLocalIndex,
                                                                 localIndex& sibLocalIndex)
{
  const OneToOneRelation& parentNodeIndex = domain.m_feNodeManager.GetOneToOneMap("parentIndex");
  const localIndex& numTotNodesDomain = domain.m_feNodeManager.m_numNodes;

  localIndex kf1, kf2;

//  if(dim==2 or dim==3)
//  {
//    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");
//
//    if (!(childFaceIndex[iC].empty()))
//    {
////      kf1 = childFaceIndex[iC][0];
////      kf2 = childFaceIndex[iC][1];
//      if(dim==2)
//      {
//        kf1 = childFaceIndex[iC][0];
//        kf2 = childFaceIndex[iC][1];
//      }
//      if(dim==3)
//      {
//        kf1 = iC;
//        kf2 = childFaceIndex[iC][0];
//      }
//    }
//  }

  if(dim==2)
  {
    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");

    if (!(childFaceIndex[iC].empty()))
    {
      kf1 = childFaceIndex[iC][0];
      kf2 = childFaceIndex[iC][1];
    }
  }
  else
  {
    const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
    const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
    bool fe;
    kf1 = domain.m_externalFaces.FaceIndex(f1[iC], fe);
    kf2 = domain.m_externalFaces.FaceIndex(f2[iC], fe);
  }

//  const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
//  const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
//  bool fe;
//  const localIndex kf1 = domain.m_externalFaces.FaceIndex(f1[iC], fe);
//  const localIndex kf2 = domain.m_externalFaces.FaceIndex(f2[iC], fe);

  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();

  const localIndex aa = a == 0 ? a : numNodes - a;
  const localIndex rowLocalToGlobalFaceIndex = domain.m_feFaceManager.m_toNodesRelation[kf1][a]; // Corresponding global index of local index a and childFaceIndex[kf][0]
  const localIndex rowLocalToGlobalSibIndex = domain.m_feFaceManager.m_toNodesRelation[kf2][aa];

  localIndex parentID, sibID;
  bool sibLocalIndexNotFound = 1, edgeLocalIndexNotFound = 1;

  if (parentNodeIndex[rowLocalToGlobalFaceIndex] < numTotNodesDomain)
    parentID = parentNodeIndex[rowLocalToGlobalFaceIndex];
  else
    parentID = rowLocalToGlobalFaceIndex;
  // If parentID is itself a child, then check for the parent of parentID till you find the original parent
  while (parentNodeIndex[parentID] < numTotNodesDomain)
  {
    parentID = parentNodeIndex[parentID];
  }

  if (parentNodeIndex[rowLocalToGlobalSibIndex] < numTotNodesDomain)
    sibID = parentNodeIndex[rowLocalToGlobalSibIndex];
  else
    sibID = rowLocalToGlobalSibIndex;
  while (parentNodeIndex[sibID] < numTotNodesDomain)
  {
    sibID = parentNodeIndex[sibID];
  }

  for (localIndex i = 0; i < localParentFaceNodes.size(); ++i)
  {
    if (localParentFaceNodes[i] == parentID)
    {
      edgeLocalIndex = i;
      edgeLocalIndexNotFound = 0; // Search for the global index of the parent of localIndex1Row in localParentFaceNodes
    }
    if (localParentFaceNodes[i] == sibID)
    {
      sibLocalIndex = i;
      sibLocalIndexNotFound = 0;
    }
    if (i == (localParentFaceNodes.size() - 1))
    {
      if (sibLocalIndexNotFound == 1)
      {
        std::cout<<"For face-pairs "<<kf1<<" and "<<kf2<<std::endl;
        throw GPException("SibLocalIndexNotFound: Sibling index not found.");
      }
      if (edgeLocalIndexNotFound == 1)
        throw GPException("edgeLocalIndexNotFound: Edge index not found.");
    }
  }
}

// GetJacobianAndShapeFunctionsOnInterface: Gets Jacobian for the surface element and shape function values on the surface

void LagrangeSmallStrainLinearElastic::GetJacobianAndShapeFunctionsOnInterface( const localIndex numNodes,
                                                                                const rArray1d& xe,
                                                                                const realT& psi,
                                                                                const realT& eta,
                                                                                realT& jcob_sub,
                                                                                rArray2d& N_sub)
{
  if (numNodes == 2)
  {
    jcob_sub = 0.5 * sqrt(pow((xe[2] - xe[0]), 2) + pow((xe[3] - xe[1]), 2)); // 2-D segment - Jacobian is just the length

    N_sub(0, 0) = 0.5 * (1 - psi);
    N_sub(0, 2) = 0.5 * (1 + psi);
    N_sub(1, 1) = 0.5 * (1 - psi);
    N_sub(1, 3) = 0.5 * (1 + psi);
  }
  else if (numNodes == 3)
  {
    realT jcobXY = 0, jcobYZ = 0, jcobXZ = 0;
    realT dXdPsi = 0, dXdEta = 0, dYdPsi = 0, dYdEta = 0, dZdPsi = 0, dZdEta = 0;

    dXdPsi = xe[0] - xe[6];
    dXdEta = xe[3] - xe[6];
    dYdPsi = xe[1] - xe[7];
    dYdEta = xe[4] - xe[7];
    dZdPsi = xe[2] - xe[8];
    dZdEta = xe[5] - xe[8];

    jcobXY = dXdPsi * dYdEta - dXdEta * dYdPsi;
    jcobXZ = dXdPsi * dZdEta - dXdEta * dZdPsi;
    jcobYZ = dYdPsi * dZdEta - dYdEta * dZdPsi;

    jcob_sub = 0.5 * (1.0 / 3.0) * sqrt(jcobXY * jcobXY + jcobXZ * jcobXZ + jcobYZ * jcobYZ);

    N_sub(0, 0) = psi;
    N_sub(0, 3) = eta;
    N_sub(0, 6) = 1 - psi - eta;
    N_sub(1, 1) = psi;
    N_sub(1, 4) = eta;
    N_sub(1, 7) = 1 - psi - eta;
    N_sub(2, 2) = psi;
    N_sub(2, 5) = eta;
    N_sub(2, 8) = 1 - psi - eta;

  }
  else if (numNodes == 4)
  {
    // Calculate Jacobian for the surface element - required for evaluating the surface integral on the parametric space. Basically, change of variables from x-y to psi-eta.
    // Even though surface element is only a plane, it could have all x, y and coordinates changing.
    realT jcobXY = 0, jcobYZ = 0, jcobXZ = 0;
    realT dXdPsi = 0, dXdEta = 0, dYdPsi = 0, dYdEta = 0, dZdPsi = 0, dZdEta = 0;

    dXdPsi = 0.25 * (xe[dim] - xe[0]) * (1 - eta) + 0.25 * (xe[2 * dim] - xe[3 * dim]) * (1 + eta);
    dXdEta = 0.25 * (xe[3 * dim] - xe[0]) * (1 - psi) + 0.25 * (xe[2 * dim] - xe[dim]) * (1 + psi);

    dYdPsi = 0.25 * (xe[dim + 1] - xe[1]) * (1 - eta) + 0.25 * (xe[2 * dim + 1] - xe[3 * dim + 1]) * (1 + eta);
    dYdEta = 0.25 * (xe[3 * dim + 1] - xe[1]) * (1 - psi) + 0.25 * (xe[2 * dim + 1] - xe[dim + 1]) * (1 + psi);

    dZdPsi = 0.25 * (xe[dim + 2] - xe[2]) * (1 - eta) + 0.25 * (xe[2 * dim + 2] - xe[3 * dim + 2]) * (1 + eta);
    dZdEta = 0.25 * (xe[3 * dim + 2] - xe[2]) * (1 - psi) + 0.25 * (xe[2 * dim + 2] - xe[dim + 2]) * (1 + psi);

    jcobXY = dXdPsi * dYdEta - dXdEta * dYdPsi;
    jcobXZ = dXdPsi * dZdEta - dXdEta * dZdPsi;
    jcobYZ = dYdPsi * dZdEta - dYdEta * dZdPsi;

    jcob_sub = sqrt(jcobXY * jcobXY + jcobXZ * jcobXZ + jcobYZ * jcobYZ);

    N_sub(0, 0)  = 0.25 * (1 - psi) * (1 - eta);
    N_sub(0, 3)  = 0.25 * (1 + psi) * (1 - eta);
    N_sub(0, 6)  = 0.25 * (1 + psi) * (1 + eta);
    N_sub(0, 9)  = 0.25 * (1 - psi) * (1 + eta);
    N_sub(1, 1)  = 0.25 * (1 - psi) * (1 - eta);
    N_sub(1, 4)  = 0.25 * (1 + psi) * (1 - eta);
    N_sub(1, 7)  = 0.25 * (1 + psi) * (1 + eta);
    N_sub(1, 10) = 0.25 * (1 - psi) * (1 + eta);
    N_sub(2, 2)  = 0.25 * (1 - psi) * (1 - eta);
    N_sub(2, 5)  = 0.25 * (1 + psi) * (1 - eta);
    N_sub(2, 8)  = 0.25 * (1 + psi) * (1 + eta);
    N_sub(2, 11) = 0.25 * (1 - psi) * (1 + eta);
  }

  if (jcob_sub == 0)
  {
    throw GPException("ZeroJacob::Zero Jacobian computed!");
  }
}
// Returns the Global indices of the parent nodes that were duplicated and their nodal coordinates
// Recursively tests for parent nodes. For example, if 2 -> 10 -> 11, the following searches for 2.

void LagrangeSmallStrainLinearElastic::GetParentFaceNodesAndCoordsAndInterfaceGaussPoints( const localIndex i,
                                                                                                const PhysicalDomainT& domain,
                                                                                                const rArray1d& gauss,
                                                                                                lArray1d& localParentFaceNodes,
                                                                                                rArray1d& xe,
                                                                                                rArray1d& psi,
                                                                                                rArray1d& eta)
{
  const localIndex& numTotNodesDomain = domain.m_feNodeManager.m_numNodes;
  const OneToOneRelation& parentNodeIndex = domain.m_feNodeManager.GetOneToOneMap("parentIndex");

  localIndex kf1;
//  if(dim==2 or dim==3)
//  {
//    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
////    kf1 = childFaceIndex[i][0];
//    if(dim==2)
//    {
//      kf1 = childFaceIndex[i][0];
//    }
//    if(dim==3)
//    {
//      kf1 = i;
//    }
//  }

  if(dim==2)
  {
    const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap( "childIndices" );
    kf1 = childFaceIndex[i][0];
  }
  else
  {
    const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
    bool fe;
    kf1 = domain.m_externalFaces.FaceIndex(f1[i], fe);
  }

//  const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
//  bool fe;
//  localIndex kf1 = domain.m_externalFaces.FaceIndex(f1[i], fe);

  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();

  xe.resize(dim*numNodes);
  localParentFaceNodes.resize(numNodes);
  psi.resize(numNodes);
  eta.resize(numNodes);

  for (localIndex iNod = 0; iNod < numNodes; ++iNod)
  {
    // localParentFaceNodes contains global indices (unsorted) of the parent face that has been duplicated
    localParentFaceNodes[iNod] = domain.m_feFaceManager.m_toNodesRelation[kf1][iNod];
    //Check if localParentFaceNodes contains indices of any children - if yes, find the parent of the corresponding child
    while (parentNodeIndex[localParentFaceNodes[iNod]] < numTotNodesDomain)
    {
      localParentFaceNodes[iNod] = parentNodeIndex[localParentFaceNodes[iNod]];
    }

    // Nodal coordinates of the parent nodes
    for (int d = 0; d < dim; ++d)
    {
      xe[iNod * dim + d] = (*domain.m_feNodeManager.m_refposition)[localParentFaceNodes[iNod]][d]; // Get nodal coordinates for interface nodes
    }
  }

  if (numNodes == 2)
  {
    psi[0] = gauss[0]; psi[1] = gauss[1]; eta[0] = 0.0; eta[1] = 0.0;
  }
  else if (numNodes == 3)
  {
    psi[0] = 2.0 / 3.0; psi[1] = 1.0 / 6.0; psi[2] = 1.0 / 6.0;
    eta[0] = 1.0 / 6.0; eta[1] = 2.0 / 3.0; eta[2] = 1.0 / 6.0;
  }
  else if (numNodes == 4)
  {
    psi[0] = gauss[0]; psi[1] = gauss[1]; psi[2] = gauss[1]; psi[3] = gauss[0];
    eta[0] = gauss[0]; eta[1] = gauss[0]; eta[2] = gauss[1]; eta[3] = gauss[1];
  }
  else
  {
    throw GPException(
        "NonRecognizedSurfaceElem::Interface element is neither a triangle, quad nor a segment");
  }
}

// Multiples matrix A: (m x n) with matrix B: (n x k) to give C: (m x k)

void LagrangeSmallStrainLinearElastic::MultiplyArray( const rArray2d& A,
                                                      const rArray2d& B,
                                                      rArray2d& C)

{
  const int m = A.Dimension(0), n = A.Dimension(1), k = B.Dimension(1);
  C.resize2(m, k);

  realT Temp=0;
  for( int iRow=0 ; iRow<m; iRow++)
  {
    for( int jCol=0; jCol < k; jCol++)
    {
      Temp = 0.;
      for( int s=0; s< n; s++)
      {
        Temp += A(iRow,s)*B(s,jCol);
      }
      C(iRow,jCol) = Temp;
    }
  }
}

// Multiples matrix A: (m x n) with vector B: (n x 1) to give vector C: (n x 1)

void LagrangeSmallStrainLinearElastic::MultiplyArray( const rArray2d& A,
                                                      const rArray1d& B,
                                                      rArray1d& C)

{
  const int m = A.Dimension(0), n = A.Dimension(1);
  C.resize(m);

  realT Temp=0;
  for( int iRow=0 ; iRow<m; iRow++)
  {
    Temp = 0.;
    for( int s=0; s< n; s++)
    {
      Temp += A(iRow,s)*B(s);
    }
    C(iRow) = Temp;
  }
}


void LagrangeSmallStrainLinearElastic::GetFaceConnToElemConnMap( const PhysicalDomainT& domain,
                                                                      const localIndex kf,
                                                                      const bool sibFlag,
                                                                      iArray1d& FaceConnToElemConn)
{
  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf].size();
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[kf][0].first;

  const localIndex numNodesParentEle = elemRegion->m_numNodesPerElem;
  const localIndex EleID = domain.m_feFaceManager.m_toElementsRelation[kf][0].second;

  lArray1d locConn(numNodesParentEle);
  for (localIndex iNod=0; iNod<numNodesParentEle; ++iNod)
  {
    locConn(iNod) = elemRegion->m_toNodesRelation[EleID][iNod];
  }

  for (localIndex iNod=0; iNod<numNodes; ++iNod)
  {
    localIndex nodIndex;
    if(sibFlag)
    {
      nodIndex = iNod == 0 ? iNod : numNodes - iNod;
    }
    else
    {
      nodIndex = iNod;
    }

    const localIndex globalIDFaceiNod = domain.m_feFaceManager.m_toNodesRelation[kf][nodIndex];

    for (localIndex jNod=0; jNod<numNodesParentEle; ++jNod)
    {
      if(globalIDFaceiNod == locConn(jNod))
      {
        FaceConnToElemConn(iNod) = jNod;
      }
    }
  }
}


void LagrangeSmallStrainLinearElastic::GetNormalVoigt( const PhysicalDomainT& domain,
                                                            const localIndex kf,
                                                            rArray2d& normalVoigt)
{
  const R1Tensor normal =  domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, kf, true);

  if(dim==2)
  {
    normalVoigt(0,0) = normal(0); normalVoigt(0,2) = normal(1); normalVoigt(1,1) = normal(1); normalVoigt(1,2) = normal(0);
  }

  if(dim==3)
  {
    normalVoigt(0,0) = normal(0); normalVoigt(0,4) = normal(2); normalVoigt(0,5) = normal(1);
    normalVoigt(1,1) = normal(1); normalVoigt(1,3) = normal(2); normalVoigt(1,5) = normal(0);
    normalVoigt(2,2) = normal(2); normalVoigt(2,3) = normal(1); normalVoigt(2,4) = normal(0);
  }
}


void LagrangeSmallStrainLinearElastic::GetElasticityTensorVoigt( const PhysicalDomainT& domain,
                                                                      const localIndex kf,
                                                                      rArray2d& D)
{
  const localIndex EleID = domain.m_feFaceManager.m_toElementsRelation[kf][0].second;
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[kf][0].first;

  const localIndex paramIndex = elemRegion->m_mat->NumParameterIndex0() > 1 ? EleID : 0 ;
  const MaterialBaseParameterData& matParams = *(elemRegion->m_mat->ParameterData(paramIndex));

  const realT G = matParams.init_shearModulus;
  realT lambda = matParams.Lame;

  if( LagrangeSolverBase::m_2dOption==LagrangeSolverBase::PlaneStress )
  {
    lambda = 2*lambda*G / ( lambda + 2*G );
  }
  if(dim==2)
  {
    D(0,0) = lambda+2*G; D(0,1) = lambda; D(1,0) = lambda; D(1,1) = lambda+2*G; D(2,2) = G;
  }
  if(dim==3)
  {
    D(0,0) = lambda+2*G; D(0,1) = lambda;     D(0,2) = lambda;
    D(1,0) = lambda;     D(1,1) = lambda+2*G; D(1,2) = lambda;
    D(2,0) = lambda;     D(2,1) = lambda;     D(2,2) = lambda+2*G;
    D(3,3) = G;          D(4,4) = G;          D(5,5) = G;
  }

}


void LagrangeSmallStrainLinearElastic::GetShapeFunctionDerivativeMatrixConstStr( const PhysicalDomainT& domain,
                                                                                      const iArray1d& FaceConnToElemConn,
                                                                                      const localIndex kf,
                                                                                      rArray2d& B)
{
  const localIndex EleID = domain.m_feFaceManager.m_toElementsRelation[kf][0].second;
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[kf][0].first;

  const localIndex numNodesParentEle = elemRegion->m_numNodesPerElem;

  // If element type is STRI (constant-strain): derivatives are constant
  for(int iRow=0;iRow<0.5*dim*(dim+1);++iRow)
  {
    for(localIndex iCol=0;iCol<numNodesParentEle;++iCol)
    {
      if(iRow<dim)
      {
        B(iRow,dim*iCol+iRow) = elemRegion->m_dNdX[EleID](0,iCol)[iRow];
      }
      else if(iRow>=dim)
      {
        if(dim==2)
        {
          for(int iDim=0;iDim<dim;++iDim)
          {
            B(iRow,dim*iCol+iDim) = elemRegion->m_dNdX[EleID](0,iCol)[(dim-1)-iDim];
          }
        }
        else if(dim==3)
        {
          if(iRow==3)
          {
            B(iRow,dim*iCol+1) = elemRegion->m_dNdX[EleID](0,iCol)[2];
            B(iRow,dim*iCol+2) = elemRegion->m_dNdX[EleID](0,iCol)[1];
          }
          else if(iRow==4)
          {
            B(iRow,dim*iCol) = elemRegion->m_dNdX[EleID](0,iCol)[2];
            B(iRow,dim*iCol+2) = elemRegion->m_dNdX[EleID](0,iCol)[0];
          }
          else if(iRow==5)
          {
            B(iRow,dim*iCol) = elemRegion->m_dNdX[EleID](0,iCol)[1];
            B(iRow,dim*iCol+1) = elemRegion->m_dNdX[EleID](0,iCol)[0];
          }
        }
      }
    }
  }
}


void LagrangeSmallStrainLinearElastic::GetShapeFunctionDerivativeMatrixQuadHex( const PhysicalDomainT& domain,
                                                                                     const iArray1d& FaceConnToElemConn,
                                                                                     const localIndex kf,
                                                                                     realT psi,
                                                                                     realT eta,
                                                                                     rArray2d& B)
{
  const localIndex EleID = domain.m_feFaceManager.m_toElementsRelation[kf][0].second;
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[kf][0].first;

  const localIndex numNodesParentEle = elemRegion->m_numNodesPerElem;

  rArray1d xEle(dim*numNodesParentEle);
  for (localIndex iNod = 0; iNod < numNodesParentEle; ++iNod)
  {
    for( int d=0 ; d<dim ; ++d )
    {
      xEle[iNod*dim+d] = (*domain.m_feNodeManager.m_refposition)[elemRegion->m_toNodesRelation[EleID] [iNod]][d];
    }
  }

  iArray1d FaceConnToElemConnSort = FaceConnToElemConn;
  std::sort(FaceConnToElemConnSort.begin(), FaceConnToElemConnSort.end());

  if(elemRegion->m_elementGeometryID.compare(0,4,"CPE4") ==0)
  {
    realT dN1dPsi = 0, dN2dPsi = 0, dN3dPsi = 0, dN4dPsi = 0;
    realT dN1dEta = 0, dN2dEta = 0, dN3dEta = 0, dN4dEta = 0;

    realT dXdPsi = 0, dXdEta = 0, dYdPsi = 0, dYdEta = 0, jcob = 0;

    realT dN1dX = 0, dN2dX = 0, dN3dX = 0, dN4dX = 0;
    realT dN1dY = 0, dN2dY = 0, dN3dY = 0, dN4dY = 0;

    if(FaceConnToElemConnSort(0)==0 and FaceConnToElemConnSort(1)==1)
    {
      eta = -1;
    }
    else if(FaceConnToElemConnSort(0)==2 and FaceConnToElemConnSort(1)==3)
    {
      eta = 1;
    }
    else if(FaceConnToElemConnSort(0)==0 and FaceConnToElemConnSort(1)==2)
    {
      eta = psi;
      psi = -1;
    }
    else if(FaceConnToElemConnSort(0)==1 and FaceConnToElemConnSort(1)==3)
    {
      eta = psi;
      psi = 1;
    }

    dN1dPsi = -0.25*(1-eta); dN2dPsi =  0.25*(1-eta); dN3dPsi = -0.25*(1+eta); dN4dPsi = 0.25*(1+eta);
    dN1dEta = -0.25*(1-psi); dN2dEta = -0.25*(1+psi); dN3dEta =  0.25*(1-psi); dN4dEta = 0.25*(1+psi);

    dXdPsi = 0.25*(1-eta)*(xEle(2)-xEle(0)) + 0.25*(1+eta)*(xEle(6)-xEle(4));
    dXdEta = 0.25*(1-psi)*(xEle(4)-xEle(0)) + 0.25*(1+psi)*(xEle(6)-xEle(2));

    dYdPsi = 0.25*(1-eta)*(xEle(3)-xEle(1)) + 0.25*(1+eta)*(xEle(7)-xEle(5));
    dYdEta = 0.25*(1-psi)*(xEle(5)-xEle(1)) + 0.25*(1+psi)*(xEle(7)-xEle(3));

    jcob = dXdPsi*dYdEta - dXdEta*dYdPsi;

    dN1dX = (dN1dPsi*dYdEta - dN1dEta*dYdPsi)/jcob, dN1dY = -(dN1dPsi*dXdEta - dN1dEta*dXdPsi)/jcob;
    dN2dX = (dN2dPsi*dYdEta - dN2dEta*dYdPsi)/jcob, dN2dY = -(dN2dPsi*dXdEta - dN2dEta*dXdPsi)/jcob;
    dN3dX = (dN3dPsi*dYdEta - dN3dEta*dYdPsi)/jcob, dN3dY = -(dN3dPsi*dXdEta - dN3dEta*dXdPsi)/jcob;
    dN4dX = (dN4dPsi*dYdEta - dN4dEta*dYdPsi)/jcob, dN4dY = -(dN4dPsi*dXdEta - dN4dEta*dXdPsi)/jcob;

    B(0,0) = dN1dX; B(0,2) = dN2dX; B(0,4) = dN3dX; B(0,6) = dN4dX;
    B(1,1) = dN1dY; B(1,3) = dN2dY; B(1,5) = dN3dY; B(1,7) = dN4dY;
    B(2,0) = dN1dY; B(2,1) = dN1dX; B(2,2) = dN2dY; B(2,3) = dN2dX;
    B(2,4) = dN3dY; B(2,5) = dN3dX; B(2,6) = dN4dY; B(2,7) = dN4dX;
  }
  else if (elemRegion->m_elementGeometryID.compare(0,4,"C3D8") ==0 )
  {
    realT chi = 0;
    if(FaceConnToElemConnSort(0) == 0 and FaceConnToElemConnSort(1) == 1 and FaceConnToElemConnSort(2) == 2 and FaceConnToElemConnSort(3) == 3)
    {
      chi = -1;
    }
    else if(FaceConnToElemConnSort(0) == 4 and FaceConnToElemConnSort(1) == 5 and FaceConnToElemConnSort(2) == 6 and FaceConnToElemConnSort(3) == 7)
    {
      chi = 1;
    }
    else if(FaceConnToElemConnSort(0) == 0 and FaceConnToElemConnSort(1) == 1 and FaceConnToElemConnSort(2) == 4 and FaceConnToElemConnSort(3) == 5)
    {
      chi = eta;
      eta = -1;
    }
    else if(FaceConnToElemConnSort(0) == 2 and FaceConnToElemConnSort(1) == 3 and FaceConnToElemConnSort(2) == 6 and FaceConnToElemConnSort(3) == 7)
    {
      chi = eta;
      eta = 1;
    }
    else if(FaceConnToElemConnSort(0) == 0 and FaceConnToElemConnSort(1) == 2 and FaceConnToElemConnSort(2) == 4 and FaceConnToElemConnSort(3) == 6)
    {
      chi = eta;
      eta = psi;
      psi = -1;
    }
    else if(FaceConnToElemConnSort(0)== 1 and FaceConnToElemConnSort(1)== 3 and FaceConnToElemConnSort(2) == 5 and FaceConnToElemConnSort(3) == 7)
    {
      chi = eta;
      eta = psi;
      psi = 1;
    }

    rArray1d dNdPsi(numNodesParentEle), dNdEta(numNodesParentEle), dNdChi(numNodesParentEle);
    rArray1d psi_a(numNodesParentEle), eta_a(numNodesParentEle), chi_a(numNodesParentEle);
    psi_a(0) = -1; psi_a(1) = 1; psi_a(2) = -1; psi_a(3) = 1; psi_a(4) = -1; psi_a(5) = 1; psi_a(6) = -1; psi_a(7) = 1;
    eta_a(0) = -1; eta_a(1) = -1; eta_a(2) = 1; eta_a(3) = 1; eta_a(4) = -1; eta_a(5) = -1; eta_a(6) = 1; eta_a(7) = 1;
    chi_a(0) = -1; chi_a(1) = -1; chi_a(2) = -1; chi_a(3) = -1; chi_a(4) = 1; chi_a(5) = 1; chi_a(6) = 1; chi_a(7) = 1;

    for (localIndex iNod = 0; iNod< numNodesParentEle; ++iNod)
    {
      dNdPsi(iNod) = 0.125*psi_a(iNod)*(1+eta_a(iNod)*eta)*(1+chi_a(iNod)*chi);
      dNdEta(iNod) = 0.125*eta_a(iNod)*(1+psi_a(iNod)*psi)*(1+chi_a(iNod)*chi);
      dNdChi(iNod) = 0.125*chi_a(iNod)*(1+psi_a(iNod)*psi)*(1+eta_a(iNod)*eta);
    }

    realT dXdPsi = 0, dXdEta = 0, dXdChi = 0;
    realT dYdPsi = 0, dYdEta = 0, dYdChi = 0;
    realT dZdPsi = 0, dZdEta = 0, dZdChi = 0;

    for (localIndex iNod = 0; iNod<numNodesParentEle; ++iNod)
    {
      dXdPsi += 0.125*psi_a(iNod)*(1+eta_a(iNod)*eta)*(1+chi_a(iNod)*chi)*xEle(dim*iNod);
      dXdEta += 0.125*eta_a(iNod)*(1+psi_a(iNod)*psi)*(1+chi_a(iNod)*chi)*xEle(dim*iNod);
      dXdChi += 0.125*chi_a(iNod)*(1+psi_a(iNod)*psi)*(1+eta_a(iNod)*eta)*xEle(dim*iNod);

      dYdPsi += 0.125*psi_a(iNod)*(1+eta_a(iNod)*eta)*(1+chi_a(iNod)*chi)*xEle(dim*iNod + 1);
      dYdEta += 0.125*eta_a(iNod)*(1+psi_a(iNod)*psi)*(1+chi_a(iNod)*chi)*xEle(dim*iNod + 1);
      dYdChi += 0.125*chi_a(iNod)*(1+psi_a(iNod)*psi)*(1+eta_a(iNod)*eta)*xEle(dim*iNod + 1);

      dZdPsi += 0.125*psi_a(iNod)*(1+eta_a(iNod)*eta)*(1+chi_a(iNod)*chi)*xEle(dim*iNod + 2);
      dZdEta += 0.125*eta_a(iNod)*(1+psi_a(iNod)*psi)*(1+chi_a(iNod)*chi)*xEle(dim*iNod + 2);
      dZdChi += 0.125*chi_a(iNod)*(1+psi_a(iNod)*psi)*(1+eta_a(iNod)*eta)*xEle(dim*iNod + 2);
    }

    realT cof11 = 0, cof12 = 0, cof13 = 0, cof21 = 0, cof22 = 0, cof23 = 0, cof31 = 0, cof32 = 0, cof33 = 0, jcob = 0;;

    cof11 = dYdEta*dZdChi - dYdChi*dZdEta;
    cof12 = dYdChi*dZdPsi - dYdPsi*dZdChi;
    cof13 = dYdPsi*dZdEta - dYdEta*dZdPsi;

    cof21 = dZdEta*dXdChi - dZdChi*dXdEta;
    cof22 = dZdChi*dXdPsi - dZdPsi*dXdChi;
    cof23 = dZdPsi*dXdEta - dZdEta*dXdPsi;

    cof31 = dXdEta*dYdChi - dXdChi*dYdEta;
    cof32 = dXdChi*dYdPsi - dXdPsi*dYdChi;
    cof33 = dXdPsi*dYdEta - dXdEta*dYdPsi;

    jcob = dXdPsi*cof11 + dXdEta*cof12 + dXdChi*cof13;

    rArray1d dNdX(numNodesParentEle), dNdY(numNodesParentEle), dNdZ(numNodesParentEle);
    for (localIndex jNod = 0; jNod< numNodesParentEle; ++jNod)
    {
      dNdX(jNod) = (dNdPsi(jNod)*cof11 + dNdEta(jNod)*cof12 + dNdChi(jNod)*cof13)/jcob;
      dNdY(jNod) = (dNdPsi(jNod)*cof21 + dNdEta(jNod)*cof22 + dNdChi(jNod)*cof23)/jcob;
      dNdZ(jNod) = (dNdPsi(jNod)*cof31 + dNdEta(jNod)*cof32 + dNdChi(jNod)*cof33)/jcob;
    }

    for (localIndex iCol = 0; iCol<numNodesParentEle; ++iCol)
    {
      B(0, dim*iCol)     = dNdX(iCol);
      B(1, dim*iCol + 1) = dNdY(iCol);
      B(2, dim*iCol + 2) = dNdZ(iCol);

      B(3, dim*iCol + 1) = dNdZ(iCol);
      B(3, dim*iCol + 2) = dNdY(iCol);

      B(4, dim*iCol)     = dNdZ(iCol);
      B(4, dim*iCol + 2) = dNdX(iCol);

      B(5, dim*iCol)     = dNdY(iCol);
      B(5, dim*iCol + 1) = dNdX(iCol);
    }
  }
}


void LagrangeSmallStrainLinearElastic::GetDBDotN( const rArray2d& D,
                                                       const rArray2d& B,
                                                       const rArray2d& normalVoigt,
                                                       const localIndex numNodesParentEle,
                                                       const realT gamma,
                                                       const PhysicalDomainT& domain,
                                                       const rArray2d& P,
                                                       rArray2d& nDotDB)
{
  rArray2d DB(0.5*dim*(dim+1), dim*numNodesParentEle);
  MultiplyArray(D, B, DB);
  MultiplyArray(normalVoigt, DB, nDotDB);

  for (int iRow = 0; iRow<dim; ++iRow)
  {
    for (localIndex iCol = 0; iCol<dim*numNodesParentEle; ++iCol)
    {
      nDotDB(iRow, iCol) = gamma*nDotDB(iRow, iCol);
    }
  }

  rArray2d nDotDB_NT(dim, dim*numNodesParentEle);
  MultiplyArray(P, nDotDB, nDotDB_NT);

  if(domain.m_contactManager.m_nitsche_tau1_active == 0)
  {
    for (localIndex iCol=0; iCol<dim*numNodesParentEle; ++iCol)
    {
      nDotDB_NT(1,iCol) = 0.0;
    }
  }
  if(dim==3 and domain.m_contactManager.m_nitsche_tau2_active ==  0)
  {
    for (localIndex iCol=0; iCol<dim*numNodesParentEle; ++iCol)
    {
      nDotDB_NT(2,iCol) = 0.0;
    }
  }

  rArray2d PTranspose(dim,dim);
  GetMatrixTranspose (P, PTranspose);

  nDotDB = 0;
  MultiplyArray(PTranspose, nDotDB_NT, nDotDB);
}


void LagrangeSmallStrainLinearElastic::GetPermutedDBDotN( const PhysicalDomainT& domain,
                                                               const rArray2d& nDotDB,
                                                               const iArray1d& FaceConnToElemConn,
                                                               const localIndex kf,
                                                               rArray2d& nDotDBPermute)
{
  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf].size();
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[kf][0].first;

  const localIndex numNodesParentEle = elemRegion->m_numNodesPerElem;

  iArray1d NodIDNotOnFace(numNodesParentEle-numNodes);
  GetNodIndicesNotOnFace(domain, kf, false, NodIDNotOnFace);

  for (int iRow=0; iRow<dim; ++iRow)
  {
    for (localIndex iCol=0;iCol<numNodesParentEle;++iCol)
    {
      for (int iDim=0; iDim<dim; ++iDim)
      {
        if(iCol<numNodes)
        {
          nDotDBPermute(iRow, dim*iCol+iDim) = nDotDB(iRow,dim*FaceConnToElemConn(iCol)+iDim);
        }
        else
        {
          int kNod = iCol-numNodes;
          nDotDBPermute(iRow, dim*iCol+iDim) = nDotDB(iRow, dim*NodIDNotOnFace(kNod)+iDim);
        }
      }
    }
  }
}


void LagrangeSmallStrainLinearElastic::GetInitialAlphaTensor( const PhysicalDomainT& domain,
                                                              const localIndex kf,
                                                              const rArray2d& P,
                                                              rArray2d& alpha )
{
  realT alpha_n, alpha_t1, alpha_t2;
  if(domain.m_contactManager.m_nitsche_active)
  {
    const rArray1d& nitscheStabNormal = domain.m_externalFaces.GetFieldData<realT>("nitscheStab_n");
    const rArray1d& nitscheStabT1 = domain.m_externalFaces.GetFieldData<realT>("nitscheStab_t1");
    const rArray1d& nitscheStabT2 = domain.m_externalFaces.GetFieldData<realT>("nitscheStab_t2");
    const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

    if(domain.m_contactManager.m_nitsche_normal_active)
      alpha_n = nitscheStabNormal(faceToExternalFaceMap(kf));
    else
      alpha_n = domain.m_contactManager.m_glob_penalty_n;
    if(domain.m_contactManager.m_nitsche_tau1_active)
      alpha_t1 = nitscheStabT1(faceToExternalFaceMap(kf));
    else
      alpha_t1 = domain.m_contactManager.m_glob_penalty_t1;
    if(dim==3)
    {
      if(domain.m_contactManager.m_nitsche_tau2_active)
        alpha_t2 = nitscheStabT2(faceToExternalFaceMap(kf));
      else
        alpha_t2 = domain.m_contactManager.m_glob_penalty_t2;
    }
  }
  else
  {
    alpha_n = domain.m_contactManager.m_glob_penalty_n;
    alpha_t1 = domain.m_contactManager.m_glob_penalty_t1;
    if(dim==3)
      alpha_t2 = domain.m_contactManager.m_glob_penalty_t2;
  }

  rArray2d alpha_NT(dim,dim);
  alpha_NT(0,0) = alpha_n; alpha_NT(1,1) = alpha_t1;
  if(dim==3)
    alpha_NT(2,2) = alpha_t2;

  //Rotate alpha_NT to alpha_XYZ
  RotateMatrix(P, alpha_NT, alpha);
}
//@annavarapusr1: Function to transform matrix from one coordinate system to another
//@brief: param[in] - Transformation tensor Q (dimension nxn)
//        param[in] - Matrix A (dimension nxn)
//        param[out] - Transformed Matrix A' = Q'*A*Q

void LagrangeSmallStrainLinearElastic::RotateMatrix( const rArray2d& Q,
                                                     const rArray2d& A,
                                                     rArray2d& APrime)
{
  const int n = Q.Dimension(0);

  rArray2d AQ(n,n);
  MultiplyArray(A, Q, AQ);

  rArray2d QTranspose(n,n);
  GetMatrixTranspose (Q, QTranspose);

  MultiplyArray(QTranspose, AQ, APrime);
}

void LagrangeSmallStrainLinearElastic::GetTransformationTensor( const PhysicalDomainT& domain,
                                                                const localIndex kf,
                                                                rArray2d& P)
{
  const R1Tensor normal =  domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, kf, true);
  if(dim==2)
  {
    R1Tensor ex, ey, ez, tau;
    ex(0) = 1; ey(1) = 1; ez(2) = 1;
    tau.Cross(ez,normal);
    P(0,0) = Dot(normal,ex); P(0,1) = Dot(normal,ey);
    P(1,0) = Dot(tau,ex);   P(1,1) = Dot(tau,ey);
  }

  else if(dim==3)
  {
    realT maxp1 = fabs(normal(0))+ fabs(normal(1));
    realT maxp2 = fabs(normal(1))+ fabs(normal(2));
    realT maxp3 = fabs(normal(0))+ fabs(normal(2));
    realT maxpair1 = (maxp1<maxp2)?maxp2:maxp1;
    realT maxpair = (maxpair1<maxp3)?maxp3:maxpair1;
    R1Tensor tau1, tau2;

    R1Tensor ex, ey, ez;
    ex(0) = 1; ey(1) = 1; ez(2) = 1;

    if(maxpair==maxp1)
    {
      tau1(0) = normal(1);
      tau1(1) = -normal(0);
      tau1(2) = 0;
    }
    else if(maxpair==maxp2)
    {
      tau1(1) = normal(2);
      tau1(2) = -normal(1);
      tau1(0) = 0;
    }
    else if(maxpair==maxp3)
    {
     tau1(2) = normal(0);
     tau1(0) = -normal(2);
     tau1(1) = 0;
    }

    tau1.Normalize();

    tau2.Cross(normal,tau1);
    tau2.Normalize();

    P(0,0) = Dot(normal,ex); P(0,1) = Dot(normal,ey); P(0,2) = Dot(normal,ez);
    P(1,0) = Dot(tau1,ex);   P(1,1) = Dot(tau1,ey);   P(1,2) = Dot(tau1,ez);
    P(2,0) = Dot(tau2,ex);   P(2,1) = Dot(tau2,ey);   P(2,2) = Dot(tau2,ez);
  }
}

void LagrangeSmallStrainLinearElastic::GetDisplacementEleNodes( const PhysicalDomainT& domain,
                                                                     const localIndex kf,
                                                                     rArray1d& uEle)
{
  const localIndex EleID = domain.m_feFaceManager.m_toElementsRelation[kf][0].second;
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[kf][0].first;

  const localIndex* const localNodeIndices = elemRegion->m_toNodesRelation[EleID];
  const Array1dT<R1Tensor>& disp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();

  for(unsigned i=0; i<elemRegion->m_numNodesPerElem; ++i)
  {
    const localIndex localNodeIndex = localNodeIndices[i];
    for( int d=0 ; d<dim ; ++d )
    {
      uEle(i*dim+d) = disp[localNodeIndex][d];
    }
  }
}

void LagrangeSmallStrainLinearElastic::GetDisplacementJumpInterface( const PhysicalDomainT& domain,
                                                                          const ElementRegionT& elemRegion,
                                                                          const localIndex kf1,
                                                                          const localIndex kf2,
                                                                          const rArray2d& N_sub,
                                                                          const localIndex iContFace,
                                                                          const lArray1d& localParentFaceNodes,
                                                                          rArray1d& uJump)
{
  rArray1d uFaceElm(dim), uFaceSib(dim);
  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();

  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegionPar = ftoe[kf1][0].first;
  const ElementRegionT* elemRegionSib = ftoe[kf2][0].first;

  const localIndex numNodesParentEle = elemRegionPar->m_numNodesPerElem;
  const localIndex numNodesSiblingEle = elemRegionSib->m_numNodesPerElem;

  rArray1d uElm(dim*numNodesParentEle), uSib(dim*numNodesSiblingEle);
  GetDisplacementEleNodes(domain, kf1, uElm);
  GetDisplacementEleNodes(domain, kf2, uSib);

  iArray1d FaceConnToElemConnPar(numNodes), FaceConnToElemConnSib(numNodes);
  GetFaceConnToElemConnMap(domain, kf1, false, FaceConnToElemConnPar);
  GetFaceConnToElemConnMap(domain, kf2, true, FaceConnToElemConnSib);

  for (localIndex iFaceNod=0; iFaceNod<numNodes; ++iFaceNod)
  {
    localIndex ParLocalIndex, SibLocalIndex;
    GetLocalIndexOnInterface(iContFace, iFaceNod, domain, localParentFaceNodes, ParLocalIndex, SibLocalIndex);
    for(int iDim=0; iDim<dim; ++iDim)
    {
      uFaceElm(iDim) += N_sub(iDim, dim*ParLocalIndex + iDim)*uElm(dim*FaceConnToElemConnPar(iFaceNod)+iDim);
      uFaceSib(iDim) += N_sub(iDim, dim*SibLocalIndex + iDim)*uSib(dim*FaceConnToElemConnSib(iFaceNod)+iDim);
    }
  }

  for (localIndex iSize=0; iSize<uFaceElm.size(); ++iSize)
  {
    uJump(iSize) = uFaceElm(iSize) - uFaceSib(iSize);
  }
}

void LagrangeSmallStrainLinearElastic::GetTractionInterface( const PhysicalDomainT& domain,
                                                             const rArray2d& N_sub,
                                                             const rArray2d& alphaPar,
                                                             const rArray2d& alphaSib,
                                                             const rArray2d& nParDotDParBParPermute,
                                                             const rArray2d& nSibDotDSibBSibPermute,
                                                             const lArray1d& localParentFaceNodes,
                                                             const localIndex kf1,
                                                             const localIndex kf2,
                                                             const localIndex iContFace,
                                                             rArray1d& tracPar,
                                                             rArray1d& tracSib)
{
  rArray1d uFaceElm(dim), uFaceSib(dim), uJump(dim), uJumpSib(dim), tracParPenalty(dim), tracSibPenalty(dim), tracParNitsche(dim), tracSibNitsche(dim);
  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegionPar = ftoe[kf1][0].first;
  const ElementRegionT* elemRegionSib = ftoe[kf2][0].first;

  const localIndex numNodesParentEle = elemRegionPar->m_numNodesPerElem;
  const localIndex numNodesSiblingEle = elemRegionSib->m_numNodesPerElem;

  rArray1d uElm(dim*numNodesParentEle), uSib(dim*numNodesSiblingEle), uElmPermute(dim*numNodesParentEle), uSibPermute(dim*numNodesSiblingEle);
  GetDisplacementEleNodes(domain, kf1, uElm);
  GetDisplacementEleNodes(domain, kf2, uSib);

  iArray1d FaceConnToElemConnPar(numNodes), FaceConnToElemConnSib(numNodes);
  GetFaceConnToElemConnMap(domain, kf1, false, FaceConnToElemConnPar);
  GetFaceConnToElemConnMap(domain, kf2, true, FaceConnToElemConnSib);

  for (localIndex iFaceNod=0; iFaceNod<numNodes; ++iFaceNod)
  {
    localIndex ParLocalIndex, SibLocalIndex;
    GetLocalIndexOnInterface(iContFace, iFaceNod, domain, localParentFaceNodes, ParLocalIndex, SibLocalIndex);
    for(int iDim=0; iDim<dim; ++iDim)
    {
      uFaceElm(iDim) += N_sub(iDim, dim*ParLocalIndex + iDim)*uElm(dim*FaceConnToElemConnPar(iFaceNod)+iDim);
      uFaceSib(iDim) += N_sub(iDim, dim*SibLocalIndex + iDim)*uSib(dim*FaceConnToElemConnSib(iFaceNod)+iDim);
    }
  }

  for (localIndex iSize=0; iSize<uFaceElm.size(); ++iSize)
  {
    uJump(iSize) = uFaceElm(iSize) - uFaceSib(iSize);
    uJumpSib(iSize) = uFaceSib(iSize) - uFaceElm(iSize);
  }

  MultiplyArray(alphaPar,uJump,tracParPenalty);
  MultiplyArray(alphaSib,uJumpSib,tracSibPenalty);

  if(domain.m_contactManager.m_nitsche_active)
  {
    GetPermutedNodalVector(domain, uElm, FaceConnToElemConnPar, kf1, uElmPermute);
    GetPermutedNodalVector(domain, uSib, FaceConnToElemConnSib, kf2, uSibPermute);

    MultiplyArray(nParDotDParBParPermute,uElmPermute,tracParNitsche);
    MultiplyArray(nSibDotDSibBSibPermute,uSibPermute,tracSibNitsche);
  }

  for(localIndex iDim=0; iDim<tracPar.size(); ++iDim)
  {
    tracPar(iDim) = tracParNitsche(iDim) - tracSibNitsche(iDim) - tracParPenalty(iDim);
    tracSib(iDim) = -tracParNitsche(iDim) + tracSibNitsche(iDim) - tracSibPenalty(iDim);
  }

/*
   * Check if crack faces are opening, if yes, set traction to zero
   * Check performed: Sign of tracNT and magnitude of uJump.
   * Opening is detected if sign of tracN is positive (tension) and if magnitude of uJump is greater than 1e-06. Smaller value of uJump essentially means numerical noise
*/
//  realT tol = domain.m_contactManager.m_traction_n_tol;
//    GetUpdatedTractionsOpeningMode(P_par, iGp, uJump, alphaPar, tol, openingGP_par, tracPar);
//    GetUpdatedTractionsOpeningMode(P_sib, iGp, uJumpSib, alphaSib, tol, openingGP_sib, tracSib);
}

void LagrangeSmallStrainLinearElastic::GetTractionInterfaceSlidingTrial( const PhysicalDomainT& domain,
                                                                         const rArray2d& N_sub,
                                                                         const rArray2d& alphaPar,
                                                                         const rArray2d& alphaSib,
                                                                         const rArray2d& P_par,
                                                                         const rArray2d& P_sib,
                                                                         const rArray2d& nParDotDParBParPermute,
                                                                         const rArray2d& nSibDotDSibBSibPermute,
                                                                         const lArray1d& localParentFaceNodes,
                                                                         const localIndex kf1,
                                                                         const localIndex kf2,
                                                                         const localIndex iContFace,
                                                                         const std::string FieldName,
                                                                         rArray1d& tracPar,
                                                                         rArray1d& tracSib)
{
  rArray1d uFaceElm(dim), uFaceSib(dim), uJump(dim), uJumpSib(dim), uJumpPlPar(dim), uJumpPlSib(dim), uJumpPlParNT(dim), uJumpPlSibNT(dim);
  rArray1d uJumpElPar(dim), uJumpElSib(dim), tracParPenalty(dim), tracSibPenalty(dim), tracParNitsche(dim), tracSibNitsche(dim);
  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();

  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegionPar = ftoe[kf1][0].first;
  const ElementRegionT* elemRegionSib = ftoe[kf2][0].first;

  const localIndex numNodesParentEle = elemRegionPar->m_numNodesPerElem;
  const localIndex numNodesSiblingEle = elemRegionSib->m_numNodesPerElem;
  const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

  const Array1dT<R1Tensor>& uJumpPl_Previous = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldName);
  for (localIndex iSize=0; iSize<uJumpPlParNT.size(); ++iSize)
  {
    uJumpPlParNT(iSize) = uJumpPl_Previous(faceToExternalFaceMap(kf1))(iSize);
    uJumpPlSibNT(iSize) = uJumpPl_Previous(faceToExternalFaceMap(kf2))(iSize);
  }
  // Get PTranspose
  rArray2d P_ParTranspose(dim,dim), P_SibTranspose(dim,dim);
  GetMatrixTranspose (P_par, P_ParTranspose);
  GetMatrixTranspose (P_sib, P_SibTranspose);

  MultiplyArray(P_ParTranspose,uJumpPlParNT,uJumpPlPar);
  MultiplyArray(P_SibTranspose,uJumpPlSibNT,uJumpPlSib);

  rArray1d uElm(dim*numNodesParentEle), uSib(dim*numNodesSiblingEle), uElmPermute(dim*numNodesParentEle), uSibPermute(dim*numNodesSiblingEle);
  GetDisplacementEleNodes(domain, kf1, uElm);
  GetDisplacementEleNodes(domain, kf2, uSib);

  iArray1d FaceConnToElemConnPar(numNodes), FaceConnToElemConnSib(numNodes);
  GetFaceConnToElemConnMap(domain, kf1, false, FaceConnToElemConnPar);
  GetFaceConnToElemConnMap(domain, kf2, true, FaceConnToElemConnSib);

  for (localIndex iFaceNod=0; iFaceNod<numNodes; ++iFaceNod)
  {
    localIndex ParLocalIndex, SibLocalIndex;
    GetLocalIndexOnInterface(iContFace, iFaceNod, domain, localParentFaceNodes, ParLocalIndex, SibLocalIndex);
    for(int iDim=0; iDim<dim; ++iDim)
    {
      uFaceElm(iDim) += N_sub(iDim, dim*ParLocalIndex + iDim)*uElm(dim*FaceConnToElemConnPar(iFaceNod)+iDim);
      uFaceSib(iDim) += N_sub(iDim, dim*SibLocalIndex + iDim)*uSib(dim*FaceConnToElemConnSib(iFaceNod)+iDim);
    }
  }

  for (localIndex iSize=0; iSize<uFaceElm.size(); ++iSize)
  {
    uJump(iSize) = uFaceElm(iSize) - uFaceSib(iSize);
    uJumpSib(iSize) = uFaceSib(iSize) - uFaceElm(iSize);
  }

  for (localIndex iSize=0; iSize<uJump.size(); ++iSize)
  {
    uJumpElPar(iSize) = uJump(iSize) - uJumpPlPar(iSize);
    uJumpElSib(iSize) = uJumpSib(iSize) - uJumpPlSib(iSize);
  }

  MultiplyArray(alphaPar,uJumpElPar,tracParPenalty);
  MultiplyArray(alphaSib,uJumpElSib,tracSibPenalty);

  if(domain.m_contactManager.m_nitsche_active)
  {
    GetPermutedNodalVector(domain, uElm, FaceConnToElemConnPar, kf1, uElmPermute);
    GetPermutedNodalVector(domain, uSib, FaceConnToElemConnSib, kf2, uSibPermute);

    MultiplyArray(nParDotDParBParPermute,uElmPermute,tracParNitsche);
    MultiplyArray(nSibDotDSibBSibPermute,uSibPermute,tracSibNitsche);
  }

  for(localIndex iDim=0; iDim<tracPar.size(); ++iDim)
  {
    tracPar(iDim) = tracParNitsche(iDim) - tracSibNitsche(iDim) - tracParPenalty(iDim);
    tracSib(iDim) = -tracParNitsche(iDim) + tracSibNitsche(iDim) - tracSibPenalty(iDim);
  }

  /*
   * Check if crack faces are opening, if yes, set traction to zero
   * Check performed: Sign of tracNT and magnitude of uJump.
   * Opening is detected if sign of tracN is positive (tension) and if magnitude of uJump is greater than 1e-06. Smaller value of uJump essentially means numerical noise
   */
//  realT tol = domain.m_contactManager.m_traction_n_tol;
//  GetUpdatedTractionsOpeningMode(P_par, iGp, uJump, alphaPar,  tol, openingGP_par, tracPar);
//  GetUpdatedTractionsOpeningMode(P_sib, iGp, uJumpSib, alphaSib, tol, openingGP_sib, tracSib);
}

void LagrangeSmallStrainLinearElastic::GetTractionInterfaceSliding( const rArray2d& N_sub,
                                                                    const rArray2d& alphaPar,
                                                                    const rArray2d& alphaSib,
                                                                    const rArray2d& P_par,
                                                                    const rArray2d& P_sib,
                                                                    const rArray2d& nParDotDParBParPermute,
                                                                    const rArray2d& nSibDotDSibBSibPermute,
                                                                    const lArray1d& localParentFaceNodes,
                                                                    const localIndex kf1,
                                                                    const localIndex kf2,
                                                                    const localIndex iContFace,
                                                                    const int iGp,
                                                                    const std::string s_previous_base,
                                                                    const std::string s_current_base,
                                                                    PhysicalDomainT& domain,
                                                                    iArray1d& stickGP_par,
                                                                    iArray1d& stickGP_sib,
                                                                    iArray1d& openingGP_par,
                                                                    iArray1d& openingGP_sib,
                                                                    rArray1d& tracPar,
                                                                    rArray1d& tracSib)
{
  // Gets the name of the field associated with Gauss point iGp for previous and current load steps
  std::stringstream ss; ss << iGp;
  std::string FieldNamePrevious = s_previous_base + ss.str(), FieldNameCurrent  = s_current_base + ss.str();

  // Gets trial traction assuming crack faces are sticking
  GetTractionInterfaceSlidingTrial(domain, N_sub, alphaPar, alphaSib, P_par, P_sib, nParDotDParBParPermute, nSibDotDSibBSibPermute,
                                   localParentFaceNodes, kf1, kf2, iContFace, FieldNamePrevious, tracPar, tracSib);

  realT tol = domain.m_contactManager.m_traction_n_tol;
  // Check if crack faces are opening, if yes, set traction to zero
  GetUpdatedTractionsOpeningMode(P_par, iGp, tol, openingGP_par, tracPar);
  GetUpdatedTractionsOpeningMode(P_sib, iGp, tol, openingGP_sib, tracSib);

  if(openingGP_par(iGp)==1)
  {
    const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

    Array1dT<R1Tensor>& uJumpPl_Previous = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNamePrevious);
    Array1dT<R1Tensor>& uJumpPl_Current  = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNameCurrent);

    uJumpPl_Current(faceToExternalFaceMap(kf1)) = uJumpPl_Previous(faceToExternalFaceMap(kf1));
  }
  // If crack faces are not opening, check for sliding and update tractions
  if(openingGP_par(iGp)==0)
  {
    // Gets the trial value value of yield function trial and a boolean flag indicating if the yield criterion is violated
    realT phiTrialPar = GetSlipStickState(domain, tracPar, P_par, stickGP_par(iGp));
    // Pulls back the trial traction to lie on the yield surface, updates stiffness and plastic slip variables (Simo and Hughes (1997)))
    GetUpdatedTractionsAndPlasticSlip(kf1, P_par, phiTrialPar, stickGP_par(iGp), FieldNamePrevious, FieldNameCurrent, alphaPar, domain, tracPar);
  }
  if(openingGP_sib(iGp)==1)
  {
    const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

    Array1dT<R1Tensor>& uJumpPl_Previous = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNamePrevious);
    Array1dT<R1Tensor>& uJumpPl_Current  = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNameCurrent);

    uJumpPl_Current(faceToExternalFaceMap(kf2)) = uJumpPl_Previous(faceToExternalFaceMap(kf2));
  }
  if(openingGP_sib(iGp)==0)
  {
    // Gets the trial value value of yield function trial and a boolean flag indicating if the yield criterion is violated
    realT phiTrialSib = GetSlipStickState(domain, tracSib, P_sib, stickGP_sib(iGp));
    // Pulls back the trial traction to lie on the yield surface, updates stiffness and plastic slip variables (Simo and Hughes (1997)))
    GetUpdatedTractionsAndPlasticSlip(kf2, P_sib, phiTrialSib, stickGP_sib(iGp), FieldNamePrevious, FieldNameCurrent, alphaSib, domain, tracSib);
  }
}

void LagrangeSmallStrainLinearElastic::GetTrialTangentModulus( const PhysicalDomainT& domain,
                                                                 const localIndex kf,
                                                                 const rArray2d& P,
                                                                 const iArray1d& FaceConnToElemConn,
                                                                 const realT psi,
                                                                 const realT eta,
                                                                 const rArray2d& D,
                                                                 const rArray2d& normalVoigt,
                                                                 const realT gamma,
                                                                 rArray2d& alpha,
                                                                 rArray2d& nDotDB)
{
  GetInitialAlphaTensor(domain, kf, P, alpha);
  if(domain.m_contactManager.m_nitsche_active)
  {
    const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
    const ElementRegionT* elemRegion = ftoe[kf][0].first;

    const localIndex numNodesParentEle = elemRegion->m_numNodesPerElem;

    rArray2d B(0.5*dim*(dim+1),dim*numNodesParentEle);
    if(elemRegion->m_elementGeometryID.compare(0,4,"CPE4") ==0 or elemRegion->m_elementGeometryID.compare(0,4,"C3D8") ==0 )
    {
      GetShapeFunctionDerivativeMatrixQuadHex(domain, FaceConnToElemConn, kf, psi, eta, B);
    }
    if(elemRegion->m_elementGeometryID.compare(0,4,"STRI") ==0 or elemRegion->m_elementGeometryID.compare(0,4,"C3D4") == 0)
    {
      GetShapeFunctionDerivativeMatrixConstStr(domain, FaceConnToElemConn, kf, B);
    }

    // GetDBDotN includes the Nitsche weights so should be renamed to GetGammaDBDotN?
    GetDBDotN(D, B, normalVoigt, numNodesParentEle, gamma, domain, P,  nDotDB);
  }
}

void LagrangeSmallStrainLinearElastic::GetPermutedNodalVector( const PhysicalDomainT& domain,
                                                                    const rArray1d& Vector,
                                                                    const iArray1d& FaceConnToElemConn,
                                                                    const localIndex kf,
                                                                    rArray1d& PermutedVector)
{
  const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
  const ElementRegionT* elemRegion = ftoe[kf][0].first;

  const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf].size();
  const localIndex numNodesParentEle = elemRegion->m_numNodesPerElem;

  iArray1d NodIDNotOnFace(numNodesParentEle-numNodes);
  GetNodIndicesNotOnFace(domain, kf, false, NodIDNotOnFace);

  for (localIndex iCol=0;iCol<numNodesParentEle;++iCol)
  {
    for (int iDim=0; iDim<dim; ++iDim)
    {
      if(iCol<numNodes)
      {
        PermutedVector(dim*iCol+iDim) = Vector(dim*FaceConnToElemConn(iCol)+iDim);
      }
      else
      {
        int kNod = iCol-numNodes;
        PermutedVector(dim*iCol+iDim) = Vector(dim*NodIDNotOnFace(kNod)+iDim);
      }
    }
  }

}

realT LagrangeSmallStrainLinearElastic::GetSlipStickState ( const PhysicalDomainT& domain,
                                                            const rArray1d& trac,
                                                            const rArray2d& P,
                                                            int& stick)
{
  rArray1d tracNT_trial(dim);

  realT tracYield, phiTrial;
  if(domain.m_contactManager.m_sliding_law == 1)
  {
    tracYield = domain.m_contactManager.m_yield_traction;
  }

  MultiplyArray(P,trac,tracNT_trial);

  if(domain.m_contactManager.m_sliding_law == 2)
  {
    tracYield = std::fabs(domain.m_contactManager.m_coulomb_coefficient*tracNT_trial(0));
  }

  if(dim==2)
    phiTrial = std::fabs(tracNT_trial(1)) - tracYield;
  else if(dim==3)
    phiTrial = sqrt(pow(tracNT_trial(1),2) + pow(tracNT_trial(2),2)) - tracYield;

   if(phiTrial <= 1e-12) // phiTrial <=0 ===> stick
   {
     stick = 1;
   }
   else                 // else ===> slip
   {
     stick = 0;
   }
   return phiTrial;
}

void LagrangeSmallStrainLinearElastic::GetUpdatedTractionsOpeningMode( const rArray2d& P,
                                                                       const localIndex iGp,
                                                                       const realT tol,
                                                                       iArray1d& openingGP,
                                                                       rArray1d& trac)
{
  rArray1d tracNT_trial(dim);
  MultiplyArray(P,trac,tracNT_trial);

  int signTracN = 0;
  if(std::fabs(tracNT_trial(0))!=0)
  {
    signTracN = tracNT_trial(0)/(std::fabs(tracNT_trial(0)));
  }

  if(signTracN>0 && std::fabs(tracNT_trial(0))>tol)
  {
    openingGP(iGp) = 1;
    trac = 0;
  }
}

void LagrangeSmallStrainLinearElastic::GetUpdatedTractionsAndPlasticSlip ( const localIndex kf,
                                                                           const rArray2d& P,
                                                                           const realT& phiTrial,
                                                                           const int stick,
                                                                           const std::string FieldNamePrevious,
                                                                           const std::string FieldNameCurrent,
                                                                           const rArray2d& alphaXYZ,
                                                                           PhysicalDomainT& domain,
                                                                           rArray1d& tracXYZ)
{
  const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

  Array1dT<R1Tensor>& uJumpPl_Previous = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNamePrevious);
  Array1dT<R1Tensor>& uJumpPl_Current  = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNameCurrent);

  if(stick == 0)
  {
    // Get PTranspose
    rArray2d PTranspose(dim,dim);
    GetMatrixTranspose (P, PTranspose);

    // alpha_NT = P*alpha_XYZ*P'
    // RotateMatrix returns alpha_NT = (P')'*alphaXYZ*P' = P*alphaPar*P'
    rArray2d  alpha_NT(dim,dim);
    RotateMatrix(PTranspose, alphaXYZ, alpha_NT);


    realT delGam;
    // Since penalty/stabilization is chosen identically in both tangential directions in 3-D
    delGam = phiTrial/alpha_NT(1,1);

    // Negative of tracXYZ is used for updating plastic slips as tracXYZ is a penalization force that is in the opposite direction of the tendency to slip
//    if(dim==2)
      GetNegativeArray(tracXYZ);

    rArray1d tracNT(dim), tracNT_trial(dim);
    MultiplyArray(P,tracXYZ,tracNT_trial);

    // Normal direction, the traction is the same as the trial traction
    tracNT(0) = tracNT_trial(0);

    if(dim==2)
    {
      tracNT(1) = tracNT_trial(1) - delGam*alpha_NT(1,1)*(tracNT_trial(1)/std::fabs(tracNT_trial(1)));
      uJumpPl_Current(faceToExternalFaceMap(kf))(1) = uJumpPl_Previous(faceToExternalFaceMap(kf))(1) + delGam*(tracNT_trial(1)/std::fabs(tracNT_trial(1)));
    }
    else if(dim==3)
    {
      realT tracTMag = sqrt(pow(tracNT_trial(1),2) + pow(tracNT_trial(2),2));
      tracNT(1) = tracNT_trial(1)*(1 - (delGam*alpha_NT(1,1)/tracTMag));
      tracNT(2) = tracNT_trial(2)*(1 - (delGam*alpha_NT(1,1)/tracTMag));

      uJumpPl_Current(faceToExternalFaceMap(kf))(1) = uJumpPl_Previous(faceToExternalFaceMap(kf))(1) + delGam*(tracNT_trial(1)/tracTMag);
      uJumpPl_Current(faceToExternalFaceMap(kf))(2) = uJumpPl_Previous(faceToExternalFaceMap(kf))(2) + delGam*(tracNT_trial(2)/tracTMag);
    }

    // Rotate tracNT to get the traction after slip
    tracXYZ = 0;
    MultiplyArray(PTranspose,tracNT,tracXYZ);

    // Change the direction of obtained traction to be opposite to that of slip
//    if(dim==2)
      GetNegativeArray(tracXYZ);
  }
  else
  {
    uJumpPl_Current(faceToExternalFaceMap(kf)) = uJumpPl_Previous(faceToExternalFaceMap(kf));
  }

}

void LagrangeSmallStrainLinearElastic::GetNegativeArray (rArray1d& myVec)
{
  for (localIndex iSize = 0; iSize<myVec.size(); ++iSize)
  {
    myVec(iSize) = -myVec(iSize);
  }
}

void LagrangeSmallStrainLinearElastic::GetMatrixTranspose (const rArray2d& myMat,
                                                           rArray2d& myMatTransposed)
{
  for (localIndex iRow = 0; iRow < myMat.Dimension(0); ++iRow)
  {
    for (localIndex iCol = 0; iCol < myMat.Dimension(1); ++iCol)
    {
      myMatTransposed(iRow, iCol) = myMat(iCol, iRow);
    }
  }
}

void LagrangeSmallStrainLinearElastic::GetUpdatedTangentModulus (const PhysicalDomainT& domain,
                                                                 const rArray2d& P,
                                                                 const iArray1d& FaceConnToElemConn,
                                                                 const localIndex kf,
                                                                 rArray1d& trac,
                                                                 rArray2d& alpha,
                                                                 rArray2d& nDotDB,
                                                                 const bool parFlag)
{
  GetUpdatedStiffnessPenalty(domain, P, trac, alpha, parFlag);
  if(domain.m_contactManager.m_nitsche_active)
  {
    if((dim==2 and domain.m_contactManager.m_nitsche_tau1_active) or (dim==3 and domain.m_contactManager.m_nitsche_tau1_active and domain.m_contactManager.m_nitsche_tau2_active))
    {
      GetUpdatedStiffnessNitsche(domain, P, trac, nDotDB);
    }
  }
}

void LagrangeSmallStrainLinearElastic::UpdateTangentModulusForOpeningAndSliding(const PhysicalDomainT& domain,
                                                                                const rArray2d& P,
                                                                                const iArray1d& FaceConnToElemConn,
                                                                                const localIndex kf,
                                                                                rArray1d& trac,
                                                                                const int openingGP,
                                                                                const int stickGP,
                                                                                const bool parFlag,
                                                                                rArray2d& nDotDB,
                                                                                rArray2d& alpha,
                                                                                rArray2d& nDotDBPermute)
{
  if(openingGP==1)
  {
    alpha = 0; nDotDBPermute = 0;
  }
  else
  {
    if(domain.m_contactManager.m_sliding_law == 1 or domain.m_contactManager.m_sliding_law == 2)
    {

      if(stickGP==0)
      {
        GetUpdatedTangentModulus(domain, P, FaceConnToElemConn, kf, trac, alpha, nDotDB, parFlag);
        nDotDBPermute = 0;
        GetPermutedDBDotN(domain, nDotDB, FaceConnToElemConn, kf, nDotDBPermute);
      }
    }
  }
}
void LagrangeSmallStrainLinearElastic::GetUpdatedStiffnessPenalty ( const PhysicalDomainT& domain,
                                                                    const rArray2d& P,
                                                                    rArray1d& trac,
                                                                    rArray2d& alphaXYZ,
                                                                    const bool parFlag)
{
  // Get PTranspose
  rArray2d PTranspose(dim,dim);
  GetMatrixTranspose (P, PTranspose);

  // alpha_NT = P*alpha_XYZ*P'
  // RotateMatrix returns alpha_NT = (P')'*alphaXYZ*P' = P*alphaPar*P'
  rArray2d  alpha_NT(dim,dim);
  RotateMatrix(PTranspose, alphaXYZ, alpha_NT);

  GetNegativeArray(trac);

  rArray1d tracNT_trial;
  MultiplyArray(P,trac,tracNT_trial);
  int signtracTauTrial = tracNT_trial(1)/(std::fabs(tracNT_trial(1)));

  if(domain.m_contactManager.m_sliding_law == 1)
  {
    if(dim==2)
    {
      alpha_NT(1,1) = 0; // For Tresca friction
    }
    if(dim==3)
    {
      //     For Tresca friction
      realT m1, m2, tracTMag;
      tracTMag = sqrt(pow(tracNT_trial(1),2)+pow(tracNT_trial(2),2));
      m1 = tracNT_trial(1)/(tracTMag); m2 = tracNT_trial(2)/tracTMag;

      realT tracYield = domain.m_contactManager.m_yield_traction;
      realT kTStar = (tracYield/tracTMag);

      realT alpha_T = alpha_NT(1,1);

      alpha_NT(1,2) = kTStar*(alpha_NT(1,2)-alpha_T*m1*m2);
      alpha_NT(2,1) = kTStar*(alpha_NT(2,1)-alpha_T*m1*m2);
      alpha_NT(1,1) = kTStar*(alpha_NT(1,1) - alpha_T*m1*m1);
      alpha_NT(2,2) = kTStar*(alpha_NT(2,2) - alpha_T*m2*m2);
    }
  }
  if(domain.m_contactManager.m_sliding_law == 2)
  {
    realT mu = domain.m_contactManager.m_coulomb_coefficient;

    if(dim==2)
    {
        alpha_NT(1,0) = mu*signtracTauTrial*alpha_NT(0,0);
        alpha_NT(1,1) = 0;
    }
    if(dim==3)
    {
      realT m1, m2, tracTMag;
      tracTMag = sqrt(pow(tracNT_trial(1),2)+pow(tracNT_trial(2),2));
      m1 = tracNT_trial(1)/(tracTMag); m2 = tracNT_trial(2)/tracTMag;

      realT alpha_N = alpha_NT(0,0), alpha_T = alpha_NT(1,1);
      realT kTStar = mu*tracNT_trial(0)/(tracTMag);

      alpha_NT = 0;

      alpha_NT(0,0) = alpha_N;
      alpha_NT(1,0) = mu*m1*alpha_N;
      alpha_NT(1,1) = kTStar*(1-m1*m1)*alpha_T;
      alpha_NT(1,2) = -kTStar*(m1*m2)*alpha_T;
      alpha_NT(2,0) = mu*m2*alpha_N;
      alpha_NT(2,1) = -kTStar*(m1*m2)*alpha_T;
      alpha_NT(2,2) = kTStar*(1-m2*m2)*alpha_T;
    }
  }

  // Rotate alpha_NT back to alphaXYZ
  alphaXYZ = 0;
  RotateMatrix(P, alpha_NT, alphaXYZ);

  GetNegativeArray(trac);

//       TODO: Add code to recompute gradXYZ_elm and gradXYZ_sib by rotating gradNT_elm and gradNT_sib
//      alpha_t = 0;
//      gradNT_elm(2,:) = 0; gradNT_sib(2,:) = 0;
//
//
//      tracPar = P'*tracNT; tracSib = P'*tracNT; gradXY_elm = P'*gradNT_elm; gradXY_sib = P'*gradNT_sib;
//      alpha_e = P'*[alpha_n 0; 0 alpha_t]*P;
}

void LagrangeSmallStrainLinearElastic::GetUpdatedStiffnessNitsche ( const PhysicalDomainT& domain,
                                                                    const rArray2d& P,
                                                                    rArray1d& trac,
                                                                    rArray2d& nDotDB)
{
  rArray2d nDotDB_NT(nDotDB.Dimension(0), nDotDB.Dimension(1));
//  MultiplyArray(P, nDotDB, nDotDB_NT);

  GetNegativeArray(trac);

  rArray1d tracNT_trial;
  MultiplyArray(P,trac,tracNT_trial);

  int signTracTauTrial = tracNT_trial(1)/(std::fabs(tracNT_trial(1)));

  if(domain.m_contactManager.m_sliding_law == 1)
  {
    if(dim==2)
    {
      MultiplyArray(P, nDotDB, nDotDB_NT);

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(1,iCol) = 0.0;
      }
    }

    if(dim==3)
    {
      rArray2d nDotDB_NT_Ref(nDotDB.Dimension(0), nDotDB.Dimension(1));
      MultiplyArray(P, nDotDB, nDotDB_NT_Ref);

      realT m1, m2, tracTMag;
      tracTMag = sqrt(pow(tracNT_trial(1),2)+pow(tracNT_trial(2),2));
      m1 = tracNT_trial(1)/(tracTMag); m2 = tracNT_trial(2)/tracTMag;

      realT tracYield = domain.m_contactManager.m_yield_traction;
      realT kTStar = (tracYield/tracTMag);

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(0,iCol) = nDotDB_NT_Ref(0,iCol);
      }

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(1,iCol) = kTStar*(nDotDB_NT_Ref(1,iCol)*(1-m1*m1)-m1*m2*nDotDB_NT_Ref(2,iCol));
      }

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(2,iCol) = kTStar*(nDotDB_NT_Ref(2,iCol)*(1-m2*m2)-m1*m2*nDotDB_NT_Ref(1,iCol));
      }

      GetNegativeArray(trac);
    }
  }

  if(domain.m_contactManager.m_sliding_law == 2)
  {
    realT mu = domain.m_contactManager.m_coulomb_coefficient;

    if(dim==2)
    {
      rArray2d nDotDB_NT_Ref(nDotDB.Dimension(0), nDotDB.Dimension(1));
      MultiplyArray(P, nDotDB, nDotDB_NT_Ref);

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(0,iCol) = nDotDB_NT_Ref(0,iCol);
      }

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(1,iCol) = mu*signTracTauTrial*nDotDB_NT_Ref(0,iCol);
      }
    }
    if(dim==3)
    {
      rArray2d nDotDB_NT_Ref(nDotDB.Dimension(0), nDotDB.Dimension(1));
      MultiplyArray(P, nDotDB, nDotDB_NT_Ref);

      realT m1, m2, tracTMag;
      tracTMag = sqrt(pow(tracNT_trial(1),2)+pow(tracNT_trial(2),2));
      m1 = tracNT_trial(1)/(tracTMag); m2 = tracNT_trial(2)/tracTMag;

      realT kTStar = mu*tracNT_trial(0)/(tracTMag);

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(0,iCol) = nDotDB_NT_Ref(0,iCol);
      }

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(1,iCol) = mu*m1*nDotDB_NT_Ref(0,iCol) + kTStar*(nDotDB_NT_Ref(1,iCol)*(1-m1*m1)-m1*m2*nDotDB_NT_Ref(2,iCol));
      }

      for (localIndex iCol=0; iCol<nDotDB.Dimension(1); ++iCol)
      {
        nDotDB_NT(2,iCol) = mu*m2*nDotDB_NT_Ref(0,iCol) + kTStar*(nDotDB_NT_Ref(2,iCol)*(1-m2*m2)-m1*m2*nDotDB_NT_Ref(1,iCol));
      }

      GetNegativeArray(trac);
    }
  }

  rArray2d PTranspose(dim,dim);
  GetMatrixTranspose (P, PTranspose);

  nDotDB = 0;
  MultiplyArray(PTranspose, nDotDB_NT, nDotDB);
}

void LagrangeSmallStrainLinearElastic::UpdatePlasticStrainsForConvergedIteration( PhysicalDomainT& domain)
{
  rArray1d gauss(2);
  gauss[0] = -1 / sqrt(3);
  gauss[1] = 1 / sqrt(3);

  const bool nitsche_active = domain.m_contactManager.m_nitsche_active;

  unsigned int totContFaces;

//  if(dim==2 or dim==3)
//  {
//    totContFaces = domain.m_feFaceManager.m_numFaces;
//  }

  if(dim==2)
  {
    totContFaces = domain.m_feFaceManager.m_numFaces;
  }
  else
  {
    totContFaces = domain.m_contactManager.DataLengths();
  }

//        unsigned int totContFaces = domain.m_contactManager.DataLengths();

  for (localIndex iContFace =  0; iContFace < totContFaces; iContFace++)
  {
    bool contActiv = false;
    localIndex kf1, kf2;

//    if(dim==2 or dim==3)
//    {
//      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");
//
//      if (!(childFaceIndex[iContFace].empty()))
//      {
//        contActiv = true;
//        //          kf1 = childFaceIndex[iContFace][0];
//        //          kf2 = childFaceIndex[iContFace][1];
//        if(dim==2)
//        {
//          kf1 = childFaceIndex[iContFace][0];
//          kf2 = childFaceIndex[iContFace][1];
//        }
//        if(dim==3)
//        {
//          kf1 = iContFace;
//          kf2 = childFaceIndex[iContFace][0];
//        }
//      }
//    }

    if(dim==2)
    {
      const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");

      if (!(childFaceIndex[iContFace].empty()))
      {
        contActiv = true;
        kf1 = childFaceIndex[iContFace][0];
        kf2 = childFaceIndex[iContFace][1];
      }
    }
    else
    {
      const iArray1d& active = domain.m_contactManager.GetFieldData<int>("activeInit");
      const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
      const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");

      if(active[iContFace]!=0)
      {
        contActiv = true;
        bool fe;
        kf1 = domain.m_externalFaces.FaceIndex(f1[iContFace], fe);
        kf2 = domain.m_externalFaces.FaceIndex(f2[iContFace], fe);
      }
    }

//    const iArray1d& active = domain.m_contactManager.GetFieldData<int>("active");
//    const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
//    const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
//
//    if(active[iContFace]!=0)
//    {
//      contActiv = true;
//      bool fe;
//      kf1 = domain.m_externalFaces.FaceIndex(f1[iContFace], fe);
//      kf2 = domain.m_externalFaces.FaceIndex(f2[iContFace], fe);
//    }

    if(contActiv)
    {
      std::string s_previous_base = "uJumpPl_Previous_gp", s_current_base  = "uJumpPl_Current_gp";

      const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();
      const Array1dT< Array1dT< std::pair< ElementRegionT*, localIndex > > >& ftoe = domain.m_feFaceManager.m_toElementsRelation;
      const ElementRegionT* elemRegionPar = ftoe[kf1][0].first;
      const ElementRegionT* elemRegionSib = ftoe[kf2][0].first;

      const localIndex numNodesParentEle = elemRegionPar->m_numNodesPerElem;
      const localIndex numNodesSiblingEle = elemRegionSib->m_numNodesPerElem;

      const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

      rArray1d xe, psi, eta;
      lArray1d localParentFaceNodes;
      // Remark: Seems like this logic assumes faces are aligned. Should probably be changed for non-conforming meshes
      GetParentFaceNodesAndCoordsAndInterfaceGaussPoints(iContFace, domain, gauss, localParentFaceNodes, xe, psi, eta);

      rArray2d P_par(dim,dim), P_sib(dim,dim);
      GetTransformationTensor(domain, kf1, P_par);
      GetTransformationTensor(domain, kf2, P_sib);

      rArray2d alphaPar(dim,dim), alphaSib(dim,dim);

      GetInitialAlphaTensor(domain, kf1, P_par, alphaPar);
      GetInitialAlphaTensor(domain, kf2, P_sib, alphaSib);

      iArray1d stickGP_par(numNodes), stickGP_sib(numNodes), openingGP_par(numNodes), openingGP_sib(numNodes);
      iArray1d FaceConnToElemConnPar(numNodes), FaceConnToElemConnSib(numNodes);
      rArray2d normalVoigtPar(dim,0.5*dim*(dim+1)), normalVoigtSib(dim,0.5*dim*(dim+1));
      rArray2d DPar(0.5*dim*(dim+1),0.5*dim*(dim+1)), DSib(0.5*dim*(dim+1),0.5*dim*(dim+1));
      rArray2d BPar(0.5*dim*(dim+1),dim*numNodesParentEle), BSib(0.5*dim*(dim+1),dim*numNodesSiblingEle);
      rArray2d nParDotDParBPar(dim,dim*numNodesParentEle), nSibDotDSibBSib(dim,dim*numNodesSiblingEle);
      rArray2d nParDotDParBParPermute(dim,dim*numNodesParentEle), nSibDotDSibBSibPermute(dim,dim*numNodesSiblingEle);
      realT gamPar, gamSib;

      if(nitsche_active)
      {
        GetFaceConnToElemConnMap(domain, kf1, false, FaceConnToElemConnPar);
        GetFaceConnToElemConnMap(domain, kf2, true, FaceConnToElemConnSib);

        GetNormalVoigt(domain, kf1, normalVoigtPar);
        GetNormalVoigt(domain, kf2, normalVoigtSib);

        GetElasticityTensorVoigt(domain, kf1, DPar);
        GetElasticityTensorVoigt(domain, kf2, DSib);

        const rArray1d& nitscheGamma = domain.m_externalFaces.GetFieldData<realT>("nitscheGamma");

        gamPar = nitscheGamma(faceToExternalFaceMap(kf1));
        gamSib = nitscheGamma(faceToExternalFaceMap(kf2));
        if(fabs(((gamPar+gamSib)-1))>1e-12)
        {
          std::cout<<"WARNING:: Nitsche weights don't sum to unity!"<<std::endl;
        }

        if(elemRegionPar->m_elementGeometryID.compare(0,4,"STRI") ==0 or elemRegionPar->m_elementGeometryID.compare(0,4,"C3D4") == 0)
        {
          GetShapeFunctionDerivativeMatrixConstStr(domain, FaceConnToElemConnPar, kf1, BPar);
          GetShapeFunctionDerivativeMatrixConstStr(domain, FaceConnToElemConnSib, kf2, BSib);

          // GetDBDotN includes the Nitsche weights so should be renamed to GetGammaDBDotN?
          GetDBDotN(DPar, BPar, normalVoigtPar, numNodesParentEle, gamPar, domain, P_par, nParDotDParBPar);
          GetDBDotN(DSib, BSib, normalVoigtSib, numNodesParentEle, gamSib, domain, P_sib, nSibDotDSibBSib);

          GetPermutedDBDotN(domain, nParDotDParBPar, FaceConnToElemConnPar, kf1, nParDotDParBParPermute);
          GetPermutedDBDotN(domain, nSibDotDSibBSib, FaceConnToElemConnSib, kf2, nSibDotDSibBSibPermute);
        }
      }

      //RHS contribution from penalty and Nitsche terms
      for (localIndex a = 0; a < numNodes; ++a)
      {
        localIndex edgeLocalIndexRow, sibLocalIndexRow;

        GetLocalIndexOnInterface(iContFace, a, domain, localParentFaceNodes, edgeLocalIndexRow, sibLocalIndexRow);

        const int aDof = 2 * a * dim;

        for (int i = 0; i < dim; ++i)
        {
          for (localIndex iGp = 0; iGp < numNodes; ++iGp)
          {
            realT jcob_sub = 0;
            rArray2d N_sub(dim,dim*numNodes);
            rArray1d tracPar(dim), tracSib(dim);

            stickGP_par(iGp) = 1; stickGP_sib(iGp) = 1;
            openingGP_par(iGp) = 0; openingGP_sib(iGp) = 0;

            GetJacobianAndShapeFunctionsOnInterface(numNodes, xe, psi(iGp), eta(iGp), jcob_sub, N_sub);

            if(nitsche_active)
            {
              if(elemRegionPar->m_elementGeometryID.compare(0,4,"CPE4") ==0 or elemRegionPar->m_elementGeometryID.compare(0,4,"C3D8") ==0 )
              {
                GetShapeFunctionDerivativeMatrixQuadHex(domain, FaceConnToElemConnPar, kf1, psi(iGp), eta(iGp), BPar);
                GetShapeFunctionDerivativeMatrixQuadHex(domain, FaceConnToElemConnSib, kf2, psi(iGp), eta(iGp), BSib);

                // GetDBDotN includes the Nitsche weights so should be renamed to GetGammaDBDotN?
                GetDBDotN(DPar, BPar, normalVoigtPar, numNodesParentEle, gamPar, domain, P_par,  nParDotDParBPar);
                GetDBDotN(DSib, BSib, normalVoigtSib, numNodesParentEle, gamSib, domain, P_sib, nSibDotDSibBSib);

                GetPermutedDBDotN(domain, nParDotDParBPar, FaceConnToElemConnPar, kf1, nParDotDParBParPermute);
                GetPermutedDBDotN(domain, nSibDotDSibBSib, FaceConnToElemConnSib, kf2, nSibDotDSibBSibPermute);

              }
            }
            GetTractionInterfaceSliding(N_sub, alphaPar, alphaSib, P_par, P_sib, nParDotDParBParPermute, nSibDotDSibBSibPermute, localParentFaceNodes, kf1, kf2, iContFace,
                                        iGp, s_previous_base, s_current_base, domain, stickGP_par, stickGP_sib, openingGP_par, openingGP_sib, tracPar, tracSib);
//            // Gets the name of the field associated with Gauss point iGp for previous and current load steps
//            std::stringstream ss; ss << iGp;
//            std::string FieldNamePrevious = s_previous_base + ss.str(), FieldNameCurrent  = s_current_base + ss.str();
//
//            // Gets trial traction assuming crack faces are sticking
//            GetTractionInterfaceSlidingTrial(domain, elemRegion, N_sub, alphaPar, alphaSib, P_par, P_sib, nParDotDParBParPermute, nSibDotDSibBSibPermute, localParentFaceNodes,
//                                             kf1, kf2, iContFace, FieldNamePrevious, iGp, openingGP_par, openingGP_sib, tracPar, tracSib);
//
//            realT tol = domain.m_contactManager.m_traction_n_tol;
//            // Check if crack faces are opening, if yes, set traction to zero
//            GetUpdatedTractionsOpeningMode(P_par, iGp, tol, openingGP_par, tracPar);
//            GetUpdatedTractionsOpeningMode(P_sib, iGp, tol, openingGP_sib, tracSib);
//
//            if(openingGP_par(iGp)==1)
//            {
//              Array1dT<R1Tensor>& uJumpPl_Previous = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNamePrevious);
//              Array1dT<R1Tensor>& uJumpPl_Current  = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNameCurrent);
//
//              uJumpPl_Current(faceToExternalFaceMap(kf1)) = uJumpPl_Previous(faceToExternalFaceMap(kf1));
//            }
//            if(openingGP_par(iGp)==0)
//            {
//              // Gets the trial value value of yield function trial and a boolean flag indicating if the yield criterion is violated
//              realT phiTrialPar = GetSlipStickState(domain, tracPar, P_par, stickGP_par(iGp));
//              // Pulls back the trial traction to lie on the yield surface, updates stiffness and plastic slip variables (Simo and Hughes (1997)))
//              GetUpdatedTractionsAndPlasticSlip(kf1, P_par, phiTrialPar, stickGP_par(iGp), FieldNamePrevious, FieldNameCurrent, alphaPar, domain, tracPar);
//            }
//            if(openingGP_sib(iGp)==1)
//            {
//              Array1dT<R1Tensor>& uJumpPl_Previous = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNamePrevious);
//              Array1dT<R1Tensor>& uJumpPl_Current  = domain.m_externalFaces.GetFieldData<R1Tensor>(FieldNameCurrent);
//
//              uJumpPl_Current(faceToExternalFaceMap(kf2)) = uJumpPl_Previous(faceToExternalFaceMap(kf2));
//            }
//            if(openingGP_sib(iGp)==0)
//            {
//              realT phiTrialSib = GetSlipStickState(domain, tracSib, P_sib, stickGP_sib(iGp));
//              GetUpdatedTractionsAndPlasticSlip(kf2, P_sib, phiTrialSib, stickGP_sib(iGp), FieldNamePrevious, FieldNameCurrent, alphaSib, domain, tracSib);
//            }
          }
        }
      }
    }
  }

}

void LagrangeSmallStrainLinearElastic::StoreHistoryVariablesForCurrentLoadStepAndResetTheField( PhysicalDomainT& domain)
{
  Array1dT<R1Tensor>& uJumpPl_Previous_gp0 = domain.m_externalFaces.GetFieldData<R1Tensor>("uJumpPl_Previous_gp0");
  Array1dT<R1Tensor>& uJumpPl_Current_gp0  = domain.m_externalFaces.GetFieldData<R1Tensor>("uJumpPl_Current_gp0");

  Array1dT<R1Tensor>& uJumpPl_Previous_gp1 = domain.m_externalFaces.GetFieldData<R1Tensor>("uJumpPl_Previous_gp1");
  Array1dT<R1Tensor>& uJumpPl_Current_gp1  = domain.m_externalFaces.GetFieldData<R1Tensor>("uJumpPl_Current_gp1");

  Array1dT<R1Tensor>& uJumpPl_Previous_gp2 = domain.m_externalFaces.GetFieldData<R1Tensor>("uJumpPl_Previous_gp2");
  Array1dT<R1Tensor>& uJumpPl_Current_gp2  = domain.m_externalFaces.GetFieldData<R1Tensor>("uJumpPl_Current_gp2");

  Array1dT<R1Tensor>& uJumpPl_Previous_gp3 = domain.m_externalFaces.GetFieldData<R1Tensor>("uJumpPl_Previous_gp3");
  Array1dT<R1Tensor>& uJumpPl_Current_gp3  = domain.m_externalFaces.GetFieldData<R1Tensor>("uJumpPl_Current_gp3");

  uJumpPl_Previous_gp0 = uJumpPl_Current_gp0; uJumpPl_Previous_gp1 = uJumpPl_Current_gp1;
  uJumpPl_Previous_gp2 = uJumpPl_Current_gp2; uJumpPl_Previous_gp3 = uJumpPl_Current_gp3;

  uJumpPl_Current_gp0 = 0; uJumpPl_Current_gp1 = 0;
  uJumpPl_Current_gp2 = 0; uJumpPl_Current_gp3 = 0;
}

/*

void LagrangeSmallStrainLinearElastic::PostProcess (PhysicalDomainT& domain,
                                                    SpatialPartition& partition,
                                                    const sArray1d& namesOfSolverRegions)
{


  RegionMap::const_iterator
  region     = domain.m_feElementManager.m_ElementRegions.begin(),
  end_region = domain.m_feElementManager.m_ElementRegions.end();

  for(; region != end_region; ++region)
  {

    const ElementRegionT& elemRegion = region->second;

    rArray1d gauss(2);
    gauss[0] = -1 / sqrt(3);
    gauss[1] = 1 / sqrt(3);

    const bool nitsche_active = domain.m_contactManager.m_nitsche_active;
    unsigned int totContFaces;

    Array1dT<R1Tensor>& tracInt = domain.m_feFaceManager.GetFieldData<R1Tensor>("tracInt");
    Array1dT<R1Tensor>& interpenetration = domain.m_feFaceManager.GetFieldData<R1Tensor>("interpenetration");

    if(dim==2)
    {
      totContFaces = domain.m_feFaceManager.m_numFaces;
    }
    else
    {
      totContFaces = domain.m_contactManager.DataLengths();
    }

//      unsigned int totContFaces = domain.m_contactManager.DataLengths();

    for (localIndex iContFace =  0; iContFace < totContFaces; iContFace++)
    {
      bool contActiv = false;
      localIndex kf1, kf2;

      R1Tensor tracAvg, interpenetrationAvg;

      if(dim==2)
      {
        const OrderedVariableOneToManyRelation& childFaceIndex = domain.m_feFaceManager.GetVariableOneToManyMap("childIndices");

        if (!(childFaceIndex[iContFace].empty()))
        {
          contActiv = true;
          kf1 = childFaceIndex[iContFace][0];
          kf2 = childFaceIndex[iContFace][1];
        }
      }
      else
      {
        const iArray1d& active = domain.m_contactManager.GetFieldData<int>("active");
        const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
        const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");

        if(active[iContFace]!=0)
        {
          contActiv = true;
          bool fe;
          kf1 = domain.m_externalFaces.FaceIndex(f1[iContFace], fe);
          kf2 = domain.m_externalFaces.FaceIndex(f2[iContFace], fe);
        }
      }

//          const iArray1d& active = domain.m_contactManager.GetFieldData<int>("active");
//          const lArray1d& f1 = domain.m_contactManager.GetFieldData<localIndex>("face1");
//          const lArray1d& f2 = domain.m_contactManager.GetFieldData<localIndex>("face2");
//
//          if(active[iContFace]!=0)
//          {
//            contActiv = true;
//            bool fe;
//            kf1 = domain.m_externalFaces.FaceIndex(f1[iContFace], fe);
//            kf2 = domain.m_externalFaces.FaceIndex(f2[iContFace], fe);
//          }

      if(contActiv)
      {
        const localIndex numNodes = domain.m_feFaceManager.m_toNodesRelation[kf1].size();
        const localIndex numNodesParentEle = elemRegion.m_numNodesPerElem;
        const lArray1d& faceToExternalFaceMap = domain.m_feFaceManager.GetFieldData<localIndex>("externalFaceIndex");

        rArray1d xe, psi, eta;
        lArray1d localParentFaceNodes;
        // Remark: Seems like this logic assumes faces are aligned. Should probably be changed for non-conforming meshes
        GetParentFaceNodesAndCoordsAndInterfaceGaussPoints(iContFace, domain, gauss, localParentFaceNodes, xe, psi, eta);

        rArray2d P_par(dim,dim), P_sib(dim,dim);
        GetTransformationTensor(domain, kf1, P_par);
        GetTransformationTensor(domain, kf2, P_sib);

        rArray2d alphaPar(dim,dim), alphaSib(dim,dim);
        GetInitialAlphaTensor(domain, kf1, P_par, alphaPar);
        GetInitialAlphaTensor(domain, kf2, P_sib, alphaSib);

        iArray1d FaceConnToElemConnPar(numNodes), FaceConnToElemConnSib(numNodes);
        rArray2d normalVoigtPar(dim,0.5*dim*(dim+1)), normalVoigtSib(dim,0.5*dim*(dim+1));
        rArray2d DPar(0.5*dim*(dim+1),0.5*dim*(dim+1)), DSib(0.5*dim*(dim+1),0.5*dim*(dim+1));
        rArray2d BPar(0.5*dim*(dim+1),dim*numNodesParentEle), BSib(0.5*dim*(dim+1),dim*numNodesParentEle);
        rArray2d nParDotDParBPar(dim,dim*numNodesParentEle), nSibDotDSibBSib(dim,dim*numNodesParentEle);
        rArray2d nParDotDParBParPermute(dim,dim*numNodesParentEle), nSibDotDSibBSibPermute(dim,dim*numNodesParentEle);
        realT gamPar, gamSib;

        if(nitsche_active)
        {
          GetFaceConnToElemConnMap(domain, elemRegion, kf1, false, FaceConnToElemConnPar);
          GetFaceConnToElemConnMap(domain, elemRegion, kf2, true, FaceConnToElemConnSib);

          GetNormalVoigt(domain, kf1, normalVoigtPar);
          GetNormalVoigt(domain, kf2, normalVoigtSib);

          GetElasticityTensorVoigt(domain, elemRegion, kf1, DPar);
          GetElasticityTensorVoigt(domain, elemRegion, kf2, DSib);

          const rArray1d& nitscheGamma = domain.m_externalFaces.GetFieldData<realT>("nitscheGamma");

          gamPar = nitscheGamma(faceToExternalFaceMap(kf1));
          gamSib = nitscheGamma(faceToExternalFaceMap(kf2));
          if(fabs(((gamPar+gamSib)-1))>1e-12)
          {
            std::cout<<"WARNING:: Nitsche weights don't sum to unity!"<<std::endl;
          }

          if(elemRegion.m_elementGeometryID.compare(0,4,"STRI") ==0 or elemRegion.m_elementGeometryID.compare(0,4,"C3D4") == 0)
          {
            GetShapeFunctionDerivativeMatrixConstStr(domain, elemRegion, FaceConnToElemConnPar, kf1, BPar);
            GetShapeFunctionDerivativeMatrixConstStr(domain, elemRegion, FaceConnToElemConnSib, kf2, BSib);

            // GetDBDotN includes the Nitsche weights so should be renamed to GetGammaDBDotN?
            GetDBDotN(DPar, BPar, normalVoigtPar, numNodesParentEle, gamPar, domain, P_par, nParDotDParBPar);
            GetDBDotN(DSib, BSib, normalVoigtSib, numNodesParentEle, gamSib, domain, P_sib, nSibDotDSibBSib);

            GetPermutedDBDotN(domain, elemRegion, nParDotDParBPar, FaceConnToElemConnPar, kf1, nParDotDParBParPermute);
            GetPermutedDBDotN(domain, elemRegion, nSibDotDSibBSib, FaceConnToElemConnSib, kf2, nSibDotDSibBSibPermute);
          }
        }

        for (int i = 0; i < dim; ++i)
        {
          for (localIndex iGp = 0; iGp < numNodes; ++iGp)
          {
            realT jcob_sub = 0;
            rArray2d N_sub(dim,dim*numNodes);
            rArray1d tracPar(dim), tracSib(dim);

            GetJacobianAndShapeFunctionsOnInterface(numNodes, xe, psi(iGp), eta(iGp), jcob_sub, N_sub);

            if(nitsche_active)
            {
              if(elemRegion.m_elementGeometryID.compare(0,4,"CPE4") ==0 or elemRegion.m_elementGeometryID.compare(0,4,"C3D8") ==0 )
              {
                GetShapeFunctionDerivativeMatrixQuadHex(domain, elemRegion, FaceConnToElemConnPar, kf1, psi(iGp), eta(iGp), BPar);
                GetShapeFunctionDerivativeMatrixQuadHex(domain, elemRegion, FaceConnToElemConnSib, kf2, psi(iGp), eta(iGp), BSib);

                // GetDBDotN includes the Nitsche weights so should be renamed to GetGammaDBDotN?
                GetDBDotN(DPar, BPar, normalVoigtPar, numNodesParentEle, gamPar, domain, P_par,  nParDotDParBPar);
                GetDBDotN(DSib, BSib, normalVoigtSib, numNodesParentEle, gamSib, domain, P_sib, nSibDotDSibBSib);

                GetPermutedDBDotN(domain, elemRegion, nParDotDParBPar, FaceConnToElemConnPar, kf1, nParDotDParBParPermute);
                GetPermutedDBDotN(domain, elemRegion, nSibDotDSibBSib, FaceConnToElemConnSib, kf2, nSibDotDSibBSibPermute);

              }
            }

            GetTractionInterface(domain, elemRegion, N_sub, alphaPar, alphaSib, nParDotDParBParPermute, nSibDotDSibBSibPermute, localParentFaceNodes, kf1, kf2, iContFace, tracPar, tracSib);

            rArray1d uJump(dim);
            GetDisplacementJumpInterface(domain, elemRegion, kf1, kf2, N_sub, iContFace, localParentFaceNodes, uJump);

            tracAvg[i] += (1./numNodes)*tracPar(i);
            interpenetrationAvg[i] += (1./numNodes)*uJump(i);
          }
        }
        tracInt[kf1] = tracAvg;
        interpenetration[kf1] = interpenetrationAvg;
      }
    }
  }
}

*/


/*

void LagrangeSmallStrainLinearElastic::ProcessElementRegion( NodeManagerT& nodeManager,
                                                                  ElementRegionT& elemRegion,
                                                                  const realT dt )
{


  static Array1dT< R1Tensor > u_local;
  static Array1dT< R1Tensor > uhat_local;
  static Array1dT<R1Tensor> f_local;
  static Array1dT< R1Tensor > x;
  static Array1dT< R1Tensor > v;
  static Array1dT<R1Tensor> dNdx;
  static Array1dT<R1Tensor> f_zemc;
  static Array1dT<R1Tensor> Q;

  Array1dT< Array2dT<R1Tensor> >& m_dNdX = elemRegion.m_dNdX;
  EnergyT& m_energy = elemRegion.m_energy;
  FiniteElement<nsdof>*& m_finiteElement = elemRegion.m_finiteElement;
  Array1dT<MaterialBaseStateDataT*>& m_matStates = elemRegion.m_matStates;
  MaterialBaseParameterDataT*& m_matParameters = elemRegion.m_matParameters;
  MaterialBaseT*& m_matComputations = elemRegion.m_matComputations;

  const unsigned int numNodesPerElem = elemRegion.m_numNodesPerElem;

  if( u_local.size() != numNodesPerElem )
  {
    u_local.resize(numNodesPerElem);
    uhat_local.resize(numNodesPerElem);
    f_local.resize(numNodesPerElem);
    x.resize(numNodesPerElem);
    v.resize(numNodesPerElem);
    dNdx.resize(numNodesPerElem);
    f_zemc.resize(numNodesPerElem);
    Q.resize(elemRegion.m_finiteElement->zero_energy_modes());
  }


  const Array1dT<R1Tensor>& incrementalDisplacement = nodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
  const Array1dT<R1Tensor>& totalDisplacement = nodeManager.GetFieldData<FieldInfo::displacement>();

//  const iArray1d& attachedToSendingGhostNode = GetFieldData<int>("attachedToSendingGhostNode");


  rArray1d& volume = elemRegion.GetFieldData<FieldInfo::volume>();
  rArray1d& volume_n = elemRegion.GetFieldData<realT>("volume_n");

  volume_n = volume;




  const Array1dT<R1Tensor>& referencePosition = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& velocity = nodeManager.GetFieldData<FieldInfo::velocity>();
  Array1dT<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();
  Array1dT<R1Tensor>& hgforce = nodeManager.GetFieldData<FieldInfo::hgforce> ();

//  hgforce = 0.0;
  Array1dT<Array1dT<R1Tensor>*> Qstiffness(elemRegion.m_finiteElement->zero_energy_modes(),NULL);
  if( elemRegion.m_finiteElement->zero_energy_modes() >= 1 )
  {
    Qstiffness[0] = &(elemRegion.GetFieldData<R1Tensor>("Qhg1"));
  }
  if( elemRegion.m_finiteElement->zero_energy_modes() >= 2 )
  {
    Qstiffness[1] = &(elemRegion.GetFieldData<R1Tensor>("Qhg2"));
  }
  if( elemRegion.m_finiteElement->zero_energy_modes() >= 3 )
  {
    Qstiffness[2] = &(elemRegion.GetFieldData<R1Tensor>("Qhg3"));
  }
  if( elemRegion.m_finiteElement->zero_energy_modes() >= 4 )
  {
    Qstiffness[3] = &(elemRegion.GetFieldData<R1Tensor>("Qhg4"));
  }

  const iArray1d& ghostRank = elemRegion.GetFieldData<FieldInfo::ghostRank>();



  m_energy.Zero();

  if( elemRegion.m_numIntegrationPointsPerElem == 1 )
  {
    CalculateElementStresses( nodeManager,
                              elementManager );
  }


  for( localIndex k=0 ; k<elemRegion.m_numElems ; ++k )
  {
//    if( ghostRank[k] < 0 )
    {

      const localIndex* const elemToNodeMap = elemRegion.m_toNodesRelation[k];


      CopyGlobalToLocal( elemToNodeMap,
                         incrementalDisplacement, totalDisplacement,
                         uhat_local, u_local );


      if( elemRegion.m_finiteElement->zero_energy_modes() )
      {
        CopyGlobalToLocal( elemToNodeMap,
                           referencePosition, velocity,
                           x, v );

        x += u_local;
      }


      volume[k] = 0.0;

      for( unsigned int a=0 ; a<elemRegion.m_numIntegrationPointsPerElem ; ++a )
      {

        R1Tensor* const dNdX = elemRegion.m_dNdX(k)[a];

        // Velocity Gradient

        volume[k] += elemRegion.m_detJ(k,a);

        Finv.Inverse(F);

        // Calculate Rate of deformation tensor and rotationAxis tensor


        // nodal force calculation


        totalStress = m_matStates[k][a].devStress;
        totalStress.PlusIdentity(m_matStates[k][a].pressure);

        f_local = 0.0;


        FiniteElementUtilities::Integrate( totalStress,
                                           m_dNdX(k)[a],
                                           m_detJ(k,a),
                                           detF,
                                           Finv,
                                           f_local.size(),
                                           f_local.data() );
        realT BB = 0.0;
        for( unsigned int b=0 ; b<numNodesPerElem ; ++b )
        {
          dNdx(b).AijBi(Finv,dNdX[b]);
          BB += Dot( dNdx(b), dNdx(b) );
        }
        realT thisdt = m_finiteElement->CalculateMaxStableExplicitTimestep( m_matParameters->init_density / fabs(detF),
                                                                            m_matParameters->Lame + 2*m_matParameters->ShearModulus,
                                                                            BB );

        if( elemRegion.m_ElementDimension == 3 )
        {
          thisdt /= sqrt(2.0);
        }

        if( thisdt < SolverBase::m_stabledt.m_maxdt )
        {
  //        timeStep.m_region = this->m_regionName;
  //        timeStep.m_index = k;
          SolverBase::m_stabledt.m_maxdt = thisdt;
        }


      }

      if( m_finiteElement->zero_energy_modes() )
      {

        for( int m=0 ; m<m_finiteElement->zero_energy_modes() ; ++m )
        {
          Q[m] = (*(Qstiffness[m]))[k];
        }
        m_finiteElement->zero_energy_mode_control( dNdx, volume[k], x, v,
                                                   elemRegion.m_hgDamp,
                                                   elemRegion.m_hgStiff*dt,
                                                   m_matParameters->init_density,
                                                   m_matParameters->Lame + 2*m_matParameters->ShearModulus ,
                                                   Q, f_zemc );

        AddLocalToGlobal( elemToNodeMap,
                          f_zemc, f_zemc,
                          force, hgforce);

      }

      AddLocalToGlobal( elemToNodeMap, f_local, force);
    }
  }

}
*/




/* Explicit Instantiations */



/* Register solver in the solver factory */

SolverRegistrator<LagrangeSmallStrainLinearElastic > reg_LagrangeSmallStrainLinearElastic;
