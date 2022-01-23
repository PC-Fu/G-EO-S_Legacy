/*
 * LagrangeLargeStrain.cpp
 *
 *  Created on: Nov 8, 2012
 *      Author: settgast1
 */

#include "PhysicsSolvers/Lagrange/LagrangeHelperFunctions.h"
#include "LagrangeLargeStrain.h"
#include "PhysicsSolvers/SolverFactory.h"
//#include "PhysicsSolvers/PhysicsSolverStrings.h"
#include "DataStructures/VectorFields/ElementRegionT.h"
#include "ElementLibrary/FiniteElement.h"
#include "ElementLibrary/FiniteElementUtilities.h"

#include "Utilities/Kinematics.h"


LagrangeLargeStrain::LagrangeLargeStrain(  const std::string& name,
                                                ProblemManagerT* const pm ):
LagrangeSolverBase(name,pm)
{
  // TODO Auto-generated constructor stub

}


LagrangeLargeStrain::~LagrangeLargeStrain()
{
  // TODO Auto-generated destructor stub
}




void LagrangeLargeStrain::ProcessElementRegion( NodeManagerT& nodeManager,
                                                     ElementRegionT& elemRegion,
                                                     const realT dt )
{

  R2Tensor A;
  R2Tensor F;
  R2Tensor Finv;
  R2Tensor dUhatdX;
  R2SymTensor totalStress;
  R2SymTensor Dadt;
  R2Tensor Rot;


  static Array1dT< R1Tensor > u_local;
  static Array1dT< R1Tensor > uhat_local;
  static Array1dT<R1Tensor> f_local;
  static Array1dT< R1Tensor > x;
  static Array1dT< R1Tensor > v;
  static Array1dT<R1Tensor> s_dNdx;
  static Array1dT<R1Tensor> f_zemc;
  static Array1dT<R1Tensor> Q;

  rArray2d& detJ = elemRegion.m_detJ;
  rArray2d& detJ_n = elemRegion.m_detJ_n;
  rArray2d& detJ_np1 = elemRegion.m_detJ_np1;

  Array2dT< R2Tensor >&  dUdX = elemRegion.m_dUdX;
  FiniteElementBase*& finiteElement = elemRegion.m_finiteElement;

  //  MaterialBaseT*& m_materialComputations = elemRegion.m_materialComputations;

  const unsigned int numNodesPerElem = elemRegion.m_numNodesPerElem;

  if( u_local.size() != numNodesPerElem )
  {
    u_local.resize(numNodesPerElem);
    uhat_local.resize(numNodesPerElem);
    f_local.resize(numNodesPerElem);
    x.resize(numNodesPerElem);
    v.resize(numNodesPerElem);
    s_dNdx.resize(numNodesPerElem);
    f_zemc.resize(numNodesPerElem);
    Q.resize(finiteElement->zero_energy_modes());
  }


  const Array1dT<R1Tensor>& incrementalDisplacement = nodeManager.GetFieldData<FieldInfo::incrementalDisplacement>();
  const Array1dT<R1Tensor>& totalDisplacement = nodeManager.GetFieldData<FieldInfo::displacement>();

//  const iArray1d& attachedToSendingGhostNode = GetFieldData<int>("attachedToSendingGhostNode");


  rArray1d& volume = elemRegion.GetFieldData<FieldInfo::volume>();
  rArray1d& volume_n = elemRegion.GetFieldData<realT>("volume_n");

  volume_n = volume;
  detJ_n = detJ_np1;








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

//  const iArray1d& ghostRank = elemRegion.GetFieldData<FieldInfo::ghostRank>();



  elemRegion.m_energy.Zero();




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

      f_local = 0.0;

//      if(0)

      const localIndex paramIndex = elemRegion.m_mat->NumParameterIndex0() > 1 ? k : 0 ;
      const MaterialBaseParameterData& param = *( elemRegion.m_mat->ParameterData(paramIndex) );
      for( unsigned int a=0 ; a<elemRegion.m_numIntegrationPointsPerElem ; ++a )
      {
        MaterialBaseStateData& state = *(elemRegion.m_mat->StateData(k,a));
        R1Tensor* const dNdX = elemRegion.m_dNdX(k)[a];

        // Velocity Gradient

        // calculate dUhat/dX at beginning of step
        CalculateGradient( dUhatdX ,uhat_local, dNdX );

        // calculate velocity gradient (mid-step)
        R2Tensor L;
	if(dt > 0)
        {
          // calculate dv/dX
          R2Tensor dvdX = dUhatdX;
          dvdX *= 1.0 / dt;

          // calculate du/dX
          F = dUhatdX;
          F *= 0.5;
          F += dUdX(k,a);
          F.PlusIdentity(1.0);

          // calculate dX/du
          Finv.Inverse(F);

          // chain rule: calculate dv/du = dv/dX * dX/du
          L.AijBjk(dvdX, Finv);
        }

        // calculate gradient (end of step)
        dUdX(k,a) += dUhatdX;
        F = dUdX(k,a);
        F.PlusIdentity(1.0);
        realT detF = F.Det();

        // calculate element volume
        detJ_np1(k,a) = detJ(k,a) * detF;
        volume[k] += detJ_np1(k,a);

        Finv.Inverse(F);

        // Calculate Rate of deformation tensor and rotationAxis tensor
        A.AijBjk(dUhatdX,Finv);
        IncrementalKinematics(A,Dadt,Rot);

        const realT rho = param.init_density / fabs(detF);

        // update state before exercising material model
        if(elemRegion.m_mat->NeedsDensity() ){
          state.SetDensity(rho);
        }
        if(elemRegion.m_mat->NeedsSpecificInternalEnergy() ){

          state.TotalStress(totalStress);
          realT ie = state.GetSpecificInternalEnergy();
          realT die = L(0,0) * totalStress(0,0)
                    + (L(0,1) + L(1,0)) * totalStress(0,1)
                    + (L(0,2) + L(2,0)) * totalStress(0,2)
                    + L(1,1)*totalStress(1,1)
                    + (L(1,2) + L(2,1)) * totalStress(1,2)
                    + L(2,2)*totalStress(2,2);
           die *= -dt/(rho+1e-64); // negative because stress is +ve in compression
           ie += die;
          state.SetSpecificInternalEnergy(ie);
        }

        // Material Update

        elemRegion.m_mat->StrainDrivenUpdateMember( k, a,
                                                    Dadt, L, Rot,
                                                    detJ_n[k][a],
                                                    detJ_np1[k][a],
                                                    dt);

        // nodal force calculation
        state.TotalStress(totalStress);

        //const realT rho = param.init_density / fabs(detF);
        const realT soundSpeed = sqrt( state.BulkModulus / rho );
        const realT trD = dt > 0 ? Dadt.Trace() / dt : 0.0;


        realT bulkQ = LagrangeHelperFunctions::BulkQ( rho,
                                                soundSpeed,
                                                this->m_bulkQLinear,
                                                this->m_bulkQQuadratic,
                                                trD,
                                                cbrt( volume[k] ) );

        totalStress.PlusIdentity( bulkQ );
        FiniteElementUtilities::Integrate( totalStress,
                                           elemRegion.m_dNdX(k)[a],
                                           detJ(k,a),
                                           detF,
                                           Finv,
                                           f_local.size(),
                                           f_local.data() );
        realT BB = 0.0;
        for( unsigned int b=0 ; b<numNodesPerElem ; ++b )
        {
          s_dNdx(b).AijBi(Finv,dNdX[b]);
          BB += Dot( s_dNdx(b), s_dNdx(b) );
        }
        realT thisdt = LagrangeHelperFunctions::CalculateMaxStableExplicitTimestep( param.init_density / fabs(detF),
                                                                            param.Lame + 2*param.init_shearModulus,
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

      if( finiteElement->zero_energy_modes() )
      {

        for( int m=0 ; m<finiteElement->zero_energy_modes() ; ++m )
        {
          Q[m] = (*(Qstiffness[m]))[k];
        }

        finiteElement->zero_energy_mode_control( s_dNdx, volume[k], x, v,
                                                   elemRegion.m_hgDamp,
                                                   elemRegion.m_hgStiff*dt,
                                                   param.init_density,
                                                   param.Lame + 2*param.init_shearModulus ,
                                                   dt, Q, f_zemc );
        for( int m=0 ; m<finiteElement->zero_energy_modes() ; ++m )
        {
          (*(Qstiffness[m]))[k] = Q[m];
        }

        AddLocalToGlobal( elemToNodeMap,
                          f_zemc, f_zemc,
                          force, hgforce);

      }

      AddLocalToGlobal( elemToNodeMap, f_local, force);
    }
  }

}

// Rui Wang, apply gravity on elements
void LagrangeLargeStrain::ApplyGravity( NodeManagerT& nodeManager,
                                        ElementRegionT& elemRegion,
                                        const realT dt )
{
  Array1dT<R1Tensor>& force = nodeManager.GetFieldData<FieldInfo::force>();
  Array1dT<realT>& mass = elemRegion.GetFieldData<FieldInfo::mass> ();
  Array1dT<R1Tensor> gravityForce;
  gravityForce.resize(elemRegion.m_numNodesPerElem);
  gravityForce *= 0.0;
  R1Tensor gravityForceElem;

  for( localIndex k=0 ; k < elemRegion.m_numElems ; ++k )
  {
    const localIndex* const elemToNodeMap = elemRegion.m_toNodesRelation[k];
    gravityForceElem = m_gravityVector * mass[k];
    for (localIndex i=0 ; i < elemRegion.m_numNodesPerElem ; ++i)
    {
      gravityForce[i] = gravityForceElem / elemRegion.m_numNodesPerElem;
    }
    AddLocalToGlobal( elemToNodeMap, gravityForce, force);
  }

}



/* Register solver in the solver factory */

SolverRegistrator<LagrangeLargeStrain> reg_LagrangeLargeStrain;
