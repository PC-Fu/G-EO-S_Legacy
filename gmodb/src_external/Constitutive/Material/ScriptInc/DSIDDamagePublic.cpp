//FUNCTION_BEGIN_PARSE
void
DSIDDamage::InitializeStates( const localIndex index )
{

//  for( localIndex a=0 ; a<m_stateData.Dimension(1) ; ++a )
//  {
//    const localIndex paramIndex = m_parameterData.size() > 1 ? index : 0 ;
//
//    m_stateData[index][a].density =  m_parameterData[paramIndex].init_density;
//    m_stateData[index][a].BulkModulus = m_parameterData[paramIndex].init_bulkModulus;
//    m_stateData[index][a].ShearModulus = m_parameterData[paramIndex].init_shearModulus;
//    m_stateData[index][a].ElasticBulkModulus = m_parameterData[paramIndex].init_bulkModulus;
//    m_stateData[index][a].ElasticShearModulus = m_parameterData[paramIndex].init_shearModulus;
//  }

}

//FUNCTION_BEGIN_PARSE
void
DSIDDamage::StrainDrivenUpdateMember( const localIndex index0,
                            const localIndex index1,
                            const R2SymTensorT < 3 >& Ddt,
                            const R2TensorT < 3 >& L,
                            const R2Tensor& Rot,
                            const realT& volume_n,
                            const realT& volume_np1,
                            const realT dt )
{
  StrainDrivenUpdateMember( index0, index1, Ddt, L, Rot, dt );
}

//FUNCTION_BEGIN_PARSE
void
DSIDDamage::StrainDrivenUpdateMember( const localIndex index0,
                                      const localIndex index1,
                                      const R2SymTensorT < 3 >& Ddt,
                                      const R2TensorT < 3 >& L,
                                      const R2Tensor& Rot,
                                      const realT dt )
{
  const localIndex paramIndex = m_parameterData.size() > 1 ? index0 : 0 ;
  const DSIDDamageParameterData& matParams = m_parameterData[ paramIndex ];
  DSIDDamageStateData& matState = m_stateData(index0,index1);


  realT& pressure = matState.pressure;
  realT& fd0 = matState.fd0; // to add
  R2SymTensorT < 3 >& devStress = matState.devStress;
  R2SymTensorT < 3 >& epsilon = matState.epsilon;
  R2SymTensorT < 3 >& omega = matState.omega;
  R2SymTensorT < 3 >& epsid = matState.epsid; // to add
  R2SymTensorT < 3 >& epsE = matState.epsE; // to add
  R2SymTensorT < 3 >& epsel = matState.epsel; // to add


  // const realT& G  = matParams.init_shearModulus;
  const realT& E0 = matParams.E;
  const realT& Nu0 = matParams.Nu;
  const realT& a1 = matParams.a1;
  const realT& a2 = matParams.a2;
  const realT& a3 = matParams.a3;
  const realT& a4 = matParams.a4;
  const realT& c0 = matParams.c0;
  const realT& c1 = matParams.c1;
  const realT& alpha = matParams.alpha;
  const R2SymTensorT < 3 >& omega0 = matParams.omega0;
  R1TensorT < 9 > Parameters;
  //Parametrs[] = // = {E0, Nu0, a1, a2, a3, a4, c0, c1, alpha};
  const realT FTOL = 10e-6;
  const int ITER = 500;


  int INC;
  realT lambda, lambda1 ,lambda2;
  realT lambda_best, lambda_best1, lambda_best2;
  realT fd_1, fd_2, fd_3, fd_best1,fd_best2,fd_best;

  R2SymTensorT < 3 > dstress,depsilon;
  R2SymTensorT < 3 > depsid,depsid1,depsid2;
  R2SymTensorT < 3 > depsE,depsel,domega;
  R2SymTensorT < 3 > stressT;
  R2SymTensorT < 3 > stressTT1, stressTT2, stressTT;
  R2SymTensorT < 3 > omegaT, omegaT1, omegaT2;
  R2SymTensorT < 3 > stress = 0.;
  R2SymTensorT < 3 > E = 0.;
  R4minSymTensorT < 3 > De,De0,Se,Se0;

  realT fd_trial;



  E.PlusIdentity(1.);
  stress.PlusIdentity(pressure);
  stress += devStress;
  depsilon = Ddt;

  EffectiveElasticStiffness(Nu0, E0, a1, a2, a3, a4, E, omega, De,
                            De0, Se, Se0);

  // Elastic Trial
  dstress.AijklBkl(De,depsilon);

  stressT = stress;
  stressT += dstress;
  fd_trial = DamageFunction(a1, a2, a3, a4, c0, c1, alpha,
                            E, stressT, omega);


  if (fd_trial>FTOL)
  {
    INC = 1;
    lambda1 = lambda_best = 0;
    lambda2 = 0.001;
    fd_3 = fd_trial;
    fd_best = 1e10;

    while (fabs(fd_3)>FTOL && INC<ITER)
    {
      Others(Nu0, E0, a1, a2, a3, a4, c0, c1, alpha,
             E, epsilon, epsid, stress, omega,lambda1,fd_1,depsilon, depsid1,omegaT1, stressTT1);

      Others(Nu0, E0, a1, a2, a3, a4, c0, c1, alpha,
             E, epsilon, epsid, stress, omega,lambda2,fd_2,depsilon, depsid2,omegaT2, stressTT2);

      if(fabs(fd_2-fd_1)<1e-8)
     {
        lambda = lambda2-fd_2*(lambda2-lambda1)/(1e-8);
     }
     else
     {
       lambda = lambda2-fd_2*(lambda2-lambda1)/(fd_2-fd_1);
     }

      Others(Nu0, E0, a1, a2, a3, a4, c0, c1, alpha,
             E, epsilon, epsid, stress, omega,lambda,fd_3,depsilon, depsid,omegaT, stressTT);


      if (fabs(fd_1)>fabs(fd_2))
      {
        fd_best1 = fd_2;
        lambda_best1 = lambda2;
      }
      else
      {
        fd_best1 = fd_1;
        lambda_best1 = lambda1;
      }

      if (fabs(fd_best1)>fabs(fd_3))
      {
        fd_best2 = fd_3;
        lambda_best2 = lambda;
      }
      else
      {
        fd_best2 = fd_best1;
        lambda_best2 = lambda_best1;
      }
      if (fabs(fd_best)>fabs(fd_best2))
      {
        fd_best = fd_best2;
        lambda_best = lambda_best2;
      }
      INC++;
      lambda1 = lambda2;
      lambda2 = lambda;
    }


    Others(Nu0, E0, a1, a2, a3, a4, c0, c1, alpha,
           E, epsilon, epsid, stress, omega,lambda_best,fd_3,depsilon, depsid,omegaT, stressTT);


    for (int ii=0; ii<3; ii++)
    {
      for (int jj=0; jj<=ii; jj++)
      {
        dstress(ii,jj) = stressTT(ii,jj)-stress(ii,jj);
        domega(ii,jj) = omegaT(ii,jj)-omega(ii,jj);
        omega(ii,jj) += domega(ii,jj);
        depsE(ii,jj) = depsilon(ii,jj)-depsid(ii,jj);
      }
    }

    fd_3 = DamageFunction(a1, a2, a3, a4, c0, c1, alpha,
                                E, stressTT, omega);

  }
  else
  {
    depsid = 0;
    depsE = depsilon;

  }

  depsel.AijklBkl(Se0,dstress);
  stress += dstress;
  epsel += depsel;
  epsE += depsE;
  epsid += depsid;
  fd0 = DamageFunction(a1, a2, a3, a4, c0, c1, alpha, E, stress, omega);
  pressure = stress.Trace();
  pressure /= 3.;
  stress.PlusIdentity(-pressure);
  devStress = stress;
  epsilon+=depsilon;
  matState.RotateState( Rot );
  return;
}

//FUNCTION_BEGIN_PARSE
void
DSIDDamage::EffectiveElasticStiffness(const realT Nu0, const realT E0,
                                      const realT a1, const realT a2,
                                      const realT a3,const realT a4,
                                      const R2SymTensorT < 3 > E,
                                      R2SymTensorT < 3 > omega,
                                      R4minSymTensorT < 3 >& De,
                                      R4minSymTensorT < 3 >& De0,
                                      R4minSymTensorT < 3 >& Se,
                                      R4minSymTensorT < 3 >& Se0)
{
    const realT B1 = (1.0+Nu0)/E0/2.0;
    const realT B2 = Nu0/E0;
    realT tromega;
    tromega = omega.Trace();
    for(int i=0;i<3;i++)
    {
      for(int j=0;j<=i;j++)
      {
        for(int k=0;k<3;k++)
        {
          for(int l=0;l<=k;l++)
          {
            Se(i,j,k,l) =
                B1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))-B2*E(i,j)*E(k,l)+2.0*a1*tromega*E(i,j)*E(k,l)+
                0.5*a2*(E(i,k)*omega(j,l)+E(i,l)*omega(j,k)+omega(i,k)*E(j,l)+omega(i,l)*E(j,k))+
                a3*(E(i,j)*omega(k,l)+omega(i,j)*E(k,l))+a4*tromega*(E(i,k)*E(j,l)+E(i,l)*E(j,k));
            Se0(i,j,k,l) = B1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))-B2*E(i,j)*E(k,l);
          }
        }
      }
    }

    De.Inverse_4(Se);
    De0.Inverse_4(Se0);
}

//FUNCTION_BEGIN_PARSE
void
DSIDDamage::DamageDrivingForce(const realT a1, const realT a2,
                               const realT a3,const realT a4,
                               const R2SymTensorT < 3 > E,
                               R2SymTensorT < 3 > stresstotal,
                               R2SymTensorT < 3 >& Yd)
  {
    R2SymTensorT < 3 > sigsig;
    realT trstresstotal,trsigsig;
    trstresstotal = stresstotal.Trace();
    sigsig.AijAjk(stresstotal);
    trsigsig = sigsig.Trace();
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<=i; j++)
      {
        Yd(i,j) = a1*pow(trstresstotal,2)*E(i,j)+a2*sigsig(i,j)
              +a3*trstresstotal*stresstotal(i,j)+a4*trsigsig*E(i,j);
      }
    }
  }


//FUNCTION_BEGIN_PARSE
void
DSIDDamage::ProjectionTensorP1(R2SymTensorT < 3 > stresstotal, R4minSymTensorT < 3 >& P1)
{
    R1TensorT < 3 > S;
    realT Eigenval[3];
    R1TensorT < 3 > EigenVector[3];
    stresstotal.EigenVals(Eigenval, 0.);
    stresstotal.EigenVecs(Eigenval, EigenVector);
    //for (int ii=0; ii<3; ii++)
    //{
    //  if(abs(Eigenval[ii]<=10e-6))
    //     Eigenval[ii] = 0.;
    //}
    for (int jj=0; jj<3; jj++)
    {
      if (Eigenval[jj]>0)
      {S(jj)=1.;}
      else
      {S(jj)=-1;}
    }
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<=i; j++)
      {
        for (int k=0; k<3; k++)
        {
          for (int l=0; l<=k; l++)
          {
            P1(i,j,k,l) = S(0)*EigenVector[0](i)*EigenVector[0](j)*EigenVector[0](k)*EigenVector[0](l)+
                          S(1)*EigenVector[1](i)*EigenVector[1](j)*EigenVector[1](k)*EigenVector[1](l)+
                          S(2)*EigenVector[2](i)*EigenVector[2](j)*EigenVector[2](k)*EigenVector[2](l);
          }
        }
      }
    }
}


//FUNCTION_BEGIN_PARSE
inline realT
DSIDDamage::DamageFunction(const realT a1, const realT a2,
                           const realT a3, const realT a4,
                           const realT c0, const realT c1, const realT alpha,
                           const R2SymTensorT < 3 > E,
                           R2SymTensorT < 3 > stresstotal,
                           R2SymTensorT < 3 > omega)
{
  R2SymTensorT < 3 > Yd, P1Y, Sij;
  realT trP1Y, SS, tromega, fd;
  R4minSymTensorT < 3 > P1;
  tromega = omega.Trace();
  DamageDrivingForce(a1, a2, a3, a4, E, stresstotal, Yd);
  ProjectionTensorP1(stresstotal, P1);
  P1Y.AijklBkl(P1,Yd);
  trP1Y = P1Y.Trace();
  Sij = P1Y;
  Sij.PlusIdentity( -trP1Y / 3.0 );
  SS = Sij.AijAij();
  fd = sqrt(0.5*SS)+alpha*trP1Y-c0-c1*tromega;  // soil convention -alpha
  return {fd};

}

//FUNCTION_BEGIN_PARSE
void
DSIDDamage::ProjectionTensorP2(R2SymTensorT < 3 > stresstotal, R4minSymTensorT < 3 >& P2)
{
  R1TensorT < 3 > S;
  realT RES;
  realT Eigenval[3];
  R1TensorT < 3 > EigenVector[3];
  stresstotal.EigenVals(Eigenval, 0.);
  stresstotal.EigenVecs(Eigenval, EigenVector);
  //      for (int ii=0; i<3; i++)
  //      {
  //        if(abs(Eigenval[i]<=10e-6))
  //           Eigenval[i] = 0.;
  //      }
  for (int jj=0; jj<3; jj++)
  {
    RES = Eigenval[jj]-Eigenval[2];
    if (RES>0)
    {S(jj)=1.;}
    else
    {S(jj)=0;}
  }
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<=i; j++)
    {
      for (int k=0; k<3; k++)
      {
        for (int l=0; l<=k; l++)
        {
          P2(i,j,k,l) = S(0)*EigenVector[0](i)*EigenVector[0](j)*EigenVector[0](k)*EigenVector[0](l)+
                       S(1)*EigenVector[1](i)*EigenVector[1](j)*EigenVector[1](k)*EigenVector[1](l)+
                       S(2)*EigenVector[2](i)*EigenVector[2](j)*EigenVector[2](k)*EigenVector[2](l);
        }
      }
    }
  }
 }

//FUNCTION_BEGIN_PARSE
void
DSIDDamage::Others(const realT Nu0, const realT E0, const realT a1, const realT a2, const realT a3, const realT a4,
                   const realT c0, const realT c1, const realT alpha, const R2SymTensorT < 3 > E,
                   R2SymTensorT < 3 > epsilon, R2SymTensorT < 3 > epsid,
                   R2SymTensorT < 3 > stress, R2SymTensorT < 3 > omega, realT& lambda, realT& fd,R2SymTensorT < 3 > depsilon,
                   R2SymTensorT < 3 >& depsid, R2SymTensorT < 3 >& omegaT, R2SymTensorT < 3 >& stressTT)
{
  R2SymTensorT < 3 > dF_domega,dF_dstress,dF_dY,dg_dY,domegaT,EP1,epsT,F1ij,F1P3,F2ij,F2P2,P1Yd,temp,Yd;
  realT F1F1,F2F2,P1YdE, trstress;
  R4minSymTensorT < 3 > dDdO_dgdY, De, De0, DT, dY_dstress, EEP1,P1,P2,P3,Se, Se0, SOMSOM,SSINV;
  R6minSymTensorT < 3 > dD_dO,dS_dO;
  trstress = stress.Trace();
  for(int i=0; i<3; i++)
    {
      for(int j=0; j<=i; j++)
      {
        for(int k=0; k<3; k++)
        {
          for(int l=0; l<=k; l++)
          {
            dY_dstress(i,j,k,l)=2*a1*trstress*E(i,j)*E(k,l)+0.5*a2*
                (E(i,k)*stress(l,j)+E(i,l)*stress(j,k)+E(j,l)*stress(i,k)
                +E(j,k)*stress(i,l))+a3*(E(k,l)*stress(i,j)+0.5*trstress*(E(i,k)*E(j,l)
                    +E(i,l)*E(j,k)))+2*a4*stress(k,l)*E(i,j);
          }
        }
      }
    }
  DamageDrivingForce(a1, a2, a3, a4, E, stress, Yd);
  ProjectionTensorP1(stress, P1);
  P1Yd.AijklBkl(P1,Yd);
  P1YdE = P1Yd.AijBij(P1Yd,E);
  ProjectionTensorP2(stress, P2);
  F2ij.AijklBkl(P2,Yd);
  F2P2.AijklBij(P2,F2ij);
  F2F2 = F2ij.AijAij();
  EP1.AijklBij(P1,E);
  EEP1.AijBkl(E,EP1);
  for(int i=0; i<3; i++)
  {
    for(int j=0; j<=i; j++)
    {
      for(int k=0; k<3; k++)
      {
        for(int l=0; l<=k; l++)
        {
          P3(i,j,k,l) = P1(i,j,k,l)-1./3.*EEP1(i,j,k,l);
        }
      }
    }
  }

  temp = P1Yd;
  temp.PlusIdentity(-1./3.*P1YdE);
  F1ij = temp;

  for(int i=0; i<3; i++)
  {
    for (int j=0; j<=i; j++)
    {
      dF_domega(i,j) = -c1*E(i,j);
      if(F2F2 == 0)
      {dg_dY(i,j) = 0;}
      else
      {dg_dY(i,j)=F2P2(i,j)/sqrt(2.*F2F2);}
     }
   }

   F1P3.AijklBij(P3,F1ij);
   F1F1 = F1ij.AijAij();

   for(int i=0; i<3; i++)
   {
     for (int j=0; j<=i; j++)
     {
       if(F1F1==0)
       {
         dF_dY(i,j)=0;
       }
       else
       {
         dF_dY(i,j)=F1P3(i,j)/sqrt(2.*F1F1)+alpha*EP1(i,j);  // soil convention uses -alpha,mechanical convention uses +alpha
       }
     }
   }

   dF_dstress.AijklBij(dY_dstress,dF_dY);

   for(int i=0; i<3; i++)
   {
     for (int j=0; j<=i; j++)
     {
       depsid(i,j) = lambda*dF_dstress(i,j);
       domegaT(i,j) = lambda*dg_dY(i,j);
       omegaT(i,j) = omega(i,j)+domegaT(i,j);
     }
   }


   for(int i=0; i<3; i++)
   {
     for(int j=0; j<=i; j++)
     {
       for(int k=0; k<3; k++)
       {
         for(int l=0; l<=k; l++)
         {
           for(int m=0; m<3; m++)
           {
             for(int n=0; n<=m; n++)
             {
           dS_dO(i,j,k,l,m,n)=2.*a1*E(i,j)*E(k,l)*E(m,n)+0.25*a2*
    (E(i,k)*(E(m,j)*E(n,l)+E(m,l)*E(n,j))+E(i,l)*(E(m,j)*E(n,k)+
    E(m,k)*E(n,j))+E(j,l)*(E(i,m)*E(n,k)+E(i,n)*E(m,k))+E(j,k)*
    (E(i,m)*E(n,l)+E(i,n)*E(m,l)))+0.5*a3*(E(i,j)*(E(k,m)*E(l,n)+
    E(k,n)*E(l,m))+E(k,l)*(E(i,m)*E(j,n)+E(i,n)*E(j,m)))+a4*
    (E(i,k)*E(j,l)+E(i,l)*E(j,k))*E(m,n);
             }
           }
         }
       }
     }
   }

   EffectiveElasticStiffness(Nu0, E0, a1, a2, a3, a4, E, omega, De,
                               De0, Se, Se0);
   SOMSOM.AijmnBmnkl(Se,Se);
   SSINV.Inverse_4(SOMSOM);


   for(int i=0; i<3; i++)
   {
     for(int j=0; j<=i; j++)
     {
       for(int k=0; k<3; k++)
       {
         for(int l=0; l<=k; l++)
         {
           for(int m=0; m<3; m++)
           {
             for(int n=0; n<=m; n++)
             {
               dD_dO(i,j,k,l,m,n)=0.;
               for(int im=0; im<3; im++)
               {
                 for(int jm=0; jm<=im; jm++)
                 {
                   dD_dO(i,j,k,l,m,n)-= SSINV(i,j,im,jm)*dS_dO(im,jm,k,l,m,n);
                 }
               }
             }
           }
         }
       }
     }
   }

   dDdO_dgdY.AijklmnBmn(dD_dO,dg_dY);


   for(int i=0; i<3; i++)
   {
     for(int j=0; j<=i; j++)
     {
       epsT(i,j) = epsilon(i,j)+depsilon(i,j)-epsid(i,j)-depsid(i,j);
       for(int k=0; k<3; k++)
       {
         for(int l=0; l<=k; l++)
         {
           DT(i,j,k,l) = De(i,j,k,l)+lambda*dDdO_dgdY(i,j,k,l);
         }
       }
     }
   }

   stressTT.AijklBkl(DT,epsT);
   fd = DamageFunction(a1, a2, a3, a4, c0, c1, alpha,
                  E, stressTT, omegaT);

}

