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
 * @file OneDReactionFront.cpp
 * @author walsh24
 * @date July 27, 2011
 **/
 
#include <cmath>
#include <vector>

#include "PhysicsSolvers/SolverFactory.h"
#include "OneDReactionFront.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"
#include "ObjectManagers/ChemistryManager.h"
#include "ObjectManagers/FunctionManager.h"
#include "ObjectManagers/UnitManager.h"
#include "ObjectManagers/LogManager.h"

//#include "Utilities/ReducedRowEchelonForm.h"
#include "Utilities/Functions.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"

using namespace GPChemistry;
using namespace ODRF;
using namespace PS_STR;

namespace{
  
  template <class T>
  void GetVectorAttribute(std::vector<T>& vals,TICPP::HierarchicalDataNode* hdn, std::string attribute,std::string sep = ",", bool doTrim = true){
    std::string vStr = hdn->GetAttributeString(attribute);
    if( !vStr.empty() ){
      sArray1d vStrs = Tokenize(vStr,sep);
      vals.resize(vStrs.size());
      for(unsigned i =0; i < vStrs.size(); ++i){
        if(doTrim) Trim(vStrs[i]," \t\n");
        vals[i] = fromString<T>(vStrs[i]);
      }
    }
  }

  inline realT dot(const rArray1d& x,const rArray1d& y ){
    if(x.size() != y.size() ) throw GPException( "Dot product - vector sizes must agree." );
    realT rv(0.0);
    for(rArray1d::size_type i = 0; i < x.size(); ++i) rv += x[i]*y[i];
    return rv;
  }
  
}
  
/*
  F(i=0:n-1) =  X[i] - exp(X[i+n]-logGamma_i)
  F(i=n:n+nA-1)  = A*X(0:n-1) - b
  F(i=n+nA:2*n-1) = K*X(n+nA:2*n-1) - p
  F(2*n)      = X[2*n]-IAPCoeff*X[0:n-1] 
*/
inline
void ODRF::EquilibriumSystem::GetFunction(rArray1d& X,rArray1d& F){
    unsigned nA = m_b.size();
    unsigned nK = m_p.size();
    unsigned n = m_K.Dimension(1); // number of species
    realT ff;


    // sanity check - adjust X when log molality or residual is too high
    for(unsigned i = 0; i < n; ++i){
        bool adjustX = log( std::max(X[i],1e-64) ) > 5  || X[i+n] > 5 ;
    	if(adjustX){
            realT logGamma = m_activityFunction.LogIAC(X, i);
    		X[i] *=0.5;
    		X[i+n] =   log(  std::max(X[i],1e-64) ) + logGamma;
    	}
    }

    F.resize(2*n+1);
    for(unsigned i = 0; i < n; ++i){
      realT logGamma = m_activityFunction.LogIAC(X, i);
      F[i] = X[i] - exp(X[i+n]-logGamma);
    }
    // switching K,A row order
    for(unsigned i(0); i < nK; ++i){
      ff = -m_p[i];
      for(unsigned j = 0; j < n; ++j) ff += m_K(i,j)*X[j+n];
      F[i+n] = ff;
    }
    for(unsigned i(0); i < nA; ++i){
      ff = -m_b[i];
      for(unsigned j = 0; j < n; ++j) ff += m_A(i,j)*X[j];
      F[i+n+nK] = ff;
    }

    ff = X[2*n];
    for(unsigned j = 0; j < n; ++j) ff -= m_ionicStrengthCoefficients(j)*X[j];
    F[2*n] = ff;


    if(m_doScaling){
      for(unsigned i =0; i < m_isolatedDofs.size(); ++i){
        unsigned dof = m_isolatedDofs[i];
        F[dof] = 0;
      }
    }
   
}
/*
  GradF = [I   -Diag(exp(X[i+n])-log(gamma_i))  exp( X[i+n]- log(gamma_i)) d(log gamma_i)/du 
           0         K                     ged    0
           A         0                         0
        IAPcoeff     0                         1    ];
*/
void ODRF::EquilibriumSystem::BuildJacobian(){
    unsigned nA = m_b.size();
    unsigned nK = m_p.size();
    unsigned n = m_K.Dimension(1); // number of species

    m_gradF = rArray2d(2*n+1,2*n+1);
    for(unsigned i = 0; i < n; ++i){
      m_gradF(i,i) = 1;
     // m_gradF(i,i+n) = -exp(X[i+n]-log gamma_i); // Added when X is supplied - ie when  jacobian requested.
     // m_gradF(i,2*n) = exp(X[i+n]-log gamma_i)*d(log gamma_i)/du ; // Added when X is supplied - ie when jacobian requested.
    }
    /*
    for(int i(0),ii(n); i < nA; ++i,++ii){
      for(int j = 0; j < n; ++j){
        m_gradF(ii,j) = m_A(i,j);        
      }
    }
    for(int i(0),ii(n+nA); i < nK; ++i,++ii){
      for(int j = 0; j < n; ++j){
        m_gradF(ii,j+n) += m_K(i,j);
      }
    }
    */
    // switching K and A row order
    for(unsigned i(0),ii(n); i < nK; ++i,++ii){
      for(unsigned j = 0; j < n; ++j){
        m_gradF(ii,j+n) += m_K(i,j);
      }
    }
    for(unsigned i(0),ii(n+nK); i < nA; ++i,++ii){
      for(unsigned j = 0; j < n; ++j){
        m_gradF(ii,j) = m_A(i,j);        
      }
    }

    m_gradF(2*n,2*n) = 1.0;
    for(unsigned j = 0; j < n; ++j){
        m_gradF(2*n,j) = -m_ionicStrengthCoefficients(j);        
    }

    if(m_doScaling) CheckActivities();
    
/*
    std::cout << "Equilibrium system" << std::endl; 
    for(int i = 0; i < 2*n+1; ++i){
    for(int j = 0; j < 2*n+1; ++j){
      std::cout << m_gradF(i,j) << " "; 
    }
      std::cout << std::endl; 
    }  
    std::cout << std::endl; 

    std::cout << "b" << std::endl;  
    for(int i = 0; i < nA; ++i){
      std::cout << m_b[i] << " "; 
    }
    std::cout << std::endl; 
*/
    
}

rArray2d& ODRF::EquilibriumSystem::GetJacobian(const rArray1d& X){
//    const unsigned nA = m_b.size();
//    const unsigned nK = m_p.size();
    const unsigned n = m_K.Dimension(1); // number of species

    for(unsigned i = 0; i < n; ++i){
      realT logGamma = m_activityFunction.LogIAC(X, i);
      realT dLogGammadU = m_activityFunction.dLogIACdmu(X, i);
      realT c = exp(X[i+n]-logGamma);
      m_gradF(i,i+n) = -c;
      m_gradF(i,2*n) = c*dLogGammadU;
    }

    if(m_doScaling){
      for(unsigned i =0; i < m_isolatedDofs.size(); ++i){
        unsigned dof = m_isolatedDofs[i];
        m_gradF(dof,dof+n) = -1;
        m_gradF(dof,2*n) = 0;  
      }
    }

    return m_gradF;
}

 
void ODRF::EquilibriumSystem::CheckActivities(void){

  const unsigned n = m_K.Dimension(1); // number of species

  m_isolatedDofs =  rArray1d();

  //column scale
  for(unsigned j = n; j < 2*n; ++j){
    bool isIsolatedDof = true;
    for(unsigned i = n; i < 2*n  && isIsolatedDof; ++i){
      isIsolatedDof =  ( fabs( m_gradF[i][j] ) == 0.0 );
    }
    
    if(isIsolatedDof) m_isolatedDofs.push_back(j-n);
  }
  
}

//////////////////////////////////////////////////////////////////////////

// Reaction front system
////////////////////////
  
inline
void ODRF::ReactionFrontSystem::GetFunction(rArray1d& X,rArray1d& F){
    unsigned nA = m_b.size();
    unsigned nK = m_p.size();
    unsigned nC = m_K.Dimension(1); // number of species per front
    unsigned n = nC*m_numFronts; 

    
    realT FF = 0.0;

    // sanity check - adjust X when log molality or residual is too high
    for(unsigned f = 0; f < m_numFronts; ++f){
      for(unsigned i(0), ii(f*nC); i < nC; ++i,++ii){

        if(X[ii] <= 0){
    	   const realT* XX(&(X[f*nC]));
           realT logGamma = m_activityFunction.LogIAC(XX, i,nC);
           if(X[ii+n] <= 1){
             X[ii] = exp(-logGamma + X[ii+n] );
           } else {
             X[ii] = 1e-64;
           } 
            //std::cout << "negative X" << std::endl;
        } 

        bool adjustX = log( std::max(X[ii],1e-64) ) > 1  || X[ii+n] > 1 ;

    	if(adjustX){
    	    const realT* XX(&(X[f*nC]));
            realT logGamma = m_activityFunction.LogIAC(XX, i,nC);
            while(log(X[ii])> 1) X[ii] *=0.5;
            X[ii] = std::max(X[ii],1e-64); 
            X[ii+n] =   log( X[ii] ) + logGamma;
            //std::cout << "adjustX" << std::endl;
    	}
      }
    }

    F.resize( (2*nC+1)*m_numFronts);
    // activity/molality
    for(unsigned f = 0; f < m_numFronts; ++f){
      const realT* XX(&(X[f*nC]));
      for(unsigned i(0), ii(f*nC); i < nC; ++i,++ii){
        realT logGamma = m_activityFunction.LogIAC(XX, i,nC);
        F[ii] = X[ii] -exp(X[ii+n] - logGamma);
        //if(fabs(F[ii]) > 1e10){
        //  std::cout << "F[ii] " << F[ii] << std::endl;
        //}
      }
    }
    // flow conditions
    for(unsigned i(0); i < nA; ++i){
      FF = -m_b[i];
      for(unsigned j = 0; j < n; ++j) FF += m_A(i,j)*X[j];
      F[i+n] = FF;
    }
    
    // equilibrium
    for(unsigned f = 0; f < m_numFronts; ++f){
      unsigned roffset = f*(nK+1)+n+nA;
      unsigned cOffset = f*nC+n;   
      // solid phase equilibrium
      FF = -m_pf[f];  
      for(unsigned j = 0; j < nC; ++j)  FF += m_Kf(f,j)*X[j+cOffset];
      F[roffset] = FF;
      // aqueous equilibrium
      for(unsigned i(0); i < nK; ++i){
        FF = -m_p[i];
        for(unsigned j = 0; j < nC; ++j) FF += m_K(i,j)*X[j+cOffset];
        F[i+roffset+1] = FF;
      }
    }

    // ionic strength coefficients
    for(unsigned f = 0; f < m_numFronts; ++f){
      realT ff = X[2*n+f];
      unsigned cOffset = f*nC;
      for(unsigned j = 0; j < nC; ++j) ff -= m_ionicStrengthCoefficients(j)*X[j+cOffset];
      F[2*n+f] = ff;
    }

    // 
    if(m_doScaling){
      for(unsigned i =0; i < m_isolatedDofs.size(); ++i){
        unsigned dof = m_isolatedDofs[i];
        F[dof] = 0;
      }
    }
   
}

/*
  GradF = [I   -Diag(exp(X[i+n]))
           A         0                 0
           0         KK                0
          -ISC        0                 1]
  where A contains information concerning flux, and charge equilibrium conditions at each front
  and KK is the matrix holding the natural log of the equilibrium equations
  KK = [ K 0 0 ... 0
         0 K 0 ... 0
              ...
         0 ...  0  K]
  (one K block per front)
*/
void ODRF::ReactionFrontSystem::BuildJacobian(){
    unsigned nA = m_b.size(); // number of A rows
    unsigned nK = m_p.size(); // number of K rows
    unsigned nC = m_K.Dimension(1); // number of species per front
    unsigned n = nC*m_numFronts; 

    m_gradF = rArray2d( (2*nC+1)*m_numFronts,(2*nC+1)*m_numFronts);
    for(unsigned i = 0; i < n; ++i){
      m_gradF(i,i) = 1;
      // m_gradF(i,i+n) = -exp(X[i+n]-log gamma_i); // Added when X is supplied - ie when  jacobian requested.
       // m_gradF(i,2*n) = exp(X[i+n]-log gamma_i)*d(log gamma_i)/du ; // Added when X is supplied - ie when jacobian requested.
    }
    
    for(unsigned i(0),ii(n); i < nA; ++i,++ii){
      for(unsigned j = 0; j < n; ++j){
        m_gradF(ii,j) = m_A(i,j);        
      }
    }
    
    for(unsigned f = 0; f < m_numFronts; ++f){
      unsigned roffset = f*(nK+1)+n+nA;
      unsigned cOffset = f*nC+n; 
      // solid phase equilibrium  
      for(unsigned j = 0; j < nC; ++j)  m_gradF(roffset,j+cOffset) += m_Kf(f,j);
      // aqueous equilibrium
      for(unsigned i(0),ii(roffset+1); i < nK; ++i,++ii){
        for(unsigned j = 0; j < nC; ++j){
          m_gradF(ii,j+cOffset) += m_K(i,j);
        }
      }  
    }

    // ionic strength coefficients
    for(unsigned f = 0; f < m_numFronts; ++f){
      unsigned row = 2*n+f;
      unsigned cOffset = f*nC;
      m_gradF(row,row) = 1.0;
      for(unsigned j = 0; j < nC; ++j){
        m_gradF(row,j+cOffset) = -m_ionicStrengthCoefficients(j);
      }
    }

    /*
    std::cout << "Reaction fronts" << std::endl; 
    for(unsigned i = 0; i < 2*n+m_numFronts; ++i){
    for(unsigned j = 0; j < 2*n+m_numFronts; ++j){
      std::cout << m_gradF(i,j) << " ";
    }
      std::cout << std::endl; 
    }
    exit(0);
    */

    
    if(m_doScaling) CheckActivities();


}

void ODRF::ReactionFrontSystem::UpdateJacobian(){
    unsigned nA = m_b.size();

    unsigned nC = m_K.Dimension(1); // number of species per front
    unsigned n = nC*m_numFronts; 

    for(unsigned i(0),ii(n); i < nA; ++i,++ii){
      for(unsigned j = 0; j < n; ++j){
        m_gradF(ii,j) = m_A(i,j);        
      }
    }


    /*
    std::cout << "Reaction front Jacobian- after update" << std::endl;
    for(unsigned i = 0; i < 2*n+m_numFronts; ++i){
    for(unsigned j = 0; j < 2*n+m_numFronts; ++j){
      std::cout << m_gradF(i,j) << " ";
    }
      std::cout << std::endl;
    }
    exit(0);
    */
}

rArray2d& ODRF::ReactionFrontSystem::GetJacobian(const rArray1d& X){

    unsigned nC = m_K.Dimension(1); // number of species per front
    unsigned n = nC*m_numFronts;    // offset for activities (Y)

    for(unsigned f = 0; f < m_numFronts; ++f){
      const realT* XX( &(X[f*nC]) );
      for(unsigned i(0), ii(f*nC); i < nC; ++i,++ii){
        realT logGamma = m_activityFunction.LogIAC(XX, i,nC);
        realT dLogGammadU = m_activityFunction.dLogIACdmu(XX, i,nC);
        //realT c = exp(XX[i+n]-logGamma);
        realT c = std::max(exp(XX[i+n]-logGamma), 1e-64 );  // can get into trouble if c too small
        m_gradF(ii,ii+n) = -c;
        m_gradF(ii,2*n+f) = c*dLogGammadU;
      }
    }

    /*
    std::cout << "Reaction front Jacobian- with IAC" << std::endl;
    for(unsigned i = 0; i < 2*n+m_numFronts; ++i){
    for(unsigned j = 0; j < 2*n+m_numFronts; ++j){
      std::cout << m_gradF(i,j) << " ";
    }
      std::cout << std::endl;
    }
    exit(0);*/


    if(m_doScaling){
      for(unsigned i =0; i < m_isolatedDofs.size(); ++i){
        unsigned dof = m_isolatedDofs[i];
        unsigned f = dof/nC;
        m_gradF(dof,dof+n) = -1;
        m_gradF(dof,2*n+f) = 0;  
      }
    }

    return m_gradF;
}

void ODRF::ReactionFrontSystem::RescaleSolution(rArray1d& dX, rArray1d& X,rArray1d& F){

    unsigned nC = m_K.Dimension(1); // number of species per front
    //unsigned n = nC*m_numFronts;


    // Under Relaxation

    realT irf =1.0;  // inverse relaxation factor

    for(unsigned f = 0; f < m_numFronts; ++f){
      for(unsigned i(0), ii(f*nC); i < nC; ++i,++ii){
        realT v = 2*dX[ii]/X[ii];  // step size at most X/2
        if(v > irf) irf = v;
        if(X[ii]< 0){
        	std::cout <<"Warning: X<0 " << X[ii] <<std::endl;
        }
      }
    }

    // under relaxation ensures positive concentrations
    realT urf = 1.0/irf;
    dX *= urf;
    
}
  
void ODRF::ReactionFrontSystem::CheckActivities(void){

  unsigned nC = m_K.Dimension(1); // number of species per front
  unsigned n = nC*m_numFronts;  // Y offset

  m_isolatedDofs =  rArray1d();

  //column scale
  for(unsigned j = n; j < 2*n; ++j){
    bool isIsolatedDof = true;
    for(unsigned i = n; i < 2*n  && isIsolatedDof; ++i){
      isIsolatedDof =  ( fabs( m_gradF[i][j] ) == 0.0 );
    }
    
    if(isIsolatedDof) m_isolatedDofs.push_back(j-n);
  }
  
}
//////////////////////////////////////////////////////////////////

  
void NewtonRaphson( NonLinearSystemBase& nls, 
                    rArray1d& x, unsigned maxNumItrs, realT tolx, bool& rootFound, unsigned& numItrs, realT& f){

    localIndex n = x.size();
	const realT TINY = 1e-64;
   
  
    rArray1d& F = nls.m_F;
    rArray1d& dx = nls.m_dx;
    rArray1d& g = nls.m_g;
    F.resize(n);
    dx.resize(n);
    g.resize(n);
    //rArray1d F(n),dx(n),g(n);

    // return values
    numItrs = 0;
    f =0.0;
    rootFound = false;

    for(;numItrs < maxNumItrs && !rootFound; ++numItrs){

      nls.GetFunction(x,F);
      f = dot(F,F); 
      realT fold = f; rArray1d xold = x;
      rArray2d& gradF = nls.GetJacobian(x);

      // dx
      int errFlag;
      bool doEquilibrateByHand = false;
      if(doEquilibrateByHand){
        // Having trouble with Teuchos equilibration - the following has been added to see if equilibration will help
        rArray1d F_scaled(n);
        rArray2d gradF_scaled(n,n);
        rArray1d Rscale(n);
        rArray1d Cscale(n);

        // row scale
        for(unsigned i =0; i < n; ++i){
          for(unsigned j =0; j < n; ++j){
            realT v = fabs(gradF[i][j]);
            if( v > Rscale[i] ) Rscale[i] = v;
          }
        }

        for(unsigned i =0; i < n; ++i){
          for(unsigned j =0; j < n; ++j){
            gradF_scaled[i][j] = gradF[i][j]/Rscale[i];
          }
          F_scaled[i] = F[i]/Rscale[i];
        }

        //column scale
        for(unsigned i =0; i < n; ++i){
          for(unsigned j =0; j < n; ++j){
            realT v = fabs(gradF_scaled[i][j]);
            if( v > Cscale[j] ) Cscale[j] = v;
          }
        }

        for(unsigned i =0; i < n; ++i){
          for(unsigned j =0; j < n; ++j){
            if(Cscale[j] > 0.0) gradF_scaled[i][j] /= Cscale[j];
          }
        }
  //debug
/*
    if(n > 50){  
        std::cout << "Grad F" << std::endl;
        for(unsigned i = 0; i < n; ++i){
        for(unsigned j = 0; j < n; ++j){
          std::cout << gradF(i,j) << " ";
        }
          std::cout << std::endl;
        }

    std::cout << "Scaled grad F" << std::endl;
    for(unsigned i = 0; i < n; ++i){
    for(unsigned j = 0; j < n; ++j){
      std::cout << gradF_scaled(i,j) << " ";
    }
      std::cout << std::endl;
    }
        std::cout << "row scale" << std::endl;
        for(unsigned i =0; i < n; ++i){
            std::cout << Rscale(i) << " ";
        }
        std::cout << std::endl;

        std::cout << "col scale" << std::endl;
        for(unsigned i =0; i < n; ++i){
            std::cout << Cscale(i) << " ";
        }
        std::cout << std::endl;
        exit(0); 


    }
     */

        errFlag = LinSolve_Local(gradF_scaled,dx,F_scaled);

        // rescale result
        for(unsigned ii =0; ii < n; ++ii){
          if(Cscale[ii] > 0.0) dx[ii] /=  Cscale[ii];

          // debug
          /*
          if(false && ii ==73 && fabs(dx[ii]) > 100) {

            // debug
            std::cout << "Scaled grad F" << std::endl;
            for(unsigned i = 0; i < n; ++i){
            for(unsigned j = 0; j < n; ++j){
              std::cout << gradF_scaled(i,j) << " ";
            }
              std::cout << std::endl;
            }

            
            std::cout << "F_scaled" << std::endl;
            for(unsigned i = 0; i < n; ++i){
              std::cout << F_scaled[i] << " ";
            }
            std::cout << std::endl;
            exit(0);

          }*/

        }

      }else {
        errFlag = LinSolve_Local(gradF,dx,F);
      }

      if (errFlag > 0){
        std::cout << "Warning: LinSolve_Local failed to converge." << std::endl;
        break;
      }

      // Account for under relaxation and/or Jacobian scaling
      nls.RescaleSolution(dx,x,F);

      dx *= -1.0;

      // check size dx
      // max step size
      realT stpmx = 10.0 * std::max( sqrt(dot(xold,xold)),double(n));
  
      // step length
      realT sl = sqrt( dot(dx,dx) );
      if(sl > stpmx)  dx *= stpmx/sl;
    
      realT tst = TINY*tolx;
      for(localIndex i =0; i< n; ++i){
        realT t = fabs(dx[i])/ std::max(fabs(xold[i]),1.0e-7);
        if(t > tst) tst = t;
      }
      realT lmin = tolx / tst;
      //std::cout << lmin << " " << tolx << " "<< tst << " " <<std::endl;

      // slope 
      //  g.AijBi(gradF,F);
      for(localIndex j =0; j < n; ++j){
        g[j] =0.0;
        for(localIndex i =0; i < n; ++i) g[j] += F[i]*gradF(i,j);
      }
      
      realT slope = dot(g,dx); 
  
      realT l = 1.0;    
      bool minFound = false;
      while(!minFound){

        // x = xold +l*dx;
        x = dx; x *= l; x += xold; 

        nls.GetFunction(x,F);

        f = dot(F,F); 

        if(l < lmin){
          x = xold;
          minFound = true;
          rootFound = true;
        } else if(f<fold+1e-4*l*slope){
          minFound = true;
        } else {
          l *= 0.5;
          if(l < 1e-64) minFound = true;   
        }
      }

      if (rootFound) break;
    }
    
    if(false){
      std::cout << "rootFound " << rootFound << std::endl;
      std::cout << "numItrs " << numItrs << std::endl;
      std::cout << "f " << f << std::endl;
    } else {
      if(f > 1.0 || numItrs > 15 || !rootFound ){
        std::cout << "rootFound " << rootFound << "\n";
        std::cout << "numItrs " << numItrs << "\n";
        std::cout << "f " << f << "\n";
      }

      if ( false && f > 1000 ){

        rArray2d& gradF = nls.GetJacobian(x);

// debug
            std::cout << "grad F" << std::endl;
            for(unsigned i = 0; i < n; ++i){
            for(unsigned j = 0; j < n; ++j){
              std::cout << gradF(i,j) << " ";
            }
              std::cout << std::endl;
            }

    std::cout << "x" << std::endl;
    for(localIndex i =0; i < n; ++i){
      std::cout << x[i] << " ";
    }
        std::cout << std::endl;

    std::cout << "dx" << std::endl;
    for(localIndex i =0; i < n; ++i){
      std::cout << dx[i] << " ";
    }
        std::cout << std::endl;

    std::cout << std::endl;
    std::cout << "F" << std::endl;
    for(localIndex i =0; i < n; ++i){
        std::cout << F[i] << " ";
    }
    exit(0);

      }
    }

    /*
    std::cout << "x" << std::endl;
    for(localIndex i =0; i < n; ++i){
      std::cout << x[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "F" << std::endl;
    for(localIndex i =0; i < n; ++i){
        std::cout << F[i] << " ";
    }
    std::cout << std::endl;
    */
}
  




//////////////////////////////////////////////////////////////////////////////////////////

// Reaction front solver

ReactionFrontSolver::ReactionFrontSolver(  const std::string& name,
                                           ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

ReactionFrontSolver::~ReactionFrontSolver()
{
  // TODO Auto-generated destructor stub
}

/**
 * <ReactionFrontSolver name="reactionFrontUpdate">    
 *        <SolidFronts
 *            names = "CaCO3_1, CaCO3_2, Portlandite"
 *            equilibrium_constants = "K_CaCO3, K_CaCO3, K_ch"
 *            initial_positions = "0.001mm, 0.002mm, 0.003mm"
 *            front_conditions = "q_Ca+2 + q_CaHCO3+ + q_CaCO3 = $:QK*q_CaHCO3+ + $:QK*q_CaCO3 + $:QK*q_CO3-2 + $:QK*q_HCO3- + $:QK*q_H2CO3 + $:QK*q_CO2;
 *                                q_Ca+2 + q_CaHCO3+ + q_CaCO3 = q_CaHCO3+ + q_CaCO3 + q_CO3-2 + q_HCO3- + q_H2CO3 + q_CO2;
 *                                q_CaHCO3+ + q_CaCO3 + q_CO3-2 + q_HCO3- + q_H2CO3 + q_CO2 = 0"
 *            dissolution_rates = "q_Ca+2 + q_CaHCO3+ + q_CaCO3;
 *                                -q_Ca+2 - q_CaHCO3+ - q_CaCO3;
 *                                 q_Ca+2 + q_CaHCO3+ + q_CaCO3"/>
 *        <Regions  
 *            names = "SiO2, CaCO3, chDepleated, Portlandite" 
 *            solid_phases = "SiO2; CaCO3, csh; csh; Portlandite, csh"/>
 *        <SolidPhases
 *            names = "SiO2, CaCO3, csh, Portlandite" 
 *            volume_fractions = "$:vf_SiO2, $:vf_CaCO3, $:vf_csh, $:vf_ch"  
 *            molar_masses="$:MM_SiO2, $:MM_CaCO3, $:MM_csh, $:MM_ch "  
 *            densities =  "$:rho_CaCO3, $:rho_SiO2, $:rho_csh, $:rho_ch" />
 *        <Species names = "H+, OH-, Ca+2, CaHCO3+, CaCO3, CO3-2, HCO3-, H2CO3, CO2" 
 *                 equilibrium_constants = "Kw, K2, K3, K4, K5, K6"               
 *                 diffusivity = "1e-10m^2/s" />
 *        <Elements names = "Ca, C" 
 *                  species = "Ca+2, CaHCO3+, CaCO3;  
 *                             CaHCO3+, CaCO3, CO3-2, HCO3-, H2CO3, CO2" 
 *                  brine_concentrations = "1e-4 mol/l,1e-4 mol/l"/>
 *        <EquilibriumData 
 *                 equilibrium_equations= 
 *                           "Kw = [H+][OH-]; 
 *                            K2 = [H+][CO3-2]/[HCO3-] ; 
 *                            K3 = [Ca+2][HCO3-]/[CaHCO3+] ; 
 *                            K4 = [Ca+2][CO3-2]/[CaCO3]; 
 *                            K5 = [H+][HCO3-]/[H2CO3] ; 
 *                            K6 = [H+][HCO3-]/[CO2][H2O]; 
 *                            K_ch = [Ca+2]/[H+]^2 ;
 *                            K_CaCO3 = [HCO3-][Ca+2]/[H+] ; 
 *                            K_csh = [Ca+2]^1.666666666/[H+]^3.33333333" 
 *                 logKdata= "Kw, 22.801,      -0.010365,    -4787.3,   -7.1321,    0.0; 
 *                            K2, -107.8871,  -0.03252849,   5151.79,   38.92561, -563713.9;  
 *                            K3, -1209.120,     -0.31294,      34765.05,  478.782,    0.0 ;   
 *                            K4, 1228.732,      0.299444,    -35512.75,  485.818,    0.0  ; 
 *                            K5, 1.000393128530148,  0.0,        0.0,      0.0,   0.0   ; 
 *                            K6, -356.3094, -0.060919964,  21834.37,  126.8339, -1684915; 
 *                            K_ch,-8.3848e+001, -1.8373e-002, 9.3154e+003, 3.2584e+001, -1.4538e+002; 
 *                            K_CaCO3,-1.4978e+002, -4.8370e-002, 4.8974e+003, 6.0458e+001, 7.6464e+001; 
 *                            K_csh,   272.602665,  0.058756104, -903.9172332, -104.1248858, 0.0 " />
 * 
 **/
void ReactionFrontSolver::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{  

  SolverBase::ReadXML(hdn);

  int rank(0);
  #if GPAC_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  using namespace GPUnits;
  
  
  // Temperature
  ///////////////
  m_Tc = hdn->GetAttributeOrDefault("temperature", "25");
  
  // Species
  //////////////
  if (rank == 0)  std::cout << "Species" << "\n";
  HierarchicalDataNode* speciesNode = hdn->GetChild("Species");
  if(!speciesNode)  
      throw GPException("ReactionFrontSolver: Must have Species defined in the input file");
  {
    GetVectorAttribute(m_species.names,speciesNode, "names");
    m_species.size = m_species.names.size();
    if (rank == 0)  std::cout << "    ";
    for(unsigned i =0; i < m_species.size; ++i){
      //if (rank == 0)  std::cout << m_species.names[i] << " ";
      m_species.indices[m_species.names[i]] = i;
    }
    if (rank == 0)  std::cout << "\n";
    
    // valences
    /*
    GetVectorAttribute(m_species.valences,speciesNode, "valences");
    if(m_species.valences.size() != m_species.size)
      throw GPException("ReactionFrontSolver: Valence list does not match species list");
    */
    
    m_species.valences.resize(m_species.size);
    m_species.ionicStrengthCoefficients.resize(m_species.size);
    for(unsigned spIndx=0; spIndx < m_species.size; ++spIndx){
      m_species.valences[spIndx] = GPChemistry::GetValence(m_species.names[spIndx]);
      if (rank == 0)   std::cout<< "    " << m_species.names[spIndx] << ": " << m_species.valences[spIndx] << std::endl ;
      m_species.ionicStrengthCoefficients[spIndx] = 0.5*m_species.valences[spIndx]*m_species.valences[spIndx];
    }
    
    GetVectorAttribute(m_species.eqConstants,speciesNode, "equilibrium_constants");
    
    m_species.brineConcentrations.resize(m_species.size);
  }
  
  // Elements
  /////////////////////
  if (rank == 0)  std::cout << "Elements" << std::endl;
  HierarchicalDataNode* elementNode = hdn->GetChild("Elements");
  if(!elementNode)  
      throw GPException("ReactionFrontSolver: Must have Elements defined in the input file");
  {
     GetVectorAttribute(m_elements.names,elementNode, "names");
     m_elements.size = m_elements.names.size();

     for(unsigned i = 0; i < m_elements.size; ++i)
       m_elements.indices[m_elements.names[i]] = i;
     
     GetVectorAttribute(m_elements.brineConcentrations, elementNode,"brine_concentrations");
     if(m_elements.brineConcentrations.size() != m_elements.size )
         throw GPException("ReactionFrontSolver: Number of brine concentrations must match number of elements. ");
  

     // get element count in derived species
     m_elements.species.resize(m_elements.size,rArray1d(m_species.size,0.0));
     
     for(unsigned i =0; i < m_species.size; ++i){
       
       std::string formula = m_species.names[i];
       std::map<std::string, realT> elMap;
       GPChemistry::GetElementCountFromFormula(formula, elMap);

       std::map<std::string, realT>::iterator
         itr = elMap.begin(),
         iend = elMap.end();
       for(;itr != iend; ++itr){
         const std::string& el = itr->first;
         realT& count = itr->second;
         if(isMember(el,m_elements.indices) ){  // ignores H,O
           unsigned elIndex = m_elements.indices[el];
           m_elements.species[elIndex][i] = count;
         }
       } 
     }

      // print data
     for(unsigned i =0; i < m_elements.size; ++i){
       if (rank == 0){
         std::cout << "    " << m_elements.names[i] << ":\t";
         for(unsigned j = 0; j < m_species.size; ++j) 
             std::cout << " " << m_elements.species[i][j];
         std::cout << std::endl;
       }
     }

     /*
     svector spStrs;
     GetVectorAttribute(spStrs,elementNode, "species",";");
     if(spStrs.size() != m_elements.size)
       throw GPException("ReactionFrontSolver: Number of elements must match list of derived species.");
     
     m_elements.species.resize(spStrs.size()); 
     for(int i =0; i < spStrs.size(); ++i){
       rArray1d eqVector(m_species.size,0.0);
       svector subStrs = Tokenize(spStrs[i],",");
       for(int j =0; j < subStrs.size(); ++j){
         realT coef = 1.0;
         
         svector coef_sp = Tokenize(subStrs[j],"*");
         std::string sp;
         
         if(coef_sp.size() == 1 ){
           sp = coef_sp[0];
         } else if(coef_sp.size() == 2){
           sp = coef_sp[1];
           coef = fromString<realT>(coef_sp[0]); // if element occurs more than once. 
         }
        
         Trim(sp);
         if(isMember(sp,m_species.indices)){
           int indx = m_species.indices[sp];
           eqVector[indx] += coef;
         } else {
           throw GPException("Element component " + sp + " not found in species list");
         }
       }
       
       m_elements.species[i] = eqVector;
      
       // print data
       if (rank == 0){
         std::cout << "    " << m_elements.names[i] << ":\t";
         for(int j =0; j < eqVector.size(); ++j) std::cout << " "<< eqVector[j];
         std::cout << std::endl;
       }
       
     }
     */
    
  }
  
  // get solid phases
  ///////////////////
  if (rank == 0)  std::cout << "Solid Phases" << std::endl;
  HierarchicalDataNode* phasesNode = hdn->GetChild("SolidPhases");
  if(!phasesNode)  
      throw GPException("ReactionFrontSolver: Must have SolidPhases defined in the input file");
  {
     GetVectorAttribute(m_solidphases.names,phasesNode, "names");
     GetVectorAttribute(m_solidphases.vfs,phasesNode, "volume_fractions");
     GetVectorAttribute(m_solidphases.MMs,phasesNode, "molar_masses"); // not used
     GetVectorAttribute(m_solidphases.rhos,phasesNode, "densities");  // not used
     m_solidphases.size = m_solidphases.names.size();
     for(unsigned i =0; i < m_solidphases.size; ++i)
       m_solidphases.indices[m_solidphases.names[i]] = i;
  }
    
  // get equilibrium data
  ///////////////////////
  if (rank == 0)  std::cout << "Equilibrium Data" << std::endl;
  HierarchicalDataNode* eqNode = hdn->GetChild("EquilibriumData");
  if(!eqNode)  
      throw GPException("ReactionFrontSolver: Must have EquilibriumData defined in the input file");
  {
     sArray1d analyticStrs, eqStrs;
     GetVectorAttribute(eqStrs,eqNode, "equilibrium_equations",";");
     GetVectorAttribute(analyticStrs,eqNode, "logKdata",";");
     
     
     for(unsigned j = 0; j < eqStrs.size(); ++j){
       sArray1d subStrs = Tokenize(eqStrs[j],"=");
       Trim(subStrs[0]);
       std::string kstr =  subStrs[0];
       sArray1d numDenom = Tokenize(subStrs[1],"/");
       
        rArray1d eqVector(m_species.size,0.0);
        
        // numerator
        sArray1d nums = Tokenize(numDenom[0],"[");
        for(unsigned i =0; i < nums.size(); ++i){
          realT p = 1.0;
          if(nums[i].size() > 0){
            sArray1d sp_pow = Tokenize(nums[i],"]"); // check for powers
            
            Trim(sp_pow[0]);
            if(sp_pow[0].size() > 0  && !streq(sp_pow[0],"H2O") ){
              
              if(!isMember(sp_pow[0],m_species.indices))
                 throw GPException("Equilibrium equation component " + sp_pow[0] + " not found in species list");
                 
              int indx = m_species.indices[sp_pow[0]];
            
              // get power
              if(sp_pow.size() > 1){
                Trim(sp_pow[1]);
                if (sp_pow[1].size() > 0){
                std::replace(sp_pow[1].begin(),sp_pow[1].end(),'^',' '); // remove "^"
                 p = fromString<realT>(sp_pow[1]); 
                }
              }
              eqVector[indx] += p;
            }
          }
        }
        
        // denominator
       if(numDenom.size() > 1){
          sArray1d denoms = Tokenize(numDenom[1],"[");
          for(unsigned i =0; i < denoms.size(); ++i){
            realT p = 1.0;
            if(denoms[i].size() > 0){
              sArray1d sp_pow = Tokenize(denoms[i],"]"); // check for powers
            
              Trim(sp_pow[0]);
              if(sp_pow[0].size() > 0 && !streq(sp_pow[0],"H2O") ){
               
                if(!isMember(sp_pow[0],m_species.indices) )
                   throw GPException("Equilibrium equation component " + sp_pow[0] + " not found in species list");
                 
                int indx = m_species.indices[sp_pow[0]];
            
                // get power
                if(sp_pow.size() > 1){
                  Trim(sp_pow[1]);
                  if (sp_pow[1].size() > 0){
                  std::replace(sp_pow[1].begin(),sp_pow[1].end(),'^',' '); // remove "^"
                  p = fromString<realT>(sp_pow[1]); 
                  }
                }
                eqVector[indx] -= p;
              }
            }
          }
        }
       
        m_equilibriumData.equations[kstr] = eqVector;
        
       // print data
        if (rank == 0){
          std::cout << "    " << kstr << ":\t";
          for(unsigned i =0; i < eqVector.size(); ++i) std::cout << " "<< eqVector[i];
          std::cout << std::endl;
        }
     }
     
     if (rank == 0)  std::cout << "Log K data" << std::endl;
     for(unsigned j = 0; j < analyticStrs.size(); ++j){
       sArray1d aSubStr = Tokenize(analyticStrs[j],",");
       std::string kstr =  aSubStr[0];
       Trim(kstr);
        rArray1d analyticVector(5,0.0);
        unsigned iend = (6<=aSubStr.size())? 6: aSubStr.size();
        for(unsigned i = 1; i < iend ; ++i){
          analyticVector[i-1] = fromString<realT>(aSubStr[i]);
        }
        m_equilibriumData.log10Kdata[kstr] = analyticVector;
        
       // print data
        if (rank == 0){
          realT logK =  Log10K(analyticVector, m_Tc); // calculate value at default temperature
          
          std::cout << "    " << kstr << ":\t" << logK << "\t(";
          for(unsigned i =0; i < analyticVector.size(); ++i) std::cout << " "<< analyticVector[i];
          std::cout << ")" << std::endl;
        }
     }
    
  }
  
  // get solid fronts
  ///////////////////
  if (rank == 0)  std::cout << "Solid fronts" << std::endl;
  HierarchicalDataNode* frontsNode = hdn->GetChild("SolidFronts");
  if(!frontsNode)  
      throw GPException("ReactionFrontSolver: Must have SolidFronts defined in the input file");
  {
     GetVectorAttribute(m_fronts.names,frontsNode, "names");
     GetVectorAttribute(m_fronts.eqConstants,frontsNode, "equilibrium_constants");
     GetVectorAttribute(m_fronts.initial_positions,frontsNode, "initial_positions");
     GetVectorAttribute(m_fronts.molar_densities,frontsNode, "molar_densities");
     
     m_fronts.size = m_fronts.names.size();
     m_fronts.current_positions = m_fronts.initial_positions;
     m_fronts.velocities.resize(m_fronts.size);
   
     m_fronts.minimum_width =frontsNode->GetAttributeOrDefault("minimum_width", "1e-4 mm");
  
     // front conditions
     sArray1d frontCondStrs;
     GetVectorAttribute(frontCondStrs,frontsNode, "front_conditions",";");
     m_fronts.front_conditions.resize(frontCondStrs.size());
     for(unsigned i =0; i < frontCondStrs.size(); ++i){
       sArray1d frontCondStrsB = Tokenize(frontCondStrs[i],",");
       m_fronts.front_conditions[i].resize(frontCondStrsB.size());
       for(unsigned ii =0; ii < frontCondStrsB.size(); ++ii){
         sArray1d lhs_rhs = Tokenize(frontCondStrsB[ii],"="); // lhs,rhs
         
         rArray1d& eqVector = m_fronts.front_conditions[i][ii];
         eqVector =  rArray1d(m_species.size,0.0);
         
         // lhs
         sArray1d lhs = TokenizeSeq(lhs_rhs[0]," +"); // lhs
         for(unsigned iii = 0; iii < lhs.size(); ++iii){
           realT coef = 1.0;
           std::string spStr;
           sArray1d coef_sp = Tokenize(lhs[iii],"*");
           if(coef_sp.size() > 1){
             size_t indx = coef_sp[1].find("_");
             spStr = coef_sp[1].substr(indx+1);
             Trim(coef_sp[0]);
             coef =  fromString<realT>(coef_sp[0]);
             if (rank == 0)  std::cout << coef <<std::endl; 
           } else {
             size_t indx = coef_sp[0].find("_");
             spStr = coef_sp[0].substr(indx+1);
           }
           Trim(spStr);
           if(isMember(spStr,m_species.indices)){
             int indx = m_species.indices[spStr];
             eqVector[indx] += coef;
           } else if( !streq(spStr,"0") ) {
                   throw GPException("Front condition component " + spStr + " not found in species list");
           }
         }
         
         // rhs
         sArray1d rhs = TokenizeSeq(lhs_rhs[1]," +"); // rhs
         for(unsigned iii = 0; iii < rhs.size(); ++iii){
           realT coef = 1.0;
           std::string spStr;
           sArray1d coef_sp = Tokenize(rhs[iii],"*");
           if(coef_sp.size() > 1){
             size_t indx = coef_sp[1].find("_");
             spStr = coef_sp[1].substr(indx+1);
             Trim(coef_sp[0]);
             coef =  fromString<realT>(coef_sp[0]);
           } else {
             size_t indx = coef_sp[0].find("_");
             spStr = coef_sp[0].substr(indx+1);
           }
           Trim(spStr);
           if(spStr.size() > 0 && isMember(spStr,m_species.indices)){
             int indx = m_species.indices[spStr];
             eqVector[indx] -= coef;
           } else if( !streq(spStr,"0") ) {
                   throw GPException("Front condition component " + spStr + " not found in species list");
           }
         }
         
        // print data
        if (rank == 0){  std::cout << "    f" << i << " " << ii << ":\t";
          for(unsigned k =0; k < eqVector.size(); ++k) std::cout << " "<< eqVector[k];
          std::cout << std::endl;
        }
         
       }
     }
     for(unsigned i =0; i < m_fronts.front_conditions.size(); ++i){
       if(m_fronts.front_conditions[i].size() != m_elements.size-1){
          throw GPException("Error: Number of front conditions at each front must be one less than the number of elements");
       }
     }
     
     // dissolution rates
     sArray1d dissRateStrs;
     GetVectorAttribute(dissRateStrs,frontsNode, "dissolution_rates",";");
     m_fronts.dissolution_rates.resize(dissRateStrs.size());
     for(unsigned i =0; i < dissRateStrs.size(); ++i){
       
       rArray1d& eqVector = m_fronts.dissolution_rates[i];
       eqVector =  rArray1d(m_species.size,0.0);
         
       // lhs
       sArray1d lhs = TokenizeSeq(dissRateStrs[i]," +"); // lhs
       for(unsigned iii = 0; iii < lhs.size(); ++iii){
         realT coef = 1.0;
         std::string spStr;
         sArray1d coef_sp = Tokenize(lhs[iii],"*");
         if(coef_sp.size() > 1){
           size_t indx = coef_sp[1].find("_");
           spStr = coef_sp[1].substr(indx+1);
           Trim(coef_sp[0]);
           coef =  fromString<realT>(coef_sp[0]);
         } else {
           size_t indx = coef_sp[0].find("_");
           spStr = coef_sp[0].substr(indx+1);
         }
         Trim(spStr);
         if(isMember(spStr,m_species.indices)){
           int indx = m_species.indices[spStr];
           eqVector[indx] += coef;
         } else if( !streq(spStr,"0") ) {
                  throw GPException("Front condition component " + spStr + " not found in species list");
         }
       }
     }
    
  }
  
  // get regions
  ///////////////////
  HierarchicalDataNode* regionsNode = hdn->GetChild("Regions");
  if(!regionsNode)  
      throw GPException("ReactionFrontSolver: Must have Regions defined in the input file");
  {
     sArray1d solidStrs;
     GetVectorAttribute(m_regions.names,regionsNode, "names");
     GetVectorAttribute(solidStrs,regionsNode, "solid_phases",";");

     GetVectorAttribute(m_regions.tortuosities,regionsNode, "tortuosities");
     
     m_regions.solidPhases.resize(solidStrs.size() );
     m_regions.phis.resize(solidStrs.size());
     if(m_regions.tortuosities.size() == 0) 
        m_regions.tortuosities.resize(solidStrs.size(),1.0);
     
     if(m_regions.solidPhases.size() != m_regions.names.size())
       throw GPException("ReactionFrontSolver: Region names and solid phases are not the same size");
     
     if(m_regions.solidPhases.size()-1 != m_fronts.names.size())
       throw GPException("ReactionFrontSolver: number of reaction fronts should be one less than number of regions");
     
     // set region solid phases, set porosities
     for(unsigned i =0; i < solidStrs.size(); ++i){
       m_regions.solidPhases[i] = Tokenize(solidStrs[i],",");
       realT phi = 1.0;
       for(unsigned j =0; j < m_regions.solidPhases[i].size(); ++j){
          Trim(m_regions.solidPhases[i][j]);
          std::string& solidPhaseStr = m_regions.solidPhases[i][j];
          int indx = m_solidphases.indices[solidPhaseStr];
          phi -= m_solidphases.vfs[indx];
       }
       if(phi < 0.0)
         throw GPException("ReactionFrontSolver: Solid phase porosity < 0 for region "+ m_regions.names[i]);
       m_regions.phis[i] = phi;   
     }
     
     // print data
     std::cout << "Region Porosities: ";
     for(unsigned i =0; i < solidStrs.size(); ++i){
       std::cout << " " << m_regions.phis[i];
     }
     std::cout << std::endl;

     // region sizes
     m_regions.lengths.resize(m_fronts.size); // last region considered infinite in extent
     SetRegionLengths();
  }
  
  // Effective diffusivity
  /////////////////////////
  m_EffectiveDiffusivityFunction = hdn->GetAttributeString("EffectiveDiffusivityFunction");
    
  FunctionManager& functionManager = FunctionManager::Instance();
  if(!m_EffectiveDiffusivityFunction.empty() ){
      m_effDiffFuncPtr = &(functionManager.GetFunction(m_EffectiveDiffusivityFunction));
  } else {
      m_effDiffFuncPtr = NULL;
  }
  
  // default diffusivity for solid phase - scaled by permeability - only used if eff diff func not supplied
  m_D = speciesNode->GetAttributeOrDefault("diffusivity", "1e-10m^2/s");// "1e-10m^2/s"
  
  // Face set
  ////////////
  m_faceSetName = hdn->GetAttributeString("faceset");
  if(m_faceSetName.empty()){
    throw GPException("ReactionFrontSolver: faceset attribute was not defined.");
  }
  
  m_recordFrontSpecies = hdn->GetAttributeOrDefault("recordFrontSpecies", false);
  m_recordBrineSpecies = hdn->GetAttributeOrDefault("recordBrineSpecies", false);
  
  m_activityFunction = ExtendedDebyeHuckelActivityFunction(m_species.names,m_Tc);
  

  // Cylindrical growth
  //m_useCylindricalGrowth = hdn->GetAttributeOrDefault("useCylindricalGrowthModel", false);
  m_Ro = hdn->GetAttributeOrDefault("cylindricalGrowthRadius", -1.0);
  if(m_Ro > 0.0){
    m_useCylindricalGrowth = true;
    std::cout << "Using cylindrical growth model, inner radius = " << m_Ro << std::endl;
  } else {
    m_useCylindricalGrowth = false;
  }

  // Timestep

  m_min_dt = hdn->GetAttributeOrDefault<realT>("initial_dt",0.001);
  m_max_dt = hdn->GetAttributeOrDefault<realT>("max_dt",1);


  // verbosity
  //unsigned verbosityLevel = hdn->GetAttributeOrDefault< unsigned >("verbosity",unsigned(LOG_NOTHING));
  //m_verbosity = (logLevelT) verbosityLevel;
  m_verbosity = hdn->GetAttributeOrDefault< logLevelT >("verbosity",LOG_NOTHING);

}

void ReactionFrontSolver::RegisterFields( PhysicalDomainT& domain )
{
  // register elements
  for(sArray1d::size_type i =0; i < m_elements.size; ++i){
    domain.m_feFaceManager.AddKeylessDataField( FieldInfo::realField,  m_elements.names[i], true, true );   
        domain.m_feFaceManager.AddKeylessDataField( FieldInfo::realField,  m_elements.names[i]+"_ReactionRate", true, true ); 
  }

  domain.m_feFaceManager.AddKeylessDataField( FieldInfo::realField, "pH", true, true );
  
  // register reaction front positions
  for(sArray1d::size_type i =0; i < m_fronts.size; ++i)
    domain.m_feFaceManager.AddKeylessDataField( FieldInfo::realField,  "x_rf_"+toString(i), true, true ); 
  
  // record data pointers
  for(sArray1d::size_type i =0; i < m_elements.size; ++i){
    Array1dT<realT>* ptr = &domain.m_feFaceManager.GetFieldData<realT>(m_elements.names[i]);
    m_elementConcPtrs.push_back(ptr);

        Array1dT<realT>* rrPtr = &domain.m_feFaceManager.GetFieldData<realT>(m_elements.names[i]+"_ReactionRate");
    m_elementReactionRatePtrs.push_back(rrPtr);
  }
  
  for(sArray1d::size_type i =0; i < m_fronts.size; ++i){
    Array1dT<realT>* ptr = &domain.m_feFaceManager.GetFieldData<realT>("x_rf_"+toString(i));
    m_frontPositionPtrs.push_back(ptr);
  }
  
  if(m_recordFrontSpecies){
    m_frontSpeciesConcPtrs.resize(m_species.size*m_fronts.size);
    for(sArray1d::size_type j =0; j < m_species.size; ++j){  
      std::string sp = m_species.names[j];
      std::replace( sp.begin(), sp.end(), '+', 'p'); //+ -> p
      std::replace( sp.begin(), sp.end(), '-', 'm'); //- -> m 
      std::replace( sp.begin(), sp.end(), '(', '_'); //( -> _ 
      std::replace( sp.begin(), sp.end(), ')', '_'); //) -> _ 
      for(sArray1d::size_type i =0; i < m_fronts.size; ++i) {
        std::string spFrtStr = sp+"_" + toString(i);
        domain.m_feFaceManager.AddKeylessDataField( FieldInfo::realField,  spFrtStr, true, true ); 
        m_frontSpeciesConcPtrs[i*m_species.size + j] 
            = &domain.m_feFaceManager.GetFieldData<realT>(spFrtStr);
      }
    }
  }
     
  if(m_recordBrineSpecies){
    m_brineSpeciesConcPtrs.resize(m_species.size);
    for(sArray1d::size_type j =0; j < m_species.size; ++j){  
      std::string sp = m_species.names[j];
      std::replace( sp.begin(), sp.end(), '+', 'p'); //+ -> p
      std::replace( sp.begin(), sp.end(), '-', 'm'); //- -> m 
      std::replace( sp.begin(), sp.end(), '(', '_'); //( -> _ 
      std::replace( sp.begin(), sp.end(), ')', '_'); //) -> _ 
      
      std::string spFrtStr = sp+"_br";
      domain.m_feFaceManager.AddKeylessDataField( FieldInfo::realField,  spFrtStr, true, true ); 
      m_brineSpeciesConcPtrs[j] 
            = &domain.m_feFaceManager.GetFieldData<realT>(spFrtStr);
      
    }
  }
    
 
}

void ReactionFrontSolver::Initialize( PhysicalDomainT& domain, SpatialPartition& partition ){
  
  // Unit scale
  /////////////
  {
    using namespace GPUnits;
    UnitManager& um = UnitManager::Instance(); 
    m_convertToMolPerL = um.ConvertTo("mol/l",1);
    m_convertFromMolPerL=  1.0/m_convertToMolPerL;
  }
   
  // Face set
  ////////////
  m_faceSet = &(domain.m_feFaceManager.GetSet(m_faceSetName));
  m_numFaces = m_faceSet->size();

  // Reaction Fronts
  ///////////////////

  unsigned n = m_species.size;
  unsigned na = m_elements.size + 1; // +1 = H
  unsigned nk = n - na;
  
  // equilibrium data for fronts
  m_KK = rArray2d(nk,n);
  m_pp= rArray1d(nk);
  for(unsigned i =0; i < nk; ++i){
    std::string& kstr = m_species.eqConstants[i];
    for(unsigned j=0; j < n; ++j) m_KK(i,j) = m_equilibriumData.equations[kstr][j];  
    rArray1d& log10Kdata = m_equilibriumData.log10Kdata[kstr];
    m_pp[i] = log(10)* Log10K(log10Kdata, m_Tc); // convert log_10 data to log_e

    // rescale
    realT kMax = fabs(m_KK(i,0));
    for(unsigned j=0; j < n; ++j){
    	if(kMax < fabs(m_KK(i,j)) ) kMax = fabs(m_KK(i,j));
    }
    for(unsigned j=0; j < n; ++j) m_KK(i,j) /= kMax;
    m_pp[i] /= kMax;
  }

  // DEBUG 
/*
  rArray2d temp =  m_KK;
  rArray1d pp = m_pp;
  std::vector<int> pivot_cols;

  std::cout << "K" << std::endl;
  for(int i =0; i < nk; ++i){
    for(int j=0; j < n; ++j) std::cout << temp(i,j)<< " ";
    std::cout << std::endl;
  }

  for(int j=0; j < nk; ++j) std::cout << pp(j)<< " ";
  std::cout << std::endl;

  reducedRowEchelonForm(temp, pp, pivot_cols);
  std::cout << "rref" << std::endl;
  for(int i =0; i < nk; ++i){
    for(int j=0; j < n; ++j) std::cout << temp(i,j)<< " ";
    std::cout << std::endl;
  }

  for(int j=0; j < nk; ++j) std::cout << pp(j)<< " ";
  std::cout << std::endl;

  exit(0);
*/
  //
 
  // equilibrium data for solid phases at fronts 
  m_Kf = rArray2d(m_fronts.size,n);
  m_pf= rArray1d(m_fronts.size);
  for(unsigned i =0; i < m_fronts.size; ++i){
    std::string& kstr = m_fronts.eqConstants[i];
    for(unsigned j=0; j < n; ++j) m_Kf(i,j) = m_equilibriumData.equations[kstr][j];  
    rArray1d& log10Kdata = m_equilibriumData.log10Kdata[kstr];
    m_pf[i] = log(10)* Log10K(log10Kdata, m_Tc); // convert log_10 data to log_e

    // rescale
    realT kMax = fabs(m_Kf(i,0));
    for(unsigned j=0; j < n; ++j){
    	if(kMax < fabs(m_Kf(i,j)) ) kMax = fabs(m_Kf(i,j));
    }
    for(unsigned j=0; j < n; ++j) m_Kf(i,j) /= kMax;
    m_pf[i] /= kMax;
  }
  
  rArray2d AA(na,n);
  rArray1d bb(na);
  for(unsigned i =0; i < na-1; ++i){
    for(unsigned j=0; j < n; ++j){
       AA(i,j) = m_elements.species[i][j];  
    }
    
    m_elements.brineConcentrations[i] *= m_convertToMolPerL;
    bb[i] = m_elements.brineConcentrations[i];
  }
  
  // charge balance
  for(unsigned j=0; j < n; ++j) AA(na-1,j) = m_species.valences[j];
  bb[na-1] = 0;
  realT maxValence = 1.0;  // rescale
  for(unsigned j=0; j < n; ++j){
	  if( fabs( AA(na-1,j) ) > maxValence) maxValence = fabs(AA(na-1,j));
  }
  for(unsigned j=0; j < n; ++j) AA(na-1,j) /= maxValence;
  bb[na-1] /= maxValence;  // redundant but included in case alkalinity added later.

    
  // Brine system
  // m_brineSystem = EquilibriumSystem(m_KK, m_pp, AA, bb, m_species.names,m_Tc);
  m_brineSystem.Initialize(m_KK, m_pp, AA, bb, m_species.names,m_Tc);
  
  // Initial guess for brine species
  //////////////////////////////////
  
  
  m_species.brineConcentrations = rArray1d (2*n+1,1e-7);  
  // ~ mean of elements
  for(unsigned i=0; i < m_elements.size; ++i){
    realT net = 0.0;
    for(unsigned ii=0; ii < n; ++ii) net += m_elements.species[i][ii];
    
    for(unsigned ii=0; ii < n; ++ii){
      if(m_elements.species[i][ii] > 0) m_species.brineConcentrations(ii) += m_elements.brineConcentrations(i)*m_elements.species[i][ii]/net;
    }
  }
  
  

  //m_species.brineConcentrations = rArray1d (2*n+1,1e-7);  
  for(unsigned i =0; i < n; ++i) m_species.brineConcentrations(i+n) = log(std::max(m_species.brineConcentrations(i),1e-64));
  
  // Build reaction front equations
  /////////////////////////////////
  rArray2d AA_rf;
  rArray1d bb_rf;
  BuildReactionFrontEquations(AA_rf,bb_rf);
  m_rfSystem.Initialize(m_KK,m_pp,m_Kf,m_pf,AA_rf,bb_rf,m_fronts.size, m_species.names,m_Tc);
  
  // Set initial front positions
  for(lSet::const_iterator fitr=m_faceSet->begin() ; fitr!=m_faceSet->end() ; ++fitr){
    localIndex fc = *fitr;
    for(sArray1d::size_type i =0; i < m_fronts.size; ++i){
      Array1dT<realT>* ptr = m_frontPositionPtrs[i];
      (*ptr)[fc] = std::max((*ptr)[fc],m_fronts.initial_positions[i]);
    }
  }
  
  // Solve for default concentrations (used as initial guess). 
  FindReactionFrontEquilibrium(false);
  m_X_br = m_species.brineConcentrations;
  m_X_s = m_species.frontConcentrations;
  
 
  if(m_recordBrineSpecies)
  {
    for(lSet::const_iterator fitr=m_faceSet->begin() ; fitr!=m_faceSet->end() ; ++fitr){
      localIndex fc = *fitr;
      for(unsigned i =0; i < m_species.size; ++i)
      { 
        (*m_brineSpeciesConcPtrs[i])[fc] = m_species.brineConcentrations[i];
      }
    }
  }

  if(m_recordFrontSpecies)
  {
    for(lSet::const_iterator fitr=m_faceSet->begin() ; fitr!=m_faceSet->end() ; ++fitr){
      localIndex fc = *fitr;
      for(unsigned i =0; i < m_species.size*m_fronts.size; ++i)
      { 
        (*m_frontSpeciesConcPtrs[i])[fc] = m_species.frontConcentrations[i];
      }
    }
  }

  m_useApertureBasedFluxLimiter = domain.m_feFaceManager.HasField<realT>(ApertureStr);

  // Aperture - only used if aperture flux limiter is invoked
  if(m_useApertureBasedFluxLimiter){
	  m_FaceAperturesPtr = &(domain.m_feFaceManager.GetFieldData<realT>(ApertureStr));
  }

}

void ReactionFrontSolver::InitializeCommunications( PartitionBase& partition )
{
  syncedFields.clear();
  sArray1d& synced_face_fields = syncedFields[PhysicalDomainT::FiniteElementFaceManager];
  for(sArray1d::size_type i =0; i < m_elements.size; ++i)
    synced_face_fields.push_back(m_elements.names[i]+"_ReactionRate");

  for(sArray1d::size_type i =0; i < m_fronts.size; ++i)
	synced_face_fields.push_back("x_rf_"+toString(i));

  partition.SetBufferSizes(syncedFields, CommRegistry::ReactionFrontSolver);
}


void ReactionFrontSolver::SetRegionLengths()
{
  m_regions.lengths[0] = m_fronts.current_positions[0];
  for(unsigned i =1; i < m_regions.lengths.size();++i)
  {
    m_regions.lengths[i] = m_fronts.current_positions[i] - m_fronts.current_positions[i-1];
  }
}



/**
 * 
 * 
**/

double ReactionFrontSolver::TimeStep( const realT& time ,
                           const realT& dt,
                           const int cycleNumber,
                           PhysicalDomainT& domain,
                           const sArray1d& namesOfSolverRegions ,
                           SpatialPartition& partition,
                           FractunatorBase* const fractunator )
{  


  m_stabledt.m_maxdt = std::numeric_limits<double>::max();
  iArray1d& faceGhostRank       = domain.m_feFaceManager.GetFieldData<FieldInfo::ghostRank>();


  int rank(0);
  //realT t1;
  #if GPAC_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif

  if(rank==0){
     std::cout << "ReactionFrontSolver - Start of timestep \n";
     //t1 = getcputime();
  }
  
  //int n = m_species.size; // number of species
  //std::cout << "FaceSet size " << m_faceSet->size() << std::endl;
  
  // pH
  rArray1d& pHField = domain.m_feFaceManager.GetFieldData<realT>("pH");
  int H_Index = m_species.indices["H+"]; 

    
  // loop over faces
  for(lSet::const_iterator fitr=m_faceSet->begin() ; fitr!=m_faceSet->end() ; ++fitr){
    
    localIndex fc = *fitr;
    if(faceGhostRank[fc] < 0){
  
      // Get elements species from face
      for(unsigned i =0; i < m_elements.size; ++i){
        m_elements.brineConcentrations(i) = (*m_elementConcPtrs[i])[fc]*m_convertToMolPerL;
      }
    
      // Get reaction front positions
      for(unsigned i =0; i < m_fronts.size; ++i){
        m_fronts.current_positions[i] =  (*m_frontPositionPtrs[i])[fc];
      }
      SetRegionLengths();
    
      // Update brine equilibrium equations
      rArray1d& bb = m_brineSystem.Get_b();
      for(unsigned i =0; i < m_elements.size; ++i){
        bb(i) = m_elements.brineConcentrations[i];
      }

      if(m_recordFrontSpecies) UpdateInitialGuess(fc);

      // Solve reaction front system
      int errFlag;
      errFlag = FindReactionFrontEquilibrium();
   
      if (errFlag == 0 || true){ // solution has converged
        // Update brine pH
        //pHField[fc] = -log10(m_species.brineConcentrations(H_Index));
        pHField[fc]  = -log10( exp(m_species.brineConcentrations(H_Index+ m_species.size) ));

        // Calculate front velocities
        /////////////////////////////
        for(unsigned fr =0; fr < m_fronts.size; ++fr){

          // Dissolving species concentrations in mol/L
          realT C = 0.0;
          realT C_m = 0.0;
          realT C_p = 0.0;

          unsigned fOff = m_species.size*fr; // front offset
          unsigned fOff_m = m_species.size*(fr-1); // -ve front offset
          unsigned fOff_p = m_species.size*(fr+1); // +ve front offset
          for(unsigned i =0; i < m_species.size; ++i){
            C += m_species.frontConcentrations[i+fOff]*m_fronts.dissolution_rates[fr][i];
            if( fr!=0 ){C_m += m_species.frontConcentrations[i+fOff_m]*m_fronts.dissolution_rates[fr][i];}
              else {C_m += m_species.brineConcentrations[i]*m_fronts.dissolution_rates[fr][i];}
            if( fr!=m_fronts.size-1 ) C_p += m_species.frontConcentrations[i+fOff_p]*m_fronts.dissolution_rates[fr][i];
          }

          // Porosities on either side of the front
          realT phi_m = m_regions.phis[fr];
          realT phi_p = m_regions.phis[fr+1];
          realT tau_m = m_regions.tortuosities[fr];
          realT tau_p = m_regions.tortuosities[fr+1];

          // Region widths on either side of the front
          realT dx_m = (fr!=0)? m_fronts.current_positions[fr] - m_fronts.current_positions[fr-1]: m_fronts.current_positions[fr];
          realT dx_p = (fr!=m_fronts.size-1)?  m_fronts.current_positions[fr+1] - m_fronts.current_positions[fr]: 0.0;

          if(m_useCylindricalGrowth){

            // Effective radius of next front current front and previous front
            realT Rp = (fr!=m_fronts.size-1)? m_fronts.current_positions[fr+1] + m_Ro: 1e-64;
            realT Rx = m_fronts.current_positions[fr] + m_Ro;
            realT Rm = (fr!=0)? m_fronts.current_positions[fr-1] + m_Ro: m_Ro;

            // correction for cylindrical growth model
            dx_m = Rx*log(Rx/Rm);
            dx_p = (fr!=m_fronts.size-1)?  Rx*log(Rp/Rx) :0.0;
          }


          realT MolarDensity = m_fronts.molar_densities[fr];

          // convert to simulation units
          C   *= m_convertFromMolPerL;
          C_p *= m_convertFromMolPerL;
          C_m *= m_convertFromMolPerL;

          m_fronts.velocities[fr] = ReactionFrontVelocity(C,C_m, C_p, dx_m, dx_p,
                                                        phi_m, phi_p,
                                                        tau_m, tau_p,
                                                        MolarDensity);

          //std::cout << "m_fronts.velocities[fr] " << m_fronts.velocities[fr] << std::endl;
        }
    

        // Front growth limiter
        ////////////////////////
        realT invFrontGrowthLimiter = 1.0;
        for(unsigned fr =0; fr < m_fronts.size; ++fr){
          realT v = 0.0;
          if(fr == 0 ){
             v = 10.0*dt*fabs( m_fronts.velocities[0] )/m_fronts.current_positions[0];// Limit front growth to 1/10 current length
          } else {
             v = 10.0*dt*fabs( m_fronts.velocities[fr]-m_fronts.velocities[fr-1] )/fabs(m_fronts.current_positions[fr] - m_fronts.current_positions[fr-1]);
          }
          if(invFrontGrowthLimiter < v) invFrontGrowthLimiter = v;
        }
        realT frontGrowthLimiter = 1.0/invFrontGrowthLimiter;
        realT scaledDt = dt*frontGrowthLimiter;


        // Determine rate of diffusion of elements to brine
        realT dx = m_fronts.current_positions[0];
        if(m_useCylindricalGrowth){
            // correction for cylindrical growth model
          dx = m_Ro*log(1.0+m_fronts.current_positions[0]/m_Ro);
        }

        if(dx > 0.0){

          const realT& phi = m_regions.phis[0];
          const realT& tau = m_regions.tortuosities[0];
          realT Deff = EffectiveDiffusivity(phi,tau);

          // Prevent concentration oscillations due to excessive flux
          if(m_useApertureBasedFluxLimiter ){
             realT aperture = (*m_FaceAperturesPtr)[fc];
             realT dtMax = fabs(aperture)*dx/Deff;
             if(dtMax < scaledDt){
               scaledDt = dtMax;
          	   const realT TINY = 1e-64;
          	   frontGrowthLimiter = scaledDt/(dt+TINY);
             }
      	     m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt,dtMax);
          }


          // diffusion into brine
          for(unsigned i =0; i < m_elements.size; ++i){
            realT gradC = -m_elements.brineConcentrations[i];
            for(unsigned j =0; j < m_species.size;++j){
              gradC += m_elements.species[i][j]*m_species.frontConcentrations(j);
            }
         
            gradC /= dx;

            //dRdT(i) = Deff*gradC;
            (*m_elementReactionRatePtrs[i])[fc] = Deff*gradC*m_convertFromMolPerL*frontGrowthLimiter;
          }


        }
        /*else {
          for(unsigned i =0; i < m_elements.size; ++i){
            (*m_elementReactionRatePtrs[i])[fc] = 0.0;
          }
        }*/
    
        // Move fronts
        ///////////////
   
        // NB if precipitation is to occur - it will need to be included here
        for(unsigned fr =0; fr < m_fronts.size; ++fr){
          if(fr == 0 ){
            m_fronts.current_positions[0] = m_fronts.current_positions[0]+scaledDt*m_fronts.velocities[0];

            // check for precipitation into channel
            if(m_fronts.current_positions[0] <= m_fronts.minimum_width){
              realT shift = m_fronts.minimum_width - m_fronts.current_positions[0];
              m_fronts.current_positions[0] = m_fronts.minimum_width;
              (*m_frontPositionPtrs[0])[fc] = m_fronts.minimum_width;

              if(false){
                // move other fronts outward - seems to lead to instabilities
                for(unsigned ff =1; ff < m_fronts.size; ++ff){
                   m_fronts.current_positions[ff] += shift;
                }

              } else {
            	// freeze fronts
                for(unsigned i =0; i < m_elements.size; ++i) (*m_elementReactionRatePtrs[i])[fc] = 0.0;
                break;
              }
            }
          } else {
            m_fronts.current_positions[fr] = std::max( m_fronts.current_positions[fr]+scaledDt*m_fronts.velocities[fr],
                                                  ((fr > 0)?  m_fronts.current_positions[fr-1]+ m_fronts.minimum_width : m_fronts.minimum_width) );
            m_fronts.current_positions[fr] = std::min( m_fronts.current_positions[fr], 10*m_fronts.current_positions[0]+ fr*m_fronts.minimum_width ); // escaping fronts.
          }

          // std::cout << fr << " " << m_fronts.current_positions[fr] << " "
          //                        << (*m_frontPositionPtrs[fr])[fc] << " " << dt*m_fronts.velocities[fr] << std::endl;

          (*m_frontPositionPtrs[fr])[fc] = m_fronts.current_positions[fr];


          // realT maxdt = 0.1*m_regions.lengths[fr]/fabs(m_fronts.velocities[fr]); // old criteria
          realT maxdt;
          if(fr == 0){
      	    maxdt = 0.1*m_regions.lengths[fr]/fabs(m_fronts.velocities[fr]);
          } else {
    	    maxdt = 0.1*m_regions.lengths[fr]/fabs(m_fronts.velocities[fr]-m_fronts.velocities[fr-1]+1e-64);
          }

          m_stabledt.m_maxdt = std::min(m_stabledt.m_maxdt,maxdt);
                                        
        }
 
        // Record Species
        /////////////////
        if(m_recordBrineSpecies)
        {
          for(unsigned i =0; i < m_species.size; ++i)
          {
            (*m_brineSpeciesConcPtrs[i])[fc] = m_species.brineConcentrations[i];
          }
        }

        if(m_recordFrontSpecies)
        {
          for(unsigned i =0; i < m_species.size*m_fronts.size; ++i)
          {
            (*m_frontSpeciesConcPtrs[i])[fc] = m_species.frontConcentrations[i];
          }
        }
      } else {
    	  // did not converge
         //  for(unsigned i =0; i < m_elements.size; ++i){
         //   (*m_elementReactionRatePtrs[i])[fc] = 0.0;
         // }
      }
    }// ghost rank
  } // loop over faces
  
  partition.SynchronizeFields( syncedFields, CommRegistry::ReactionFrontSolver);

  if(m_stabledt.m_maxdt < m_min_dt || m_stabledt.m_maxdt >= std::numeric_limits<double>::max()){
    m_stabledt.m_maxdt = m_min_dt;
  }

  // limit maximum increase in timestep
  if(m_stabledt.m_maxdt > 2*dt && dt > 0.0){
    m_stabledt.m_maxdt = 2*dt;
  }

  if(m_stabledt.m_maxdt > m_max_dt){
    m_stabledt.m_maxdt = m_max_dt;
  }

  return dt;
}

/**
 * 
 * Find brine concentrations and use result to solve reaction front
 * 
 **/
int ReactionFrontSolver::FindReactionFrontEquilibrium(bool haveCalculatedInitialGuess)
{
  
  const unsigned maxNumItrs = 1000;
  const realT tolx = 2*std::numeric_limits<realT>::epsilon();
  //std::cout << tolx << std::endl;
  //const int n = m_species.size;
  bool rootFound(false); unsigned numItrs(0); realT ff(0.0);
  int err(0);
  
  // Brine system
  /////////////////
  //std::cout << "Brine System"  << std::endl;
  NewtonRaphson( m_brineSystem, m_species.brineConcentrations, maxNumItrs, tolx, rootFound, numItrs, ff);
  
  // Reaction front system
  /////////////////////////
  
  // Update equilibrium equations
  rArray2d& AA_rf = m_rfSystem.Get_A();
  rArray1d& bb_rf = m_rfSystem.Get_b();
  BuildReactionFrontEquations(AA_rf,bb_rf);
  m_rfSystem.UpdateJacobian();
  
  // Initial guess (if not supplied)
  if(haveCalculatedInitialGuess){
    //m_species.frontConcentrations = m_X_s; // use brine equ as initial guess
  } else {
    // initial guess = brine eq at each front
    m_species.frontConcentrations.resize( (2*m_species.size+1)*m_fronts.size);

    unsigned off = m_species.size*m_fronts.size; 
    for(unsigned f =0; f < m_fronts.size; ++f){
      for(unsigned i =0; i < m_species.size; ++i){
        m_species.frontConcentrations(i+f*m_species.size) 
          = m_species.brineConcentrations(i);

        m_species.frontConcentrations(i+f*m_species.size+off) 
          = m_species.brineConcentrations(i+m_species.size);
      }
    
      m_species.frontConcentrations(2*off+f) 
         = m_species.brineConcentrations(2*m_species.size);
    }
  } 
  
    
  // Solve
  //std::cout << "Reaction Front System" << std::endl;
  NewtonRaphson( m_rfSystem, m_species.frontConcentrations, maxNumItrs, tolx,rootFound, numItrs, ff);


  // Try again with another start (or record start if gave reasonable answer)
  if(haveCalculatedInitialGuess && false){
    if(ff > 1.0 ) {
      // initial guess = 0.9*brine eq at each front
      std::cout << "Restarted Reaction Front System" << std::endl;
      m_species.frontConcentrations = m_X_s; // use last sucessful solution
      NewtonRaphson( m_rfSystem, m_species.frontConcentrations, maxNumItrs, tolx,rootFound, numItrs, ff);
      std::cout << rootFound << " " << numItrs << " " << ff << std::endl;
    } else if(haveCalculatedInitialGuess) {
      m_X_s = m_species.frontConcentrations ;
    }
  } 

  if(ff > 1.0){
	  err = 1;
  }

  return err;
}

/**
*
*  Use solution from previous timestep as initial guess for next timestep
*
**/
void ReactionFrontSolver::UpdateInitialGuess(unsigned fc)
{
  // recover  solution from previous timestep
  unsigned nn = m_species.size*m_fronts.size;
  unsigned nC = m_species.size; // number of species per front
  for(unsigned i =0; i < nn; ++i)
  { 
    m_species.frontConcentrations[i] = (*m_frontSpeciesConcPtrs[i])[fc];
  }
       
  // activity/molality
  for(unsigned f = 0; f < m_fronts.size; ++f){
    const realT* XX(&(m_species.frontConcentrations[f*nC]));
    for(unsigned i(0), ii(f*nC); i < nC; ++i,++ii){
      realT logGamma = m_activityFunction.LogIAC(XX, i,nC);
      m_species.frontConcentrations[ii+nn] = log(m_species.frontConcentrations[ii]) + logGamma;
    }
  }
  // ionic strength coefficients
  for(unsigned f = 0; f < m_fronts.size; ++f){
    realT ff = 0.0;
    unsigned cOffset = f*nC;
    for(unsigned j = 0; j < nC; ++j) ff += m_rfSystem.GetIonicStrengthCoefficient(j)*m_species.frontConcentrations[j+cOffset];
    m_species.frontConcentrations[2*nn+f] = ff;
  }
}  



/**
 * 
 * Build governing equations
 * 
**/
void ReactionFrontSolver::BuildReactionFrontEquations(rArray2d& AA, rArray1d& bb)
{
  const realT TINY = 1e-64;
  const unsigned n = m_species.size;
  const unsigned na = m_elements.size; // (1 valence eq +(ne-1) flux eqs)
                                  // missing flux equation -> equilibrium with front
  const int nf = m_fronts.size;
  
  if(bb.size() != na*nf ){
    AA = rArray2d ( na*nf,n*nf);
    bb = rArray1d ( na*nf);
  }
  
  realT maxValence = 1.0;
  for(unsigned ii =0; ii < n; ++ii){
	  if (maxValence < fabs(m_species.valences[ii])) maxValence = fabs(m_species.valences[ii]);
  }

  // loop over fronts
  for(int f = 0; f < nf; ++f){
    int rOff = na*f; // row offset
    int cOff = n*f; // col offset
    int cOff_m = n*(f-1); // col offset minus (previous front concentrations)
    int cOff_p = n*(f+1); // col offset plus (next front concentrations)

    // valence
    bb(rOff) = 0/maxValence;  // redundant but included in case alkalinity added later
    for(unsigned ii =0; ii < n; ++ii){
      AA(rOff,ii+cOff) = m_species.valences[ii]/maxValence;
    }

/*
    // front flux conditions // fixme currently assumes D_eff = phi D

    realT lp = (f<nf-1)? m_regions.lengths[f+1]: TINY;
    realT phi_p = (f<nf-1)? m_regions.phis[f+1]: 0.0;
    realT lm =    m_regions.lengths[f];
    realT phi_m = m_regions.phis[f];

    realT lsum = lp+lm;  
    
    realT scale = -lsum*(phi_m/lm + phi_p/lp);
    realT scale_m = lsum*phi_m/lm;
    realT scale_p = lsum*phi_p/lp;  
    
    for(unsigned i =0; i < na-1; ++i){ // equilibrium with front removes one condition
      bb(rOff+i+1) = 0;
          
      for(unsigned ii =0; ii < n; ++ii){    
        // this front    
        AA(rOff+i+1,ii+cOff) = scale*m_fronts.front_conditions[f][i][ii];
        
        // previous front
        if(f>0){
          AA(rOff+i+1,ii+cOff_m) = scale_m*m_fronts.front_conditions[f][i][ii];
        } else {
          // brine boundary
          bb(rOff+i+1) -= scale_m*m_fronts.front_conditions[f][i][ii]*m_species.brineConcentrations[ii];
        }
        
        // next front
        if(f<nf-1){
          AA(rOff+i+1,ii+cOff_p) = scale_p*m_fronts.front_conditions[f][i][ii];
        }
      }
    }
*/
    // front flux conditions
    realT lp = (f<nf-1)? m_regions.lengths[f+1]: TINY;
    realT Deff_p = (f<nf-1)? 
                      EffectiveDiffusivity(m_regions.phis[f+1],m_regions.tortuosities[f+1]):
                      0.0;
    realT lm =    m_regions.lengths[f];
    realT Deff_m = EffectiveDiffusivity(m_regions.phis[f],m_regions.tortuosities[f]);

    
    realT denom = -(Deff_m/lm + Deff_p/lp);
    realT scale_m = Deff_m/(lm*denom);
    realT scale_p = Deff_p/(lp*denom);  
    
    // cylindrical growth model
    if(m_useCylindricalGrowth){
      realT Rp = (f<nf-1)? m_fronts.current_positions[f+1] + m_Ro: m_fronts.current_positions[f] + m_Ro+TINY;
      realT Rx = m_fronts.current_positions[f] + m_Ro;
      realT Rm = (f>0)? m_fronts.current_positions[f-1] + m_Ro: m_Ro;

      realT Dm = Deff_m/(Rx*log(Rx/Rm));
      realT Dp = Deff_p/(Rx*log(Rp/Rx));

      denom = -(Dm + Dp);


      scale_m = Dm/denom;
      scale_p = Dp/denom;

    }

    for(unsigned i =0; i < na-1; ++i){ // equilibrium with front mineral removes one condition
      bb(rOff+i+1) = 0;

      realT maxFC = TINY;
      for(unsigned ii =0; ii < n; ++ii){
        if(maxFC < fabs(m_fronts.front_conditions[f][i][ii])) maxFC = fabs(m_fronts.front_conditions[f][i][ii]);
      }
          
      for(unsigned ii =0; ii < n; ++ii){    
        // this front    
        AA(rOff+i+1,ii+cOff) = m_fronts.front_conditions[f][i][ii]/maxFC;
        
        // previous front
        if(f>0){
          AA(rOff+i+1,ii+cOff_m) = scale_m*m_fronts.front_conditions[f][i][ii]/maxFC;
        } else {
          // brine boundary
          bb(rOff+i+1) -= scale_m*m_fronts.front_conditions[f][i][ii]*m_species.brineConcentrations[ii]/maxFC;
        }
        
        // next front
        if(f<nf-1){
          AA(rOff+i+1,ii+cOff_p) = scale_p*m_fronts.front_conditions[f][i][ii]/maxFC;
        }
      }
    }
  }
  
  
}


/**
 * 
 * SetMaxStableTimeStep
 *
**/
void ReactionFrontSolver::SetMaxStableTimeStep( PhysicalDomainT& ,
                                                           const sArray1d& namesOfSolverRegions,
                                                           SpatialPartition& partition __attribute__((unused)) )
{
  m_stabledt.m_maxdt =  m_min_dt;

  m_stabledt.m_maxdt *= this->m_courant;

}


/**
 *
 * Effective diffusivity
 * 
**/
inline realT ReactionFrontSolver::EffectiveDiffusivity(realT phi, realT tau){
    
    if(m_effDiffFuncPtr){
      realT x[2] = {phi,tau};
      return (*m_effDiffFuncPtr)(x[0]);
    } else {
      // Simplest expression: Effective diffusivity = porosity * FluidDiffusivity/tortuosity
      // Ulm Lemarchard Heukamp, Elements of chemomechanics of calcium leaching of cement-based materials at different scales. Eng. Frac. Mech. 70 (2003) 871-889
       return phi*m_D/tau;
     }
}

/**
 * 
 *  Calculate the velocity of the reaction front
 * 
**/
inline realT ReactionFrontSolver::ReactionFrontVelocity(realT C,realT C_m,realT C_p,
                                     realT dx_m,realT dx_p,
                                     realT phi_m,realT phi_p,
                                     realT tau_m,realT tau_p,realT MolarDensity){
    // C = concentration of the control species at the front
    // C_m = concentration of the control species on the -ve side of the front
    // C_p = concentration of the control species on the +ve side of the front
    // dx_m = distance to the control species measurement on the -ve side of the front
    // dx_p = distance to the control species measurement on the +ve side of the front
    // phi_m = porosity on the -ve side of the front
    // phi_p = porosity on the +ve side of the front
    // D = diffusivity of the control species (in open water)
    // Molar density = molar density per unit of the solid phase of the species controlling the dissolution rate

    realT dCdx_m = (dx_m > 0.0 )?(C - C_m)/dx_m : 0.0;
    realT dCdx_p = (dx_p > 0.0 )?(C_p - C)/dx_p : 0.0;
    realT D_m = EffectiveDiffusivity(phi_m,tau_m);
    realT D_p = EffectiveDiffusivity(phi_p,tau_p);
    return (D_m*dCdx_m - D_p*dCdx_p )/MolarDensity; // dX_dt 
}

REGISTER_SOLVER( ReactionFrontSolver )
