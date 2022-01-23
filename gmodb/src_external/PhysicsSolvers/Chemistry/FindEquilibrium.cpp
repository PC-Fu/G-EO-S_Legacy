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
 * @file FindEquilibrium.cpp
 * @author walsh24
 * @date July 25, 2011
 */
 
#include <cmath>
#include <vector>

#include "PhysicsSolvers/SolverFactory.h"
#include "FindEquilibrium.h"
#include "Common/Common.h"
#include "Common/intrinsic_typedefs.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/TrilinosUtilities.h"
#include "ObjectManagers/ChemistryManager.h"

using namespace GPChemistry;

namespace{

realT dot(const rArray1d& x,const rArray1d& y ){
  if(x.size() != y.size() ) throw GPException( "Dot product - vector sizes must agree." );
  realT rv(0.0);
  for(rArray1d::size_type i = 0; i < x.size(); ++i) rv += x[i]*y[i];
  return rv;
}
}

   
  void EquilibriumSystem::GetFunction(const rArray1d& X,rArray1d& F){
    int nA = m_b.size();
    int nK = m_p.size();
    int n = nA+nK;

    F.resize(2*n);
    for(int i = 0; i < n; ++i){
      F[i] = X[i] -exp(X[i+n]);
    }
    for(int i(0),ii(n); i < nA; ++i,++ii){
      F[ii] = -m_b[i];
      for(int j = 0; j < n; ++j){
        F[ii] += m_A(i,j)*X[j];
      }
    }
    for(int i(0),ii(n+nA); i < nK; ++i,++ii){
      F[ii] = -m_p[i];
      for(int j = 0; j < n; ++j){
        F[ii] += m_K(i,j)*X[j+n];
      }
    }
   
  }


  void EquilibriumSystem::GetJacobian(const rArray1d& X,rArray2d& gradF){
    int nA = m_b.size();
    int nK = m_p.size();
    int n = nA+nK;

    gradF = rArray2d(2*n,2*n);
    for(int i = 0; i < n; ++i){
      gradF(i,i) = 1;
      gradF(i,i+n) = -exp(X[i+n]);
    }
    for(int i(0),ii(n); i < nA; ++i,++ii){
      for(int j = 0; j < n; ++j){
        gradF(ii,j) = m_A(i,j);
      }
    }
    for(int i(0),ii(n+nA); i < nK; ++i,++ii){
      for(int j = 0; j < n; ++j){
        gradF(ii,j+n) += m_K(i,j);
      }
    }
    
    /*
    std::cout << "find equilibrium" << std::endl; 
    for(int i = 0; i < 2*n; ++i){
    for(int j = 0; j < 2*n; ++j){
    	std::cout << gradF(i,j) << " "; 
    }
    	std::cout << std::endl; 
    }*/
  }
  
  void NewtonRaphson( NonLinearSystemBase& nls, 
                    rArray1d& x, int maxNumItrs, realT tolx){

    
    localIndex n = x.size();
    realT f(0.0);
  
    rArray1d F(n),dx(n),g(n);
    rArray2d gradF(n,n);

    int numItrs(0);
    bool rootFound = false;
    for(;numItrs < maxNumItrs && !rootFound; ++numItrs){

      nls.GetFunction(x,F);
      f = dot(F,F);  
      realT fold = f; rArray1d xold = x;
      nls.GetJacobian(x, gradF);

      // dx
      LinSolve_Local(gradF,dx,F); dx *= -1.0;

      // check size dx
      // max step size
      realT stpmx = 10.0 * std::max( sqrt(dot(xold,xold)),double(n));
  
      // step length
      realT sl = sqrt( dot(dx,dx) );
      if(sl > stpmx)  dx *= stpmx/sl;
    
      realT tst = 1e-64;
      for(localIndex i =0; i< n; ++i){
        realT t = fabs(dx[i])/ std::max(fabs(xold[i]),1.0);
        if(t > tst) tst = t;
      }
      realT lmin = tolx / tst;
     // std::cout << lmin << " " << tolx << " "<< tst << " " <<std::endl;

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
          if(false && l < 1e-64) minFound = true;   
        }
      }

      if (rootFound) break;
    }
    
    std::cout << "numItrs " << numItrs << std::endl;
    std::cout << "f " << f << std::endl;
  
  }
  
namespace {
  /// Data
  realT A[3][9]= { {1,-1,2,1,0,-2,-1,0,0},  // Valence charge balance coeff
                   {0,0,1,1,1,0,0,0,0},     // Ca count
                   {0,0,0,1,1,1,1,1,1} };   // C count
  
  realT K[6][9]= {{ 1,   1,   0,  0,  0,  0,  0,  0,  0},  // [H+][OH-] = Kw
                  { 1,   0,   0,  0,  0,  1, -1,  0,  0},  // [H+][CO3-2] = K2 [HCO3-]
                  { 0,   0,   1, -1,  0,  0,  1,  0,  0},  // [Ca+2][HCO3-] = K3 [CaHCO3+] 
                  { 0,   0,   1,  0, -1,  1,  0,  0,  0},  // {Ca+2]{CO3-2] = K4 {CaCO3]
                  { 1,   0,   0,  0,  0,  0,  1, -1,  0},  // [H+][HCO3-] = K5 [H2CO3]
                  { 1,   0,   0,  0,  0,  0,  1,  0, -1} };   // [H+][HCO3-] = K6 [CO2][H2O]
                  
  realT log10Kdata[6][5] = 
    {{22.801,      -0.010365,    -4787.3,   -7.1321,    0.0      }, // [H+][OH-] = Kw
     {-107.8871,  -0.03252849,   5151.79,   38.92561, -563713.9  }, // [H+][CO3-2] = K2 [HCO3-]
     {-1209.120,     -0.31294,      34765.05,  478.782,    0.0   }, // [Ca+2][HCO3-] = K3 [CaHCO3+]  
     {1228.732,      0.299444,    -35512.75,  485.818,    0.0    }, // [Ca+2][CO3-2] = K4 [CaCO3]
     {1.000393128530148,  0.0,        0.0,      0.0,   0.0      },  // [H+][HCO3-] = K5 [H2CO3]
     {-356.3094, -0.060919964,  21834.37,  126.8339, -1684915 }};   // [H+][HCO3-] = K6 [CO2][H2O]
     
 
  realT CBrine = 0.0001;
  realT CaBrine = 0.0001;
  realT Tc = 25.0;
}
  
//}



//////////////////////////////////////////////////////////////////////////////////////////

// Upate field with function

FindEquilibrium::FindEquilibrium( const std::string& name,
                                  ProblemManagerT* const pm ):
SolverBase(name,pm)
{
}

FindEquilibrium::~FindEquilibrium()
{
  // TODO Auto-generated destructor stub
}

/*
 * <FindEquilibrium name="CaUpdate"       * name of the solver
 *            objecttype="Face"                  * location of field (Default=Node)
 *            species="Scalar Scalar" />   * variableTypes (assumed scalar for all if omitted) 
 */
void FindEquilibrium::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{  
  SolverBase::ReadXML(hdn);

  std::string objectTypeStr = hdn->GetAttributeStringOrDefault("objecttype", PhysicalDomainT::FiniteElementNodeManagerStr() );
  m_objectKey = PhysicalDomainT::GetObjectDataStructureKey(objectTypeStr);
  m_regionName = hdn->GetAttributeStringOrDefault("regionname",""); // used for element regions only
  
  std::string speciesStr = hdn->GetAttributeString("species");
  m_species = Tokenize(speciesStr," ");
    
}

void FindEquilibrium::RegisterFields( PhysicalDomainT& domain )
{

  ObjectDataStructureBaseT& objectManager = domain.GetObjectDataStructure(m_objectKey,m_regionName);

  // register variables
  for(sArray1d::size_type i =0; i < m_species.size(); ++i)
  	objectManager.AddKeylessDataField( FieldInfo::realField,  m_species[i], true, true );     
}

void FindEquilibrium::Initialize( PhysicalDomainT& domain , SpatialPartition& partition ){


  rArray2d KK(6,9);
  for(int i =0; i < 6; ++i){
    for(int j=0; j < 9; ++j) KK(i,j) = K[i][j];	
  }
  rArray2d AA(3,9);
  for(int i =0; i < 3; ++i){
    for(int j=0; j < 9; ++j) AA(i,j) = A[i][j];	
  }
  rArray1d pp(6);
  for(int i =0; i < 6; ++i) pp[i] = log(10)* Log10K(log10Kdata[i], Tc);
  
  rArray1d bb(3);
  bb[0] = 0.0;
  bb[1] = CaBrine;
  bb[2] = CBrine;
  
  m_equilibriumSystem = EquilibriumSystem(KK, pp, AA, bb);
}



/**
 * 
 * 
**/

double FindEquilibrium::TimeStep( const realT& time ,
                                const realT& dt ,
                                const int cycleNumber,
                                PhysicalDomainT& domain ,
                                const sArray1d& namesOfSolverRegions ,
                                SpatialPartition& partition ,
                                FractunatorBase* const fractunator )
{  

  int n = 9; // number of species
  rArray1d conc(2*n);
  int maxNumItrs = 300;
  realT tolx = 1e-64;
  
  // initial guess
  realT Hi = 1e-7;
  realT OHi = 1e-7;
  realT Cai = 0.9*CaBrine;
  realT CaHCO3i = 0.09*CaBrine;
  realT CaCO3i = 0.01*CaBrine;
  realT CO3i = std::max(0.5*(CBrine-CaCO3i),1e-16);
  realT HCO3i = std::max(0.4*(CBrine-CaCO3i),1e-16);
  realT H2CO3i  = std::max(0.05*(CBrine-CaCO3i),1e-16);
  realT CO2i = std::max(0.05*(CBrine-CaCO3i),1e-16);
  
  realT xi[9] = {Hi,OHi,Cai,CaHCO3i,CaCO3i,CO3i,HCO3i,H2CO3i,CO2i};
  
  for(int i =0; i < n; ++i) conc(i) = xi[i];
  for(int i =0; i < n; ++i) conc(i+n) = log(std::max(conc(i),1e-64));
  
  // solve
  	std::cout << " Newton Raphson started " << std::endl;
  NewtonRaphson( m_equilibriumSystem, conc, maxNumItrs, tolx);
  
  // report result
  for(int i =0; i < 9; ++i){
  	std::cout << i << " " << conc(i) << std::endl;
  }
  
  exit(0);
  return dt;
}


REGISTER_SOLVER( FindEquilibrium )
