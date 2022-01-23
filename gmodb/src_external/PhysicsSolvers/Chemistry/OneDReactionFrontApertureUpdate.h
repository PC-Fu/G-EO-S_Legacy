//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2013, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)     Stuart Walsh(walsh24@llnl.gov)
//  Scott Johnson (johnson346@llnl.gov)        Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//
//  LLNL-CODE-6182322
//  GPAC, Version 2.0
//
//  All rights reserved.
//  This file is part of GPAC.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file SteadyStateParallelPlateFlowSolver_TwoD.h
 * @author walsh24
 * @date June 1, 2011
 */

#ifndef REACTIONFRONTAPERTUREUPDATE_H_
#define REACTIONFRONTAPERTUREUPDATE_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"
#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include <set>

#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif


class ReactionFrontApertureUpdate : public SolverBase
{
public:
  ReactionFrontApertureUpdate( const std::string& name, ProblemManagerT* const problemManager);
  virtual ~ReactionFrontApertureUpdate();
  
  void ReadXML( TICPP::HierarchicalDataNode* const hdn) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);
  void InitializeCommunications( PartitionBase & partition  );

  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 
  /// name of solver class
  /// NB: Parent class's Name() function returns name of specific solver objects 
  static const char* SolverName(){return "ReactionFrontApertureUpdate";};

private:
  const lSet* m_faceSet;
  std::string m_faceSetName;
  
  realT Ku,Kam,Kd; // elastic moduli
  realT Kamp,Kdp; // plastic moduli
  realT sigmaYam, sigmaYd; // yield stresses
  realT pillarRadius;
  realT rr; // pillarRadius**2
  realT bmin; // min value of fracture aperture

  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;

protected:

  realT CalculateNewAperture(realT bu, realT b,realT sigmaEff,realT amSiExtent, realT depletedExtent);

  void UpdateSpringSliderEtaEtaP(realT K, realT Kp, realT sigma,realT sigmay,realT& eta,realT& etaP){
	  realT dSigma = sigma - Kp*etaP- K*eta;
	  if(sigma - Kp*etaP <= sigmay){
	    eta = (sigma - Kp*etaP)/K;
	  } else {
	    if( K*eta < sigmay ){
	      eta = sigmay/K;
	    }
	    dSigma = sigma - sigmay - Kp*etaP;
	    etaP += dSigma/Kp;
	    eta += (1.0/K + 1.0/Kp)*dSigma;
	  }
	};

  realT FindInitialEtaP(realT K,realT sigma,realT sigmay,realT eta){
	  realT etaP = 0;
    if (sigma > sigmay){
      etaP = eta - sigma/K;
    }
    return etaP;
  };

  realT FindSpringSliderK(realT K,realT Kp,realT sigma,realT sigmay){
	  realT KK = K;
    if(sigma > sigmay){
      KK = 1.0/(1.0/K + 1.0/Kp);
    }
    return KK;
  };

  realT FindSpringSliderYieldStress(realT K,realT Kp,realT sigmay,realT eta){
	realT sigma = K*eta;
    if(sigma > sigmay){
      realT  KK = 1.0/(1.0/K + 1.0/Kp);
      sigma = KK*(eta-sigmay/K) + sigmay;
    }
    return sigma;
   };


  void CalculateRelativeAreas(realT amSiExtent, realT depletedExtent,
		                       realT& aratio,realT& dratio, realT& uratio){

	if(amSiExtent > pillarRadius) amSiExtent = pillarRadius;
    if(depletedExtent > pillarRadius) depletedExtent = pillarRadius;

    realT aa = rr-std::pow(pillarRadius-amSiExtent,2);
    realT da = rr-aa-std::pow(pillarRadius-depletedExtent,2);

	aratio = aa/rr;
	dratio = da/rr;
	uratio = 1.0 - aratio - dratio;
   };

};


#endif /* ReactionFrontApertureUpdate_H_ */
