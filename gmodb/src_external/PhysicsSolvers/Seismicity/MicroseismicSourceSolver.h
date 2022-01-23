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
 * @file MicroseismicSourceSolver.h
 * @author Scott Johnson
 * @date created on June 22, 2012
 */

#ifndef MICROSEISMICSOURCESOLVER_H_
#define MICROSEISMICSOURCESOLVER_H_

#include "PhysicsSolvers/SolverBase.h"
#include "FaultRuptureBEMSolver.h"
#include "FaultRupture.h"
#include "MicroseismicSourceSolver.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"
#include "SurfaceGeneration/JointSetT.h"

class MicroseismicSourceSolver: public SolverBase
{
public:
  MicroseismicSourceSolver(  const std::string& name,
                             ProblemManagerT* const pm );

  ~MicroseismicSourceSolver();

  double
  TimeStep(const realT& time, const realT& dt,
           const int cycleNumber,
           PhysicalDomainT& domain,
           const sArray1d& namesOfSolverRegions, SpatialPartition& partition,
           FractunatorBase* const fractunator);

  virtual void SetMaxStableTimeStep( PhysicalDomainT& domain,
                                     const sArray1d& namesOfSolverRegions,
                                     SpatialPartition& partition);

  void Initialize(PhysicalDomainT& domain, SpatialPartition& partition  );

  void InitializeCommunications( PartitionBase& partition );

  virtual void
  RegisterFields(PhysicalDomainT& domain);

  virtual void
  ReadXML( TICPP::HierarchicalDataNode* const hdn );

  static const char*
  SolverName() {return "MicroseismicSourceSolver";};

  virtual void WriteSilo( SiloFile& siloFile ) const;

  virtual void ReadSilo( const SiloFile& siloFile );

private:
  R2Tensor m_Twpp;

  StatisticalDistributionBaseT m_cohesionDistribution;
  StatisticalDistributionBaseT m_frictionDistribution;
  StatisticalDistributionBaseT m_ADistribution, m_ABRatioDistribution, m_DcDistribution, m_vstarDistribution, m_thetaDistribution;

  JointPopulator m_jointSets;

  //state variable: time
  realT m_timePrevious, m_M0MtotalRatio, m_criticalStateLimit;
  Table<3, realT>* m_tocTable;

  void RandomPointInElement(const localIndex a,
                            const NodeManagerT& nodeManager,
                            const FaceManagerT& faceManager,
                            const ElementRegionT& elementRegion,
                            R1Tensor& point);

  static R1Tensor OrderedPrincipalMoments(const R2Tensor& moment);

  static void TOCToRSParameters(realT tocPlusClayWt, realT A, realT& mu, realT& B, realT& Dc);

  static void HudsonParameters(const R1Tensor& principalMoments,
                               R1Tensor& T_k_rho);

  static void TapeParameters(const R1Tensor& principalMoments,
                               R1Tensor& gamma_delta_rho);

  realT CentroidAndVolumeOfElement(const localIndex a,
                                   const NodeManagerT& nodeManager,
                                   const FaceManagerT& faceManager,
                                   const ElementRegionT& elementRegion,
                                   R1Tensor& centroid);

  static void PrintRupture(const localIndex i,
                           const realT time,
                           const Array1dT<R1Tensor>& ref,
                           const Array1dT<R1Tensor>& disp,
                           const realT mmag,
                           const R2Tensor& T,
                           const R2Tensor& momentTensor)
  {
    std::cout << "rupture " << time;
    {
      for (localIndex ii = 0; ii < nsdof; ii++)
        std::cout << " " << (ref[i](ii) + disp[i](ii));
      std::cout << " " << mmag;
      for (localIndex ii = 0; ii < nsdof; ii++)
        for (localIndex jj = 0; jj < nsdof; jj++)
          std::cout << " " << T(ii,jj);
      for (localIndex ii = 0; ii < nsdof; ii++)
        for (localIndex jj = 0; jj < nsdof; jj++)
          std::cout << " " << momentTensor(ii,jj);
    }
    std::cout << "\n";
  }

  static void PrintWPP(const realT time,
                       const R1Tensor& ref,
                       const R1Tensor& disp,
                       const realT mmag,
                       const R2Tensor& Twpp,
                       const R2Tensor& momentTensor)
  {
    R1Tensor location(ref);
    location += disp;

    R1Tensor wpploc;
    wpploc.AijBi(Twpp, location);

    R2Tensor twpp;
    GeometryUtilities::TransformTensorFrameTranspose(Twpp, momentTensor, twpp);

    std::cout << "source type=Dirac t0=" << time
        << " x=" << wpploc(0)
        << " y=" << wpploc(1)
        << " z=" << wpploc(2)
        << " m0=" << mmag
        << " mxx=" << twpp(0,0)
        << " myy=" << twpp(1,1)
        << " mzz=" << twpp(2,2)
        << " mxy=" << twpp(0,1)
        << " mxz=" << twpp(0,2)
        << " myz=" << twpp(1,2)
        << std::endl;
  }
};

#endif /* MICROSEISMICSOURCESOLVER_H_ */
