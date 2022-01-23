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
 * @file MicroseismicSourceSolver.cpp
 * @author Scott Johnson
 * @date created on June 22, 2012
 */

#include "MicroseismicSourceSolver.h"
#include "SurfaceGeneration/FractalVolume.h"
#include "PhysicsSolvers/SolverFactory.h"
#include "Utilities/Utilities.h"
#include "Utilities/Kinematics.h"
#include "Utilities/StringUtilities.h"

#include "Common/Common.h"

// Boundary Conditions
#include "BoundaryConditions/BoundaryConditions.h"
#include "BoundaryConditions/ApplyBoundaryConditions.h"

#include "ObjectManagers/TableManager.h"
#include <algorithm>
#include <math.h>

#if GPAC_MPI
#include <mpi.h>
#endif

#include "SurfaceGeneration/JointPopulator.h"

using namespace BoundaryConditionFunctions;

MicroseismicSourceSolver::MicroseismicSourceSolver(  const std::string& name,
                                                     ProblemManagerT* const pm ) :
    SolverBase(name, pm),
    m_cohesionDistribution(),
    m_frictionDistribution(),
    m_ADistribution(), m_ABRatioDistribution(), m_DcDistribution(), m_vstarDistribution(), m_thetaDistribution(),
    m_jointSets(),
    m_timePrevious(0.0), m_M0MtotalRatio(1.0), m_criticalStateLimit(0.0), m_tocTable(0)
{
}

MicroseismicSourceSolver::~MicroseismicSourceSolver()
{
}

void MicroseismicSourceSolver::ReadXML( TICPP::HierarchicalDataNode* const hdn )
{
  SolverBase::ReadXML(hdn);

  m_M0MtotalRatio = hdn->GetAttributeOrDefault<realT>("M0_Mtotal_ratio", 1e-6);
  m_criticalStateLimit = hdn->GetAttributeOrDefault<realT>("criticalStateLimit", 0.0);

  //Look for Total Organic Solids percentage table
  {
    const std::string tname = hdn->GetAttributeString("tocTableName");
    if(!tname.empty())
    {
      const std::map<std::string, Table<3, realT> >::iterator it = TableManager::Instance().Tables<3>().find(tname);
      if(it !=  TableManager::Instance().Tables<3>().end())
      {
        m_tocTable = &it->second;
      }
    }
  }

  //Cohesion distribution
  {
    realT tmp1, tmp2;

    tmp1 = hdn->GetAttributeOrDefault<realT>("meanCohesion", 1e6);
    //m_failureDistribution.AddParameter(StatisticalDistributionBaseT::MEAN, tmp1);

    tmp2 = hdn->GetAttributeOrDefault<realT>("stdevCohesion", 1e4);
    //m_failureDistribution.AddParameter(StatisticalDistributionBaseT::STANDARD_DEVIATION, tmp2);
    m_cohesionDistribution.AddWeibullParameters(tmp1, tmp2);

    tmp1 = hdn->GetAttributeOrDefault<realT>("minCohesion", 1e3);
    m_cohesionDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp1);

    tmp2 = hdn->GetAttributeOrDefault<realT>("maxCohesion", 1e7);
    m_cohesionDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp2);
  }

  //Friction distribution
  {
    realT tmp;

    tmp = hdn->GetAttributeOrDefault<realT>("meanFrictionCoefficient", 0.6);
    m_frictionDistribution.AddParameter(StatisticalDistributionBaseT::MEAN, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("stdevFrictionCoefficient", 0.1);
    m_frictionDistribution.AddParameter(StatisticalDistributionBaseT::STANDARD_DEVIATION, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minFrictionCoefficient", 0.3);
    m_frictionDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxFrictionCoefficient", 0.9);
    m_frictionDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);
  }

  //Rate and state
  {
    realT tmp;

    tmp = hdn->GetAttributeOrDefault<realT>("minA", 0.01);
    m_ADistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxA", 0.01);
    m_ADistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minABRatio", 0.67);
    m_ABRatioDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxABRatio", 0.67);
    m_ABRatioDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minDc", 1e-5);
    m_DcDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxDc", 1e-5);
    m_DcDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minVstar", 1e-6);
    m_vstarDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxVstar", 1e-6);
    m_vstarDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("minTheta", 1e8);
    m_thetaDistribution.AddParameter(StatisticalDistributionBaseT::MINIMUM_VALUE, tmp);

    tmp = hdn->GetAttributeOrDefault<realT>("maxTheta", 1e8);
    m_thetaDistribution.AddParameter(StatisticalDistributionBaseT::MAXIMUM_VALUE, tmp);
  }

  //Joint set
  m_jointSets.ReadXML(hdn->Next(true));
  if(m_jointSets.m_elementRegionName.length() == 0)
    throw GPException("You must associate the joint set with an element region!");

  //set WPP transform
  {
    R1Tensor n = m_jointSets.North(), u = m_jointSets.Up();
    GeometryUtilities::TransformWPPFrame(n, u, m_Twpp);
  }
}

void MicroseismicSourceSolver::WriteSilo(SiloFile& siloFile) const
{
  SolverBase::WriteSilo( siloFile );
  siloFile.DBWriteWrapper("timePrevious", m_timePrevious);
}

void MicroseismicSourceSolver::ReadSilo(const SiloFile& siloFile)
{
  SolverBase::ReadSilo( siloFile );
  siloFile.DBReadWrapper("timePrevious", m_timePrevious);
}

void MicroseismicSourceSolver::RegisterFields(PhysicalDomainT& domain)
{
  domain.m_microseismicSourceNodes.AddKeyedDataField<FieldInfo::displacement>();
  domain.m_microseismicSourceNodes.AddKeyedDataField<FieldInfo::referencePosition>();

  domain.m_microseismicSourceNodes.AddKeylessDataField<int>("failed", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("timeFailure", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("magnitude", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<R2Tensor>("momentTensor", false, false);

  //element mechanics
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("shearStress", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("cohesion", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("frictionCoefficient", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<localIndex>("elementIndex", true, true);
  //rate-and-state specific
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("A", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("B", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("Dc", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("vstar", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("theta", true, true);


  //element geometry
  domain.m_microseismicSourceNodes.AddKeylessDataField<R1Tensor>("normal", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<R1Tensor>("strike", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<R1Tensor>("dip", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("strikeDimension", true, true);
  domain.m_microseismicSourceNodes.AddKeylessDataField<realT>("area", true, true);
}

void MicroseismicSourceSolver::Initialize(PhysicalDomainT& domain, SpatialPartition& partition)
{

  Array1dT<R1Tensor>& positions = domain.m_microseismicSourceNodes.GetFieldData<
      FieldInfo::referencePosition>();
  Array1dT<R1Tensor>& normals = domain.m_microseismicSourceNodes.GetFieldData<R1Tensor>("normal");
  Array1dT<R1Tensor>& strikes = domain.m_microseismicSourceNodes.GetFieldData<R1Tensor>("strike");
  Array1dT<R1Tensor>& dips = domain.m_microseismicSourceNodes.GetFieldData<R1Tensor>("dip");
  rArray1d& areas = domain.m_microseismicSourceNodes.GetFieldData<realT>("area");
  rArray1d& strikeDims = domain.m_microseismicSourceNodes.GetFieldData<realT>("strikeDimension");
  rArray1d& cohesion = domain.m_microseismicSourceNodes.GetFieldData<realT>("cohesion");
  rArray1d& mu = domain.m_microseismicSourceNodes.GetFieldData<realT>("frictionCoefficient");

  rArray1d& A = domain.m_microseismicSourceNodes.GetFieldData<realT>("A");
  rArray1d& B = domain.m_microseismicSourceNodes.GetFieldData<realT>("B");
  rArray1d& Dc = domain.m_microseismicSourceNodes.GetFieldData<realT>("Dc");
  rArray1d& vstar = domain.m_microseismicSourceNodes.GetFieldData<realT>("vstar");
  rArray1d& theta = domain.m_microseismicSourceNodes.GetFieldData<realT>("theta");
  rArray1d& magnitude = domain.m_microseismicSourceNodes.GetFieldData<realT>("magnitude");

  lArray1d& elements = domain.m_microseismicSourceNodes.GetFieldData<localIndex>("elementIndex");

  ElementRegionT* elementRegion = stlMapLookupPointer(domain.m_feElementManager.m_ElementRegions,
                                                      m_jointSets.m_elementRegionName);
  int size = 1;
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if(!elementRegion)
  {
    if(rank == 0)
      throw GPException("MicroseismicSourceSolver::Initialize -> Could not find element region: " +  m_jointSets.m_elementRegionName);
    else
      return;
  }
  const iArray1d& elementGhost = elementRegion->GetFieldData<FieldInfo::ghostRank>();

  const unsigned int seed = rank + 54732645;

  const Array1dT<R1Tensor>& nodesRef = domain.m_feNodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& nodesDisp = domain.m_feNodeManager.GetFieldData<FieldInfo::displacement>();

  //handle the case of read in joint sets
  if(!m_jointSets.m_fileName.empty())
  {
    StatisticalDistributionBaseT::InitializeRandom(seed);

    Array1dT<R1Tensor> tpositions, tnormals, tstrikes, tdips;
    const bool success = m_jointSets.Populate(tpositions, tnormals, tstrikes, tdips);
    if(!success)
      throw GPException("Cannot read the JointPopulator file in its m_fileName field!");

    //remove any out-of-bound joints
    {
      Array1dT<R1Tensor>::iterator itp = tpositions.begin(),
          itn = tnormals.begin(),
          its = tstrikes.begin(),
          itd = tdips.begin();

      R1Tensor min, max;
      partition.getSizes(min, max);
      while (itp != tpositions.end())
      {
        R1Tensor& pos = *itp;
        if (pos(0) > max(0) || pos(1) > max(1) || pos(2) > max(2) ||
            pos(0) < min(0) || pos(1) < min(1) || pos(2) < min(2))
        {
          itp = tpositions.erase(itp);
          itn = tnormals.erase(itn);
          its = tstrikes.erase(its);
          itd = tdips.erase(itd);
        }
        else
        {
          ++itp;
          ++itn;
          ++its;
          ++itd;
        }
      }
    }

    //add remaining to the elements
    {
      localIndex icurr = 0;
      Array1dT<R1Tensor>::iterator itp = tpositions.begin(),
          itn = tnormals.begin(),
          its = tstrikes.begin(),
          itd = tdips.begin();
      domain.m_microseismicSourceNodes.resize(tpositions.size(), true);
      while (itp != tpositions.end())
      {
        positions[icurr] = *itp;
        normals[icurr] = *itn;
        const realT sdim = (*its).L2_Norm();
        strikeDims[icurr] = sdim;
        areas[icurr] = sdim * (*itd).L2_Norm();

        //Cohesion
        //TODO: change to correlative Weibull relationship
        cohesion[icurr] = m_cohesionDistribution.WeibullSample();

        A[icurr] = m_ADistribution.UniformSample();
        vstar[icurr] = m_vstarDistribution.UniformSample();
        theta[icurr] = m_thetaDistribution.UniformSample();

        //If a wt fraction of TOC + Clay exists, then use that table
        if(m_tocTable)
        {
          const realT tocPlusClayWt = m_tocTable->Lookup(*itp,  TableInterpolation::linear);
          TOCToRSParameters(tocPlusClayWt, A[icurr], mu[icurr], B[icurr], Dc[icurr]);
        }
        else
        {
          mu[icurr] = m_frictionDistribution.NormalSample();
          B[icurr] = A[icurr] / m_ABRatioDistribution.UniformSample();
          Dc[icurr] = m_DcDistribution.UniformSample();
        }

        //Element index
        elements[icurr] = std::numeric_limits<localIndex>::max();
        for (localIndex ii = 0; ii < elementRegion->m_localToGlobalMap.size(); ++ii)
        {
          if(elementGhost[ii] < 0 && GeometryUtilities::PointInPolyhedron(*itp, ii,
                                                 elementRegion->m_toFacesRelation,
                                                 elementRegion->m_toNodesRelation,
                                                 domain.m_feFaceManager.m_toNodesRelation,
                                                 *domain.m_feNodeManager.m_refposition,
                                                 *domain.m_feNodeManager.m_displacement))
          {
            elements[icurr] = ii;
            break;
          }
        }

        ++itp;
        ++itn;
        ++its;
        ++itd;
        ++icurr;
      }
    }
    return;
  }


  ///////////////////////////////////////////////////////////////////////////////////
  // USE STATISTICAL SPATIAL DISTRIBUTION IF NO FILE
  ///////////////////////////////////////////////////////////////////////////////////

  //populate the joints
  Array1dT<rArray1d> frequencies;
  {
    R1Tensor min, max;
    partition.getSizes(min, max);
    //std::cout << "min: " << min(0) << " " << min(1) << " " << min(2) << std::endl;
    //std::cout << "max: " << max(0) << " " << max(1) << " " << max(2) << std::endl;
    StatisticalDistributionBaseT::InitializeRandom(seed+2319);
    m_jointSets.Populate(*elementRegion, nodesRef, nodesDisp, min, max, frequencies);
    //note: frequencies will have the same spatial distribution across processes
    //if they all share the same random seed
  }

  //Set the hypocenters, normals, etc in my local domain -- now, everyone has the frequencies
  //const unsigned sret = StatisticalDistributionBaseT::InitializeRandom(seed);
  StatisticalDistributionBaseT::InitializeRandom(seed);

  for (localIndex ii = 0; ii < elementRegion->m_localToGlobalMap.size(); ++ii)
  {

    //Get the element centroid -> you only want one process to deal with any particular element
    R1Tensor centroid;
    const realT volume = CentroidAndVolumeOfElement(ii, domain.m_feNodeManager,
                                                    domain.m_feFaceManager, *elementRegion,
                                                    centroid);

    //NOTE: the following will have different realizations based solely on the spatial
    //discretization; we could guarantee the same results for any partitioning (given the same seed)
    //but this would require some thought to make sure the same number of offsets in the random
    //number generator is handled
    if(elementGhost[ii] < 0)//partition.IsCoordInPartition(centroid))
    {
      for (Array1dT<rArray1d>::iterator itf = frequencies.begin(); itf != frequencies.end(); ++itf)
      {
        const realT vf = (*itf)[ii] * volume;
        const int npoints = StatisticalDistributionBaseT::FrequencyToNumber(vf);
        if(npoints <= 0)
          continue;

//#if GPAC_MPI
//        for(int irank = 0; irank < size; irank++)
//        {
//          MPI_Barrier(MPI_COMM_WORLD);
//          if(irank == rank && npoints > 0)
//          {
//            std::cout << "xf " << rank << " " << centroid(0) << " " << centroid(1) << " " << centroid(2) << " " << npoints << " " << vf << std::endl;
//            //std::cout << "xf " << centroid(0) << " " << centroid(1) << " " << centroid(2) << " " << vf << std::endl;
//          }
//        }
//#endif

        localIndex icurr = domain.m_microseismicSourceNodes.DataLengths();
        const localIndex ncurr = icurr + npoints;
        domain.m_microseismicSourceNodes.resize(ncurr, true);
        for (; icurr < ncurr; ++icurr)
        {
          //Position
          RandomPointInElement(ii, domain.m_feNodeManager, domain.m_feFaceManager, *elementRegion,
                               positions[icurr]);
          //Normal, strike, dip
          m_jointSets.Next(strikes[icurr], dips[icurr], normals[icurr], strikeDims[icurr],
                           areas[icurr]);
          areas[icurr] *= strikeDims[icurr];
          //Cohesion
          //TODO: change to correlative Weibull relationship
          cohesion[icurr] = m_cohesionDistribution.WeibullSample();

          A[icurr] = m_ADistribution.UniformSample();
          vstar[icurr] = m_vstarDistribution.UniformSample();
          theta[icurr] = m_thetaDistribution.UniformSample();

          //If a wt fraction of TOC + Clay exists, then use that table
          if(m_tocTable)
          {
            const realT tocPlusClayWt = m_tocTable->Lookup(positions[icurr],  TableInterpolation::linear);
            TOCToRSParameters(tocPlusClayWt, A[icurr], mu[icurr], B[icurr], Dc[icurr]);
          }
          else
          {
            mu[icurr] = m_frictionDistribution.NormalSample();
            B[icurr] = A[icurr] / m_ABRatioDistribution.UniformSample();
            Dc[icurr] = m_DcDistribution.UniformSample();
          }

          //Element index
          elements[icurr] = ii;
        }
      }//for each joint set
    }//only do this for the elements in the partition
  }//for each element in the region

  magnitude = -20;
}

void MicroseismicSourceSolver::TOCToRSParameters(realT tocPlusClayWt,
                                                 realT A, realT& mu, realT& B, realT& Dc)
{
  //following derived from Kohli and Zoback (2013) "Frictional Properties of Shale Reservoir Rocks"
  //values are all non-dimensional
  //note: Figure 5 indicates that the mu is relatively invariant with confinement
  //note: only applies in TOC [0.092 : 0.52], so thresheld to that range
  const realT toc1 = tocPlusClayWt < 0.092 ? 0.092 : (tocPlusClayWt > 0.52 ? 0.52 : tocPlusClayWt);
  const realT toc2 = toc1 * toc1;

  //FRICTION COEFFICIENT RELATIONSHIP
  //(1) R^2 = 0.98862
  mu = 0.7916*toc2 - 1.3725*toc1 + 0.8934;

  //A-B RELATIONSHIP
  {
    //A-B (MIN) RELATIONSHIP
    //(2) R^2 = 0.99349
    const realT minVal = -2.2917*toc2*toc2 + 2.4657*toc2*toc1 - 0.8114*toc2 + 0.122*toc1 - 0.0093;
    //A-B (MAX) RELATIONSHIP
    //(2) R^2 = 0.99852
    const realT maxVal = -3.0635*toc2*toc2 + 3.554*toc2*toc1 - 1.3227*toc2 + 0.2223*toc1 - 0.0124;

    //The data looks pretty well distributed within this region, so assume uniform
    const realT AminB = StatisticalDistributionBaseT::UniformSample(minVal, maxVal);

    B = A - AminB;
  }

  //Dc RELATIONSHIP
  {
    //NOTE: this should be decreasing with shear strain! This is not captured here (note: max corresponds to strain of 1um and min to 2.5um)
    //ALSO, THE EQN ARE IN MICRONS, SO CONVERT AT END

    //Dc (MIN) RELATIONSHIP : 79138x5 - 121931x4 + 68982x3 - 17451x2 + 1886.2x - 53.048
    //(2) R^2 = 0.94692
    const realT minVal = 79138*toc2*toc2*toc1 - 121931*toc2*toc2 + 68982*toc2*toc1 - 17451*toc2 + 1886.2*toc1 - 53.048;
    //Dc (MAX) RELATIONSHIP : 13024x5 - 11702x4 - 230.37x3 + 2633x2 - 747.94x + 94.53
    //(2) R^2 = 0.94249
    const realT maxVal = 13024*toc2*toc2*toc1 - 11702*toc2*toc2 - 230.37*toc2*toc1 + 2633*toc2 - 747.94*toc1 + 94.53;
    Dc = (1e-6) * StatisticalDistributionBaseT::UniformSample(minVal, maxVal);
  }
}

void MicroseismicSourceSolver::InitializeCommunications(
    PartitionBase& partition )
{
}

R1Tensor
MicroseismicSourceSolver::OrderedPrincipalMoments(const R2Tensor& moment)
{
  //symmetrize the moment tensor ... should already be symmetric
  R1Tensor principalMoments;
  {
    R2SymTensor smom(0.0);
    for(localIndex i = 0; i < 3; i++)
    {
      smom(i,i) += 0.5 * moment(i,i);
      for(localIndex j = 0; j < 3; j++)
        smom(i,j) += 0.5 * moment(i,j);
    }
    smom.EigenVals(principalMoments.Data(), 0.0);
  }

  //order the moments in decreasing order
  R1Tensor orderedPrincipalMoments(principalMoments);
  {
    if(orderedPrincipalMoments(2) > orderedPrincipalMoments(1))
    {
      const realT t = orderedPrincipalMoments(2);
      orderedPrincipalMoments(2) = orderedPrincipalMoments(1);
      orderedPrincipalMoments(1) = t;
    }
    if(orderedPrincipalMoments(0) < orderedPrincipalMoments(1))
    {
      if(orderedPrincipalMoments(0) < orderedPrincipalMoments(2))
      {
        const realT t = orderedPrincipalMoments(2);
        orderedPrincipalMoments(2) = orderedPrincipalMoments(0);
        orderedPrincipalMoments(0) = orderedPrincipalMoments(1);
        orderedPrincipalMoments(1) = t;
      }
      else
      {
        const realT t = orderedPrincipalMoments(0);
        orderedPrincipalMoments(0) = orderedPrincipalMoments(1);
        orderedPrincipalMoments(1) = t;
      }
    }
  }
  return orderedPrincipalMoments;
}

void MicroseismicSourceSolver::HudsonParameters(const R1Tensor& orderedPrincipalMoments,
                                                R1Tensor& T_k_rho)
{
  //get the total moment
  T_k_rho(2) = orderedPrincipalMoments.L2_Norm();

  //decompose into isotropic and deviatoric components
  R1Tensor Mdeviatoric(orderedPrincipalMoments);
  R1Tensor Miso(1.0/3.0);
  const realT M = Dot(Miso,Mdeviatoric);
  Mdeviatoric -= M;

  //go through cases for least principal moment
  if(!isZero(Mdeviatoric(2)))
  {
    const realT fct0 = Mdeviatoric(2) > 0 ?
        -2.0 / Mdeviatoric(1) :
        2.0 / Mdeviatoric(0);
    T_k_rho(0) = fct0*Mdeviatoric(2);
    const realT Mbar = fct0*M;
    T_k_rho(1) = Mbar / (2.0 + fabs(Mbar));
  }
  else
  {
    T_k_rho(0) = 0.0;
    if(isZero(Mdeviatoric(0)))
    {
      T_k_rho(1) = M < 0 ? -1.0 : (isZero(M) ? 0.0 : 1.0);
    }
    else
    {
      const realT fct0 = 2.0 / Mdeviatoric(0);
      const realT Mbar = fct0*M;
      T_k_rho(1) = Mbar / (2.0 + fabs(Mbar));
    }
  }
}

void MicroseismicSourceSolver::TapeParameters(const R1Tensor& orderedPrincipalMoments,
                                              R1Tensor& gamma_delta_rho)
{
  //From Tape (2012a) "A geometric setting for moment tensors"
  const realT invsqrt3 = 1.0 / sqrt(3.0);
  const realT pi_over_2 = (2.0 * atan(1));

  //get the total moment
  gamma_delta_rho(2) = orderedPrincipalMoments.L2_Norm();

  //get tan gamma
  const realT tan_gamma = (-orderedPrincipalMoments(0) + 2*orderedPrincipalMoments(1) - orderedPrincipalMoments(2)) * invsqrt3 / (orderedPrincipalMoments(0) - orderedPrincipalMoments(2));
  const realT cos_beta = (orderedPrincipalMoments(0) + orderedPrincipalMoments(1) + orderedPrincipalMoments(2)) * invsqrt3 / gamma_delta_rho(2);

  gamma_delta_rho(0) = atan(tan_gamma);
  gamma_delta_rho(1) = pi_over_2 - acos(cos_beta);
}

double MicroseismicSourceSolver::TimeStep(
    const realT& time, const realT& dt,
    const int cycleNumber,
    PhysicalDomainT& domain,
    const sArray1d& namesOfSolverRegions,
    SpatialPartition& partition,
    FractunatorBase* const fractunator)
{
  iArray1d& failed = domain.m_microseismicSourceNodes.GetFieldData<int>("failed");
  rArray1d& timeFailure = domain.m_microseismicSourceNodes.GetFieldData<realT>("timeFailure");
  rArray1d& shearStress = domain.m_microseismicSourceNodes.GetFieldData<realT>("shearStress");
  rArray1d& magnitude = domain.m_microseismicSourceNodes.GetFieldData<realT>("magnitude");
  Array1dT<R2Tensor>& momentTensor = domain.m_microseismicSourceNodes.GetFieldData<R2Tensor>(
      "momentTensor");

  //immutable
  const rArray1d& cohesion = domain.m_microseismicSourceNodes.GetFieldData<realT>(
      "cohesion");
  const rArray1d& mu = domain.m_microseismicSourceNodes.GetFieldData<realT>(
      "frictionCoefficient");

  const rArray1d& A = domain.m_microseismicSourceNodes.GetFieldData<realT>("A");
  const rArray1d& B = domain.m_microseismicSourceNodes.GetFieldData<realT>("B");
  const rArray1d& Dc = domain.m_microseismicSourceNodes.GetFieldData<realT>("Dc");
  const rArray1d& vstar = domain.m_microseismicSourceNodes.GetFieldData<realT>("vstar");
  rArray1d& theta = domain.m_microseismicSourceNodes.GetFieldData<realT>("theta");

//  const lArray1d& elementIndex = domain.m_microseismicSourceNodes.GetFieldData<localIndex>("elementIndex");
  ElementRegionT& elementRegion = stlMapLookup(domain.m_feElementManager.m_ElementRegions,
                                               m_jointSets.m_elementRegionName);
  const rArray1d& volume = elementRegion.GetFieldData<FieldInfo::volume>();

  //element geometry
  const lArray1d& elements = domain.m_microseismicSourceNodes.GetFieldData<localIndex>("elementIndex");
  const Array1dT<R1Tensor>& normal = domain.m_microseismicSourceNodes.GetFieldData<R1Tensor>(
      "normal");
  const Array1dT<R1Tensor>& ref = domain.m_microseismicSourceNodes.GetFieldData<
      FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& disp = domain.m_microseismicSourceNodes.GetFieldData<
      FieldInfo::displacement>();

  //const rArray1d& strikeDimension = domain.m_microseismicSourceNodes.GetFieldData<realT>( "strikeDimension" );
  const rArray1d& area = domain.m_microseismicSourceNodes.GetFieldData<realT>("area");

  const R1Tensor up = m_jointSets.Up();
  const unsigned int inormal = 1, istrike = 0, idip = 2;

  //retrieve parent element's material state
  R2Tensor deviatorPlaneFrame, T, sigma, dadt, dadtPlaneFrame;
  R1Tensor orderedPrincipalMoments, params;
  for (localIndex i = 0; i < domain.m_microseismicSourceNodes.DataLengths(); i++)
  {
    if (failed[i] > 0)
      continue;

    //Rotate the deviatoric strain to the frame of reference of the discontinuity
    const realT p = elementRegion.m_mat->StateData(elements[i],0)->pressure;
    {
      const R2SymTensor& S = elementRegion.m_mat->StateData(elements[i],0)->devStress;
      const R2SymTensor& Dadt = elementRegion.m_Dadt[elements[i]][0];
      for (localIndex ii = 0; ii < nsdof; ii++)
      {
        for (localIndex jj = 0; jj < nsdof; jj++)
        {
          dadt(ii, jj) = Dadt(ii, jj); // + (ii == jj ? p : 0.0);
          sigma(ii, jj) = S(ii, jj); // + (ii == jj ? p : 0.0);
        }
      }
      //Get the transform matrix from plane frame to space frame
      GeometryUtilities::TransformPlaneFrame(normal[i], up, T);
      //Transform the stress to the plane frame
      GeometryUtilities::TransformTensorFrameTranspose(T, sigma, deviatorPlaneFrame);
      GeometryUtilities::TransformTensorFrameTranspose(T, dadt, dadtPlaneFrame);
    }

    //Determine the current friction coefficient via rate-and-state friction
    realT muCurr = mu[i];
    {
      //Get an estimate of the long-range slip rate
      const realT ddot = sqrt( dadtPlaneFrame(istrike, inormal) * dadtPlaneFrame(istrike, inormal) +
                               dadtPlaneFrame(idip, inormal) * dadtPlaneFrame(idip, inormal)) *
                               pow(volume[elements[i]], 1.0/3.0);
      //Evaluate the Dieterich-Ruina rate-and-state friction law
      muCurr += A[i] * log((ddot < vstar[i] ? vstar[i] : ddot)/vstar[i]) + B[i] * log( theta[i] * vstar[i] / Dc[i]);
      theta[i] += (1-theta[i]*ddot / Dc[i]) * dt;
    }

    //get the shear stress on the element (as the full stress magnitude orthogonal to the normal)
    const realT nstress = fabs(deviatorPlaneFrame(inormal,inormal) + p);
    const realT sfail = cohesion[i] + nstress * muCurr;

    //sstress = |sigma_ij*nj - sigma_n*ni|
    const realT sstress = sqrt( deviatorPlaneFrame(istrike, inormal) * deviatorPlaneFrame(istrike, inormal) +
                                deviatorPlaneFrame(idip, inormal) * deviatorPlaneFrame(idip, inormal));

    if (fabs(sstress) > sfail)
    {
      failed[i] = 1;
      const realT dss = fabs(sstress) - fabs(shearStress[i]);
      timeFailure[i] = m_timePrevious + (
          dss > 0 ? (time - m_timePrevious) * (fabs(sstress) - sfail) / dss : 0.0);

      //Currently, just assumes energy components orthogonal to normal are lost ... may want to revisit
      for (localIndex ii = 0; ii < nsdof; ii++)
      {
        deviatorPlaneFrame(inormal, ii) = ii == inormal ? (nstress > m_criticalStateLimit ? -1.0 : 1.0) *
            StatisticalDistributionBaseT::UniformSample(0, 0.2*nstress) : 0.0;
        deviatorPlaneFrame(ii, inormal) = ii == inormal ? deviatorPlaneFrame(ii, inormal) : 0.0;
      }
      if (domain.m_feFaceManager.m_toNodesRelation[0].size() == 2)  // This is 2D.  We create a sphere because the real tensor will be a single line segment.
      {
        deviatorPlaneFrame(istrike, istrike) = fabs(deviatorPlaneFrame(0, 0));
        deviatorPlaneFrame(inormal, inormal) = deviatorPlaneFrame(0, 0);
        deviatorPlaneFrame(idip, idip) = deviatorPlaneFrame(0, 0);
      }
      //Transform the stress back to space frame (moment tensor * (1/(strain*area))
      GeometryUtilities::TransformTensorFrame(T, deviatorPlaneFrame, momentTensor[i]);

      //note: according to Wells and Coppersmith (1994) and Wyss (1979), there is good correlation between
      //magnitude and rupture area in larger magnitude events, such that:
      //  Wells and Coppersmith: M = 4.07 + 0.98 log(A)
      //  Wyss:                  M = 4.15 + log(A)
      //
      //where M is given by Hanks and Kanamori (1979) as M = (2/3)*(log(M0) - 10.7)
      //Using the simpler form in Wyss:
      //        log(M0) = 1.5*(4.15 + log(A)) + 10.7
      //                = 16.925 + 1.5*log(A)
      //     -> M0 = (10^16.925) * A^(3/2) = 8.414e16 Pa * A^(3/2)
      // NOTE: the functional form is dependent on A^(3/2), so ...
      const realT area_strain = area[i] * sqrt(area[i]);

      //Multiply by (strain*area) to resolve the moment tensor
      //ALSO, adjust by the seismic transfer coefficient M0/Mtotal (Lapusta, 2009, JGR)
      const realT fct = fabs(sfail / sstress) * m_M0MtotalRatio;
      momentTensor[i] *= fct * area_strain;

      realT mmag = -20;
      {
        //Print out source parameters for plotting
        orderedPrincipalMoments = OrderedPrincipalMoments(momentTensor[i]);
        TapeParameters(orderedPrincipalMoments, params);
        params *= 180.0 / (4 * atan(1));

        if(A[i]-B[i] > 0)
          std::cout << "AseismicTape " << params(0) << " " << params(1) << " " << params(2);
        else
          std::cout << "Tape " << params(0) << " " << params(1) << " " << params(2);

        HudsonParameters(orderedPrincipalMoments, params);

        std::cout << " Hudson " << params(0) << " " << params(1) << " " << params(2) << " nstress " << nstress << " time " << timeFailure[i] << std::endl;

        mmag = params(2);

        //get the format expected by WPP
        //PrintRupture(i, time, ref, disp, mmag, T, momentTensor);
        PrintWPP(timeFailure[i], ref[i], disp[i], mmag, m_Twpp, momentTensor[i]);
      }
      magnitude[i] = EarthquakeSimulation::FaultRupture::Magnitude(mmag);//note: expects SI
    }
    shearStress[i] = sstress;
  }

  //Calculate next stable time step
  SetMaxStableTimeStep(domain, namesOfSolverRegions, partition);
  m_timePrevious = time;
  return dt;
}

void MicroseismicSourceSolver::SetMaxStableTimeStep(
    PhysicalDomainT& domain ,
    const sArray1d& namesOfSolverRegions ,
    SpatialPartition& partition )
{
  m_stabledt.m_maxdt = 1e30; //0.9*std::numeric_limits<double>::max();
}

realT MicroseismicSourceSolver::CentroidAndVolumeOfElement(const localIndex a,
                                                           const NodeManagerT& nodeManager,
                                                           const FaceManagerT& faceManager,
                                                           const ElementRegionT& elementRegion,
                                                           R1Tensor& centroid)
{
  const Array1dT<R1Tensor>& nodePos = nodeManager.GetFieldData<FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& nodeDisp = nodeManager.GetFieldData<FieldInfo::displacement>();

  //get the centroid
  centroid = 0.0;
  {
    for (localIndex i = 0; i < elementRegion.m_toNodesRelation.Dimension(1); i++)
    {
      centroid += nodePos[elementRegion.m_toNodesRelation(a, i)];
      centroid += nodeDisp[elementRegion.m_toNodesRelation(a, i)];
    }
    centroid *= 1.0 / elementRegion.m_toNodesRelation.Dimension(1);
  }

  // We need to figure out whether this is 2D or 3D.
  localIndex nNdPerFace;
  {
    localIndex iFace = elementRegion.m_toFacesRelation(a, 0);
    nNdPerFace = faceManager.m_toNodesRelation[iFace].size();
  }

  realT volume = 0.0;
  if (nNdPerFace > 2)
  {
    //decompose into tets to get the volume
    for (localIndex i = 0; i < elementRegion.m_toFacesRelation.Dimension(1); i++)
    {
      const localIndex faceIndex = elementRegion.m_toFacesRelation(a, i);
      const localIndex nnodes = faceManager.m_toNodesRelation[faceIndex].size();
      R1Tensor x0(nodePos[faceManager.m_toNodesRelation[faceIndex][0]]);
      x0 += nodeDisp[faceManager.m_toNodesRelation[faceIndex][0]];
      for (localIndex j = 2; j < nnodes; j++)
      {
        R1Tensor x1(nodePos[faceManager.m_toNodesRelation[faceIndex][j - 1]]);
        x1 += nodeDisp[faceManager.m_toNodesRelation[faceIndex][j - 1]];
        R1Tensor x2(nodePos[faceManager.m_toNodesRelation[faceIndex][j]]);
        x2 += nodeDisp[faceManager.m_toNodesRelation[faceIndex][j]];
        R1Tensor xx;
        volume += GeometryUtilities::CentroidAndVolume_3DTetrahedron(centroid, x0, x1, x2, xx);
      }
    }
  }
  else if (nNdPerFace == 2)
  {
    const localIndex* elemToNodeMap = elementRegion.m_toNodesRelation[a];

    lArray1d nodeList;
    for (localIndex i=0; i<elementRegion.m_toNodesRelation.Dimension(1); ++i)
    {
      nodeList.push_back(elemToNodeMap[i]);
    }

    R1Tensor junk1, junk2;

    nodeManager.SortNodeOnPlane(nodeList);

    volume = GeometryUtilities::Centroid_3DPolygon(nodeList,
                                                 *(nodeManager.m_refposition),
                                                 *(nodeManager.m_displacement),
                                                 junk1,
                                                 junk2 );
  }
  else
  {
    throw GPException("MicroseismicSourceSolver: The number of nodes per face is not right!");
  }
  return volume;
}

void MicroseismicSourceSolver::RandomPointInElement(const localIndex a,
                                                    const NodeManagerT& nodeManager,
                                                    const FaceManagerT& faceManager,
                                                    const ElementRegionT& elementRegion,
                                                    R1Tensor& point)
{
  //Get extrema of the element
  R1Tensor xmin(std::numeric_limits<realT>::max());
  R1Tensor xmax(-std::numeric_limits<realT>::max());
  R1Tensor xel(0.0), xfc(0.0), d(0.0), position(0.0), center(0.0), normal(0.0);
  {
    const Array1dT<R1Tensor>& nodePos = nodeManager.GetFieldData<FieldInfo::referencePosition>();
    const Array1dT<R1Tensor>& nodeDisp = nodeManager.GetFieldData<FieldInfo::displacement>();
    for (localIndex i = 0; i < elementRegion.m_toNodesRelation.Dimension(1); i++)
    {
      position = nodePos[elementRegion.m_toNodesRelation(a, i)];
      position += nodeDisp[elementRegion.m_toNodesRelation(a, i)];
      xel += position;
      xmin.SetMin(position);
      xmax.SetMax(position);
    }
    xel *= 1.0 / elementRegion.m_toNodesRelation.Dimension(1);
  }

  //Get a valid random point inside the convex hull of the element
  while (true)
  {
    bool ok = true;

    //get a trial point
    for (localIndex i = 0; i < nsdof; i++)
      point(i) = StatisticalDistributionBaseT::UniformSample(xmin(i), xmax(i));

    //determine whether to reject the point
    for (localIndex i = 0; i < elementRegion.m_toFacesRelation.Dimension(1); i++)
    {
      const localIndex faceIndex = elementRegion.m_toFacesRelation(a, i);
      faceManager.FaceCenterAndNormal(nodeManager, faceIndex, center, normal);
      //make sure the normal is facing out
      {
        xfc = center;
        xfc -= xel;
        if (Dot(normal, xfc) < 0)
          normal *= -1.0;
      }
      //make sure the point does not lay outside of the planar approximation of the face
      d = point;
      d -= center;
      if (Dot(normal, d) > 0)
      {
        ok = false;
        break;
      }
    }
    if (ok)
    {
      break;
    }
  }
}

/// Register solver in the solver factory
REGISTER_SOLVER( MicroseismicSourceSolver)
