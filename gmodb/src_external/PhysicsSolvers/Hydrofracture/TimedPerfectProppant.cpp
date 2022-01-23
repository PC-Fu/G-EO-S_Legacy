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
//  LLNL-CODE-618232
//  GPAC, Version 2.0
//
//  All rights reserved.
//
//  This file is part of GPAC. For details, please contact Scott Johnson or Randolph Settgast. Please also read "Additional BSD Notice" below.
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
 * @file TimedPerfectProppant.h
 * @author fu4
 * @date September 24, 2014
 * This solver will prop fracture open at pre-specified time by adding an offset (=aperture size at the time) to the contact model.
 */

#include "TimedPerfectProppant.h"

#include "PhysicsSolvers/SolverFactory.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"



TimedPerfectProppant::TimedPerfectProppant( const std::string& name,
                                                  ProblemManagerT* const pm ):
    SolverBase(name,pm)
{}

TimedPerfectProppant::~TimedPerfectProppant()
{
}

void TimedPerfectProppant::ReadXML( TICPP::HierarchicalDataNode* const hdn  )
{
  SolverBase::ReadXML( hdn );

  m_fractureSetNames = hdn->GetStringVector("fractureSetNames");
  m_proppingTimes = hdn->GetAttributeVector<realT>("proppingTimes");
  m_efficiencyCoefficients = hdn->GetAttributeVector<realT>("efficiencyCoefficients");

  if (m_fractureSetNames.size() != m_proppingTimes.size() || m_fractureSetNames.size() == 0)
  {
    throw GPException("Error! fractureSetNames and proppingTimes must be provided and must be of the same size.");
  }

  if (m_efficiencyCoefficients.size() == 0)
  {
    m_efficiencyCoefficients.resize(m_fractureSetNames.size());
    m_efficiencyCoefficients = 1.0;
  }
  else if (m_efficiencyCoefficients.size() != m_fractureSetNames.size())
  {
    throw GPException("Error! efficiencyCoefficients must be of the same size as fractureSetNames.");
  }


}


void TimedPerfectProppant::RegisterFields( PhysicalDomainT& domain )
{
  domain.m_feFaceManager.AddKeylessDataField<realT>("proppedWidth",true, true);
}


void TimedPerfectProppant::Initialize( PhysicalDomainT& domain, SpatialPartition& partition )
{}





void TimedPerfectProppant::InitializeCommunications( PartitionBase& partition )
{

}




double TimedPerfectProppant::TimeStep( const realT& time,
                                        const realT& dt,
                                        const int cycleNumber,
                                        PhysicalDomainT& domain,
                                        const sArray1d& namesOfSolverRegions ,
                                        SpatialPartition& partition,
                                        FractunatorBase* const fractunator )
{

  FaceManagerT& faceManager = domain.m_feFaceManager;
  rArray1d& proppedWidth = faceManager.GetFieldData<realT>("proppedWidth");
  const iArray1d& flowFaceType = faceManager.GetFieldData<int>("flowFaceType");

  for( localIndex i = 0; i < m_fractureSetNames.size(); ++i )
  {
    if ( time < m_proppingTimes[i] && time + dt >= m_proppingTimes[i])
    {
      std::map< std::string, lSet >::const_iterator setMap = faceManager.m_Sets.find( m_fractureSetNames[i] );
      if( setMap != faceManager.m_Sets.end() )
      {
        for( lSet::const_iterator kf=setMap->second.begin() ; kf!=setMap->second.end() ; ++kf )
        {
          if (flowFaceType[*kf] == 0)
          {
            //We calculate the current gap
            R1Tensor gap;
            R1Tensor N;

            localIndex numChildren = domain.m_feFaceManager.m_childIndices[*kf].size();
            if (numChildren <= 1)
            {
              N = faceManager.FaceNormal( domain.m_feNodeManager, *kf );
            }
            else
            {
              N = faceManager.FaceNormal( domain.m_feNodeManager, domain.m_feFaceManager.m_childIndices[*kf][0] );
            }

            gap = faceManager.CalculateGapVector( domain.m_feNodeManager, *kf );

            realT width = Dot(gap,N);

            if (width > 0)
            {
              proppedWidth[*kf] = width * m_efficiencyCoefficients[i];
            }
          }
        }
      }
      else
      {
        throw GPException("Error! Cannot find a face set where to apply the TimedPerfectProppant solver.");
      }
    }
  }
  m_stabledt.m_maxdt = 10000;
  return dt;
}





/// Register solver in the solver factory
REGISTER_SOLVER( TimedPerfectProppant )
