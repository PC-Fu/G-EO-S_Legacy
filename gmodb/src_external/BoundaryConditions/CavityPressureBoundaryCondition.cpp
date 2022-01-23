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
 * @file CavityPressureBoundaryCondition.cpp
 * @author johnson346
 * @date June 24, 2014
 */

#include "CavityPressureBoundaryCondition.h"

////////////////////////////
//
/// Cavity Pressure Boundary condition
/**
 * @author johnson346
 * @brief Hydraulic pressure BC
 *
 **/

CavityPressureBoundaryCondition::CavityPressureBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm):
  RadialHydraulicPressureBoundaryCondition(hdn,pm, false), m_setNamesHydro(), m_radialDistance(-1.0), m_mass(0.0),
  m_volume(-1.0), m_K0(0.0), m_rho0(0.0)
{
  m_radialDistance = hdn->GetAttributeOrDefault<realT>("radius", -1.0);
  m_K0 = hdn->GetAttributeOrDefault<realT>("bulkModulus", 2e9);
  m_rho0 = hdn->GetAttributeOrDefault<realT>("referenceDensity", 1000);
  m_pressureCap = hdn->GetAttributeOrDefault<realT>("pressureCap", std::numeric_limits<realT>::max());

  const realT lEff = hdn->GetAttributeOrDefault<realT>("effectiveBoreholeLength", m_axialDistance);
  m_volume = 2.0 * (lEff > 0 ? lEff : 0);
  m_volume *= m_radialDistance > 0 ? (m_radialDistance * m_radialDistance) : 0;
  m_volume *= 3.14159265358979323846;

  m_volume = hdn->GetAttributeOrDefault<realT>("volume", m_volume);


  m_mass = hdn->GetAttributeOrDefault<realT>("initialMass", m_rho0 * m_volume);

  m_setNamesHydro = hdn->GetStringVector("setnamesHydro");
}


// Individual solvers should use this function rather than calculate the traction directly in the solver.
// This allows the traction bc to be developed separately from the individual solvers
R1Tensor CavityPressureBoundaryCondition::GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,  realT& time)
{
  R1Tensor traction(0);
  if(m_axialDistance <= 0 || m_volume <= 0)
    return traction;

  R1Tensor position;
  {
    domain.m_feFaceManager.FaceCenter( domain.m_feNodeManager, *fc, position );
    position -= m_origin;

    const realT distance = Dot(position, m_axis);
    if(fabs(distance) > m_axialDistance)
      return traction;

    R1Tensor n(m_axis);
    n *= distance;
    position -= n;
  }
  //We can assume the face set has been declared to describe the well bore wall
//  const realT radius = position.L2_Norm();
//  if(radius > m_radialDistance)
//    return traction;

  traction = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *fc );
  traction *= -GetPressure();

  return traction;
}


void CavityPressureBoundaryCondition::UpdateMass(const realT massLoss, const realT dt, const realT time)
{
  if (!(m_timeTableName.empty()))
  {
    rArray1d t(1);
    t[0] = time;
    const realT tableval = TableManager::Instance().LookupTable<1>(m_timeTableName, t);
    m_value = m_scale * tableval;
  }
  else if (!(m_functionName.empty()))
  {
    m_value = m_scale * (*m_function)(time);
  }
  m_time = time;
  m_mass += m_value * m_rho0 * dt;
  m_mass -= massLoss;
  //return GetMassVolumeModRho0();
}

realT CavityPressureBoundaryCondition::GetPressure() const
{
  const realT rho = m_mass / m_volume;
  const realT pressure = m_K0 * (rho/m_rho0 - 1);
  return pressure > 0 ? (pressure < m_pressureCap ? pressure : m_pressureCap) : 0;
}

REGISTER_BoundaryCondition( CavityPressureBoundaryCondition )
