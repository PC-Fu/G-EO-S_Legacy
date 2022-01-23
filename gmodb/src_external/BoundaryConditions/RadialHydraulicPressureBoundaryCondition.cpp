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
 * @file BoundaryConditions.cpp
 * @author walsh24
 * @date December 5, 2011
 */

#include "RadialHydraulicPressureBoundaryCondition.h"

////////////////////////////
//
/// Radial Hydraulic Pressure Boundary condition
/**
 * @author johnson346
 * @brief Hydraulic pressure BC
 *
 **/

RadialHydraulicPressureBoundaryCondition::RadialHydraulicPressureBoundaryCondition( TICPP::HierarchicalDataNode* hdn, const ProblemManagerT* const pm, bool t2Flag):
  TractionBoundaryCondition(hdn,pm), m_useNormalFlag(true),
  m_origin(0), m_axis(0), m_axialDistance(-1.0)
{
  m_useNormalFlag = true;
  m_isConstantInSpace = false;
  m_isConstantInTime = false;

  m_origin = hdn->GetAttributeTensor("origin");
  m_axis = hdn->GetAttributeTensor("axis");
  m_axialDistance = 0.5 * hdn->GetAttributeOrDefault<realT>("axialLength", -1.0);
}

// Individual solvers should use this function rather than calculate the traction directly in the solver.
// This allows the traction bc to be developed separately from the individual solvers
R1Tensor RadialHydraulicPressureBoundaryCondition::GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,  realT& time)
{
  R1Tensor traction(0);
  if(m_axialDistance <= 0)
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
  const realT radius = position.L2_Norm();

  traction = domain.m_feFaceManager.FaceNormal( domain.m_feNodeManager, *fc );
  m_value = 0;
  if (!(m_timeTableName.empty()))
  {
    rArray1d t(2);
    t[0] = time;
    t[1] = radius;
    const realT tableval = -TableManager::Instance().LookupTable<2>(m_timeTableName, t);
    m_value = m_scale * tableval;
  }
  else if (!(m_functionName.empty()))
  {
    const realT xx[] = {time, radius};
    m_value = m_scale * (*m_function)(*((realT*)xx));
  }
  traction *= -m_value;
  return traction;
}

REGISTER_BoundaryCondition( RadialHydraulicPressureBoundaryCondition )
