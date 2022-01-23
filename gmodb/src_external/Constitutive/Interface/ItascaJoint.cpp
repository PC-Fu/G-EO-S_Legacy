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


/*
 * ItascaJoint.cpp
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 

#include "Utilities/GeometryUtilities.h"
#include "Utilities/FindRoots.h"
#include "Utilities/MaterialUtilities.h"
#include "IO/ticpp/HierarchicalDataNode.h"

#include "ItascaJoint.h"
#include "IO/ticpp/HierarchicalDataNode.h"
#include "Constitutive/Interface/InterfaceFactory.h"
#include <typeinfo>
#include <assert.h>

ItascaJoint::ItascaJoint( ):
JointIntermediate( sizeof(ParameterClass), sizeof(StateClass) )
{
  // TODO Auto-generated constructor stub
}

ItascaJoint::~ItascaJoint()
{

}

void
ItascaJoint::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                             InterfaceBaseStateData& matStateBase) const
{
  //get temporary references
  const ItascaJointParameterData& matParams = static_cast < const ItascaJointParameterData& > ( matParamsBase );
  ItascaJointStateData& matState = static_cast < ItascaJointStateData& > ( matStateBase );

  if(matState.ifail == 1)
  {
    //once failed, always failed
    matState.mu = matParams.muResidual;
  }
  else
  {
    const realT xsmag = fabs(matState.xs);
    if(xsmag < matParams.xsResidual)
    {
      //not yet failed ...
      const realT slope = xsmag / (matParams.xsResidual-xsmag);
      const realT phi1 = (1.0 - slope) * (matState.mu - matParams.muResidual);
      matState.mu = matParams.muResidual + ( (phi1 < 0.0) ? 0.0 : phi1 );
    }
    else
    {
      //just failed ...
      matState.ifail = 1;
      matState.mu = matParams.muResidual;
    }
  }
}

realT
ItascaJoint::DilationalStressIncrement( const InterfaceBaseParameterData& matParamsBase,
                                       InterfaceBaseStateData& matStateBase,
                                       const realT dxs) const
{
  ItascaJointStateData& matState = static_cast<ItascaJointStateData&> (matStateBase);
  const ItascaJointParameterData& matParams = static_cast<const ItascaJointParameterData&> (matParamsBase);

  //check whether dilation is set
  if(matParams.dilationCoefficient0 <= 0)
    return 0.0;

  const realT xsmag = fabs(matState.xs);

  // Check if we have reached relative tangential displacement limit
  if (xsmag < matParams.xsDilationLimit &&
      matState.normalApproach > matParams.normalApproachDilationInit)
  {
    // We have not reached relative tangential displacement limit
    // Dilate using tangentOfDilationAngleAtZeroStress
    // Soften tangentOfDilationAngleAtZeroStress with increasing relative tangential displacement
    realT stressSoftening = matParams.xsDilationLimit > 0.0 ?
        (matParams.xsDilationLimit - xsmag) / matParams.xsDilationLimit :
        1.0;

    // Soften tangentOfDilationAngleAtZeroStress with increasing normal stress
    if (matParams.stressDilationLimit > 0.0)
    {
      stressSoftening *= (1.0 - matState.stress / matParams.stressDilationLimit);
      stressSoftening = stressSoftening < 0 ? 0 : (stressSoftening > 1.0 ? 1.0 : stressSoftening);
    }

    // Replace stiffness with an equilibrium stiffness
    realT kcurr = matParams.kCoefficientElastic * (matParams.kdilation > 0.0 ? matParams.kdilation /
        (matState.kcurrent + matParams.kdilation) : 1.0);
    kcurr *= stressSoftening * matParams.dilationCoefficient0;

    //return change in stress wrt time
    return kcurr * sqrt(xsmag) * dxs;
  }

  // Here we are either in residual or before the onset of dilation
  return 0.0;
}

/// Register class in the class factory
REGISTER_INTERFACE( ItascaJoint )
