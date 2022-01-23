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
 * @file PerforatedCasedWellboreBoundaryCondition.h
 * @author johnson346
 * @date June 24, 2014
 */

#ifndef PerforatedCasedWellboreBoundaryCondition_H_
#define PerforatedCasedWellboreBoundaryCondition_H_

#include "RadialHydraulicPressureBoundaryCondition.h"

#include "Wellbore/PerforatedCasedWellboreProblemInterface.h"
#include "Wellbore/PerforatedCasedWellboreProblemInterfaceQA.h"
#include "Wellbore/PerforatedCasedWellboreProblemInterfaceQB.h"
#include "Wellbore/PerforatedCasedWellboreProblemInterfaceQC.h"
#include "Wellbore/PerforatedCasedWellboreProblemInterfaceQD.h"


#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"

class PerforatedCasedWellboreBoundaryCondition: public RadialHydraulicPressureBoundaryCondition
{
public:
  PerforatedCasedWellboreBoundaryCondition(TICPP::HierarchicalDataNode* BoundaryConditionNode,
                                           const ProblemManagerT* const problemManager);
  virtual ~PerforatedCasedWellboreBoundaryCondition()
  {
  }

  static const char* BoundaryConditionName()
  {
    return "PerforatedCasedWellboreBoundaryCondition";
  }

  virtual const char* GetBoundaryConditionName()
  {
    return BoundaryConditionName();
  }

  virtual R1Tensor GetTractionOnFace(PhysicalDomainT& domain, const lSet::const_iterator& fc,
                                     realT& time);

  void Apply(PhysicalDomainT& domain, realT time, realT dt);

  sArray1d m_setNamesHydro;

protected:
  void SetPerforationLocations(PhysicalDomainT& domain, bool force = false);

  void NonlinearSolve(const rArray1d& fracturePressures, const rArray1d& facePerms,  realT alpha,
                      const rArray1d& beta, realT qwellbore,
                      const bool printFlag = false);

  int Fill(const rArray1d& fracturePressures,const rArray1d& facePerms, rArray1d& betas, rArray1d& xs, rArray1d& ps, rArray1d& fp,
           iArray1d& ns);

  int ExtractFaceData( const FaceManagerT& fm, rArray1d& pressures, rArray1d& fperms );

  void Solve(const realT time, const rArray1d& fracturePressures, const rArray1d& facePerms,rArray1d& fractureFluxes);
//  void LinSolve(const realT time, const rArray1d& fracturePressures, rArray1d& fractureFluxes);

#if GPAC_MPI
  void FracturePressuresToRoot(const FaceManagerT& fm, rArray1d& pressures, rArray1d& facePerms );
  void FractureFluxesToAll(FaceManagerT& fm, realT dt, const rArray1d& fluxes);
#endif

  realT m_radialDistance, m_mass, m_volume, m_K0, m_rho0, m_pressureCap, m_pipeFrictionFactor,
      m_perforationAlpha, m_perforationSnapDistance;
  rArray1d m_perforationAlphas;
  bool m_locationsSet;
  Array1dT<ValueGUID> m_perforationLocations;
  rArray1d m_initialGuess;
#if GPAC_MPI
  std::map<globalIndex, localIndex> m_globalToNeighborRank, m_globalToNeighborOrder,
      m_globalToLocal;
  iArray1d m_sizes;
#endif
};
#endif
