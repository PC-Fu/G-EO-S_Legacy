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
 * @file PerforatedCasedWellboreProblemInterface.h
 * @author johnson346
 * @date July 4, 2014
 */

#ifndef PerforatedCasedWellboreProblemInterface_H_
#define PerforatedCasedWellboreProblemInterface_H_

#include "Common/Common.h"
#include "Utilities/Utilities.h"

#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"

#if 1
class ValueGUID
{
public:
  ValueGUID() :
      m_index(0),
      m_x(0.0),
      m_y(0.0)
  {

  }

  ValueGUID(globalIndex index, realT x, realT y) :
      m_index(index),
      m_x(x),
      m_y(y)
  {
  }

  ~ValueGUID()
  {

  }

  //int rankOfOriginatingProcess;
  globalIndex m_index;
  realT m_x;
  realT m_y;
};
#else
struct ValueGUID
{
  //int rankOfOriginatingProcess;
  globalIndex index;
  realT x;
};
#endif

class PerforatedCasedWellboreProblemInterface: public NOX::Epetra::Interface::Jacobian,
                                               public NOX::Epetra::Interface::Required
{

public:

  PerforatedCasedWellboreProblemInterface(){}

  virtual ~PerforatedCasedWellboreProblemInterface() {}

  void Initialize(const rArray1d& fracturePressures, const rArray1d& facePerms, const rArray1d& alphas,
                  const rArray1d& beta, realT qwellbore);

  virtual void InitialGuess(Epetra_Vector& x) = 0;

  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& f,
                        NOX::Epetra::Interface::Required::FillType F) = 0;

  bool constructJacobian(const Epetra_Vector& x,
                         Epetra_CrsMatrix& jacobianMatrix) const;

  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& jacobianMatrix);

  bool computePrecMatrix(const Epetra_Vector& x, Epetra_RowMatrix& ePetraMatrix);

  bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& ePetraOperator);

  virtual int Size() const
  {
    return 3 * m_fracturePressures.size() - 1;
  }

  static realT dxabsxdx(const realT x)
  {
//    if(isZero(x))
//      throw GPException("dxabsxdx: Zero x");
    return 2 * fabs(x); //x * x / fabs(x) + fabs(x);
  }

  realT m_alpha;
  rArray1d m_beta, m_fracturePressures, m_faceK;
  realT m_qwellbore;

  // SW
  rArray1d m_alphas;

protected:

  virtual void Jacobian(const localIndex ieqn, const Epetra_Vector& x, rArray1d& values) const = 0;

  // x = {q_f:i=1..N,
  //      q_w:i=1..N-1,
  //      p_w:i=1..N} <-- p_w only for "A" version
  static realT GetQf(const localIndex i, const size_t nclusters, const Epetra_Vector& x)
  {
    return x[i];
  }

  static realT GetQw(const localIndex i, const size_t nclusters, const Epetra_Vector& x)
  {
    return x[i + nclusters];
  }

  static realT GetPw(const localIndex i, const size_t nclusters, const Epetra_Vector& x)
  {
    return x[i + 2 * nclusters - 1];
  }
};

#endif
