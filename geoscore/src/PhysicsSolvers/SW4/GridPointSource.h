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
//  LLNL-CODE-656616
//  GEOS-CORE, Version 1.0
//
//  All rights reserved.
//
//  This file is part of GEOS-CORE. For details, please contact Scott Johnson (scott.johnson@alum.mit.edu) or Randolph Settgast (rrsettgast@gmail.com).
//
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
#ifndef GRID_POINT_SOURCE_H
#define GRID_POINT_SOURCE_H

#include <iostream>
#include "TimeDep.h"

class GridPointSource
{
   friend std::ostream& operator<<(std::ostream& output, const GridPointSource& s);
public:

  GridPointSource(double frequency, double t0,
		  int i0, int j0, int k0, int g,
		  double Fx, double Fy, double Fz,
		  timeDep tDep, int ncyc, double* pars, int npar, int* ipars, int nipar );

 ~GridPointSource();

  int m_i0,m_j0,m_k0; // grid point index
  int m_grid;

  void getFxyz( double t, double* fxyz ) const;
  void getFxyztt( double t, double* fxyz ) const;
  void getFxyz_notime( double* fxyz ) const;

  double getTimeFunc(double t) const;
  double evalTimeFunc_t(double t) const;
  double evalTimeFunc_tt(double t) const;
  double evalTimeFunc_ttt(double t) const;
  double evalTimeFunc_tttt(double t) const;

  void limitFrequency(double max_freq);

  void print_info() const;

 private:

  GridPointSource();

  void initializeTimeFunction();
  double mForces[3];
  double mFreq, mT0;

  timeDep mTimeDependence;
   double (*mTimeFunc)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_t)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_tt)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_ttt)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_om)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_omtt)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_tttt)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_tttom)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_ttomom)(double f, double t,double* par, int npar, int* ipar, int nipar );
  double (*mTimeFunc_tom)(double f, double t,double* par, int npar, int* ipar, int nipar );
   double (*mTimeFunc_omom)(double f, double t,double* par, int npar, int* ipar, int nipar );

  double* mPar;
  int* mIpar; 
  int  mNpar, mNipar;

  int mNcyc;
};

#endif
