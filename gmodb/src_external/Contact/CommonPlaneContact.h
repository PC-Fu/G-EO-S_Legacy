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
/************************************************************
 * @file CommonPlaneContact.cpp
 * @date May 24, 2011
 * @author Scott Johnson
 *
 * @brief Common plane contact class
 * Adapted from the work of Oleg Vorobiev for semi-implicit
 * contact in a Lagrangian FE code
 ************************************************************/

#ifndef COMMONPLANECONTACT_H_
#define COMMONPLANECONTACT_H_

#ifndef USE_STL_IN_CPC
#define USE_STL_IN_CPC

#ifdef USE_STL_IN_CPC
#include <vector>
#endif

#include <math.h>

class CommonPlaneContact
{
public:
static bool FaceNodeContributions (
    //---NODES---
    const int nnodes,
    const double* const xNodes,
    const double* const dxNodes,
    double* const ddxNodes,
    const double* const mNodes,

    //---FACES---
    const int nfaces,
    const int* const faceNodeCounts,
    const int* const faceNodeOffsets,
    const int* const facePairCounts,
    const int* const facePairOffsets,
    const int* const faceNodes,
    const int* const faceIDs,
    const int* const faceGhosts,
    //stiffness and damping
    const double* const kbulk,
    const double* const kshear,
    const double* const viscosity,
    //cohesion and friction
    const double* const cohesion,
    const double* const cohesionLimit,
    const double* const frictionSlope,
    const double* const residualSlope,
    const double* const criticalShearSlip,
    //dilation
    const double* const dilation,
    double* const dilationStress,
    const double* const dilationAngle,
    const double* const dilationLimit,
    //softening
    const double* const tensileSoftening,
    const double* const softening,
    double* const tensileDamage,
    const double* const tensileCutOff,
    //pore fluid
    const double* const pf,

    //---FACEPAIRS---
    const unsigned long* const face2,
    int* const active,
    //normal
    double* const aperture,
    double* const normalApproach,
    double* const normalApproachMax,
    double* const normal,
    double* const forceNormal,
    double* const forceDilation,
    //tangential
    double* const forceTangential,
    double* const slip,

    //---INTERSECTION POLYGON DEFS---
     int& npoly,
     int* const poly_facepairs,
     int* const poly_offsets,
     int* const poly_counts,
#ifdef USE_STL_IN_CPC
     std::vector<double>& xpoly,
#else
     double* const xpoly,
#endif

    //---MISC---
    const double dt,
    /*
    const double penetrationFraction = 0.5,
    const double boxTol = 1.e-2,
    const double tol = 2.,
    const double normalTol = 0.01,
    const double cosMin = 9.e-1,
    const double big = 1e100,
    const double small = 1e-10);
    */
    const double penetrationFraction = 0.5,
    const double boxTol = 0.01,
    const double tol = 2.,
    const double normalTol = 0.01,
    const double cosMin = 1e-2,
    const double big = 1e100,
    const double small = 1e-10);


static void FaceProperties(
    const int ifc,
    const double* const xNodes,
    const double* const dxNodes,
    const double* const ddxNodes,
    const int* const faceNodes,
    const int* const faceNodeCounts,
    const int* const faceNodeOffsets,
    const double dt,
    double* const xfc,
    double* const dxfc,
    double* const nx,
#ifdef USE_STL_IN_CPC
    std::vector<double>& xs,
    std::vector<double>& dxs,
#else
    double* const xs,
    double* const dxs,
#endif
    double* const xmin,
    double* const xmax,
    double& area,
    double& length
    );


static bool CheckFaces(
    const int ifc1,
    const int ifc2,
    const double* const xNodes,
    const double* const dxNodes,
    const int* const faceNodes,
    const int* const faceNodeCounts,
    const int* const faceNodeOffsets,
    const double aper,
    const double dt,
    const double* const xfc1,
    const double* const dxfc1,
    const double* const nx1,
#ifdef USE_STL_IN_CPC
    const std::vector<double>& xs1,
    const std::vector<double>& dxs1,
#else
    const double* const xs1,
    const double* const dxs1,
#endif
    const double* const xmin1,
    const double* const xmax1,
    const double area1,
    const double length1,
#ifdef USE_STL_IN_CPC
    std::vector<double>& penetr1,
#else
    double* const penetr1,
#endif
    int& num1,
    const double* const xfc2,
    const double* const dxfc2,
    const double* const nx2,
#ifdef USE_STL_IN_CPC
    std::vector<double>& xs2,
    std::vector<double>& dxs2,
#else
    const double* const xs2,
    const double* const dxs2,
#endif
    const double* const xmin2,
    const double* const xmax2,
    const double area2,
    const double length2,
#ifdef USE_STL_IN_CPC
    std::vector<double>& penetr2,
#else
    double* const penetr2,
#endif
    int& num2,
    const double boxTol,
    const double normalTol,
    const double cosMin,
    const double big,
    const double small,
    double* const xr,
    double* const xrn,
    double& prod1_cp,
    double& prod2_cp,
    double* const e1,
    double* const e2,
    double* const penetr_sum,
    double& contactArea,
#ifdef USE_STL_IN_CPC
    std::vector<double>& xpoly,
#else
    double* const xpoly,
#endif
    unsigned int& npoly,
    const double tol = 0.5
    );

static void ProjectFaceToCommonPlane(
#ifdef USE_STL_IN_CPC
    const std::vector<double>& x,
    const std::vector<double>& dist,
#else
    const double* const x,
    const double* const dist,
#endif
    const double* const e1,
    const double* const e2,
    const double* const xr,
    const double* const xrn,
    const int num_nod,
#ifdef USE_STL_IN_CPC
    std::vector<double>& a,
#else
    double* const a,
#endif
    double* const ctr,
    int& num_p,
    const bool reverse = false);

static double PolygonalAreaOfIntersection(
#ifdef USE_STL_IN_CPC
    const std::vector<double>& a1,
    const std::vector<double>& a2,
#else
    const double* const a1,
    const double* const a2,
#endif
    const unsigned int nmax1,
    const unsigned int nmax2,
    const double big = 1.e100);

static double PolygonalAreaOfIntersection2(
#ifdef USE_STL_IN_CPC
							const std::vector<double>& a,
                                                        const std::vector<double>& b,
#else
							const double* const a,
                                                        const double* const b,
#endif
                                                        const unsigned int na,
                                                        const unsigned int nb);

static double PolygonalAreaOfIntersection3(
#ifdef USE_STL_IN_CPC
    const std::vector<double>& a1,
    const std::vector<double>& a2,
#else
    const double* const a1,
    const double* const a2,
#endif
    const unsigned int nmax1,
    const unsigned int nmax2,
#ifdef USE_STL_IN_CPC
    std::vector<double>& aa,
#else
    double* const aa,
#endif
    unsigned int& nmax);

static double PolygonalArea(
#ifdef USE_STL_IN_CPC
    const std::vector<double>& aa,
#else
    const double* const aa,
#endif 
    const unsigned int na);

private:

static unsigned int CoincidentIntersection(
#ifdef USE_STL_IN_CPC
							const std::vector<double>& a1, const unsigned int nmax1,
							const std::vector<double>& a2, const unsigned int nmax2,
							std::vector<double>& aa
#else
							const double* const a1, const unsigned int nmax1,
							const double* const a2, const unsigned int nmax2,
							double* const aa
#endif
							);

static unsigned int PolygonOfIntersection(
#ifdef USE_STL_IN_CPC
						       const std::vector<double>& a1, const std::vector<int>& in1a, const std::vector<int>& in1b, const  int nmax1,
						       const std::vector<double>& a2, const std::vector<int>& in2a, const std::vector<int>& in2b, const  int nmax2,
						       const std::vector<double>& aaa, const std::vector<unsigned int>& iaaa, const unsigned int naaa, const unsigned int nlim,
						       std::vector<double>& aa
#else
						       const double* const a1, const int* const in1a, int* const in1b, const unsigned int nmax1,
						       const double* const a2, const int* const in2a, int* const in2b, const unsigned int nmax2,
						       const double* const aaa, const unsigned int* const iaaa, const unsigned int naaa, const unsigned int nlim,
						       double* const aa
#endif
						       );

static unsigned int Intersections(
#ifdef USE_STL_IN_CPC
					       const std::vector<double>& a1, std::vector<int>& in1a, std::vector<int>& in1b, const unsigned int nmax1,
					       const std::vector<double>& a2, std::vector<int>& in2a, std::vector<int>& in2b, const unsigned int nmax2,
					       std::vector<double>& aaa, std::vector<unsigned int>& iaaa
#else
					       const double* const a1, int* const in1a, int* const in1b, const unsigned int nmax1,
					       const double* const a2, int* const in2a, int* const in2b, const unsigned int nmax2,
					       double* const aaa, unsigned int* const iaaa
#endif
					       );
public:

static int LineIntersection(const double* const uv1, const double* const uv2, double* const intr, const double small = 1e-5);

static int PointInTriangle(const double* const uv, const double* const uv0);

static void SlideContact( const double* const xrn,
                   const double massAverage,
                   const double frictionSlopeContact,
                   const double velnContact,
                   const double* const veltContact,
                   const double gapContact,
                   const double penetrationFraction,
                   const double dtst,
                   double* const force,
                   const double small = 1.e-20);

static void StickyContact( const double* const xrn,
                    const double massAverage,
                    const double velnContact,
                    const double* const veltContact,
                    const double gapContact,
                    const double penetrationFraction,
                    const double dtst,
                    double* const force,
                    const double small = 1.e-20);

static void ExplicitContact( const double* const xrn,
                   const double massAverage,
                   const double contactArea,
                   const double frictionSlopeContact,
                   const double apertureContact,
                   const double kbulkContact,
                   const double kshearContact,
                   const double viscContact,
                   const double cohesionContact,
                   const double dilationStressContact,
                   const double dilationAngleContact,
                   const double residualSlopeContact,
                   const double criticalShearSlipContact,
                   const double tensileCutOffContact,
                   const double climContact,
                   double& dnContact,
                   double& dnContactMax,
                   double& dsContact,
                   const double pfContact,
                   double* const fsContact,
                   double& fnContact,
                   const double velnContact,
                   const double* const veltContact,
                   const double gapContact,
                   const double dtst,
                   double* const force,
                   const double small = 1.e-20);

static void ImplicitContact( const double* const xrn,
                      const double massAverage,
                      const double frictionSlopeContact,
                      const double apertureContact,
                      const double kbulkContact,
                      const double kshearContact,
                      const double viscContact,
                      const double cohesionContact,
                      const double dilationStressContact,
                      const double dilationAngleContact,
                      const double residualSlopeContact,
                      const double criticalShearSlipContact,

                      const double tensileSoftening,
                      const double dilationLimit,
                      const double softening,
                      double& tensileDamage,
                      double dilationContact,
                      double& fnDilation,

                      double& dnContact,
                      double& dnContactMax,
                      double& dsContact,
                      const double pfContact,
                      double* const fsContact,
                      double& fnContact,
                      const double contactArea,
                      const double velnContact,
                      const double* const veltContact,
                      const double gapContact,
                      const double dtst,
                      double* const force,
                      const double small = 1.e-20);
};
#endif
#endif
