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
#include <math.h>
#include <assert.h>
#include <cstring>
#include <iostream>
#include "CommonPlaneContact.h"
#include "Utilities/Utilities.h"

extern "C"
{
  #include "aip.h"
}
/*******************************************************
 * @file CommonPlaneContact.cpp
 * @author Scott Johnson
 * @date May 24, 2011
 * @brief Common plane contact class
 * Adapted from the work of Oleg Vorobiev for semi-implicit
 * contact in a Lagrangian framework. I/O structures are built to be generic
 * to facilitate portability between codes
 *******************************************************/

/**
 * @brief FaceNodeContributions determines nodal accelerations due to contacts
 * @author Scott Johnson
 * Iterate through face pairs to determine nodal force contributions due to contacts
 *
 * @param[in] nnodes Number of external nodes
 * @param[in] xNodes Coordinates of the nodes (length 3*nnodes)
 * @param[in] mNodes Masses of the nodes (length nnodes)
 * @param[in,out] dxNodes Velocities of the nodes (length 3*nnodes)
 * @param[out] ddxNodes Accelerations of the nodes (length 3*nnodes)
 *
 * @param[in] nfaces Number of external faces
 * @param[in] faceNodeCounts Number of nodes in each face (length nfaces)
 * @param[in] faceNodeOffsets Index offset for each face's node list in faceNodes (length nfaces)
 * @param[in] facePairCounts Number of face pairs associated with each face (length nfaces)
 * @param[in] facePairOffsets Index offset for each face pair's state variable array (length nfaces)
 * @param[in] faceNodes Nodes associated with each face (length sum(faceNodeCounts))
 * @param[in] faceIDs ID of the owning element for each face (length nfaces)
 * @param[in] faceGhosts Flag describing whether the face is a ghost (length nfaces)
 * @param[in] kbulk Normal stiffnesses (length nfaces)
 * @param[in] kshear Shear stiffnesses (length nfaces)
 * @param[in] viscosity Viscosity (length nfaces)
 * @param[in] cohesion Cohesion (length nfaces)
 * @param[in] cohesionLimit Cohesion limit (length nfaces)
 * @param[in] dilationStress Dilation stress (length nfaces)
 * @param[in] dilationAngle Dilation angle (length nfaces)
 * @param[in] dilation Dilation (length nfaces)
 * @param[in] dilationLimit Dilation limit (length nfaces)
 * @param[in] frictionSlope Tangent of friction angle (length nfaces)
 * @param[in] residualSlope (length nfaces)
 * @param[in] criticalShearSlip Critical plastic slip in shear (length nfaces)
 * @param[in] tensileSoftening Tensile softening parameter (length nfaces)
 * @param[in] softening Softening parameter (length nfaces)
 * @param[in,out] tensileDamage Tensile damage state variable (length nfaces)
 * @param[in] tensileCutOff Cutoff for the tensile damage (length nfaces)
 * @param[in] pf Pore fluid (length nfaces)
 * @param[in] penetrationFraction Penetration fraction factor (length nfaces)
 *
 * @param[in] face2 second face in the pair (length nfacepairs = sum(facePairCounts))
 * @param[in] active Flag describing whether the contact already exists (length nfacepairs)
 * @param[in,out] aperture Aperture (length nfacepairs)
 * @param[in,out] normalApproach Normal approach of the faces (length nfacepairs)
 * @param[in,out] normalApproachMax Maximum normal approach of the faces (length nfacepairs)
 * @param[in,out] normal Normal unit vector array (length 3*nfacepairs)
 * @param[in,out] forceNormal Magnitude of the normal force (length nfacepairs)
 * @param[in,out] forceDilation Magnitude of the dilation force along the normal (length nfacepairs)
 * @param[in,out] forceTangential Tangential force vector array (length 3*nfacepairs)
 * @param[in,out] slip Cumulative shear slip (length nfacepairs)
 *
 * @param[out] npoly Number of polygons of intersection
 * @param[out] poly_offsets Offsets of coordinates in the xpoly list (length npoly)
 * @param[out] poly_count Counts of coordinates in the xpoly list (length npoly)
 * @param[out] xpoly Coordinates (in tuples of x,y,z) of the points composing each polygon of intersection
 *
 * @param[in] dt Timestep
 * @param[in] boxTol Tolerance for determining axis-aligned bounding box overlap
 * @param[in] tol Tolerance for determining distance between face centers
 * @param[in] normalTol Tolerance for relative distance between face centers along the normal
 * @param[in] cosMin Tolerance for the minimum magnitude of the dot product of the faces' normals
 * @param[in] big Large number used for initializing mins and maxes
 * @param[in] small Small number used for tolerance
 */
bool CommonPlaneContact::FaceNodeContributions (
                                                //---NODES---
                                                const int nnodes ,
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
                                                const double* const dilation ,
                                                double* const dilationStress,
                                                const double* const dilationAngle,
                                                const double* const dilationLimit ,
                                                //softening
                                                const double* const tensileSoftening ,
                                                const double* const softening ,
                                                double* const tensileDamage ,
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
                                                double* const forceDilation ,
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
                                                const double penetrationFraction ,
                                                const double boxTol,
                                                const double tol,
                                                const double normalTol,
						const double cosMin,
                                                const double big,
                                                const double small
                                                )
{
//  //reset the acceleration increments to zero
//  {
//    for (int i = 0; i < 3*nnodes; ++i)
//      ddxNodes[i] = 0.;
//  }

  npoly = 0;

  //go through each face
  for(int ifc1 = 0; ifc1 < nfaces; ++ifc1)
  {
    if(facePairCounts[ifc1]==0)
      continue;
    double xfc1[] = {0., 0., 0.};
    double dxfc1[] = {0., 0., 0.};
    double xmin1[] = {big,big,big};
    double xmax1[] = {-big,-big,-big};
#ifdef USE_STL_IN_CPC
    std::vector<double> x1;    x1.resize(3*faceNodeCounts[ifc1],0);
    std::vector<double> dx1;  dx1.resize(3*faceNodeCounts[ifc1],0);
#else
    double* x1 = new double[3*faceNodeCounts[ifc1]];
    double* dx1 = new double[3*faceNodeCounts[ifc1]];
#endif
    double area1 = 0.;
    double length1 = 0.;
    double nx1[] = {0., 0., 0.};
    FaceProperties(
        ifc1,
        xNodes,
        dxNodes,
        ddxNodes,
        faceNodes,
        faceNodeCounts,
        faceNodeOffsets,
        dt,
        xfc1,
        dxfc1,
        nx1,
        x1,
        dx1,
        xmin1,
        xmax1,
        area1,
        length1
        );

    for(int ifp1 = 0; ifp1 < facePairCounts[ifc1]; ++ifp1)
    {
      int ifp = facePairOffsets[ifc1] + ifp1;
      unsigned long ifc2 = face2[ifp];

      //Do not attempt contact of two ghosts or of two faces from the same element
      if ((faceGhosts[ifc1] && faceGhosts[ifc2]) || faceIDs[ifc1]==faceIDs[ifc2])
        continue;

      double xfc2[] = {0.,0.,0.};
      double dxfc2[] = {0.,0.,0.};
      double xmin2[] = {big,big,big};
      double xmax2[] = {-big,-big,-big};
      double nx2[] = {0., 0., 0.};
      double area2 = 0.;
      double length2 = 0.;
#ifdef USE_STL_IN_CPC
      std::vector<double> x2;   x2.resize(3*faceNodeCounts[ifc2],0);
      std::vector<double> dx2; dx2.resize(3*faceNodeCounts[ifc2],0);
#else
      double* x2 =  new double[3*faceNodeCounts[ifc2]];
      double* dx2 = new double[3*faceNodeCounts[ifc2]];
#endif
      FaceProperties(
          ifc2,
          xNodes,
          dxNodes,
          ddxNodes,
          faceNodes,
          faceNodeCounts,
          faceNodeOffsets,
          dt,
          xfc2,
          dxfc2,
          nx2,
          x2,
          dx2,
          xmin2,
          xmax2,
          area2,
          length2
          );

      double xr[] = {0., 0., 0.};
      double xrn[] = {0., 0., 0.};
      double prod1_cp = 0.;
      double prod2_cp = 0.;
      double e1[] = {0., 0., 0.};
      double e2[] = {0., 0., 0.};
      double penetr_sum[] = {0., 0.};
      double contactArea = 0.;
#ifdef USE_STL_IN_CPC
      std::vector<double> penetr1; penetr1.resize(faceNodeCounts[ifc1],0);
      std::vector<double> penetr2; penetr2.resize(faceNodeCounts[ifc2],0);
#else
      double* penetr1 = new double[faceNodeCounts[ifc1]];
      double* penetr2 = new double[faceNodeCounts[ifc2]];
#endif
      int num1 = 0, num2 = 0;

      //TODO: replace with material properties
      if(aperture[ifp]<=0)
        aperture[ifp] = 1e-4;

#ifdef USE_STL_IN_CPC
      std::vector<double> xpolypts;
#else
      //determine whether these are in contact
      //TODO: replace 60 hard-code (20 points) with user-defined number?
      double xpolypts[60];
#endif
      unsigned int npts = 0;
      active[ifp] = CheckFaces( ifc1,
                                ifc2,
                                xNodes,
                                dxNodes,
                                faceNodes,
                                faceNodeCounts,
                                faceNodeOffsets,
                                aperture[ifp],
                                dt,
                                xfc1,
                                dxfc1,
                                nx1,
                                x1,
                                dx1,
                                xmin1,
                                xmax1,
                                area1,
                                length1,
                                penetr1,
                                num1,
                                xfc2,
                                dxfc2,
                                nx2,
                                x2,
                                dx2,
                                xmin2,
                                xmax2,
                                area2,
                                length2,
                                penetr2,
                                num2,
                                boxTol,
                                normalTol,
                                cosMin,
                                big,
                                small,
                                xr,
                                xrn,
                                prod1_cp,
                                prod2_cp,
                                e1,
                                e2,
                                penetr_sum,
                                contactArea,
                                xpolypts,
                                npts,
                                tol);
      if(contactArea<=-2.0)
      {
#ifndef USE_STL_IN_CPC
        delete[] x2;
        delete[] dx2;
        delete[] x1;
        delete[] dx1;
        delete[] penetr1;
        delete[] penetr2;
#endif
        return false;
      }

      //if in contact, then resolve the contact
      if(active[ifp])
      {
        //polygon of intersection
        poly_facepairs[npoly] = ifp;
        poly_offsets[npoly] = npoly==0 ? 0 : poly_offsets[npoly-1] + poly_counts[npoly-1];
        poly_counts[npoly] = npts;
        for(unsigned int i = 0; i < npts; i++)
        {
#ifdef USE_STL_IN_CPC
          for(unsigned int j = 0; j < 3; ++j)
            xpoly.push_back(xpolypts[3*i+j]);
#else
          xpoly[3*(i+poly_offsets[npoly])] = xpolypts[3*i];
          xpoly[3*(i+poly_offsets[npoly])+1] = xpolypts[3*i+1];
          xpoly[3*(i+poly_offsets[npoly])+2] = xpolypts[3*i+2];
#endif
        }
        ++npoly;

        //contact area
        contactArea = fabs(contactArea);

        //normal velocity
        double velnContact = 0.;
        for(unsigned int i = 0; i < 3; i++)
          velnContact += xrn[i] * (dxfc1[i] - dxfc2[i]);

        //tangential velocity
        double veltContact[] = {0.,0.,0.};
        for(unsigned int i = 0; i < 3; i++)
          veltContact[i] = (dxfc1[i] - dxfc2[i]) - xrn[i] * velnContact;
        double gapContact = (penetr_sum[0] + penetr_sum[1])/(num1 + num2);

        //get the mass
        double massAverage = 0.;
        {
          double mass1 = 0.;
          for( int i = 0; i < faceNodeCounts[ifc1]; i++)
            mass1 += mNodes[faceNodes[faceNodeOffsets[ifc1]+i]];
          mass1 *= 0.5;
          double mass2 = 0.;
          for( int i = 0; i < faceNodeCounts[ifc2]; i++)
            mass2 += mNodes[faceNodes[faceNodeOffsets[ifc2]+i]];
          mass2 *= 0.5;
          {
            //distance^-1 weighted by the area/face
            double frac1 = contactArea / (fabs(prod1_cp) * area1);
            double frac2 = contactArea / (fabs(prod2_cp) * area2);
            frac1 = 1. < frac1 ? 1. : frac1;
            frac1 = small > frac1 ? small : frac1;
            frac2 = 1. < frac2 ? 1. : frac2;
            frac2 = small > frac2 ? small : frac2;
            //<--why this weighting?
            //estimated mass of impacted faces
            mass1 *= frac1;
            mass2 *= frac2;
            massAverage = mass1 * mass2 / (mass1 + mass2);
          }
        }

        //get composite properties
//        double fct1 = area1 / (area1 + area2);
        double fct1 = 0.5;
        double fct2 = 1-fct1;
        double friction = fct1 * frictionSlope[ifc1] + fct2 * frictionSlope[ifc2];
        double kb = fct1 * kbulk[ifc1] + fct2 * kbulk[ifc2];
        double ks = fct1 * kshear[ifc1] + fct2 * kshear[ifc2];
        double visc = fct1 * viscosity[ifc1] + fct2 * viscosity[ifc2];
        double coh = fct1 * cohesion[ifc1] + fct2 * cohesion[ifc2];
        double dils = fct1 * dilationStress[ifc1] + fct2 * dilationStress[ifc2];
        double dila = fct1 * dilationAngle[ifc1] + fct2 * dilationAngle[ifc2];
        double ress = fct1 * residualSlope[ifc1] + fct2 * residualSlope[ifc2];
        double critss = fct1 * criticalShearSlip[ifc1] + fct2 * criticalShearSlip[ifc2];
        double tco = fct1 * tensileCutOff[ifc1] + fct2 * tensileCutOff[ifc2];
        double clim = fct1 * cohesionLimit[ifc1] + fct2 * cohesionLimit[ifc2];

        //TODO: replace with material properties
        friction = friction <= 0 ? 0.3 : friction;
        kb = kb <= 0 ? 1e9 : kb;
        ks = ks <= 0 ? kb*0.1 : ks;
        //visc = visc <= 0 ? 0.01 : visc;
        //coh = coh <= 0 ? 1e3 : visc;
        critss = critss <= 0 ? 1e-4 : critss;
        clim = clim <= 0 ? 0.9 : clim;

        //cache the normal
        normal[3*ifp] = xrn[0];
        normal[3*ifp+1] = xrn[1];
        normal[3*ifp+2] = xrn[2];

        //set the force vector ON FACE 2 !!
        double force[] = {0., 0., 0.};

#if 0
        SlideContact(xrn,
                     massAverage,
                     friction,
                     velnContact,
                     veltContact,
                     gapContact,
                     penetrationFraction,
                     dt,
                     force);

        StickyContact(xrn,
                      massAverage,
                      velnContact,
                      veltContact,
                      gapContact,
                      penetrationFraction,
                      dt,
                      force);
#endif
        ExplicitContact(xrn,
                        massAverage,
                        contactArea,
                        friction,
                        aperture[ifp],
                        kb,
                        ks,
                        visc,
                        coh,
                        dils,
                        dila,
                        ress,
                        critss,
                        tco,
                        clim,
                        normalApproach[ifp],
                        normalApproachMax[ifp],
                        slip[ifp],
                        pf[ifp],
                        &forceTangential[3*ifp],
                        forceNormal[ifp],
                        velnContact,
                        veltContact,
                        gapContact,
                        dt,
                        force,
                        small);
#if 0
        ImplicitContact(xrn,
                        massAverage,
                        friction,
                        aperture[ifp],
                        kb,
                        ks,
                        visc,
                        coh,
                        dils,
                        dila,
                        ress,
                        critss,
                        tensileSoftening[ifp],
                        dilationLimit[ifp],
                        softening[ifp],
                        tensileDamage[ifp],
                        dilationContact[ifp],
                        fnDilation[ifp],
                        normalApproach[ifp],
                        normalApproachMax[ifp],
                        slip[ifp],
                        pf[ifp],
                        &forceTangential[3*ifp],
                        forceNormal[ifp],
                        contactArea,
                        velnContact,
                        veltContact,
                        gapContact,
                        dt,
                        force);
  #endif
        double factor = 1.; //temp factor to scale the forc
        // find contact force(accelerations) on the nodes:

        //face 1
        for ( int i = 0; i < faceNodeCounts[ifc1]; ++i)
        {
          if(penetr1[i]<0)
          {
            int iNode = faceNodes[faceNodeOffsets[ifc1] + i];
            double weight = (-penetr1[i] > small ? -penetr1[i] : small) / penetr_sum[0];
            double nodmas = mNodes[iNode];
            double coeff = -factor * (prod1_cp < 0 ? -1. : 1.) * weight / nodmas;
            int ii = 3*iNode;
            ddxNodes[ii] += coeff * force[0];
            ddxNodes[ii+1] += coeff * force[1];
            ddxNodes[ii+2] += coeff * force[2];
          }
        }

        //face 2
        for ( int i = 0; i < faceNodeCounts[ifc2]; ++i)
        {
          if(penetr2[i]<0)
          {
            int iNode = faceNodes[faceNodeOffsets[ifc2] + i];
            double weight = (-penetr2[i] > small ? -penetr2[i] : small) / penetr_sum[0];
            double nodmas = mNodes[iNode];
            double coeff = -factor * (prod2_cp < 0 ? -1. : 1.) * weight / nodmas;
            int ii = 3*iNode;
            ddxNodes[ii] += coeff * force[0];
            ddxNodes[ii+1] += coeff * force[1];
            ddxNodes[ii+2] += coeff * force[2];
          }
        }
      }//if in contact
#ifndef USE_STL_IN_CPC
      delete[] penetr1;
      delete[] penetr2;
      delete[] x2;
      delete[] dx2;
#endif
    }//for each face 2
#ifndef USE_STL_IN_CPC
    delete[] x1;
    delete[] dx1;
#endif
  }//for each face 1
  return true;
}

/**
 * @brief Get area, length, center, velocity, normal, etc of face
 * @author Scott Johnson
 * @param[in] ifc Index of the face
 * @param[in] xNodes Array of nodal coordinates
 * @param[in] dxNodes Array of nodal velocities
 * @param[in] faceNodes Array of nodal indices associated with each face
 * @param[in] faceNodeCounts Array of node counts for each face
 * @param[in] faceNodeOffsets Array of faceNodes array index offsets for each face
 * @param[in] dt Time step
 * @param[out] xfc Center of the face
 * @param[out] dxfc Velocity of the face
 * @param[out] nx Normal to the face
 * @param[out] xs Coordinate array for just the nodes associated with the face
 * @param[out] dxs Velocity array for just the nodes associated with the face
 * @param[out] xmin Lower coordinate of the AABB
 * @param[out] xmax Upper coordinate of the AABB
 * @param[out] area Area of the face
 * @param[out] length Characteristic dimension of the face
 */
void CommonPlaneContact::FaceProperties(
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
    )
{
  //get AABB, center, and velocity for each face
  for ( int i = 0; i < faceNodeCounts[ifc]; ++i)
  {
    int js = 3*i;
    int iNode = faceNodes[faceNodeOffsets[ifc] + i];
    int jNode = 3*iNode;
    for(unsigned int j = 0; j < 3; ++j, ++jNode, ++js) {
      dxs[js] = dxNodes[jNode] + 0.5* dt *ddxNodes[jNode];
      xs[js] = xNodes[jNode] + dt*dxs[js];
      xfc[j] += xs[js];
      dxfc[j] += dxs[js];
      xmax[j] = xs[js] > xmax[j] ? xs[js] : xmax[j];
      xmin[j] = xs[js] < xmin[j] ? xs[js] : xmin[j];
    }
  }
  {
    double tmp = 1.0/faceNodeCounts[ifc];
    for(unsigned int i=0;i<3;i++)
    {
      xfc[i] *= tmp;
      dxfc[i] *= tmp;
    }
  }

  //get area and normal
  area = 0.0;
  length = 0.0;
  for ( int i = 0; i < faceNodeCounts[ifc]; ++i)
  {
    double x1 = xfc[0] - xs[3*i];
    double y1 = xfc[1] - xs[3*i+1];
    double z1 = xfc[2] - xs[3*i+2];
    int i1 = i < (faceNodeCounts[ifc]-1) ? i + 1 : 0;
    double x2 = xfc[0] - xs[3*i1];
    double y2 = xfc[1] - xs[3*i1+1];
    double z2 = xfc[2] - xs[3*i1+2];

    area += sqrt((y1 * z2 - z1 * y2) * (y1 * z2 - z1 * y2) + (z1 * x2 - z2 * x1) * (z1 * x2
        - z2 * x1) + (x1 * y2 - y1 * x2) * (x1 * y2 - y1 * x2));
    if (i == 0)
    {
      nx[0] = -(y1 * z2 - z1 * y2) / area;
      nx[1] = -(z1 * x2 - z2 * x1) / area;
      nx[2] = -(x1 * y2 - y1 * x2) / area;
#if 1
      nx[0] *= -1.; nx[1] *= -1.; nx[2] *= -1.;
#endif
    }
  }

  // up until here, area is actually twice the area
  length = sqrt(area * 0.5);// square root of the face area!
}//end of calculation of face properties for the pair

/**
 * @brief Get area, length, center, velocity, normal, etc of face
 * @author Scott Johnson
 * @param[in] ifc1 Index of the first face in the pair (no checking)
 * @param[in] ifc2 Index of the second face in the pair
 * @param[in] xNodes Array of nodal coordinates
 * @param[in] dxNodes Array of nodal velocities
 * @param[in] faceNodes Array of nodal indices associated with each face
 * @param[in] faceNodeCounts Array of node counts for each face
 * @param[in] faceNodeOffsets Array of faceNodes array index offsets for each face
 * @param[in] dt Time step
 * @param[in] xfc1 Center of the first participating face
 * @param[in] dxfc1 Velocity of the first participating face
 * @param[in] nx1 Normal to the first participating face
 * @param[in] x1 Coordinate array for just the nodes associated with the first participating face
 * @param[in] dx1 Velocity array for just the nodes associated with the first participating face
 * @param[in] xmin1 Lower coordinate of the AABB
 * @param[in] xmax1 Upper coordinate of the AABB
 * @param[in] area1 Area of the first participating face
 * @param[in] length1 Characteristic dimension of the first participating face
 * @param[in] xfc2 Center of the second participating face
 * @param[in] dxfc2 Velocity of the second participating face
 * @param[in] nx2 Normal to the second participating face
 * @param[in] x2 Coordinate array for just the nodes associated with the second participating face
 * @param[in] dx2 Velocity array for just the nodes associated with the second participating face
 * @param[in] xmin2 Lower coordinate of the AABB
 * @param[in] xmax2 Upper coordinate of the AABB
 * @param[in] area2 Area of the second participating face
 * @param[in] length2 Characteristic dimension of the second participating face
 * @param[in] boxTol Tolerance for the AABB calculation
 * @param[in] normalTol Tolerance for relative distance between face centers along the normal
 * @param[in] cosMin Tolerance for the minimum magnitude of the dot product of the faces' normals
 * @param[in] big Suitably large number used for initializing mins and maxes
 * @param[in] small Small number for tolerance
 * @param[out] xr Center of the common plane
 * @param[out] xrn Normal to the common plane
 * @param[out] prod1_cp Sum of the normal penetration distances into the common plane
 * @param[out] prod1_cp Sum of the normal penetration distances into the common plane
 * @param[out] e1 Unit vector normal to the common plane normal
 * @param[out] e2 Unit vector normal to both the common plane and the e1 vector
 * @param[out] penetr_sum Sum of the penetration distances
 * @param[out] contactArea Area of contact between the contact pair
 * @param[out] xpoly Coordinates in global frame for the polygon of intersection
 * @param[out] npoly Number of coordinates in the polygon of intersection
 * @return Whether the participating faces are in contact
 */
bool CommonPlaneContact::CheckFaces(
    const int ifc1,
    const int ifc2,
    const double* const xNodes ,
    const double* const dxNodes ,
    const int* const faceNodes ,
    const int* const faceNodeCounts,
    const int* const faceNodeOffsets ,
    const double aper,
    const double dt ,
    const double* const xfc1,
    const double* const dxfc1 ,
    const double* const nx1,
#ifdef USE_STL_IN_CPC
    const std::vector<double>& xs1,
    const std::vector<double>& dxs1 ,
#else
    const double* const xs1,
    const double* const dxs1 ,
#endif
    const double* const xmin1,
    const double* const xmax1,
    const double area1 ,
    const double length1,
#ifdef USE_STL_IN_CPC
    std::vector<double>& penetr1,
#else
    double* const penetr1,
#endif
    int& num1,
    const double* const xfc2,
    const double* const dxfc2 ,
    const double* const nx2,
#ifdef USE_STL_IN_CPC
    std::vector<double>& xs2,
    std::vector<double>& dxs2 ,
#else
    const double* const xs2,
    const double* const dxs2 ,
#endif
    const double* const xmin2,
    const double* const xmax2,
    const double area2 ,
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
    const double tol
    )
{

  //--------------------------------------------------------
  // CHECK FACE PAIR
  //--------------------------------------------------------
  //+++++++++++++++++++++++++++++++
  //CRITERION 1: CHECK AABB OVERLAP
  double asize = 0.5 * (length1 + length2);
  {
    double boxTolFacePair = asize * boxTol;

    //TODO: (HYDROFRAC) make sure this accommodates del_hydro
    for(unsigned int i=0; i<3; ++i)
      if(xmin1[i] > xmax2[i]+boxTolFacePair || xmin2[i] > xmax1[i]+boxTolFacePair)
        return false;
  }

  //+++++++++++++++++++++++++++++++
  //CRITERION 2: CHECK DISTANCE
  {
    double distance = 0.;
    for(unsigned int i=0; i<3; ++i)
    {
      double tmp = xfc1[i] - xfc2[i];
      tmp *= tmp;
      distance += tmp;
    }
    if(distance > tol*tol*asize*asize)
      return false;
  }

  //+++++++++++++++++++++++++++++++
  //CRITERION 3: CHECK OPPOSITE POINTING NORMAL
  {
    double prod = 0.;
    for(unsigned int i=0;i<3;++i)
      prod += nx1[i] * nx2[i];
    if (prod > -cosMin)
      return false;
  }

  // estimate cp from face centers
  {
    double tmp = 0.;
    for(unsigned int i = 0; i < 3; ++i) {
      xr[i] = (xfc1[i] + xfc2[i]) * 0.5;
      xrn[i] = (nx1[i] - nx2[i]) * 0.5;
      tmp += xrn[i] * xrn[i];
    }
    if( !isZero(tmp) && !isEqual(tmp,1) )
    {
      tmp = 1./sqrt(tmp);
      for(unsigned int i = 0; i < 3; ++i)
        xrn[i] *= tmp;
    }
  }

  //+++++++++++++++++++++++++++++++
  //CRITERION 4: CHECK DISTANCE TO COMMON PLANE
  // skip faces for which the sum of their distances to
  // the common plane is more than a fraction of the face sizes
  {
    double normdist = 0.;
    {
      for(unsigned int i=0; i<3; i++)
        normdist += (xfc1[i]-xfc2[i])*xrn[i];
      normdist = fabs(normdist);
    }
    //joint aperture

    //TODO: (HYDROFRAC) add del_hydro into the comparison
    if (normdist > (aper > (normalTol * asize) ? aper : normalTol * asize))
      return false;
  }

  //+++++++++++++++++++++++++++++++
  //CRITERION 5: CHECK FACES VERSUS COMMON PLANE
  // faces should be on opposite sides of the common plane
  // check common plane normal relative to face normals
  prod1_cp = 0.;
  prod2_cp = 0.;
  {
    for(unsigned int i=0; i<3; ++i)
    {
      prod1_cp += (xfc1[i]-xr[i])*xrn[i];
      prod2_cp += (xfc2[i]-xr[i])*xrn[i];
    }
    if(prod1_cp * prod2_cp > 0.)
      return false;
  }

  //+++++++++++++++++++++++++++++++
  //CRITERION 6: CHECK WHETHER NODES INTERPENETRATE OPPOSING PLANE
  for(unsigned int ii = 0; ii < 3; ++ii)
  {
    e1[ii] = 0.;
    e2[ii] = 0.;
  }

  bool eSet = false;
  {
    //find overpenetration for the nodes of the faces
    //face1
    //distance to the plane for each node
    //(dist>0 no penetration, < 0 penetration)
    bool ifpenetr[] = {false,false};
    for(unsigned int ii=0;ii<2;++ii)
    {
      penetr_sum[ii] = 0.;
#ifdef USE_STL_IN_CPC
      const std::vector<double>& xs = ii ? xs2 : xs1;
      std::vector<double>& distance = ii ? penetr2 : penetr1;
#else
      const double* xs = ii ? xs2 : xs1;
      double* distance = ii ? penetr2 : penetr1;
#endif
      int count = ii ? faceNodeCounts[ifc2] : faceNodeCounts[ifc1];
      double dist_min = big;
      const double dfct = ii == 0 ? 1. : -1.;
      for ( int i = 0; i < count; ++i)
      {
        distance[i] = 0.;
        for(unsigned int j=0; j<3; ++j)
          distance[i] += dfct * xrn[j] * (xr[j] - xs[3*i + j]);
        dist_min = dist_min > distance[i] ? distance[i] : dist_min;

        //TODO: (HYDROFRAC) add del_hydro/2 into the comparison
        if (distance[i] < 0.0) {
          ifpenetr[ii] = true;
          if(!eSet)
          {
            double tdist = 0.;
            for(unsigned int j=0;j<3;++j)
            {
              e1[j] = (xs[3*i+j] - xr[j]) + distance[i] * dfct * xrn[j];
              tdist += e1[j]*e1[j];
            }
            if( isZero(tdist) )
            {
              tdist = 1.0/sqrt(tdist);
              for(unsigned int j=0;j<3;++j)
                e1[j] *= tdist;
              eSet = true;
            }
          }
        }
        //TODO: (HYDROFRAC) offset distance[i] by -del_hydro/2
        penetr_sum[ii] += -distance[i] > small ? -distance[i] : small;
      }
    }
    if ((!ifpenetr[0]) && (!ifpenetr[1])) {
      return false;
    }
  }

  // at least one of the nodes from each face is through the common plane
  // (penetrating). find projection of penetrating nodes onto cp
  //define two orthogonal directions in plane (one is taken from the first
  //node penetrating at an angle. the second one is found from the vector
  //product of the plane normal with the first direction

  //find another in-plane direction as e2= n x e1
  assert(eSet);
  {
    e2[0] = xrn[1] * e1[2] - xrn[2] * e1[1];
    e2[1] = xrn[2] * e1[0] - xrn[0] * e1[2];
    e2[2] = xrn[0] * e1[1] - xrn[1] * e1[0];
  }

  // project face 1
  int np1 = 0;
  {
    if(ifc1 == 9 && ifc2 == 16)
      np1 = 1;
    else if(ifc1 == 11 && ifc2 == 24)
      np1 = 2;
    else if(ifc1 == 5 && ifc2 == 18)
      np1 = 3;
  }
  np1 = faceNodeCounts[ifc1];
  double ctr1[3];
#ifdef USE_STL_IN_CPC
  std::vector<double> a1;
#else
  double* a1 = new double[2*(np1+2)];
#endif
  double area1p = -1.;
  num1 = -1;
  {
    //TODO: (HYDROFRAC) add del_hydro into the projection calculation
    ProjectFaceToCommonPlane(xs1, penetr1, e1, e2,
                              xr, xrn, np1,
                              a1, ctr1, num1);
    area1p = PolygonalArea(a1, num1);
    if(area1p < small) {
#ifndef USE_STL_IN_CPC
      delete[] a1;
#endif
      return false;
    }
  }

  // project face 2
  int np2 = faceNodeCounts[ifc2];
  double ctr2[3];
#ifdef USE_STL_IN_CPC
  std::vector<double> a2;
#else
  double* a2 = new double[2*(np2+2)];
#endif
  double area2p = -1.;
  num2 = -1;
  {
    double xn2[3];
    memcpy(xn2, xrn, 3*sizeof(double));
    for(unsigned int i=0;i<3;++i)
      xn2[i] *= -1.;
    //TODO: (HYDROFRAC) add del_hydro into the projection calculation
    ProjectFaceToCommonPlane(xs2, penetr2, e1, e2,
                              xr, xn2, np2,
                              a2, ctr2, num2, true);
    area2p = PolygonalArea(a2, num2);
    if(area2p < small) {
#ifndef USE_STL_IN_CPC
      delete[] a1;
      delete[] a2;
#endif
      return false;
    }
  }

  //update location of the common plane's center
  for(unsigned int j=0;j<3;++j)
    xr[j] = 0.5 * (ctr1[j] + ctr2[j]);

  // find contact area and polygonal geometry
//  if(1)
  {
#ifdef USE_STL_IN_CPC
    std::vector<double> aa;
//    unsigned int npolyMax = 2*(num1 > num2 ? num1 : num2);
    //aa.resize(2*npolyMax, 0);
#else
    double* aa = new double[4*(num1 > num2 ? num1 : num2)];
#endif
//    int numaa;
    contactArea = PolygonalAreaOfIntersection3(a1, a2, num1, num2, aa, npoly);

    //transform to global coordinates
    for(unsigned int i = 0; i < npoly; i++)
    {
#ifndef USE_STL_IN_CPC
      xpoly[3*i] = e1[0] * aa[2*i] + e2[0] * aa[2*i+1] + xr[0];
      xpoly[3*i+1] = e1[1] * aa[2*i] + e2[1] * aa[2*i+1] + xr[1];
      xpoly[3*i+2] = e1[2] * aa[2*i] + e2[2] * aa[2*i+1] + xr[2];
#else
      for(unsigned int j = 0; j < 3; ++j)
        xpoly.push_back(e1[j] * aa[2*i] + e2[j] * aa[2*i+1] + xr[j]);
#endif
    }
#ifndef USE_STL_IN_CPC
    delete[] aa;
#endif
  }
//  else
//  {
//    npoly = 0;
//    contactArea = PolygonalAreaOfIntersection2(a1, a2, num1, num2);
//  }

#ifndef USE_STL_IN_CPC
  delete[] a1;
  delete[] a2;
#endif
  return contactArea > 0.;
}

/**
 * @brief Project the given face onto the common plane
 * @author Scott Johnson
 * @param[in] x Array of nodal positions on the face
 * @param[in] dist Distance of each node into the common plane (num_nod)
 * @param[in] e1 Unit vector in the plane
 * @param[in] e2 Unit vector in the plane that is also orthogonal to e2
 * @param[in] xr Center of the common plane
 * @param[in] xrn Normal to the common plane
 * @param[in] num_nod Number of nodes composing the face
 *
 * @param[out] a Coordinates of the nodes in the plane along e1 and e2 that penetrate the common plane (pair tuples)
 * @param[out] ctr Center of the polygon of intersection with the common plane
 * @param[out] nds Nodal index offsets for the face's nodal index array of nodes that penetrate the common plane (-1 = temporary node)
 * @param[out] num_p Number of nodes penetrating the common plane
 *
 * @param[in] reverse Reverse the ordering of the nodes
 */
void CommonPlaneContact::ProjectFaceToCommonPlane(
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
    const bool reverse)
{
  num_p = 0;
  int nplus1;
  double xc = 0.0;
  double yc = 0.0;
  double zc = 0.0;
  double dxr, dyr, dzr, prod;
  for ( int iii = 0; iii < num_nod; ++iii)
  {
    int n = reverse ? (num_nod - iii - 1) : iii;
    nplus1 = reverse ? (iii == (num_nod-1) ? (num_nod-1) : n-1) : (n == (num_nod-1) ? 0 : n + 1);

    //node n penetrated find proj
    if (dist[n] < 0.0)
    {
      //get x,y,z common plane normal distance between nodal position and common plane
      dxr = (x[3*n]   - xr[0]) + dist[n] * xrn[0];
      dyr = (x[3*n+1] - xr[1]) + dist[n] * xrn[1];
      dzr = (x[3*n+2] - xr[2]) + dist[n] * xrn[2];

      // add distance to xc, yc, zc
      xc += dxr;
      yc += dyr;
      zc += dzr;

      // set the penetrating nodes' e1-e2 coordinate vector position
      // set the local nodal index
      // increment the array position
#ifdef USE_STL_IN_CPC
      double val = e1[0] * dxr + e1[1] * dyr + e1[2] * dzr;
      a.push_back(val);
      val = e2[0] * dxr + e2[1] * dyr + e2[2] * dzr;
      a.push_back(val);
#else
      a[2*num_p] = e1[0] * dxr + e1[1] * dyr + e1[2] * dzr;
      a[2*num_p+1] = e2[0] * dxr + e2[1] * dyr + e2[2] * dzr;
#endif
      ++num_p;
    }

    //find intersection of straight segment [n,n+1] with the common plane
    //only search for an intersection point if the nodal positions
    //lie on opposite sides of the common plane
    if (dist[n] * dist[nplus1] < 0.0)
    {
      //get the fraction of the distance along the line segment at which the line segment intersects the common plane
      prod = fabs(dist[n]) / (fabs(dist[n]) + fabs(dist[nplus1]));
      dxr = x[3*n] + (x[3*nplus1] - x[3*n]) * prod - xr[0];
      dyr = x[3*n+1] + (x[3*nplus1+1] - x[3*n+1]) * prod - xr[1];
      dzr = x[3*n+2] + (x[3*nplus1+2] - x[3*n+2]) * prod - xr[2];

      //dxr, dyr, and dzr are now a new point of intersection with the common plane
      xc += dxr;
      yc += dyr;
      zc += dzr;
#ifdef USE_STL_IN_CPC
      double val = e1[0] * dxr + e1[1] * dyr + e1[2] * dzr;
      a.push_back(val);
      val = e2[0] * dxr + e2[1] * dyr + e2[2] * dzr;
      a.push_back(val);
#else
      a[2*num_p] = e1[0] * dxr + e1[1] * dyr + e1[2] * dzr;
      a[2*num_p+1] = e2[0] * dxr + e2[1] * dyr + e2[2] * dzr;
#endif
      //temporary nodal positions like this must be registered but not referenced
      ++num_p;
    }
  }

  //find the center of projection in global coords
  if (num_p > -1)
  {
    ctr[0] = xc / num_p + xr[0];
    ctr[1] = yc / num_p + xr[1];
    ctr[2] = zc / num_p + xr[2];
  }
}

/**
 * @brief Determine the area of intersection of two arbitrary, coplanar polygons
 * @author Scott Johnson
 * @param[in] a Local planar coordinates of the projected face 1
 * @param[in] b Local planar coordinates of the projected face 2
 * @param[in] na Maximum number of coordinates in a1
 * @param[in] nb Maximum number of coordinates in a2
 * @param[in] big Large value to use to initialize mins
 * @return Area of intersection
 */
double CommonPlaneContact::PolygonalAreaOfIntersection2(
#ifdef USE_STL_IN_CPC
							                                          const std::vector<double>& a,
                                                        const std::vector<double>& b,
#else
                                                        const double* const a,
                                                        const double* const b,
#endif
                                                        const unsigned int na,
                                                        const unsigned int nb
                                                        )
{
#ifdef USE_STL_IN_CPC
  double area = inter((__C_point*)a.data(), (int)na, (__C_point*)b.data(), (int)nb);
#else
  double area = inter((__C_point*)a, (int)na, (__C_point*)b, (int)nb);
#endif
//  std::cout << a[0] << " " << a[1] << "\n"
//      << a[2] << " " << a[3] << "\n"
//      << a[4] << " " << a[5] << "\n"
//      << b[0] << " " << b[1] << "\n"
//      << b[2] << " " << b[3] << "\n"
//      << b[4] << " " << b[5] << "\n"
//      << na << " " << nb << " " << area << "\n";
  return area;
}

/**
 * @brief Determine the area of intersection of two coplanar, convex polygons
 * @author Scott Johnson
 * @param[in] a1 Local planar coordinates of the projected face 1
 * @param[in] a2 Local planar coordinates of the projected face 2
 * @param[in] nmax1 Maximum number of coordinates in a1
 * @param[in] nmax2 Maximum number of coordinates in a2
 * @param[in] big Large value to use to initialize mins
 * @return Area of intersection
 */
double CommonPlaneContact::PolygonalAreaOfIntersection(
#ifdef USE_STL_IN_CPC
    const std::vector<double>& a1,
    const std::vector<double>& a2,
#else
    const double* const a1,
    const double* const a2,
#endif
    const unsigned int nmax1,
    const unsigned int nmax2,
    const double big)
{
  double area = 0.;
  //if either "surface" has fewer than 3 points, it has no area
  if (nmax1 < 3 || nmax2 < 3)
    return 0.;

  // calculate origin shift to avoid roundoff errors
  double zero = 0.0;
  double uvorg[] = {big,big};//u-v min for surfaces 1 & 2
  double uvmin1[] = {big, big};//u-v min for surface 1
  double uvmax1[] = {-big, -big};//u-v max for surface 1
  double uvmin2[] = {big, big};//u-v min for surface 2
  double uvmax2[] = {-big, -big};//u-v max for surface 2

  //get min and max of u and v for the two bodies
  for (unsigned int i1 = 0; i1 < nmax1; ++i1)
  {
    for(unsigned int j = 0; j < 2; ++j)
    {
      if (a1[2*i1 + j] < uvmin1[j])
        uvmin1[j] = a1[2*i1 + j];
      if (a1[2*i1 + j] > uvmax1[j])
        uvmax1[j] = a1[2*i1 + j];
      uvorg[j] = uvorg[j] < a1[2*i1 + j] ? uvorg[j] : a1[2*i1 + j];
    }
  }
  for (unsigned int i2 = 0; i2 < nmax2; ++i2)
  {
    for(unsigned int j = 0; j < 2; ++j)
    {
      if (a2[2*i2 + j] < uvmin2[j])
        uvmin2[j] = a2[2*i2 + j];
      if (a2[2*i2 + j] > uvmax2[j])
        uvmax2[j] = a2[2*i2 + j];
      uvorg[j] = uvorg[j] < a2[2*i2 + j] ? uvorg[j] : a2[2*i2 + j];
    }
  }

  //if the min > max, this is invalid (if the same shape) or not overlapping
  for(unsigned int j = 0; j < 2; ++j) {
    if(uvmin1[j] > uvmax1[j] || uvmin1[j] > uvmax2[j] || uvmin2[j] > uvmax2[j] || uvmin2[j] > uvmax1[j])
      return 0.;
  }

  // loop over segments of polygon #1
  double duv1[] = {0.,0.};//vector difference of uv1p and uv1
  double duv2[] = {0.,0.};//vector difference of uv2p and uv2
  double uv1[] = {0.,0.};//u-v coordinates of current index for loop #1
  double uv1p[] = {0.,0.};//u-v coordinates of current index + 1 for loop #2
  double uv2[] = {0.,0.};//u-v coordinates of current index for loop #1
  double uv2p[] = {0.,0.};//u-v coordinates of current index + 1 for loop #2
  double uvl[] = {0.,0.};
  double uvr[] = {0.,0.};
  double uvm[] = {0.,0.};

  double slope1 = 0.;
  double slope2 = 0.;
  double s = 0., dslope = 0.;
  int i1p = 0, i2p = 0;

  //loop over edges of polygon 1
  for (unsigned int i1 = 0; i1 < nmax1; ++i1)
  {
    i1p = i1 == (nmax1 - 1) ? 0 : (i1 + 1);
    bool cont = false;
    for(unsigned int j=0;j<2;++j)
    {
      uv1[j] = a1[2*i1 + j] - uvorg[j];
      uv1p[j] = a1[2*i1p + j] - uvorg[j];
      duv1[j] = uv1p[j] - uv1[j];
      if( j==0 && isZero(duv1[j]) ) {
        cont = true;
        break;
      }
    }
    if(cont)
      continue;
    slope1 = duv1[1] / duv1[0];

    // loop over edges of polygon #2
    for (unsigned int i2 = 0; i2 < nmax2; ++i2)
    {
      i2p = i2 + 1;
      if (i2 == (nmax2 - 1))
        i2p = 0;
      bool cont2 = false;
      for(unsigned int j=0;j<2;++j)
      {
        uv2[j] = a2[2*i2 + j] - uvorg[j];
        uv2p[j] = a2[2*i2p + j] - uvorg[j];
        duv2[j] = uv2p[j] - uv2[j];
        if( j==0 && isEqual(duv2[j],zero) ) {
          cont2 = true;
          break;
        }
      }
      if(cont2)
        continue;
      slope2 = duv2[1] / duv2[0];

      // determine sign of volume of intersection
      s = duv1[0] * duv2[0];

      // calculate left and right coordinates of overlap
      uvl[0] = uv1[0] < uv1p[0] ? uv1[0] : uv1p[0];
      {
        double tmp = uv2[0] < uv2p[0] ? uv2[0] : uv2p[0];
        if(tmp > uvl[0])
          uvl[0] = tmp;
      }
      uvr[0] = uv1[0] > uv1p[0] ? uv1[0] : uv1p[0];
      {
        double tmp = uv2[0] > uv2p[0] ? uv2[0] : uv2p[0];
        if(tmp < uvr[0])
          uvr[0] = tmp;
      }
      if(uvl[0] > uvr[0])
        continue;

      uvl[1] = uv1[1] + (uvl[0] - uv1[0]) * slope1;
      {
        double tmp = uv1p[1] + (uvl[0] - uv1p[0]) * slope2;
        if(tmp < uvl[1])
          uvl[1] = tmp;
      }
      uvr[1] = uv1[1] + (uvr[0] - uv1[0]) * slope1;
      {
        double tmp = uv1p[1] + (uvr[0] - uv1p[0]) * slope2;
        if(tmp < uvr[1])
          uvr[1] = tmp;
      }

      // check whether lines intersect
      dslope = slope1 - slope2;
      if (!isEqual(dslope,zero) )
      {
        uvm[0] = (uv1p[1] - uv1[1] + slope1 * uv1[0] - slope2 * uv1p[0]) / dslope;
        uvm[1] = uv1[1] + slope1 * (uvm[0] - uv1[0]);
        if (uvm[0] > uvl[0] && uvm[0] < uvr[0])
        {
          // lines intersect, case ii
          area += (s < 0 ? -1. : 1.) * 0.5 * fabs((uvl[1] + uvm[1]) * (uvm[0] - uvl[0]) + (uvm[1] + uvr[1]) * (uvr[0] - uvm[0]));
          continue;
        }
      }
      // lines do not intersect, case i
      area += (s < 0 ? -1. : 1.) * 0.5 * fabs((uvr[0] - uvl[0]) * (uvr[1] + uvl[1]));
    }
  }
  return area;
}

/**
 * @brief Determine the area of intersection and polygon of intersection of two coplanar, convex polygons
 * @author Scott Johnson
 * @param[in] a1 Local planar coordinates of the projected face 1
 * @param[in] a2 Local planar coordinates of the projected face 2
 * @param[in] nmax1 Maximum number of coordinates in a1
 * @param[in] nmax2 Maximum number of coordinates in a2
 * @param[out] aa Local planar coordinates of the polygon of intersection
 * @param[out] nmax Number of entries in aa
 * @return Area of intersection
 *
 * NOTE: for this to work, both a1 and a2 have to be listed in the same loop direction (either CW or CCW)
 * THERE IS NO CHECKING FOR THIS CONDITION, SO BE CAREFUL
 */
double CommonPlaneContact::PolygonalAreaOfIntersection3(
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
    unsigned int& nmax)
{
  //if this is true, there is no contact area
  if(nmax1 < 3 || nmax2 < 3)
    return -1.;

  //there can be no more than 2 intersections per side for interfering convex polygons
  unsigned int nlim = 2*(nmax1 > nmax2 ? nmax1 : nmax2);
  double area = 0.;
  nmax = 0;

#ifdef USE_STL_IN_CPC
  std::vector<int> in1a; in1a.resize(nmax1, -1);
  std::vector<int> in1b; in1b.resize(nmax1, -1);
  std::vector<int> in2a; in2a.resize(nmax2, -1);
  std::vector<int> in2b; in2b.resize(nmax2, -1);
  std::vector<double> aaa;
  std::vector<unsigned int> iaaa;
#else
  //each line segment can have up to 2 intersections
  //for contact between strictly convex polygons
  int* in1a = new int[nmax1];
  int* in1b = new int[nmax1];
  int* in2a = new int[nmax2];
  int* in2b = new int[nmax2];
  //there can be no more than twice as many intersections as line segments
  double* aaa = new double[2*nlim];
  unsigned int* iaaa = new unsigned int[2*nlim];
  for(unsigned int ii1 = 0; ii1 < nmax1; ++ii1)
  {
    in1a[ii1] = -1;
    in1b[ii1] = -1;
  }
  for(unsigned int ii2 = 0; ii2 < nmax2; ++ii2)
  {
    in2a[ii2] = -1;
    in2b[ii2] = -1;
  }
#endif

  //get the number of intersections as well as the points and connectivity of intersection
  unsigned int naaa = Intersections(a1, in1a, in1b, nmax1,
				    a2, in2a, in2b, nmax2,
				    aaa, iaaa);

  unsigned int nlim0 = (nmax1 > nmax2 ? nmax1 : nmax2) + naaa;
  if(nlim0 > nlim)
    nlim0 = nlim;

  //SPECIAL CASE: ONE OF THE POLYGONS IS CONTAINED IN THE OTHER (IE NO INTERSECTIONS)
  if(naaa==0)
    nmax = CoincidentIntersection(a1, nmax1, a2, nmax2, aa);
  else
    nmax = PolygonOfIntersection(a1, in1a, in1b, nmax1, 
				 a2, in2a, in2b, nmax2,
				 aaa, iaaa, naaa, nlim0, aa);

  //THIS IS AN EXCEPTION ... HANDLE APPROPRIATELY
  //THIS FLAGS THAT THE POLYGONS DON'T OVERLAP
  if(nmax == 0) {
    area = -2.0;
  } else {
    area = PolygonalArea(aa, nmax);
  }

#ifndef USE_STL_IN_CPC
  //delete the arrays
  delete[] aaa;
  delete[] iaaa;
  delete[] in1a;
  delete[] in2a;
  delete[] in1b;
  delete[] in2b;
#endif
  return area;
}//PolygonalAreaOfIntersection3


/**
 * @brief Get the area of the convex polygon with the given ordered points
 * @author Scott Johnson
 * @param aa List of points in 2D
 * @return Area of the given polygon
 */
double CommonPlaneContact::PolygonalArea(
#ifdef USE_STL_IN_CPC
    const std::vector<double>& aa,
#else
    const double* const aa,
#endif 
    const unsigned int na)
{
  double area = 0.;
  if(na < 3)
    return area;

  double uvi[] = {0.,0.};
  double uvip[] = {0.,0.};
  double tmp;

  //now, we have the list of vertices in order in aa,
  //so get the area
  for(unsigned int i = 1; i < (na-1); ++i)
  {
    unsigned int ip = i+1;
    //get the cross product of uv(i+1)-uv(0) and uv(i) - uv(0)
    //uv(i)-uv(0)
    uvi[0] = aa[2*i] - aa[0];
    uvi[1] = aa[2*i+1] - aa[1];
    //uv(i+1)-uv(0)
    uvip[0] = aa[2*ip] - aa[0];
    uvip[1] = aa[2*ip+1] - aa[1];
    //add magnitude of cross
    tmp = uvip[0]*uvi[1] - uvi[0]*uvip[1];
    tmp *= tmp;
    area += sqrt(tmp);
  }
  //area is half of the stored value (1/2 |cross|)
  area *= 0.5;
  return area;
}

unsigned int CommonPlaneContact::PolygonOfIntersection(
#ifdef USE_STL_IN_CPC
						       const std::vector<double>& a1, const std::vector<int>& in1a, const std::vector<int>& in1b, const int nmax1,
						       const std::vector<double>& a2, const std::vector<int>& in2a, const std::vector<int>& in2b, const int nmax2,
						       const std::vector<double>& aaa, const std::vector<unsigned int>& iaaa, const unsigned int naaa , const unsigned int nlim,
						       std::vector<double>& aa
#else
						       const double* const a1, const int* const in1a, int* const in1b, const unsigned int nmax1,
						       const double* const a2, const int* const in2a, int* const in2b, const unsigned int nmax2,
						       const double* const aaa, const unsigned int* const iaaa, const unsigned int naaa , const unsigned int nlim,
						       double* const aa
#endif
						       )
{
  //1) INITIALIZE VARIABLES
  unsigned int nmax = 0;
  int i1 = 0, i2 = 0, next = -1, first = -1;

  //1) FIND FIRST INTERSECTION ON LOOP 1
  for (i1 = 0; i1 < nmax1; ++i1)
  {
    if (in1a[i1] > -1)
    {
      first = in1a[i1];

      //add first intersection to list of points
#ifdef USE_STL_IN_CPC
      aa.push_back(aaa[2 * first]);
      aa.push_back(aaa[2 * first + 1]);
#else
      aa[2*nmax] = aaa[2*first];
      aa[2*nmax+1] = aaa[2*first+1];
#endif
      ++nmax;

      //get the next intersection point (-1 if not another on the segment)
      next = in1b[i1] > -1 ? in1b[i1] : -1;
      if (next > -1)
      {
        i1 = iaaa[2 * next];
        i2 = iaaa[2 * next + 1];
      }
      else
      {
        i2 = iaaa[2 * first + 1];
        //check to make sure there isn't another intersection point on 2's segment after "first"
        next = in2a[i2] == first ? in2b[i2] : -1;
      }
      break;
    }
  }//end 1

  //2) TRAVERSE THE INTERSECTION LOOP STARTING FROM THE FIRST INTERSECTION
  while (next != first)
  {
    //THIS IS TANTAMOUNT TO AN EXCEPTION ... HANDLE APPROPRIATELY
    if (nmax > nlim)
      return 0;

    //2a) If the "next" in the loop is an intersection, we just go there
    if (next > -1)
    {
      //add the point to the polygon list
#ifdef USE_STL_IN_CPC
      aa.push_back(aaa[2 * next]);
      aa.push_back(aaa[2 * next + 1]);
#else
      aa[2*nmax] = aaa[2*next];
      aa[2*nmax+1] = aaa[2*next+1];
#endif
      ++nmax;

      //find the next choices
      i1 = iaaa[2 * next];
      i2 = iaaa[2 * next + 1];
      int curr = next;
      if (in1a[i1] == next)
        next = in1b[i1];
      else
        next = -1;
      if (next == -1)
        if (in2a[i2] == curr)
          next = in2b[i2];
    }
    else
    {
      //2b) There is no intersection next on the loop, so first, determine whether we have a choice
      if (i1 < 0)
      {
#ifdef USE_STL_IN_CPC
        aa.push_back(a2[2 * i2]);
        aa.push_back(a2[2 * i2 + 1]);
#else
        aa[2*nmax] = a2[2*i2];
        aa[2*nmax+1] = a2[2*i2+1];
#endif
        ++nmax;

        //go to the next point on 2
        ++i2;
        if (i2 == nmax2)
          i2 = 0;
        next = in2a[i2];

        //check to make sure we're not finished
        if (next == first)
          return nmax;

        if (in2b[i2] > -1)
        {
          //add the point to the polygon list
#ifdef USE_STL_IN_CPC
          aa.push_back(aaa[2 * next]);
          aa.push_back(aaa[2 * next + 1]);
#else
          aa[2*nmax] = aaa[2*next];
          aa[2*nmax+1] = aaa[2*next+1];
#endif
          ++nmax;
          next = in2b[i2];
        }
        i1 = (next > -1) ? static_cast<int> (iaaa[2 * next]) : -1;
      }
      else if (i2 < 0)
      {
#ifdef USE_STL_IN_CPC
        aa.push_back(a1[2 * i1]);
        aa.push_back(a1[2 * i1 + 1]);
#else
        aa[2*nmax] = a1[2*i1];
        aa[2*nmax+1] = a1[2*i1+1];
#endif
        ++nmax;

        //go to the next point on 1
        ++i1;
        if (i1 == nmax1)
          i1 = 0;
        next = in1a[i1];

        //check to make sure we're not finished
        if (next == first)
          return nmax;
        if (in1b[i1] > -1)
        {
          //add the point to the polygon list
#ifdef USE_STL_IN_CPC
          aa.push_back(aaa[2 * next]);
          aa.push_back(aaa[2 * next + 1]);
#else
          aa[2*nmax] = aaa[2*next];
          aa[2*nmax+1] = aaa[2*next+1];
#endif
          ++nmax;
          next = in1b[i1];
        }
        i2 = (next > -1) ? static_cast<int> (iaaa[2 * next + 1]) : -1;
      }
      else
      {
        //we now need to determine whether to go on to the next point on loop 1 or the next on loop 2
        double aa1[] =
        { 0., 0. };
        double aa2[] =
        { 0., 0. };
        {
          //get the square of the distance between the next point on loop 1 and the last point on the intersection loop
          double dd = (a1[2 * i1] - aa[2 * (nmax - 1)]) * (a1[2 * i1] - aa[2 * (nmax - 1)]) + (a1[2
              * i1 + 1] - aa[2 * (nmax - 1) + 1]) * (a1[2 * i1 + 1] - aa[2 * (nmax - 1) + 1]);
          if (dd > 0)
          {
            dd = 1. / sqrt(dd);
            aa1[0] = dd * (a1[2 * i1] - aa[2 * (nmax - 1)]);
            aa1[1] = dd * (a1[2 * i1 + 1] - aa[2 * (nmax - 1) + 1]);
          }
          //aa1 now contains the unit vector between the last point on the intersection loop and the next point on loop 1

          //get the square of the distance between the next point on loop 2 and the last point on the intersection loop
          dd = (a2[2 * i2] - aa[2 * (nmax - 1)]) * (a2[2 * i2] - aa[2 * (nmax - 1)]) + (a2[2 * i2
              + 1] - aa[2 * (nmax - 1) + 1]) * (a2[2 * i2 + 1] - aa[2 * (nmax - 1) + 1]);
          if (dd > 0)
          {
            dd = 1. / sqrt(dd);
            aa2[0] = dd * (a2[2 * i2] - aa[2 * (nmax - 1)]);
            aa2[1] = dd * (a2[2 * i2 + 1] - aa[2 * (nmax - 1) + 1]);
          }
          //aa2 now contains the unit vector between the last point on the intersection loop and the next point on loop 2
        }

        //get the dot product of the unit vector along loop 1 with the position of the last point on the intersection loop
        double d1 = aa1[0] * aa[2 * (nmax - 1)] + aa1[1] * aa[2 * (nmax - 1) + 1];
        //get the dot product of the unit vector along loop 2 with the position of the last point on the intersection loop
        double d2 = aa2[0] * aa[2 * (nmax - 1)] + aa2[1] * aa[2 * (nmax - 1) + 1];

        //whichever dot product is less is the loop which goes closest to the center, so take that one
        if (d1 < d2)
        {
          //follow 1
#ifdef USE_STL_IN_CPC
          aa.push_back(a1[2 * i1]);
          aa.push_back(a1[2 * i1 + 1]);
#else
          aa[2*nmax] = a1[2*i1];
          aa[2*nmax+1] = a1[2*i1+1];
#endif
          ++nmax;
          ++i1;
          if (i1 >= nmax1)
            i1 = 0;
          next = in1a[i1];
          i2 = -1;
        }
        else
        {
          //follow 2
#ifdef USE_STL_IN_CPC
          aa.push_back(a2[2 * i2]);
          aa.push_back(a2[2 * i2 + 1]);
#else
          aa[2*nmax] = a2[2*i2];
          aa[2*nmax+1] = a2[2*i2+1];
#endif
          ++nmax;
          ++i2;
          if (i2 >= nmax2)
            i2 = 0;
          next = in2a[i2];
          i1 = -1;
        }
      }
    }//end 2
  }
  return nmax;
}

unsigned int CommonPlaneContact::CoincidentIntersection(
#ifdef USE_STL_IN_CPC
							const std::vector<double>& a1, const unsigned int nmax1,
							const std::vector<double>& a2, const unsigned int nmax2,
							std::vector<double>& aa
#else
							const double* const a1, const unsigned int nmax1,
							const double* const a2, const unsigned int nmax2,
							double* const aa
#endif
							)
{
  unsigned int nmax = 0;
  double uv[6];
  double uv0[2];
  unsigned int allIn1 = 1;
  unsigned int allIn2 = 0;
//  unsigned int anyOverlap = 0;

  //FIRST, CHECK WHETHER ALL POINTS OF 2 ARE IN 1
  for (unsigned int ii = 0; ii < nmax2; ii++)
  {
    int found = 0;
    uv0[0] = a2[2 * ii];
    uv0[1] = a2[2 * ii + 1];
    for (unsigned int i = 2; i < nmax1; i++)
    {
      uv[0] = a1[2 * (i - 2)];
      uv[1] = a1[2 * (i - 2) + 1];
      uv[2] = a1[2 * (i - 1)];
      uv[3] = a1[2 * (i - 1) + 1];
      uv[4] = a1[2 * i];
      uv[5] = a1[2 * i + 1];
      if (PointInTriangle(uv, uv0) > 0)
      {
        found = 1;
        break;
      }
    }
    if (found == 0)
    {
      allIn1 = 0;
      break;
    }
  }
  if(allIn1==0)
  {
    //SECOND, CHECK WHETHER ALL POINTS OF 1 ARE IN 2
    allIn2 = 1;
//    int found = 0;
    for (unsigned int ii = 0; ii < nmax1; ii++)
    {
      int found = 0;
      uv0[0] = a1[2 * ii];
      uv0[1] = a1[2 * ii + 1];
      for (unsigned int i = 2; i < nmax2; i++)
      {
	uv[0] = a2[2 * (i - 2)];
	uv[1] = a2[2 * (i - 2) + 1];
	uv[2] = a2[2 * (i - 1)];
	uv[3] = a2[2 * (i - 1) + 1];
	uv[4] = a2[2 * i];
	uv[5] = a2[2 * i + 1];
	if (PointInTriangle(uv, uv0) > 0)
	{
	  found = 1;
	  break;
	}
      }
      if (found == 0)
      {
	allIn2 = 0;
	break;
      }
    }
  }


  //determine whether 1 is in 2 or vice versa
  if (allIn1)
  {
    //all the points of polygon 2 are contained in 1
    //so populate the list with the points in 2
    for (unsigned int i = 0; i < 2 * nmax2; i++)
#ifdef USE_STL_IN_CPC
      aa.push_back(a2[i]);
#else
    aa[i] = a2[i];
#endif
    nmax = nmax2;
  }
  else if(allIn2)
  {
    //all the points of polygon 1 are contained in 2
    //so populate the list with the points in 1
    for (unsigned int i = 0; i < 2 * nmax1; i++)
#ifdef USE_STL_IN_CPC
      aa.push_back(a1[i]);
#else
    aa[i] = a1[i];
#endif
    nmax = nmax1;
  }
  else
  {
    //assume by exclusion that the polygons don't indeed overlap!
    //so populate the list with nothing!
    nmax = 0;
  }
  return nmax;
}



/*
 * @brief Determine whether the 2D point is in the 2D triangle
 * @param[in] uv List of points on the triangle (u0,v0,u1,v1,u2,v2)
 * @param[in] uv0 Point to test
 * @return 1 if inside, else 0
 */
int CommonPlaneContact::PointInTriangle(const double* const uv, const double* const uv0)
{
  double duv0[] = {0.,0.};
  duv0[0] = uv[2] - uv[0];
  duv0[1] = uv[3] - uv[1];
  double duv1[] = {0.,0.};
  duv1[0] = uv[4] - uv[0];
  duv1[1] = uv[5] - uv[1];
  double duv[] = {0.,0.};
  duv[0] = uv0[0] - uv[0];
  duv[1] = uv0[1] - uv[1];

  //check that signs of cross products agree
  double d0 = duv[0]*duv0[1] - duv[1]*duv0[0];
  double d1 = duv1[0]*duv[1] - duv1[1]*duv[0];
  if(d0 * d1 < 0)
    return 0;

  duv0[0] *= -1;
  duv0[1] *= -1;
  duv1[0] = uv[4] - uv[2];
  duv1[1] = uv[5] - uv[3];
  duv[0] = uv0[0] - uv[2];
  duv[1] = uv0[1] - uv[3];

  //check that signs of cross products agree
  d0 = duv[0]*duv0[1] - duv[1]*duv0[0];
  d1 = duv1[0]*duv[1] - duv1[1]*duv[0];

  if(d0 * d1 < 0)
    return 0;
  else
    return 1;
}

/**
 * @date July 18, 2011
 * @author Scott Johnson
 * @brief Line intersection finds the intersection (if it exists) between two 2D line segments
 * @param[in] uv1 The pair of u-v coordinates for segment 1 (u0,v0,u1,v1)
 * @param[in] uv2 The pair of u-v coordinates for segment 2 (u0,v0,u1,v1)
 * @param[out] intr The u-v coordinates of the intersection (if it exists)
 * @return 1 if intersection found, else 0
 */
int CommonPlaneContact::LineIntersection(
    const double* const uv1, const double* const uv2,
    double* const intr, const double small)
{
  //line 1
  double uvl1[] = {0.,0.};
  uvl1[0] = uv1[2] - uv1[0];
  uvl1[1] = uv1[3] - uv1[1];
  //normal to line 1
  double nv1[] = {0.,0.};
  nv1[0] = -uvl1[1];
  nv1[1] = uvl1[0];
  //relative positions of line 2 nodes to line 1a
  double duv2[] = {0.,0.};
  duv2[0] = nv1[0] * (uv2[0] - uv1[0]) + nv1[1] * (uv2[1] - uv1[1]);
  duv2[1] = nv1[0] * (uv2[2] - uv1[0]) + nv1[1] * (uv2[3] - uv1[1]);
  if(duv2[0] * duv2[1] > 0)
    return 0;
  //get the position of the intersection along line 2
  double dd = fabs(duv2[0]/(duv2[0] - duv2[1]));
  //double intr[] = {0., 0.};
  intr[0] = uv2[0] + (uv2[2] - uv2[0])*dd;
  intr[1] = uv2[1] + (uv2[3] - uv2[1])*dd;
  //project intr onto line segment 1
  dd = (intr[0] - uv1[0]) * (uvl1[0]) + (intr[1] - uv1[1]) * (uvl1[1]);
  dd /= (uvl1[0] * uvl1[0]) + (uvl1[1] * uvl1[1]);
  //check to make sure the (normalized) distance is within tolerance
  return dd > small && dd < (1. -  small) ? 1 : 0;
}


unsigned int CommonPlaneContact::Intersections(
#ifdef USE_STL_IN_CPC
					       const std::vector<double>& a1, std::vector<int>& in1a, std::vector<int>& in1b, const unsigned int nmax1,
					       const std::vector<double>& a2, std::vector<int>& in2a, std::vector<int>& in2b, const unsigned int nmax2,
					       std::vector<double>& aaa, std::vector<unsigned int>& iaaa
#else
					       const double* const a1, int* const in1a, int* const in1b, const unsigned int nmax1,
					       const double* const a2, int* const in2a, int* const in2b, const unsigned int nmax2,
					       double* const aaa, unsigned int* const iaaa
#endif
					       )
{
  unsigned int naaa = 0;
  //identify intersections
  for (unsigned int i1 = 1; i1 <= nmax1; ++i1)
  {
    unsigned int ii1 = i1 == nmax1 ? 0 : i1;
    in1a[ii1] = -1;
    in1b[ii1] = -1;
    //get segment 1
    double uv1[] =
    { 0., 0., 0., 0. };
    uv1[0] = a1[2 * i1 - 2];
    uv1[1] = a1[2 * i1 - 1];
    uv1[2] = a1[2 * ii1];
    uv1[3] = a1[2 * ii1 + 1];
    for (unsigned int i2 = 1; i2 <= nmax2; ++i2)
    {
      unsigned int ii2 = i2 == nmax2 ? 0 : i2;
      //get segment 2
      double uv2[] =
      { 0., 0., 0., 0. };
      uv2[0] = a2[2 * i2 - 2];
      uv2[1] = a2[2 * i2 - 1];
      uv2[2] = a2[2 * ii2];
      uv2[3] = a2[2 * ii2 + 1];
      //find intersection of segments 1 and 2
      double intr[] =
      { 0., 0. };
      int ii = LineIntersection(uv1, uv2, intr);
      //if there is one, add it to the list of intersection positions and indices
      if (ii > 0)
      {
#ifdef USE_STL_IN_CPC
        aaa.push_back(intr[0]);
        aaa.push_back(intr[1]);
        iaaa.push_back(ii1);
        iaaa.push_back(ii2);
#else
        aaa[2*naaa] = intr[0];
        aaa[2*naaa+1] = intr[1];
        iaaa[2*naaa] = ii1;
        iaaa[2*naaa+1] = ii2;
#endif

        //add entry to in1a/in1b
        if (in1a[ii1] == -1)
        {
          in1a[ii1] = (int) naaa;
        }
        else
        {
          //make sure the intersection points are ordered correctly
          int ia = in1a[ii1];
          double atmp[2] =
          { 0., 0. };
          atmp[0] = aaa[2 * ia] - uv1[0];
          atmp[1] = aaa[2 * ia + 1] - uv1[1];
          double dda = atmp[0] * atmp[0] + atmp[1] * atmp[1];
          atmp[0] = intr[0] - uv1[0];
          atmp[1] = intr[1] - uv1[1];
          double ddb = atmp[0] * atmp[0] + atmp[1] * atmp[1];
          if (ddb < dda)
          {
            in1b[ii1] = in1a[ii1];
            in1a[ii1] = (int) naaa;
          }
          else
          {
            in1b[ii1] = (int) naaa;
          }
        }

        //add entry to in2a/in2b
        if (in2a[ii2] == -1)
        {
          in2a[ii2] = (int) naaa;
        }
        else
        {
          //make sure the intersection points are ordered correctly
          int ia = in2a[ii2];
          double atmp[2] =
          { 0., 0. };
          atmp[0] = aaa[2 * ia] - uv2[0];
          atmp[1] = aaa[2 * ia + 1] - uv2[1];
          double dda = atmp[0] * atmp[0] + atmp[1] * atmp[1];
          atmp[0] = intr[0] - uv2[0];
          atmp[1] = intr[1] - uv2[1];
          double ddb = atmp[0] * atmp[0] + atmp[1] * atmp[1];
          if (ddb < dda)
          {
            in2b[ii2] = in2a[ii2];
            in2a[ii2] = (int) naaa;
          }
          else
          {
            in2b[ii2] = (int) naaa;
          }
        }
        ++naaa;
      }
    }
  }//end id intersections
  return naaa;
}

/**
 * @brief Calculate forces due to a slide contact
 * @author Scott Johnson
 * @param[in] xrn Normal to the contact
 * @param[in] massAverage Average mass of the interface
 * @param[in] frictionSlopeContact Friction coefficient for the contact
 * @param[in] velnContact Relative normal velocity magnitude of the contact
 * @param[in] veltContact Relative tangential velocity vector of the contact
 * @param[in] gapContact Gap between participating faces
 * @param[in] penetrationFraction
 * @param[in] dtst Timestep size
 * @param[out] force Force calculated for the slide contact
 * @param[in] small Lower tolerance value
 */
void CommonPlaneContact::SlideContact(const double* const xrn,
                   const double massAverage,
                   const double frictionSlopeContact,
                   const double velnContact,
                   const double* const veltContact,
                   const double gapContact,
                   const double penetrationFraction,
                   const double dtst,
                   double* const force,
                   const double small)
{
//  double one = 1.0;
  double zero = 0.0;
  double fric_mu = frictionSlopeContact;
  double veln = velnContact;
  double mass_avr = massAverage;

  // This is Equation (24) of Vorobiev UCRL-TR-227085
  double pen_force = 2. * gapContact / dtst / dtst;
  //this is user controlled parameter from 0 to 1
  double frac_total_penalty = penetrationFraction;
  double fnx = mass_avr * (veln / dtst + frac_total_penalty * pen_force);
  double fny = mass_avr * (veln / dtst + frac_total_penalty * pen_force);
  double fnz = mass_avr * (veln / dtst + frac_total_penalty * pen_force);

  //no tensile force
  double xn = xrn[0];
  double yn = xrn[1];
  double zn = xrn[2];

  fnx = (fnx > zero ? fnx : zero) * xn;
  fny = (fny > zero ? fny : zero) * yn;
  fnz = (fnz > zero ? fnz : zero) * zn;
  double fnContact = fnx * xn + fny * yn + fnz * zn;

  double fsx = 0.;
  double fsy = 0.;
  double fsz = 0.;
  //add frictional force
  if (fric_mu > small)
  {
    double veltx = veltContact[0];
    double velty = veltContact[1];
    double veltz = veltContact[2];
    double veltmod = sqrt(veltx * veltx + velty * velty + veltz * veltz);
    if (fabs(veltmod) > small)
    {
      double veltxu = veltx / veltmod;
      double veltyu = velty / veltmod;
      double veltzu = veltz / veltmod;
      //unconstrained friction fource (Ft=Fn*mu)
      fsx = fnContact * fric_mu * veltxu;
      fsy = fnContact * fric_mu * veltyu;
      fsz = fnContact * fric_mu * veltzu;
    }
  }
  // total force
  force[0] = fnContact * xn + fsx;
  force[1] = fnContact * yn + fsy;
  force[2] = fnContact * zn + fsz;
}

/**
 * @brief Calculate forces due to a sticky contact
 * @author Scott Johnson
 * @param[in] xrn Normal to the contact
 * @param[in] massAverage Average mass of the interface
 * @param[in] velnContact Relative normal velocity magnitude of the contact
 * @param[in] veltContact Relative tangential velocity vector of the contact
 * @param[in] gapContact Gap between participating faces
 * @param[in] penetrationFraction
 * @param[in] dtst Timestep size
 * @param[out] force Force calculated for the sticky contact
 * @param[in] small Lower tolerance value
 */
void CommonPlaneContact::StickyContact( const double* const xrn,
                    const double massAverage,
                    const double velnContact,
                    const double* const veltContact,
                    const double gapContact,
                    const double penetrationFraction,
                    const double dtst,
                    double* const force,
                    const double small  )
{
  //small
//  double one = 1.0;
  double mass_avr = massAverage;
  double veln = velnContact;
  double veltx = veltContact[0];
  double velty = veltContact[1];
  double veltz = veltContact[2];

  // This is Equation (24) of Vorobiev UCRL-TR-227085
  double pen_force = 2. * gapContact / dtst / dtst;
  //this is user controlled parameter from 0 to 1
  double frac_total_penalty = penetrationFraction;
  double fnx = mass_avr * (veln / dtst + frac_total_penalty * pen_force);
  double fny = mass_avr * (veln / dtst + frac_total_penalty * pen_force);
  double fnz = mass_avr * (veln / dtst + frac_total_penalty * pen_force);

  double xn = xrn[0];
  double yn = xrn[1];
  double zn = xrn[2];

  force[0] = fnx * xn;
  force[1] = fny * yn;
  force[2] = fnz * zn;

  double fsx = mass_avr * veltx / dtst;
  double fsy = mass_avr * velty / dtst;
  double fsz = mass_avr * veltz / dtst;

//  double fnContact = fnx * xn + fny * yn + fnz * zn;

  // total force (to be applied to the nodes)
  force[0] = fnx + fsx;
  force[1] = fny + fsy;
  force[2] = fnz + fsz;
}

/**
 * @brief Calculate forces due to an explicit contact
 * @author Scott Johnson
 * @param[in] xrn Normal to the contact
 * @param[in] massAverage Average mass of the interface
 * @param[in] contactArea Contact area
 * @param[in] frictionSlopeContact Friction coefficient for the contact
 * @param[in] apertureContact Current aperture of the contact
 * @param[in] kbulkContact Bulk modulus of the contact
 * @param[in] kshearContact Shear modulus of the contact
 * @param[in] viscContact Viscosity of the contact
 * @param[in] cohesionContact Cohesion of the contact
 * @param[in] dilationStressContact Stress magnitude due to dilation
 * @param[in] dilationAngleContact Dilation angle
 * @param[in] residualSlopeContact Post peak internal friction coefficient
 * @param[in] criticalShearSlipContact Critical slip in shear
 * @param[in] tensileCutOffContact Tensile cut-off force
 * @param[in] climContact Cohesion limit
 * @param[in,out] dnContact Normal displacement increment
 * @param[in,out] dnContactMax Normal displacement maximum
 * @param[in,out] dsContact Shear displacement increment
 * @param[in] pfContact Pore fluid
 * @param[in,out] fsContact Magnitude of the shear force
 * @param[in,out] fnContact Magnitude of the normal force
 * @param[in] velnContact Relative normal velocity magnitude of the contact (v1 - v0)
 * @param[in] veltContact Relative tangential velocity vector of the contact (v1 - v0)
 * @param[in] gapContact Gap between participating faces
 * @param[in] dtst Timestep size
 * @param[out] force Force calculated for the explicit contact
 * @param[in] small Lower tolerance value
 */
void CommonPlaneContact::ExplicitContact( const double* const xrn,
                   const double massAverage ,
                   const double contactArea,
                   const double frictionSlopeContact,
                   const double apertureContact,
                   const double kbulkContact,
                   const double kshearContact,
                   const double viscContact,
                   const double cohesionContact,
                   const double dilationStressContact,
                   const double dilationAngleContact ,
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
                   const double gapContact ,
                   const double dtst,
                   double* const force,
                   const double small)
{
/*
  std::cout<<"  massAverage = "<<massAverage<<std::endl;
  std::cout<<"  contactArea = "<<contactArea<<std::endl;
  std::cout<<"  frictionSlopeContact = "<<frictionSlopeContact<<std::endl;
  std::cout<<"  apertureContact = "<<apertureContact<<std::endl;
  std::cout<<"  kbulkContact = "<<kbulkContact<<std::endl;
  std::cout<<"  kshearContact = "<<kshearContact<<std::endl;
  std::cout<<"  viscContact = "<<viscContact<<std::endl;
  std::cout<<"  cohesionContact = "<<cohesionContact<<std::endl;
  std::cout<<"  dilationStressContact = "<<dilationStressContact<<std::endl;
  std::cout<<"  dilationAngleContact = "<<dilationAngleContact<<std::endl;
  std::cout<<"  residualSlopeContact = "<<residualSlopeContact<<std::endl;
  std::cout<<"  criticalShearSlipContact = "<<criticalShearSlipContact<<std::endl;
  std::cout<<"  tensileCutOffContact = "<<tensileCutOffContact<<std::endl;
  std::cout<<"  climContact = "<<climContact<<std::endl;
  std::cout<<"  dnContact = "<<dnContact<<std::endl;
  std::cout<<"  dnContactMax = "<<dnContactMax<<std::endl;
  std::cout<<"  dsContact = "<<dsContact<<std::endl;
  std::cout<<"  pfContact = "<<pfContact<<std::endl;
  std::cout<<"  velnContact = "<<velnContact<<std::endl;
  std::cout<<"  gapContact = "<<gapContact<<std::endl;
  std::cout<<"  dtst = "<<dtst<<std::endl;
  std::cout<<"  small = "<<small<<std::endl;
  */
  //explicit joint with stiff hardening law
  //small
  double one = 1.0;
  double zero = 0.0;

  // contact parameters defined by the boundary condition on the faces:
  //joint aperture
  double aper = apertureContact;
  // squared aperture
  double aper2 = aper * aper;
  // cut-off how stiff the joints can be
  double clim = climContact * aper;
  // reference normal modulus
  double kbulk_0 = kbulkContact;
  // reference shear modulus
  double kshear_0 = kshearContact;
  // viscosity
  double visco = viscContact;
  // cohesion
  double cohx = cohesionContact;
  // dilation stress
  double dil_stressx = dilationStressContact;
  // dilation angle
//  double dil_ang = dilationAngleContact;
  // friction slope tan(fi)
  double slope = frictionSlopeContact;
  // residual friction slope
  double res_slope = residualSlopeContact;
  // critical plastic slip when the friction drops to the residual slope
  double sslip_crit = criticalShearSlipContact;
  // tensile cut-off
  double ten_cut_off = tensileCutOffContact;

  // history variables
  // contact area
  double area = contactArea;
  // normal disp
  double undisp_old = dnContact;
  // normal velocity
  double veln = velnContact;
  // shear slip
  double sslip_old = dsContact;
  // fluid pressure
  double p_f = pfContact;

  // copy shear forces
  double fsx = fsContact[0];
  double fsy = fsContact[1];
  double fsz = fsContact[2];

  //projected normal penetration:
  double undisp = undisp_old + veln * dtst;
  if(aper < undisp)
    undisp = aper;

  //  save new normal closure
  dnContact += veln * dtst;

  // save max closer (normal disp)
  double gap_maximal = undisp > dnContactMax ? undisp : dnContactMax;
  dnContactMax = gap_maximal;

  // stiffness calculation
  double kbulk = kbulk_0;
  {
    //find effective stiffness
    double dAper = (aper - gap_maximal) > clim ? (aper - gap_maximal) : clim;
    if (gap_maximal <= undisp)
    {
      double tmp = aper - 0.5 * (gap_maximal + undisp_old);
      dAper = tmp > clim ? tmp : clim;
    }
    dAper *= dAper;
    double fct = isZero(dAper,0.0) ? (isZero(aper2,0.0) ? 0. : 100) : aper2 / dAper;
    if(fct > one)
      kbulk *= fct > 100 ? 100 : fct;
  }
//  std::cout<<"  kbulk = "<<kbulk<<std::endl;
  //constant shear
  double kshear = kshear_0;

  // find normal force increment in the direction normal to face 2 (i.e., -xrn)
  double dfnorm = area * veln * dtst * kbulk;
  double fn = fnContact + dfnorm;

  /*
  std::cout<<"FORCE CALCULATION"<<std::endl;
  std::cout<<"  fnContact = "<<fnContact<<std::endl;
  std::cout<<"       area = "<<area<<std::endl;
  std::cout<<"       veln = "<<veln<<std::endl;
  std::cout<<"       dtst = "<<dtst<<std::endl;
  std::cout<<"      kbulk = "<<kbulk<<std::endl;
  std::cout<<"     dfnorm = "<<dfnorm<<std::endl;
*/

  /*
  std::cout.precision(15);
  std::cout<<"dfnorm = "<<dfnorm<<std::endl;
  std::cout<<"area = "<<area<<std::endl;
  std::cout<<"veln = "<<veln<<std::endl;
  std::cout<<"dtst = "<<dtst<<std::endl;
  std::cout<<"kbulk = "<<kbulk<<std::endl;
  */

  //CURRENTLY, TENSILE STRENGTH IS NOT HANDLED; IF IT OCCURS
  //NULLIFY THE CONTACT
//  if(fn < 0.)
//  {
//    // normal force to be saved
//    fnContact = 0;
//
//    // save shear forces
//    fsContact[0] = 0;
//    fsContact[1] = 0;
//    fsContact[2] = 0;
//
//    // total force (to be applied to the nodes)
//    force[0] = 0;
//    force[1] = 0;
//    force[2] = 0;
//
//    return;
//  }

  // find shear force increment
  // increment shear forces opposite to velt
  fsx += area * veltContact[0] * dtst * kshear;
  fsy += area * veltContact[1] * dtst * kshear;
  fsz += area * veltContact[2] * dtst * kshear;

  // check coulomb friction law and scale shear forces
  // and correct the normal force for dilation

  // convert dilation stress to force
  double dil_stress = dil_stressx * area;

  // convert cohesion stress to force
  double coh = cohx * area;

  // find current slope
  double cur_slope = 0.;
  {
    double tmp = 1. - sslip_old / sslip_crit;
    cur_slope = res_slope + (slope - res_slope) * (zero > tmp ? zero : tmp);
  }
  double fn_minus_p_f = fn - p_f;
  double fsmax = coh + (zero > fn_minus_p_f ? zero : fn_minus_p_f) * cur_slope;
  if(small > fsmax)
    fsmax = small;

  double fsmod = sqrt(fsx * fsx + fsy * fsy + fsz * fsz);

  if (fn < dil_stress && dil_stress > 0)
  {
//    double cof1 = kbulk / kshear * dil_ang * (one - fn / dil_stress);
//    double cof2 = cof1 * cur_slope;

    //dilation adjusted normal force;
//    fn += cof1 * ((fsmod - fsmax) < zero ? zero : (fsmod-fsmax)) / (one + cof2);
  }

  // new shear strength based on normal force corrected for dilation
  {
    double tmp = zero > fn_minus_p_f ? zero : fn_minus_p_f;
    tmp *= slope;
    tmp += coh;
    fsmax = small > tmp ? small : tmp;
  }

  double dslip = zero;
  if (fsmod > fsmax)
  {
    double scale = fsmax / fsmod;
    fsx *= scale;
    fsy *= scale;
    fsz *= scale;
    dslip = (fsmod - fsmax) / (area * kshear);
  }

  // new normal force on face 2 including viscosity
  double fnx = (fn + veln * visco * area) * xrn[0];
  double fny = (fn + veln * visco * area) * xrn[1];
  double fnz = (fn + veln * visco * area) * xrn[2];

  // apply tensile cut-off
  if (fn < -area * ten_cut_off)
  {
    fnx = 0.0;
    fny = 0.0;
    fnz = 0.0;
    fsx = 0.0;
    fsy = 0.0;
    fsz = 0.0;
  }

  //  save plastic slip
  dsContact = sslip_old + dslip;

  // normal force to be saved
  fnContact = fn;

  fsx = fsy = fsz = 0.0;

  // save shear forces
  fsContact[0] = fsx;
  fsContact[1] = fsy;
  fsContact[2] = fsz;

  // total force (to be applied to the nodes)
  force[0] = fnx + fsx;
  force[1] = fny + fsy;
  force[2] = fnz + fsz;
}



/**
 * @brief Calculate forces due to an implicit contact
 * @author Scott Johnson
 * @param[in] xrn Normal to the contact
 * @param[in] massAverage Average mass of the interface
 * @param[in] frictionSlopeContact Friction coefficient for the contact
 * @param[in] apertureContact Current aperture of the contact
 * @param[in] kbulkContact Bulk modulus of the contact
 * @param[in] kshearContact Shear modulus of the contact
 * @param[in] viscContact Viscosity of the contact
 * @param[in] cohesionContact Cohesion of the contact
 * @param[in] dilationStressContact Stress magnitude due to dilation
 * @param[in] dilationAngleContact Dilation angle
 * @param[in] residualSlopeContact Post peak internal friction coefficient
 * @param[in] criticalShearSlipContact Critical slip in shear
 *
 * @param[in] tensileSoftening Tensile softening rate parameter
 * @param[in] dilationLimit Dilation limit parameter
 * @param[in] softening Softening rate parameter
 * @param[in,out] tensileDamage Tensile damage parameter
 * @param[in] dilationContact Previous dilation
 * @param[in,out] fnDilation Change in normal force due to dilation
 *
 * @param[in,out] dnContact Normal displacement increment
 * @param[in,out] dnContactMax Normal displacement maximum
 * @param[in,out] double dsContact Shear displacement increment
 * @param[in] double pfContact Pore fluid
 * @param[in,out] fsContact Magnitude of the shear force
 * @param[in,out] fnContact Magnitude of the normal force
 * @param[in] contactArea Area of the contact
 * @param[in] faceArea1 Area of the first participating face
 * @param[in] faceArea2 Area of the second participating face
 * @param[in] massContact1 Mass associated with the first participating face
 * @param[in] massContact2 Mass associated with the second participating face
 * @param[in] velnContact Relative normal velocity magnitude of the contact
 * @param[in] veltContact Relative tangential velocity vector of the contact
 * @param[in] gapContact Gap between participating faces
 * @param[in] dtst Timestep size
 * @param[out] force Force calculated for the implicit contact
 * @param[in] small Lower tolerance value
 */
void CommonPlaneContact::ImplicitContact( const double* const xrn,
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
                      const double gapContact ,
                      const double dtst,
                      double* const force,
                      const double small)
{

  //implicit joint with stiff hardening law
  //small
  double one = 1.0;
  double zero = 0.0;
  double mass_avr = massAverage;

  // contact parameters defined by the boundary condition on the faces:
  // joint aperture
  double aper = apertureContact;
  // squared aperture
  double aper2 = aper * aper;
  // cut-off how stiff the joints can be
  double clim = aper * 1.e-8;
  // reference normal modulus
  double kbulk_0 = kbulkContact;
  // reference shear modulus
  double kshear_0 = kshearContact;
  // viscosity
  double visco = viscContact;

  // dilation stress
  double dil_stressx = dilationStressContact;
  // dilation angle
  double dil_ang = dilationAngleContact;
  // friction slope tan(fi)
  double slope = frictionSlopeContact;
  // residual friction slope
  double res_slope = residualSlopeContact;
  // critical plastic slip when the friction drops to the residual slope
  double sslip_crit = criticalShearSlipContact;

  // history variables:
  // contact area
  double area = contactArea;
  // normal disp
  double undisp_old = dnContact;
  // normal velocity
  double veln = velnContact;
  // shear slip
  double sslip_old = dsContact;
  // fluid pressure
  double p_f = pfContact;

//   RES_SLOPE=CONT_PARAM(ifr1)  //final slope
//   SLOPE=CONT_PARAM(ifr0)      //initial slope
//   SSLIP_CRIT=CONT_PARAM(icrs) //current plastic slip
//
//   CUR_SLOPE=RES_SLOPE+(SLOPE-RES_SLOPE)*MAX(ZERO,ONE-SSLIP/SSLIP_CRIT) //CURRENT SLOPE
//   SIG_T0=CONT_PARAM(icoh)/SLOPE     //initial tensile limit found from cohesion
//
//   //COF1,COF2,COF3,COF4 -are derivative for cohesion,normal stress, fluid pressure and friction slope
//   //with respect to the shear slip COF4 and COF1 are calculated below
//   //and the others later when the dilation is active and the new normal modulus is known
//
//   COF4=(-SLOPE+RES_SLOPE)/SSLIP_CRIT                   //Dslope /DU_sp
//   COF1=-CONT_PARAM(isoft)*CUR_SLOPE + SIG_T0*COF4     //dCOHX/dU_sp
//
//   TEN_SOFT=CONT_PARAM(itensoft)                       //tensile softening rate for tensile limit
//   DLIMIT=CONT_PARAM(idlim)              //dilation limit
//   APER=CONT_PARAM(iaper)             //aperture
//   VISCO=CONT_PARAM(ivis)             //viscosity
//   DIL_STRESSX=CONT_PARAM(icrf)   //critical stress for dilation
//   DIL_ANG=CONT_PARAM(idil)       //dilation slope
//
//   KSHEAR_0=CONT_PARAM(iks)//initial shear stiffness
//   KBULK_0=CONT_PARAM(ikn) //initial normal stiffness
//
//   //find current cohesion (MPa)
//   COHX=MAX(ZERO,CONT_PARAM(icoh) +COF1*SSLIP-TEN_SOFT*CONT_HIST(itdam)) //linear function of slip

  //current slope
  double cur_slope = res_slope;
  {
    double tmp = one - sslip_old / sslip_crit;
    cur_slope += (slope - res_slope) * (tmp > zero ? tmp : zero);
  }
  double sig_to = cohesionContact / slope;

  //TODO: HERE
  double ten_soft = tensileSoftening; // tensile softening rate parameter
  double dlimit = dilationLimit;// dilation limit parameter

  // current cohesion
  double cof4 = (res_slope - slope) / criticalShearSlipContact;
  double cof1 = -softening * cur_slope + sig_to * cof4;
  double dcohx = cof1 * sslip_old - ten_soft * tensileDamage;
  double cohx = zero > (cohesionContact + dcohx) ? zero : (cohesionContact - dcohx);

  //use constant shear
  double kshear = kshear_0;

  //find effective stiffness
  double gap_maximal = dnContactMax;
  double stiffness = (aper - gap_maximal) > clim ? (aper - gap_maximal) : clim;
  double kbulk = kbulk_0;
  {
    double tmp = aper2 / (stiffness * stiffness);
    if(tmp > one)
      kbulk *= tmp;
  }

  // copy shear forces
  double fsx = fsContact[0];
  double fsy = fsContact[1];
  double fsz = fsContact[2];

  //tensile cut-off
  double ten_cut_off = area * cohx / slope;

  //find trial displacement assuming that kbulk is constant
  double fn_min = -ten_cut_off;
  double fnor = fn_min;
  bool cut_off_applied = true;
  double fn_back = fnContact;

  if (fn_back > fn_min)
  {
    fnor = fn_back; //current normal force acting on the common plane
    cut_off_applied = false;
  }

  //coefficient alpha 0- explicit,0.5 half_explicit,1-implicit
  double co_a = 0.5;
  double co_1 = 0.5 * dtst * dtst / mass_avr;
  double dln_1 = veln * dtst - co_1 * fnor; //displacement if const force

  double dfnorm = zero;
  double dln = zero;

  //if joint is closing account for the change of the stiffness
  if ( (dln_1 > zero) && isEqual(undisp_old,dnContactMax) )
  {
    //nonlinear loading
    //solve quadratic equation
    double co_b = (one + co_1 * area * kbulk * co_a) * ((aper - undisp_old) > clim ? (aper - undisp_old) : clim) + dln_1;
    double co_c = dln_1 * ((aper - undisp_old) > clim ? (aper - undisp_old): clim);
    double det = co_b * co_b - 4 * co_c;
    if (det > zero)
      dln = 0.5 * (co_b - sqrt(det));
    else
      dln = dln_1 / (one + co_1 * area * kbulk * co_a); //displacement if const bulk

    //find normal force increment
    dfnorm = (veln * dtst - dln) / co_1 - fnor;
  }
  else
  {
    //linear reload/unload
    if (cut_off_applied)
      dln = dln_1; //unload
    else
      dln = dln_1 / (one + co_1 * area * kbulk * co_a); //reload with constant stiffness: Force is assumed to be changing

    // find normal force increment
    dfnorm = area * dln * kbulk;
  }

  //increment normal force and displacement
  double fn = fn_back + dfnorm;
  double undisp = undisp_old + dln;

  if (undisp > dnContactMax)
    dnContactMax = undisp;

  // save normal displacement
  dnContact = undisp;

  // find shear force increment
  // increment shear forces
  fsx += area * veltContact[0] * dtst * kshear;
  fsy += area * veltContact[1] * dtst * kshear;
  fsz += area * veltContact[2] * dtst * kshear;

  // check coulomb friction law
  double dil_stress = dil_stressx * area; //convert dilation stress to force
  double coh = cohx * area; //convert cohesion stress to force
  double fn_eff = zero > (fn - p_f) ? zero : (fn - p_f); //effective normal force
  double fsmax = small > (coh + fn_eff * cur_slope) ? small : (coh + fn_eff * cur_slope);//current yield
  double fsmod = sqrt(fsx * fsx + fsy * fsy + fsz * fsz); //elastic trial
  double dslip = zero;
  double scale = one;

  //check for yield condition
  if (fsmod > fsmax)
  {
    //yielding
    //account for dilation
    if (fn_eff < dil_stress)
    {
      //dilation law rate multipliers
      double cof = dil_ang * (zero > (one - fn_eff / dil_stress) ? zero : (one - fn_eff / dil_stress)); //limit dilation as fn increases
      double dila = dilationContact;// old dilation
      double cof2 = cof * (zero > (one - dila / aper) ? zero : (one - dila / aper)); //limit dilation as u-> a
      double cof3 = zero; //dependence of fluid pressure on shear slip
      double aa = fn_eff * cof4 + area * (cof1 + (kbulk * cof2 - cof3) * cur_slope);
      //implicit slip increment consistent with dilation and fluid pressure
      dslip = (fsmod - fsmax) / (area * kshear + aa);

      //dilation opening integration u_dila< a
      double dopen = dila;
      //new dilation:
      dilationContact = dlimit - (dlimit - dila) * exp(-cof * dslip / aper);
      dopen = dilationContact - dopen; // change in dilation opening

      fnDilation += area * kbulk * dopen; //save fn increment due to dilation

      fn = fn + area * kbulk * dopen; // increment total normal stress
      //increment cohesion
      coh = coh - cof1 * dslip * area;

      //scale shear stresses
      fsmax = fsmax + aa * dslip;

      scale = fsmax / fsmod;
      fsx *= scale;
      fsy *= scale;
      fsz *= scale;
    }
    else
    {
      //* no dilation case
      double aa = fn_eff * cof4 + area * cof1;
      //implicit slip increment consistent with dilation and fluid pressure
      dslip = (fsmod - fsmax) / (area * kshear + aa);

      //increment cohesion
      coh = coh - cof1 * dslip * area;

      //scale shear stresses
      fsmax = fsmax + aa * dslip;

      scale = fsmax / fsmod;
      fsx *= scale;
      fsy *= scale;
      fsz *= scale;
    }
  }

  ten_cut_off = coh / cur_slope; //tensile force current limit includes tensile softening

  double fnx = zero;
  double fny = zero;
  double fnz = zero;

  if (fn + veln * visco * area < -ten_cut_off)
  {
    fnx = -ten_cut_off * xrn[0];
    fny = -ten_cut_off * xrn[1];
    fnz = -ten_cut_off * xrn[2];

    //increment tensile damage
    tensileDamage += (-ten_cut_off - fn - veln * visco * area) / (area * (kbulk_0
        - ten_soft));
  }
  else
  {
    //new normal force
    fnx = (fn + veln * visco * area) * xrn[0];
    fny = (fn + veln * visco * area) * xrn[1];
    fnz = (fn + veln * visco * area) * xrn[2];
  }

  //  save plastic slip
  dsContact = sslip_old + dslip;

  // normal force to be saved
  fnContact = fn;

  // save shear forces
  fsContact[0] = fsx;
  fsContact[1] = fsy;
  fsContact[2] = fsz;

  // total force (to be applied to the nodes)
  force[0] = fnx + fsx;
  force[1] = fny + fsy;
  force[2] = fnz + fsz;
}

#ifdef USE_STL_IN_CPC
#undef USE_STL_IN_CPC
#endif
