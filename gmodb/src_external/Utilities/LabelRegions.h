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
//  All rights reserved.
//
//  This file is part of GPAC.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file LabelRegions.h
 * @author walsh24
 * @date Nov 21, 2013
 */


#ifndef LABELREGIONS_H
#define LABELREGIONS_H

#include "DisjointSet.h"
#include <map>
#include <vector>

/** LabelRegions
 * Assigns a unique identifier to contiguous regions and returns the number of regions found
 * isRegion - binary region geometry: 1 = a region node, 0 = a boundary node
 * regions - record of the region ids, -1 indicates non-region (boundary node)
 * nx - x dimension
 * ny - y dimension
 * nz - z dimension
 * idOffset - offset for the unique id assigned to each region
**/ 
template<class BooleanVector>
int LabelRegions(const BooleanVector isRegion, std::vector<int>& regions, int nx, int ny, int nz, int idOffset){
  int nxy = nx*ny;
  int nxyz = nxy*nz;
  regions.resize(nxyz);
  std::vector< DSN > sets( nxyz,  DSN(-1,0) );

  // build regions
  for(int k = 0; k < nz; ++k){
    int kp = (k==0)? 0: k-1;
    for(int j = 0; j < ny; ++j){
      int jp = (j==0)? 0: j-1;
      for(int i = 0; i < nx; ++i){
        int ip = (i==0)? 0: i-1;
  
        int indx = (k*ny+j)*nx+i;

        if(isRegion[indx]){

          DSN& thisNode = sets[indx];
      
          int prevXindx = (k *ny+ j)*nx+ip;
          int prevYindx = (k *ny+jp)*nx+ i;
          int prevZindx = (kp*ny+ j)*nx+ i;

          if(isRegion[prevXindx]){
            DSN& prevXNode = sets[prevXindx];
            thisNode.Union(prevXNode);
          }

          if(isRegion[prevYindx]){
            DSN& prevYNode = sets[prevYindx];
            thisNode.Union(prevYNode);
          }

          if(isRegion[prevZindx]){
            DSN& prevZNode = sets[prevZindx];
            thisNode.Union(prevZNode);
          }

        }
      } // i loop
    } // j loop
  } // k loop

  // build set of unique roots
  std::map<DSN*,int> rootMap;
  for(int i = 0; i < nxyz; ++i){
    if(isRegion[i]){
      DSN* root = sets[i].FindRoot();
      rootMap[root]=0;
    }
  }

  // number roots
  std::map<DSN*,int>::iterator itr_end = rootMap.end();
  int regionCount = 0;
  for(  std::map<DSN*,int>::iterator itr = rootMap.begin();itr != itr_end; ++itr){
    itr->second = regionCount + idOffset;
    ++regionCount;
  }
  
  // record regions
  for(int i = 0; i < nxyz; ++i){
    if(isRegion[i]){
      DSN* root = sets[i].FindRoot();
      regions[i] = rootMap[root];
    } else {
      regions[i] = -1;
    }
  }
  
  return regionCount;
}

#endif
