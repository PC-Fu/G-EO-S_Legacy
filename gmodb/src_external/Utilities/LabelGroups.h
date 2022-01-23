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
 * @file LabelGroups.h
 * @author walsh24
 * @date Nov 21, 2013
 */


#ifndef LABELGROUPS_H
#define LABELGROUPS_H

#include "DisjointSet.h"
#include <map>
#include <vector>

/** LabelGroups
 * Assigns a unique identifier to connected groups and returns the number of groups found
 * connectivityMap - connections between one node and its immediate neighbors
 * groups - vector of pairs linking unique node identifiers to group ids
 * idOffset - offset for the unique id assigned to each group
**/ 

template<class LabelType>
int LabelGroups(std::map<LabelType,std::vector<LabelType> >& connectivity,
                std::vector<std::pair<LabelType,int> >& groups, int idOffset){
  int nn = connectivity.size();
  groups.resize(nn);
  std::map<LabelType, DSN > sets;

  // build groups
  int count = 0;
  typename std::map<LabelType,std::vector<LabelType> >::iterator it = connectivity.begin();
  for (; it!=connectivity.end(); ++it){
  
    LabelType indx = it->first;
    
    groups[count] = std::pair<LabelType,int> (indx,-1); 
    count++;
    
    sets[indx] = DSN(-1,0); // -1: node info - not used here
    DSN& thisNode = sets[indx];
    std::vector<LabelType>& nbrs =  it->second;
    for(int i =0; i <nbrs.size(); ++i ){
      int nbrIndx = nbrs[i];
      if(indx > nbrIndx){
        DSN& prevNode = sets[nbrIndx];
        thisNode.Union(prevNode);
      }
    }
  } 

  // build set of unique roots
  std::map<DSN*,int> rootMap;
  for(int i = 0; i < nn; ++i){
      DSN* root = sets[i].FindRoot();
      rootMap[root]=0;
  }

  // number roots
  std::map<DSN*,int>::iterator itr_end = rootMap.end();
  int regionCount = 0;
  for(  std::map<DSN*,int>::iterator itr = rootMap.begin();itr != itr_end; ++itr){
    itr->second = regionCount + idOffset;
    ++regionCount;
  }
  
  // record groups
  for(int i = 0; i < nn; ++i){
    LabelType id = groups[i].first;
    DSN* root = sets[id].FindRoot();
    groups[i].second = rootMap[root];
  }
  
  return regionCount;
}

#endif
