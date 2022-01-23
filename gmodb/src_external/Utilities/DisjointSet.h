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
 * @file DisjointSet.h
 * @author walsh24
 * @date Nov 21, 2013
 */

#ifndef DISJOINTSET_H
#define DISJOINTSET_H
#include <iostream>

template<class T, class rank_T>
class DisjointSetNode{

  public:
    DisjointSetNode();
    DisjointSetNode(T data);
    DisjointSetNode(T data, const rank_T& rank);
    DisjointSetNode(const DisjointSetNode<T, rank_T>& dsn);
    DisjointSetNode<T, rank_T>& operator=(const DisjointSetNode<T, rank_T>& dsn);
    DisjointSetNode<T, rank_T>* FindRoot(void);
    void Union(DisjointSetNode<T, rank_T>& ds);
    T m_data;  // data carried on the node (does not affect the ranking of the node)
    rank_T m_rank; // data used to rank the node (must have < operator defined)
  private:
    DisjointSetNode<T, rank_T>* m_parent;
};

struct ParallelRank{
  int local_rank;
  int omp_rank;
  int mpi_rank;
  ParallelRank(int r):local_rank(r),omp_rank(r),mpi_rank(r){};
  bool operator<(const ParallelRank & other) const { 
    if (mpi_rank < other.mpi_rank){
      return true;
    } else if (mpi_rank == other.mpi_rank){
      if (omp_rank < other.omp_rank){
        return true;
      } else if (omp_rank == other.omp_rank){
        return (local_rank < other.local_rank);
      }
    } 
    return false; 
  }

  bool operator==(const ParallelRank & other) const{
    return (local_rank == other.local_rank)&&(omp_rank == other.omp_rank)&&(mpi_rank == other.mpi_rank);
  }
   
};

template<class T, class rank_T> 
DisjointSetNode<T,rank_T>::DisjointSetNode(T data):
  m_data(data),
  m_rank(0)
{
  m_parent =  this;
}

template<class T, class rank_T> 
DisjointSetNode<T,rank_T>::DisjointSetNode():
  m_data(),
  m_rank(0)
{
  m_parent =  this;
}


template<class T, class rank_T> 
DisjointSetNode<T,rank_T>::DisjointSetNode(T data, const rank_T& rank):
  m_data(data),
  m_rank(rank)
{
  m_parent =  this;

}

template<class T, class rank_T> 
DisjointSetNode<T,rank_T>::DisjointSetNode(const DisjointSetNode<T,rank_T>& rhs):
  m_data(rhs.m_data),
  m_rank(rhs.m_rank)
{
  m_parent =(rhs.m_parent == &rhs)? this : rhs.m_parent;
}

template<class T, class rank_T> 
DisjointSetNode<T,rank_T>&
DisjointSetNode<T,rank_T>::operator=(const DisjointSetNode<T, rank_T>& rhs){
  if(this != &rhs){
    m_data = rhs.m_data;
    m_rank = rhs.m_rank;
    m_parent =(rhs.m_parent == &rhs)? this : rhs.m_parent;
  }
  return *this;
}

template<class T, class rank_T> 
void DisjointSetNode<T,rank_T>::Union(DisjointSetNode<T, rank_T>& ds){
  DisjointSetNode<T, rank_T>* xRoot = this->FindRoot();
  DisjointSetNode<T, rank_T>* yRoot = ds.FindRoot();
  if( (*yRoot).m_rank <  (*xRoot).m_rank ){
    (*yRoot).m_parent = xRoot;
  } else if ((*xRoot).m_rank < (*yRoot).m_rank){
    (*xRoot).m_parent = yRoot;
  } else {
    (*yRoot).m_parent = xRoot;
    (*xRoot).m_rank++;
  }
}



template<class T, class rank_T> 
DisjointSetNode<T, rank_T>* DisjointSetNode<T,rank_T>::FindRoot(void){
  //find the root
  DisjointSetNode<T, rank_T>* root = m_parent;
  while((*root).m_parent != root) root = (*root).m_parent;

  //update parent pointers (flatten tree)
  DisjointSetNode<T, rank_T>* x = this;
  DisjointSetNode<T, rank_T>* xp;
  while((*x).m_parent != x){
    xp = (*x).m_parent;
    (*x).m_parent = root;
    x = xp;
  }
  return root;
}


// typedefs
typedef DisjointSetNode< int, int> DSN;

#endif
