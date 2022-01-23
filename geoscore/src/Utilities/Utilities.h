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
/*
 * Utilities.h
 *
 *  Created on: Sep 13, 2010
 *      Author: settgast1
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "Common/typedefs.h"
#include <sys/resource.h>
#include "Common/GPException.h"

#include <map>
#include <set>
#include <algorithm>

/////////////////////////////////////////////////
// Forward declaration of templated functions

template< class T >
void PushFieldForwardInTime( const realT& dt,
                             const Array1dT< T >& dfield,
                             Array1dT< T >& field );

template< class T >
void IntegrateFieldInTime( const realT& dt,
                           const Array1dT< T >& field,
                           Array1dT< T >& Ifield );

template< class T >
inline void IntegrateField( const realT& dt,
                            const T& dfdt,
                            T& df );

template< typename T >
void SetConstPointer( T* const& pointer,  T*  newpointer );

template < typename T1, typename T2 >
void ClearStlMapValues( std::map<T1,T2>& Map );

template <class ElementClass, class SetClass>
bool isMember(const ElementClass& x, const SetClass& aSetOrMap);

template< typename T1, typename T2 >
const T2* stlMapLookupPointer( const std::map<T1,T2>& Map, const T1& key );

template< typename T1, typename T2 >
T2* stlMapLookupPointer( std::map<T1,T2>& Map, const T1& key );


template< typename T1, typename T2 >
const T2& stlMapLookup( const std::map<T1,T2>& Map, const T1& key, const std::string& message="" );

template< typename T1, typename T2 >
T2& stlMapLookup( std::map<T1,T2>& Map, const T1& key, const std::string& message="" );


rArray1d logspace(realT start, realT stop, int count=100);
rArray1d linspace(realT start, realT stop, int count=100);



/////////////////////////////////////////////////





template< class T >
inline void PushFieldForwardInTime( const realT& dt,
                                    const Array1dT< T >& dfield,
                                    Array1dT< T >& field )
{
  T dfieldDt;

  const int N = field.size();
  for( int a=0 ; a<N ; ++a )
  {
    dfieldDt = dfield(a);
    dfieldDt *= dt;

    field(a) += dfieldDt ;
  }
}


template< class T >
inline void IntegrateFieldInTime( const realT& dt,
                                  const Array1dT< T >& field,
                                  Array1dT< T >& Ifield )
{
  const int N = field.size();
  for( int a=0 ; a<N ; ++a )
  {
    Ifield(a) = field(a) ;
    Ifield(a) *= dt;
  }
}


template< class T >
inline void IntegrateField( const realT& dt,
                            const T& dfdt,
                            T& df )
{
  df = dfdt;
  df *= dt;
}



template< typename T >
inline void CopyGlobalToLocal(const lArray1d& globalToLocalRelation,
                              const Array1dT< T >& globalField,
                              Array1dT< T >& localField)
{
  const typename Array1dT<T>::size_type N = globalToLocalRelation.size() ;

  for( typename Array1dT<T>::size_type a=0 ; a<N ; ++a )
  {
    localField[a] = globalField[ globalToLocalRelation[a] ];
  }
}



template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              const Array1dT< T >& globalField,
                              Array1dT< T >& localField)
{
  const typename Array1dT<T>::size_type N = localField.size() ;

  for( typename Array1dT<T>::size_type a=0 ; a<N ; ++a )
  {
    localField[a] = globalField[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              const Array1dT< T >& globalField1,
                              const Array1dT< T >& globalField2,
                              Array1dT< T >& localField1,
                              Array1dT< T >& localField2 )
{
  const typename Array1dT<T>::size_type N = localField1.size() ;

  for( typename Array1dT<T>::size_type a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void CopyGlobalToLocal(const localIndex* __restrict__ const globalToLocalRelation,
                              const Array1dT< T >& globalField1,
                              const Array1dT< T >& globalField2,
                              const Array1dT< T >& globalField3,
                              Array1dT< T >& localField1,
                              Array1dT< T >& localField2,
                              Array1dT< T >& localField3 )
{
  const typename Array1dT<T>::size_type N = localField1.size() ;

  for( typename Array1dT<T>::size_type a=0 ; a<N ; ++a )
  {
    localField1[a] = globalField1[ globalToLocalRelation[a] ];
    localField2[a] = globalField2[ globalToLocalRelation[a] ];
    localField3[a] = globalField3[ globalToLocalRelation[a] ];
  }
}

template< typename T >
inline void AddLocalToGlobal( const localIndex* __restrict__ const globalToLocalRelation,
                              const Array1dT< T >& localField,
                              Array1dT< T >& globalField)
{
  const typename Array1dT<T>::size_type N = localField.size() ;

  for( typename Array1dT<T>::size_type a=0 ; a<N ; ++a )
  {
    globalField[ globalToLocalRelation[a] ] += localField[a];
  }
}

template< typename T >
inline void AddLocalToGlobal( const localIndex* __restrict__ const globalToLocalRelation,
                              const Array1dT< T >& localField1,
                              const Array1dT< T >& localField2,
                              Array1dT< T >& globalField1,
                              Array1dT< T >& globalField2 )
{
  const typename Array1dT<T>::size_type  N = localField1.size() ;

  for( typename Array1dT<T>::size_type a=0 ; a<N ; ++a )
  {
    globalField1[ globalToLocalRelation[a] ] += localField1[a];
    globalField2[ globalToLocalRelation[a] ] += localField2[a];
  }
}



template< typename T >
inline bool listsHaveEqualPermutations( const T* const list1, const T* const list2, const localIndex n )
{
  localIndex count = 0;
  bool rval = false;
  for( localIndex i1=0 ; i1<n ; ++i1 )
    for( localIndex i2=0 ; i2<n ; ++i2 )
    {
      if( list1[i1] == list2[i2] )
      {
        ++count;
        break;
      }
    }

  if( count == n )
    rval = true;

  return rval;

}

inline bool isEven(int x) {
  return !(x&1);
}

inline bool isOdd(int x) {
  return (x&1);
}


/// find if object is member of vector
template <class ElementClass>
bool isMember(const ElementClass& x, const std::vector<ElementClass>& aVec) {return ( std::find(aVec.begin(), aVec.end(), x) !=  aVec.end() );}
inline bool isMember(localIndex x, const lArray1d& aVec) {return ( std::find(aVec.begin(), aVec.end(), x) !=  aVec.end() );}


/// find if object is member of set or map
template <class ElementClass, class SetClass>
bool isMember(const ElementClass& x, const SetClass& aSetOrMap) {return ( aSetOrMap.find(x) != aSetOrMap.end() );}



/// permutation tensor
inline int eijk(int i,int j,int k){ 
  return ((i-j)*(j-k)*(k-i))/2; 
}


template< typename T >
void SetConstPointer( T* const& pointer,  T*  newpointer )
{
   T** temp = const_cast< T**>(&pointer);
  *temp = newpointer;

  return;
}



template< typename T1, typename T2 >
T2* stlMapLookupPointer( std::map<T1,T2>& Map, const T1& key )
{
  T2* rval = NULL;
  typename std::map<T1,T2>::iterator MapIter = Map.find( key );
  if( MapIter!=Map.end()  )
  {
    rval = &(MapIter->second);
  }


  return rval;
}


template< typename T1, typename T2 >
const T2* stlMapLookupPointer( const std::map<T1,T2>& Map, const T1& key )
{
  const T2* rval = NULL;
  typename std::map<T1,T2>::const_iterator MapIter = Map.find( key );
  if( MapIter!=Map.end()  )
  {
    rval = &(MapIter->second);
  }


  return rval;
}


template< typename T1, typename T2 >
T2& stlMapLookup( std::map<T1,T2>& Map, const T1& key, const std::string& message )
{
  typename std::map<T1,T2>::iterator MapIter = Map.find( key );
  if( MapIter==Map.end()  )
  {
    std::cout<<std::endl;
    std::stringstream st;
    st << "Error in stlMapLookup. Key not found in map! key: " << key << " message: " << message <<"\n";
    throw GPException(st.str().c_str());
  }

  return MapIter->second;
}


template< typename T1, typename T2 >
const T2& stlMapLookup( const std::map<T1,T2>& Map, const T1& key, const std::string& message )
{
  return (stlMapLookup( const_cast<std::map<T1,T2>&>(Map), key, message ));
}


/*
 * Initialize map values.
 * std::map<int,int> aMap = CreateStlMap<int, int >(1,2)(3,4)(5,6);
 */
template <typename K, typename V>
class CreateStlMap
{
private:
    std::map<K, V> m_map;
public:
    CreateStlMap(const K& key, const V& value){ m_map[key] = value; }

    CreateStlMap<K, V>& operator()(const K& key, const V& value)
    {
        m_map[key] = value;
        return *this;
    }

    operator std::map<K, V>() { return m_map; }
};


template < typename T1, typename T2 >
void ClearStlMapValues( std::map<T1,T2>& Map )
{
  for( typename std::map<T1,T2>::iterator iter=Map.begin() ; iter!=Map.end() ; ++iter )
  {
    iter->second.clear();
  }
}

/*
 * Initialize vector values.
 * Usage:
 *   std::vector<int> aVect;
 *   aVect += 1,1,2,3,4;
 */
template <class T> class vector_inserter{
public:
    std::vector<T>& v;
    vector_inserter(std::vector<T>& vv):v(vv){}
    vector_inserter& operator,(const T& val){v.push_back(val);return *this;}
};
template <class T> vector_inserter<T>& operator+=(std::vector<T>& v,const T& x);

template <class T> vector_inserter<T>& operator+=(std::vector<T>& v,const T& x){
    return vector_inserter<T>(v),x;
}



/*
 * Initialize string vector values with chars
 * Usage:
 *   sArray1d aVect;
 *   aVect += "The","quick","brown","fox";
 */
class svector_inserter{
public:
    sArray1d& v;
    svector_inserter(sArray1d& vv):v(vv){}
    svector_inserter& operator,(const char* val){v.push_back(std::string(val));return *this;}
};
svector_inserter& operator+=(sArray1d& v,const std::string& x);

inline
svector_inserter& operator+=(sArray1d& v,const char* x){
    return svector_inserter(v),x;
}



/// Static Assert
/// Check satement at compile time
template <bool b>
struct gp_static_assert{};
// Specialization with member function
template <>
struct gp_static_assert<true>
{
  static void is_valid() {}
};
// use: gp_static_assert<TEST>is_valid();



inline bool isGTE0( const int i )
{ return (i>=0); }

/**
 * @author Randy Settgast
 * @param val1
 * @param val2
 * @param tolfac
 * @return
 */
inline bool isEqual( const realT& val1, const realT& val2, const realT& tolfac=0.0 )
{
  realT tol = 0.0;
  if( tolfac > 1.0e-15 )
    tol = fabs(tolfac) * (fabs(val1)+fabs(val2))*0.5;
  return val1<=(val2+tol) && val1>=(val2-tol);
}

inline realT Power(const realT val, const realT exponent)
{
  if(isEqual(exponent, 0.5))
    return sqrt(val);
  else if(isEqual(exponent, 1.5))
    return val*sqrt(val);
  else if(isEqual(exponent, 2.0))
    return val*val;
  else
    return pow(val, exponent);
}

inline bool isZero( const realT& val, const realT& tol=0.0 )
{
  if( val<=tol && val>=-tol )
  {
    return true;
  }
  else
  {
    return false;
  }
}


/**
 * @author Randy Settgast
 * @return cpu usage
 *
 * This function uses the rusage structure to query elapsed system time and user time, and returns
 * the result.
 */
inline realT getcputime(void)        
{ 
  struct timeval tim;        
  struct rusage ru;        
  getrusage(RUSAGE_SELF, &ru);        

  tim=ru.ru_utime;        
  realT t=(realT)tim.tv_sec + (realT)tim.tv_usec / 1.0e6;

  tim=ru.ru_stime;        
  t+=(realT)tim.tv_sec + (realT)tim.tv_usec / 1.0e6;
  return t; 
}

template< typename TYPE >
inline void Intersection( const std::set<TYPE>& set1, const std::set<TYPE>& set2, std::set<TYPE>& intersection )
{
  intersection.clear();
  typename std::set<TYPE>::const_iterator iter_1 = set1.begin();
  typename std::set<TYPE>::const_iterator iter_2 = set2.begin();

  while( iter_1!=set1.end() && iter_2!=set2.end() )
  {
    if( *iter_1 == *iter_2 )
    {
      intersection.insert(*iter_1);
      ++iter_1;
      ++iter_2;
    }
    else if( (*iter_1)<(*iter_2) )
    {
      ++iter_1;
    }
    else if( (*iter_1)>(*iter_2) )
    {
      ++iter_2;
    }
  }
}

template< typename TYPE >
inline void Intersection( const std::set<TYPE>& set, const Array1dT<TYPE>& array, std::set<TYPE>& intersection )
{
  intersection.clear();

  for( typename Array1dT<TYPE>::const_iterator iter_arr=array.begin() ; iter_arr!=array.end() ; ++iter_arr )
  {
    if( set.count( *iter_arr ) == 1 )
    {
      intersection.insert(*iter_arr);
    }
  }
}

template< typename TYPE >
inline void Intersection( const std::set<TYPE>& set, const Array1dT<TYPE>& array, Array1dT<TYPE>& intersection )
{
  intersection.clear();

  for( typename std::set<TYPE>::const_iterator iter_arr=array.begin() ; iter_arr!=array.end() ; ++iter_arr )
  {
    if( set.count( *iter_arr ) == 1 )
    {
      intersection.push_back(*iter_arr);
    }
  }
}

inline
rArray1d linspace(realT start, realT stop, int count){
  rArray1d rv(count,start);

  if(count > 1){

    realT dL = (stop-start)/double(count-1);
    realT sum = start;
    for(int i = 1;i < count-1; ++i){
      sum += dL;
      rv[i] =  sum;
    }

    rv[count -1] = stop;
  }
  return rv;
}


inline
rArray1d logspace(realT start, realT stop, int count){
  rArray1d rv(count,start);

  if(count > 1){

    realT dL = std::pow(stop/start,1.0/double(count-1));
    realT prod = start;
    for(int i = 1;i < count-1; ++i){
      prod *= dL;
      rv[i] =  prod;
    }

    rv[count -1] = stop;
  }
  return rv;
}




#endif /* UTILITIES_H_ */
