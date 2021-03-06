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
/**
 * File: FileManagerT.h
 * Class provides file IO
 * created : RRS (10/11/2001)
 */
 
#ifndef _BIN_STREAM_H_
#define _BIN_STREAM_H_

// ***** Included Headers *****************************************************
#include "Common/Common.h"
#include "Common/GPException.h"
#include <map>
#include <set>
#include <fstream>
#include <iostream>


// ****************************************************************************
// ***** BINSTREAM CLASS DECLARATION ******************************************
// ****************************************************************************
class BinStream
{
public:
  BinStream(void);
  virtual ~BinStream(void);

  virtual void open( const char* filename, const bool truncate ) = 0;
  virtual void close(void) = 0;

protected:

};




// ****************************************************************************
// ***** oBINSTREAM CLASS DECLARATION *****************************************
// ****************************************************************************
class oBinStream : public BinStream
{
public:
  oBinStream(void);
  ~oBinStream(void);

  void open( const char* filename, const bool truncate = false );
  void close(void);

  template<class TYPE>
  void write( const TYPE* const p_var , const int var_length )
  {
    output.write( reinterpret_cast<const char*>(p_var) , sizeof(TYPE)*var_length );
  }
  
  template< typename TYPE >
  void write( const TYPE& val )
  {
    this->write( &val , 1 );
  }

  void write( const std::string& str )
  {
    std::string::size_type length = str.size();
    this->write( length );
    this->write( str.c_str(), str.size() );
  }

  template< typename TYPE >
  void write( const Array1dT<TYPE>& array )
  {
    const typename Array1dT<TYPE>::size_type length = array.size();
    this->write( length );
    this->write( array.data(), array.size() );
  }

  template< typename TYPE >
  void write( const std::set<TYPE>& set )
  {
    const typename std::set<TYPE>::size_type length = set.size();
    this->write( length );
    for( typename std::set<TYPE>::const_iterator i=set.begin() ; i!=set.end() ; ++i )
    {
      this->write( *i );
    }
  }


  template< typename TYPE >
  void write( const Array2dT<TYPE>& array )
  {
    const typename Array2dT<TYPE>::size_type length = array.size();
    this->write( length );

    const typename Array2dT<TYPE>::size_type dimension[2] = { array.Dimension(0), array.Dimension(1) };
    this->write( dimension, 2 );
    this->write( array.data(), array.size() );
  }


  template< typename TYPE >
  void write( const Array1dT<Array1dT<TYPE> >& array )
  {
    const typename Array1dT<Array1dT<TYPE> >::size_type length0 = array.size();
    this->write( length0 );

    for( typename Array1dT<Array1dT<TYPE> >::const_iterator i=array.begin() ; i!=array.end() ; ++i )
    {
      this->write(*i);
    }
  }


  template< typename TYPE >
  void write( const Array1dT<std::set<TYPE> >& array )
  {

    const typename Array1dT<std::set<TYPE> >::size_type length0 = array.size();
    this->write( reinterpret_cast<const char*>(&length0), sizeof(typename Array1dT<std::set<TYPE> >::size_type) );

    for( typename Array1dT<std::set<TYPE> >::const_iterator i=array.begin() ; i!=array.end() ; ++i )
    {
      this->write( *i );
    }
  }


  template< typename T1, typename T2 >
  void write( const std::map< T1, T2 >& datamap )
  {
    this->write( datamap.size() );
    for( typename std::map<T1,T2>::const_iterator i=datamap.begin() ; i!=datamap.end() ; ++i )
    {
      this->write( i->first );
      this->write( i->second );
    }
  }

  template< typename T >
  void write( const std::map< std::string, T>& member )
  {



    // first get the number of map entries that we are writing.
    typename std::map< std::string, T >::size_type numEntries = 0;
    // Iterate over all entries in the member map
    for( typename std::map< std::string, T >::const_iterator i = member.begin() ; i!=member.end() ; ++i )
    {
      // the field name is the key to the map
      const std::string fieldName = i->first;

      // check to see if the field should be written
      std::map<std::string, FieldBase*>::const_iterator fieldAttributes = FieldInfo::AttributesByName.find(fieldName);

      if( fieldAttributes == FieldInfo::AttributesByName.end() )
      {
        ++numEntries;
      }
      else if( fieldAttributes->second->m_WriteToRestart )
      {
        ++numEntries;
      }
    }
    // write the number of entries that we are writing out.
    this->write(numEntries);



    // iterate over all entries in the member map
    for( typename std::map< std::string, T >::const_iterator i = member.begin() ; i!=member.end() ; ++i )
    {
      // the field name is the key to the map
      const std::string fieldName = i->first;

      // check to see if the field should be written
      std::map<std::string, FieldBase*>::const_iterator fieldAttributes = FieldInfo::AttributesByName.find(fieldName);

      bool writeField = false;
      if( fieldAttributes == FieldInfo::AttributesByName.end() )
      {
        writeField = true;
      }
      else if( fieldAttributes->second->m_WriteToRestart )
      {
        writeField = true;
      }

      if( writeField )
      {
//        std::cout<<"    writing "<<fieldName<<std::endl;

        // the field data is mapped value
        const T& fieldData = i->second;

        // write the field name
        this->write( fieldName );

        // write the field data
        this->write( fieldData );
      }
    }
  }

  //***** Data Member Declarations **********************************************
private:
  std::ofstream output;

public:
  std::ofstream::pos_type tellp(void)
  { return output.tellp(); }

};








// ****************************************************************************
// ***** iBINSTREAM CLASS DECLARATION *****************************************
// ****************************************************************************
class iBinStream : public BinStream
{
public:
  iBinStream(void);
  ~iBinStream(void);

  void open( const char* filename, const bool truncate = false );
  void close(void);

  template<class TYPE>
  void read( TYPE* const p_var , const int var_length )
  {
    input.read( reinterpret_cast<char*>(p_var) , sizeof(TYPE)*var_length );
  }

  template<class TYPE>
  void read( TYPE& p_var )
  {
    this->read( &p_var , 1 );
  }

  void read( std::string& str )
  {
    std::string::size_type length;
    this->read( length );

    Array1dT<char> readstring;
    readstring.resize(length);
    this->read( readstring.data(), readstring.size() );

    str.assign( readstring.begin(),readstring.end() );

  }

  template< typename TYPE >
  void read( Array1dT<TYPE>& array, const bool realloc = true )
  {
    const typename Array1dT<TYPE>::size_type length = array.size();
    typename Array1dT<TYPE>::size_type readLength;
    this->read( readLength );

    if( readLength != length )
    {
     if( realloc )
     {
       array.resize(readLength);
     }
     else
     {
      throw GPException( "BinStream::read(Array1dT<TYPE>& array): length mismatch\n");
     }
    }

    this->read( array.data(), array.size() );
  }


  template< typename TYPE >
  void read( std::set<TYPE>& set, const bool realloc = true )
  {
    typename std::set<TYPE>::size_type readLength;
    this->read( readLength );

    if( readLength != set.size() && !realloc )
    {
      throw GPException( "BinStream::read(std::set<TYPE>& set): length mismatch\n");
    }

    for( typename std::set<TYPE>::size_type i=0 ; i<readLength ; ++i )
    {
      TYPE readVal;
      this->read( readVal );
      set.insert(readVal);
    }
  }

  template< typename TYPE >
  void read( Array2dT<TYPE>& array, const bool realloc = true )
  {
    const typename Array2dT<TYPE>::size_type length = array.size();
    typename Array2dT<TYPE>::size_type readLength;
    this->read( readLength );
    if( readLength != length )
      throw GPException( "BinStream::read(Array2dT<TYPE>& array): length mismatch\n");

    size_t dimension[2] = { array.Dimension(0), array.Dimension(1) };
    size_t readDimension[2];
    this->read( readDimension, 2 );

    if( dimension[0] != readDimension[0] || dimension[1] != readDimension[1] )
    {
      if( realloc )
      {
        array.resize2( readDimension[0], readDimension[1] );
      }
      else
      {
        throw GPException( "BinStream::read(Array2dT<TYPE>& array): dimension mismatch\n");
      }
    }

    this->read( array.data(), array.size() );
  }


  template< typename TYPE >
  void read( Array1dT<Array1dT<TYPE> >& array, const bool realloc = true )
  {
    const typename Array1dT<Array1dT<TYPE> >::size_type length = array.size();
    typename Array1dT<Array1dT<TYPE> >::size_type readLength;
    this->read( readLength );

    if( readLength != length )
    {
      if( realloc )
      {
        array.resize(readLength);
      }
      else
      {
        throw GPException( "BinStream::read(Array1dT<Array1dT<TYPE> >& array): length mismatch\n");
      }
    }

    for( typename Array1dT<Array1dT<TYPE> >::iterator i=array.begin() ; i!=array.end() ; ++i )
    {
      read( *i, realloc );
    }
  }


  template< typename TYPE >
  void read( Array1dT<std::set<TYPE> >& array, const bool realloc = true )
  {
    const typename Array1dT<std::set<TYPE> >::size_type length = array.size();
    typename Array1dT<std::set<TYPE> >::size_type readLength;
    this->read( readLength );

    if( readLength != length )
    {
      if( realloc )
      {
        array.resize(readLength);
      }
      else
      {
        throw GPException( "BinStream::read(Array1dT<std::set<TYPE> >& array): length mismatch\n");
      }
    }
    for( typename Array1dT<std::set<TYPE> >::iterator i=array.begin() ; i!=array.end() ; ++i )
    {
      this->read( *i, realloc );
    }
  }


  template< typename T1, typename T2 >
  void read( std::map< T1, T2 >& member )
  {
    typename std::map< T1, T2 >::size_type length;
    this->read( length );
    for( typename std::map< T1, T2 >::size_type i=0 ; i<length ; ++i  )
    {
      T1 key;
      T2 val;
      this->read( key );
      this->read( val );
      member[key] = val;
    }
  }


  template< typename T >
  void read( std::map< std::string, T>& member, const bool realloc = true )
  {

    typename std::map< std::string, T >::size_type numEntries = 0;
    this->read(numEntries);



    for( typename std::map< std::string, T >::size_type i=0 ; i<numEntries ; ++i )
    {
      // read the field name
      std::string readFieldName;
      this->read(readFieldName);

      std::map<std::string, FieldBase*>::iterator fieldAttributes = FieldInfo::AttributesByName.find(readFieldName);

      if( fieldAttributes == FieldInfo::AttributesByName.end()  )
      {
        if( realloc )
        {
          FieldInfo::AttributesByName[readFieldName]  = new FieldBase( FieldInfo::noKey, readFieldName, true, true );
          fieldAttributes = FieldInfo::AttributesByName.find(readFieldName);
        }
        else
        {
          throw GPException("ObjectDataStructureBaseT::ReadMapFromRestart: name not found\n");
        }
      }


//      std::cout<<"    reading "<<readFieldName<<std::endl;

      T& fieldData = member[readFieldName];

      // read the field data
      this->read( fieldData );


    }

    /*
    // iterate over all entries in the member map
    for( typename std::map< std::string, T >::iterator i = member.begin() ; i!=member.end() ; ++i )
    {
      // the field name is the key to the map
      const std::string fieldName = i->first;

      std::map<std::string, FieldBase*>::const_iterator fieldAttributes = FieldInfo::AttributesByName.find(fieldName);
      bool readField = false;
      if( fieldAttributes == FieldInfo::AttributesByName.end() )
      {
        readField = true;
      }
      else if( fieldAttributes->second->m_WriteToRestart )
      {
        readField = true;
      }

      if( readField )
      {
        std::cout<<"    reading "<<fieldName<<std::endl;
        // the field data is mapped value
        T& fieldData = i->second;

        // read the field name
        std::string readFieldName;
        this->read(readFieldName);

        if( fieldName != readFieldName )
          throw GPException("ObjectDataStructureBaseT::ReadMapFromRestart: field name mismatch\n");

        // write the field data
        this->read( fieldData );
      }
    }*/
  }

private:
  //***** Data Member Declarations **********************************************
  std::ifstream input;

public:
  std::ifstream::pos_type tellg(void)
  { return input.tellg(); }

};



#endif
