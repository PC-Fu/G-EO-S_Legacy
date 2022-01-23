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
 * @file IOUtilities.h
 * @author walsh24
 * @date Oct 13, 2011
 */

#ifndef IOUTILITIES_H_
#define IOUTILITIES_H_

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <set>

#include "StringUtilities.h"
#include "Common/GPException.h"


/// Read from a character deliminated file into a vector
template<class ARRAY>
void dlmreadVector(const std::string& filename, ARRAY& values, char delim = ' ', int skipLines = 0)
{

  std::ifstream inputStream(filename.c_str());
  if (inputStream)
  {
    std::string lineString;

    for (int i = 0; i < skipLines; ++i)
      std::getline(inputStream, lineString);

    // read every line from the stream
    while (std::getline(inputStream, lineString))
    {
      std::istringstream lineStream(lineString);
      std::string cellString;
      // read every element from the line separated by the delimiter
      while (std::getline(lineStream, cellString, delim))
      {
        values.push_back(fromString<typename ARRAY::value_type>(cellString));
      }
    }
    inputStream.close();
  }
  else
  {
    throw GPException("dlmreadVector: Failed to load file:" + filename + " \n");
  }
}


///////////////////////////////////////////////////
//
///// Read from a character deliminated file into a vector of vectors
template<class ARRAY>
void dlmreadArray(const std::string& filename, ARRAY& values, char delim = ' ', int skipLines = 0)
{

  std::ifstream inputStream(filename.c_str());
  if (inputStream)
  {
    std::string lineString;

    for (int i = 0; i < skipLines; ++i)
      std::getline(inputStream, lineString);

    // read every line from the stream
    while (std::getline(inputStream, lineString))
    {
      std::istringstream lineStream(lineString);
      std::vector < typename ARRAY::value_type > row;
      std::string cellString;
      // read every element from the line separated by the delimiter
      while (std::getline(lineStream, cellString, delim))
      {
        row.push_back(fromString<typename ARRAY::value_type>(cellString));
      }
      values.push_back(row);
    }
    inputStream.close();
  }
  else
  {
    throw GPException("dlmreadArray: Failed to load file:" + filename + " \n");
  }
}

///// Read from a character deliminated file into a vector of vectors and return the transpose
///// TODO: Could be improved - not particularly fast or memory efficient
template<class ARRAY>
void dlmreadArrayTranspose(const std::string& filename, ARRAY& values, char delim = ' ', int skipLines = 0)
{

  std::ifstream inputStream(filename.c_str());
  std::string lineString;

  if (inputStream)
  {
    for (int i = 0; i < skipLines; ++i)
      std::getline(inputStream, lineString);

    std::vector < std::vector<typename ARRAY::value_type> > valuesT;

    // read every line from the stream
    while (std::getline(inputStream, lineString))
    {
      std::istringstream lineStream(lineString);
      std::vector < typename ARRAY::value_type > row;
      std::string cellString;
      // read every element from the line separated by the delimiter
      while (std::getline(lineStream, cellString, delim))
      {
        row.push_back(fromString<typename ARRAY::value_type>(cellString));
      }
      valuesT.push_back(row);
    }
    inputStream.close();

    localIndex n = valuesT.size();
    localIndex m = valuesT[0].size();
    values.resize(m);
    for (localIndex i = 0; i < m; ++i)
    {
      values[i].resize(n);
      for (localIndex j = 0; j < n; ++j)
        values[i][j] = valuesT[j][i];
    }
  }
  else
  {
    throw GPException("dlmreadArrayTranspose: Failed to load file:" + filename + " \n");
  }
}

#endif /*IOUTILITIES_H_*/
