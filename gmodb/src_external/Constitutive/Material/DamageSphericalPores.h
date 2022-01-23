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

#ifndef DAMAGESPHERICALPORES_H_
#define DAMAGESPHERICALPORES_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/Material/LinearElasticIntermediate.h"

/*
 * DamageSphericalPores.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class DamageSphericalPoresParameterData : public LinearElasticIntermediateParameterData
{

public:

  typedef LinearElasticIntermediateParameterData base;
  realT f1;
  realT f2;
  realT initialVoidRatio;


  DamageSphericalPoresParameterData():
    base(),
    f1(0),
    f2(0),
    initialVoidRatio(0)
  {}

  DamageSphericalPoresParameterData( const DamageSphericalPoresParameterData& source):
    base( source ),
    f1(source.f1),
    f2(source.f2),
    initialVoidRatio(source.initialVoidRatio)
  {}

  ~DamageSphericalPoresParameterData() {}
  friend class ConstitutiveBase;
  friend class DamageSphericalPores;

  static void GetVariableCounts( localIndex& intVarCounts,
                                 localIndex& realVarCounts,
                                 localIndex& R1TensorVarCounts,
                                 localIndex& R2TensorVarCounts,
                                 localIndex& R2SymTensorVarCounts )
  {
    base::GetVariableCounts( intVarCounts,
                             realVarCounts,
                             R1TensorVarCounts,
                             R2TensorVarCounts,
                             R2SymTensorVarCounts );
    realVarCounts += 3;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("crackInitiationYield");
    realNames.push_back("crackPropgationYield");
    realNames.push_back("initialVoidRatio");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["crackInitiationYield"] = (char*)(&f1) - (char*)this;
    realOffsets["crackPropgationYield"] = (char*)(&f2) - (char*)this;
    realOffsets["initialVoidRatio"] = (char*)(&initialVoidRatio) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["crackInitiationYield"] = f1;
    realValues["crackPropgationYield"] = f2;
    realValues["initialVoidRatio"] = initialVoidRatio;
  }

  void Serialize(const localIndex index,
                  Array1dT<iArray1d*>& intVars,
                  Array1dT<rArray1d*>& realVars,
                  Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                  Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                  Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& intVarCounts,
                  localIndex& realVarCounts,
                  localIndex& R1TensorVarCounts,
                  localIndex& R2TensorVarCounts,
                  localIndex& R2SymTensorVarCounts  ) const
  {
    base::Serialize(index, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
    (*(realVars[realVarCounts]))[index] = f1; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = f2; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = initialVoidRatio; realVarCounts++;
  }


  void  Deserialize( const localIndex index,
                     const Array1dT<iArray1d*>& intVars,
                     const Array1dT<rArray1d*>& realVars,
                     const Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                     const Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                     const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& intVarCounts,
                  localIndex& realVarCounts,
                  localIndex& R1TensorVarCounts,
                  localIndex& R2TensorVarCounts,
                  localIndex& R2SymTensorVarCounts )
  {
    base::Deserialize(index, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
    f1 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    f2 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    initialVoidRatio = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline DamageSphericalPoresParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    f1 *= factor;
    f2 *= factor;
    initialVoidRatio *= factor;
    return *this;
  }

  inline DamageSphericalPoresParameterData&
  operator=(const DamageSphericalPoresParameterData& datum)
  {
    base::operator=(datum);
    f1 = datum.f1;
    f2 = datum.f2;
    initialVoidRatio = datum.initialVoidRatio;
    return *this;
  }

  inline DamageSphericalPoresParameterData&
  operator+=(const DamageSphericalPoresParameterData& datum)
  {
    base::operator+=(datum);
    f1 += datum.f1;
    f2 += datum.f2;
    initialVoidRatio += datum.initialVoidRatio;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, DamageSphericalPoresParameterData& p0, DamageSphericalPoresParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, f1, p0.f1, p1.f1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, f2, p0.f2, p1.f2);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, initialVoidRatio, p0.initialVoidRatio, p1.initialVoidRatio);

  }

  void MapFromRegion(const DamageSphericalPoresParameterData& p0, const DamageSphericalPoresParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.f1, p1.f1, fct0, fct1, f1);
    GeometryUtilities::MapFromRegion(p0.f2, p1.f2, fct0, fct1, f2);
    GeometryUtilities::MapFromRegion(p0.initialVoidRatio, p1.initialVoidRatio, fct0, fct1, initialVoidRatio);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    LinearElasticIntermediateParameterData::ReadXML( node );
    f1 = node.GetAttributeOrDefault("crackInitiationYield", 0.0);
    f2 = node.GetAttributeOrDefault("crackPropgationYield", 0.0);
    initialVoidRatio = node.GetAttributeOrDefault("initialVoidRatio", 0.0);
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class DamageSphericalPoresStateData : public LinearElasticIntermediateStateData
{

public:

  typedef LinearElasticIntermediateStateData base;
  realT radiusRatio;
  realT fractureNumberDensity;
  realT yAffinity;
  realT gAffinity;
  realT chordLengthMean;
  realT m01;
  realT m02;
  realT m11;
  realT m12;
  realT m21;
  realT m22;


  DamageSphericalPoresStateData():
    base(),
    radiusRatio(0),
    fractureNumberDensity(0),
    yAffinity(0),
    gAffinity(0),
    chordLengthMean(0),
    m01(0),
    m02(0),
    m11(0),
    m12(0),
    m21(0),
    m22(0)
  {}

  DamageSphericalPoresStateData( const DamageSphericalPoresStateData& source):
    base( source ),
    radiusRatio(source.radiusRatio),
    fractureNumberDensity(source.fractureNumberDensity),
    yAffinity(source.yAffinity),
    gAffinity(source.gAffinity),
    chordLengthMean(source.chordLengthMean),
    m01(source.m01),
    m02(source.m02),
    m11(source.m11),
    m12(source.m12),
    m21(source.m21),
    m22(source.m22)
  {}

  ~DamageSphericalPoresStateData() {}
  friend class ConstitutiveBase;
  friend class DamageSphericalPores;

  static void GetVariableCounts( localIndex& intVarCounts,
                                 localIndex& realVarCounts,
                                 localIndex& R1TensorVarCounts,
                                 localIndex& R2TensorVarCounts,
                                 localIndex& R2SymTensorVarCounts )
  {
    base::GetVariableCounts( intVarCounts,
                             realVarCounts,
                             R1TensorVarCounts,
                             R2TensorVarCounts,
                             R2SymTensorVarCounts );
    realVarCounts += 11;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("radiusRatio");
    realNames.push_back("fractureNumberDensity");
    realNames.push_back("yAffinity");
    realNames.push_back("gAffinity");
    realNames.push_back("chordLengthVariance");
    realNames.push_back("M01");
    realNames.push_back("M02");
    realNames.push_back("M11");
    realNames.push_back("M12");
    realNames.push_back("M21");
    realNames.push_back("M22");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["radiusRatio"] = (char*)(&radiusRatio) - (char*)this;
    realOffsets["fractureNumberDensity"] = (char*)(&fractureNumberDensity) - (char*)this;
    realOffsets["yAffinity"] = (char*)(&yAffinity) - (char*)this;
    realOffsets["gAffinity"] = (char*)(&gAffinity) - (char*)this;
    realOffsets["chordLengthVariance"] = (char*)(&chordLengthMean) - (char*)this;
    realOffsets["M01"] = (char*)(&m01) - (char*)this;
    realOffsets["M02"] = (char*)(&m02) - (char*)this;
    realOffsets["M11"] = (char*)(&m11) - (char*)this;
    realOffsets["M12"] = (char*)(&m12) - (char*)this;
    realOffsets["M21"] = (char*)(&m21) - (char*)this;
    realOffsets["M22"] = (char*)(&m22) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["radiusRatio"] = radiusRatio;
    realValues["fractureNumberDensity"] = fractureNumberDensity;
    realValues["yAffinity"] = yAffinity;
    realValues["gAffinity"] = gAffinity;
    realValues["chordLengthVariance"] = chordLengthMean;
    realValues["M01"] = m01;
    realValues["M02"] = m02;
    realValues["M11"] = m11;
    realValues["M12"] = m12;
    realValues["M21"] = m21;
    realValues["M22"] = m22;
  }

  void Serialize(const localIndex index,
                  const unsigned int stride,
                  const localIndex elemNum,
                  Array1dT<iArray1d*>& intVars,
                  Array1dT<rArray1d*>& realVars,
                  Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                  Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                  Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& intVarCounts,
                  localIndex& realVarCounts,
                  localIndex& R1TensorVarCounts,
                  localIndex& R2TensorVarCounts,
                  localIndex& R2SymTensorVarCounts  ) const
  {
    base::Serialize(index, stride, elemNum, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
    (*(realVars[realVarCounts]))[elemNum] = radiusRatio; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = fractureNumberDensity; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = yAffinity; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = gAffinity; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = chordLengthMean; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = m01; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = m02; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = m11; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = m12; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = m21; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = m22; realVarCounts += stride;
  }


  void  Deserialize( const localIndex index,
                  const unsigned int stride,
                  const localIndex elemNum,
                     const Array1dT<iArray1d*>& intVars,
                     const Array1dT<rArray1d*>& realVars,
                     const Array1dT<Array1dT<R1Tensor>*>& R1Vars,
                     const Array1dT<Array1dT<R2Tensor>*>& R2Vars,
                     const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars,
                  localIndex& intVarCounts,
                  localIndex& realVarCounts,
                  localIndex& R1TensorVarCounts,
                  localIndex& R2TensorVarCounts,
                  localIndex& R2SymTensorVarCounts )
  {
    base::Deserialize(index, stride, elemNum, intVars, realVars, R1Vars, R2Vars, R2SymVars,
                    intVarCounts, realVarCounts, R1TensorVarCounts, R2TensorVarCounts, R2SymTensorVarCounts );
    radiusRatio = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    fractureNumberDensity = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    yAffinity = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    gAffinity = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    chordLengthMean = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    m01 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    m02 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    m11 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    m12 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    m21 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    m22 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
  }
  inline DamageSphericalPoresStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    radiusRatio *= factor;
    fractureNumberDensity *= factor;
    yAffinity *= factor;
    gAffinity *= factor;
    chordLengthMean *= factor;
    m01 *= factor;
    m02 *= factor;
    m11 *= factor;
    m12 *= factor;
    m21 *= factor;
    m22 *= factor;
    return *this;
  }

  inline DamageSphericalPoresStateData&
  operator=(const DamageSphericalPoresStateData& datum)
  {
    base::operator=(datum);
    radiusRatio = datum.radiusRatio;
    fractureNumberDensity = datum.fractureNumberDensity;
    yAffinity = datum.yAffinity;
    gAffinity = datum.gAffinity;
    chordLengthMean = datum.chordLengthMean;
    m01 = datum.m01;
    m02 = datum.m02;
    m11 = datum.m11;
    m12 = datum.m12;
    m21 = datum.m21;
    m22 = datum.m22;
    return *this;
  }

  inline DamageSphericalPoresStateData&
  operator+=(const DamageSphericalPoresStateData& datum)
  {
    base::operator+=(datum);
    radiusRatio += datum.radiusRatio;
    fractureNumberDensity += datum.fractureNumberDensity;
    yAffinity += datum.yAffinity;
    gAffinity += datum.gAffinity;
    chordLengthMean += datum.chordLengthMean;
    m01 += datum.m01;
    m02 += datum.m02;
    m11 += datum.m11;
    m12 += datum.m12;
    m21 += datum.m21;
    m22 += datum.m22;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, DamageSphericalPoresStateData& p0, DamageSphericalPoresStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, radiusRatio, p0.radiusRatio, p1.radiusRatio);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, fractureNumberDensity, p0.fractureNumberDensity, p1.fractureNumberDensity);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, yAffinity, p0.yAffinity, p1.yAffinity);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, gAffinity, p0.gAffinity, p1.gAffinity);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, chordLengthMean, p0.chordLengthMean, p1.chordLengthMean);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, m01, p0.m01, p1.m01);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, m02, p0.m02, p1.m02);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, m11, p0.m11, p1.m11);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, m12, p0.m12, p1.m12);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, m21, p0.m21, p1.m21);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, m22, p0.m22, p1.m22);

  }

  void MapFromRegion(const DamageSphericalPoresStateData& p0, const DamageSphericalPoresStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.radiusRatio, p1.radiusRatio, fct0, fct1, radiusRatio);
    GeometryUtilities::MapFromRegion(p0.fractureNumberDensity, p1.fractureNumberDensity, fct0, fct1, fractureNumberDensity);
    GeometryUtilities::MapFromRegion(p0.yAffinity, p1.yAffinity, fct0, fct1, yAffinity);
    GeometryUtilities::MapFromRegion(p0.gAffinity, p1.gAffinity, fct0, fct1, gAffinity);
    GeometryUtilities::MapFromRegion(p0.chordLengthMean, p1.chordLengthMean, fct0, fct1, chordLengthMean);
    GeometryUtilities::MapFromRegion(p0.m01, p1.m01, fct0, fct1, m01);
    GeometryUtilities::MapFromRegion(p0.m02, p1.m02, fct0, fct1, m02);
    GeometryUtilities::MapFromRegion(p0.m11, p1.m11, fct0, fct1, m11);
    GeometryUtilities::MapFromRegion(p0.m12, p1.m12, fct0, fct1, m12);
    GeometryUtilities::MapFromRegion(p0.m21, p1.m21, fct0, fct1, m21);
    GeometryUtilities::MapFromRegion(p0.m22, p1.m22, fct0, fct1, m22);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class DamageSphericalPores: public LinearElasticIntermediate
{
public:

  typedef DamageSphericalPoresParameterData ParameterClass;
  typedef DamageSphericalPoresStateData     StateClass;

  typedef Array1dT<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "DamageSphericalPores"; }

  DamageSphericalPores();
  virtual ~DamageSphericalPores();

  inline void MapToRegion(const realT fctNormal, const realT fct0, const realT fct1, 
                          const localIndex from0, const localIndex from1,
                          StateClass& s0, StateClass& s1)
  {
    StateData(from0, from1)->MapToRegion(fctNormal, fct0, fct1, s0, s1);
    //ParameterData(from)->MapToRegion(fctNormal, fct0, fct1, p0, p1);
  }

  inline void MapFromRegion(const realT fct0, const realT fct1, 
                          const StateClass& s0, const StateClass& s1, 
                          const localIndex to0, const localIndex to1)
  {
    StateData(to0, to1)->MapFromRegion(s0, s1, fct0, fct1);
    //ParameterData(to)->MapFromRegion(p0, p1, fct0, fct1);
  }
 
  inline void MapToRegion(const realT fctNormal, const realT fct0, const realT fct1, 
                          const localIndex from,
                          ParameterClass& p0, ParameterClass& p1)
  {
    //StateData(from0, from1)->MapToRegion(fctNormal, fct0, fct1, s0, s1);
    ParameterData(from)->MapToRegion(fctNormal, fct0, fct1, p0, p1);
  }

  inline void MapFromRegion(const realT fct0, const realT fct1, 
                          const ParameterClass& p0, const ParameterClass& p1, 
                          const localIndex to)
  {
    //StateData(to0, to1)->MapFromRegion(s0, s1, fct0, fct1);
    ParameterData(to)->MapFromRegion(p0, p1, fct0, fct1);
  }
  
  
  virtual void ZeroStates()
  {
    for(localIndex j = 0; j < m_stateData.Dimension(1); j++)
    {
      for(localIndex i = 0; i < m_stateData.Dimension(0); i++)
      {
        m_stateData(i,j) *= 0.0;
      }
    }
  }
  
  virtual void SetVariableParameters(const bool varParams, const localIndex newSize = 0)
  { SetVariableParametersFromDerived<DamageSphericalPores>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<DamageSphericalPores>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<DamageSphericalPores>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<DamageSphericalPores>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<DamageSphericalPores>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<DamageSphericalPores>( num ); }
 
  void GetVariableNames( sArray1d& intVars, sArray1d& realVars, sArray1d& R1TensorVars, sArray1d& R2TensorVars, sArray1d& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<DamageSphericalPores>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<DamageSphericalPores>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<DamageSphericalPores>(name, type ); }
  
  bool GetStateValues( const std::string& name, rArray1d& values ) const
  { return GetStateValuesFromDerived<DamageSphericalPores>(name, values); }

  bool GetParameterValues( const std::string& name, rArray1d& values ) const
  { return GetParameterValuesFromDerived<DamageSphericalPores>(name, values); }

  bool SetStateValues( const std::string& name, const rArray1d& values )
  { return SetStateValuesFromDerived<DamageSphericalPores>(name, values); }

  bool SetParameterValues( const std::string& name, const rArray1d& values )
  { return SetParameterValuesFromDerived<DamageSphericalPores>(name, values); }

  virtual void Serialize( Array1dT<iArray1d*>& intVars, Array1dT<rArray1d*>& realVars, Array1dT<Array1dT<R1Tensor>*>& R1Vars, Array1dT<Array1dT<R2Tensor>*>& R2Vars, Array1dT<Array1dT<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<DamageSphericalPores>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const Array1dT<iArray1d*>& intVars, const Array1dT<rArray1d*>& realVars, const Array1dT<Array1dT<R1Tensor>*>& R1Vars, const Array1dT<Array1dT<R2Tensor>*>& R2Vars, const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<DamageSphericalPores>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<DamageSphericalPores>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<DamageSphericalPores>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<DamageSphericalPores>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<DamageSphericalPores>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<DamageSphericalPores>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<DamageSphericalPores>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  


private:
  DamageSphericalPores(const DamageSphericalPores&);
  DamageSphericalPores& operator=(const DamageSphericalPores&);
  

};
#endif /* DAMAGESPHERICALPORES_H_ */
