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

#ifndef DSIDDAMAGE_H_
#define DSIDDAMAGE_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/Material/MaterialBase.h"

/*
 * DSIDDamage.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class DSIDDamageParameterData : public MaterialBaseParameterData
{

public:

  typedef MaterialBaseParameterData base;
  realT alpha;
  realT a1;
  realT a2;
  realT a3;
  realT a4;
  realT c0;
  realT c1;
  R2SymTensor omega0;


  DSIDDamageParameterData():
    base(),
    alpha(0),
    a1(0),
    a2(0),
    a3(0),
    a4(0),
    c0(0),
    c1(0),
    omega0(0)
  {}

  DSIDDamageParameterData( const DSIDDamageParameterData& source):
    base( source ),
    alpha(source.alpha),
    a1(source.a1),
    a2(source.a2),
    a3(source.a3),
    a4(source.a4),
    c0(source.c0),
    c1(source.c1),
    omega0(source.omega0)
  {}

  ~DSIDDamageParameterData() {}
  friend class ConstitutiveBase;
  friend class DSIDDamage;

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
    realVarCounts += 7;
    R2SymTensorVarCounts += 1;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("alpha");
    realNames.push_back("parameter1");
    realNames.push_back("parameter2");
    realNames.push_back("parameter3");
    realNames.push_back("parameter4");
    realNames.push_back("initialDamageThreshold");
    realNames.push_back("damageHardeningVariable");
    R2SymTensorNames.push_back("initialDamage");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["alpha"] = (char*)(&alpha) - (char*)this;
    realOffsets["parameter1"] = (char*)(&a1) - (char*)this;
    realOffsets["parameter2"] = (char*)(&a2) - (char*)this;
    realOffsets["parameter3"] = (char*)(&a3) - (char*)this;
    realOffsets["parameter4"] = (char*)(&a4) - (char*)this;
    realOffsets["initialDamageThreshold"] = (char*)(&c0) - (char*)this;
    realOffsets["damageHardeningVariable"] = (char*)(&c1) - (char*)this;
    R2SymTensorOffsets["initialDamage"] = (char*)(&omega0) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["alpha"] = alpha;
    realValues["parameter1"] = a1;
    realValues["parameter2"] = a2;
    realValues["parameter3"] = a3;
    realValues["parameter4"] = a4;
    realValues["initialDamageThreshold"] = c0;
    realValues["damageHardeningVariable"] = c1;
    R2SymTensorValues["initialDamage"] = omega0;
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
    (*(realVars[realVarCounts]))[index] = alpha; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = a1; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = a2; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = a3; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = a4; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c0; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c1; realVarCounts++;
    (*(R2SymVars[R2SymTensorVarCounts]))[index] = omega0; R2SymTensorVarCounts++;
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
    alpha = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    a1 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    a2 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    a3 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    a4 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c0 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c1 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    omega0 = (*(R2SymVars[R2SymTensorVarCounts]))[index]; R2SymTensorVarCounts++;
  }
  inline DSIDDamageParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    alpha *= factor;
    a1 *= factor;
    a2 *= factor;
    a3 *= factor;
    a4 *= factor;
    c0 *= factor;
    c1 *= factor;
    omega0 *= factor;
    return *this;
  }

  inline DSIDDamageParameterData&
  operator=(const DSIDDamageParameterData& datum)
  {
    base::operator=(datum);
    alpha = datum.alpha;
    a1 = datum.a1;
    a2 = datum.a2;
    a3 = datum.a3;
    a4 = datum.a4;
    c0 = datum.c0;
    c1 = datum.c1;
    omega0 = datum.omega0;
    return *this;
  }

  inline DSIDDamageParameterData&
  operator+=(const DSIDDamageParameterData& datum)
  {
    base::operator+=(datum);
    alpha += datum.alpha;
    a1 += datum.a1;
    a2 += datum.a2;
    a3 += datum.a3;
    a4 += datum.a4;
    c0 += datum.c0;
    c1 += datum.c1;
    omega0 += datum.omega0;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, DSIDDamageParameterData& p0, DSIDDamageParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, alpha, p0.alpha, p1.alpha);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, a1, p0.a1, p1.a1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, a2, p0.a2, p1.a2);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, a3, p0.a3, p1.a3);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, a4, p0.a4, p1.a4);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c0, p0.c0, p1.c0);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c1, p0.c1, p1.c1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, omega0, p0.omega0, p1.omega0);

  }

  void MapFromRegion(const DSIDDamageParameterData& p0, const DSIDDamageParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.alpha, p1.alpha, fct0, fct1, alpha);
    GeometryUtilities::MapFromRegion(p0.a1, p1.a1, fct0, fct1, a1);
    GeometryUtilities::MapFromRegion(p0.a2, p1.a2, fct0, fct1, a2);
    GeometryUtilities::MapFromRegion(p0.a3, p1.a3, fct0, fct1, a3);
    GeometryUtilities::MapFromRegion(p0.a4, p1.a4, fct0, fct1, a4);
    GeometryUtilities::MapFromRegion(p0.c0, p1.c0, fct0, fct1, c0);
    GeometryUtilities::MapFromRegion(p0.c1, p1.c1, fct0, fct1, c1);
    GeometryUtilities::MapFromRegion(p0.omega0, p1.omega0, fct0, fct1, omega0);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    MaterialBaseParameterData::ReadXML( node );
    alpha = node.GetAttributeOrDefault("alpha", 0.0);
    a1 = node.GetAttributeOrDefault("parameter1", 0.0);
    a2 = node.GetAttributeOrDefault("parameter2", 0.0);
    a3 = node.GetAttributeOrDefault("parameter3", 0.0);
    a4 = node.GetAttributeOrDefault("parameter4", 0.0);
    c0 = node.GetAttributeOrDefault("initialDamageThreshold", 0.0);
    c1 = node.GetAttributeOrDefault("damageHardeningVariable", 0.0);
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class DSIDDamageStateData : public MaterialBaseStateData
{

public:

  typedef MaterialBaseStateData base;
  realT fd0;
  R2SymTensor omega;
  R2SymTensor epsilon;
  R2SymTensor epsid;
  R2SymTensor epsE;
  R2SymTensor epsel;


  DSIDDamageStateData():
    base(),
    fd0(0),
    omega(0),
    epsilon(0),
    epsid(0),
    epsE(0),
    epsel(0)
  {}

  DSIDDamageStateData( const DSIDDamageStateData& source):
    base( source ),
    fd0(source.fd0),
    omega(source.omega),
    epsilon(source.epsilon),
    epsid(source.epsid),
    epsE(source.epsE),
    epsel(source.epsel)
  {}

  ~DSIDDamageStateData() {}
  friend class ConstitutiveBase;
  friend class DSIDDamage;

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
    realVarCounts += 1;
    R2SymTensorVarCounts += 5;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("DamageIndicator");
    R2SymTensorNames.push_back("currentDamage");
    R2SymTensorNames.push_back("strain");
    R2SymTensorNames.push_back("IrreversibleStrain");
    R2SymTensorNames.push_back("TotalElasticStrain");
    R2SymTensorNames.push_back("pureElasticStrain");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["DamageIndicator"] = (char*)(&fd0) - (char*)this;
    R2SymTensorOffsets["currentDamage"] = (char*)(&omega) - (char*)this;
    R2SymTensorOffsets["strain"] = (char*)(&epsilon) - (char*)this;
    R2SymTensorOffsets["IrreversibleStrain"] = (char*)(&epsid) - (char*)this;
    R2SymTensorOffsets["TotalElasticStrain"] = (char*)(&epsE) - (char*)this;
    R2SymTensorOffsets["pureElasticStrain"] = (char*)(&epsel) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["DamageIndicator"] = fd0;
    R2SymTensorValues["currentDamage"] = omega;
    R2SymTensorValues["strain"] = epsilon;
    R2SymTensorValues["IrreversibleStrain"] = epsid;
    R2SymTensorValues["TotalElasticStrain"] = epsE;
    R2SymTensorValues["pureElasticStrain"] = epsel;
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
    (*(realVars[realVarCounts]))[elemNum] = fd0; realVarCounts += stride;
    (*(R2SymVars[R2SymTensorVarCounts]))[elemNum] = omega; R2SymTensorVarCounts += stride;
    (*(R2SymVars[R2SymTensorVarCounts]))[elemNum] = epsilon; R2SymTensorVarCounts += stride;
    (*(R2SymVars[R2SymTensorVarCounts]))[elemNum] = epsid; R2SymTensorVarCounts += stride;
    (*(R2SymVars[R2SymTensorVarCounts]))[elemNum] = epsE; R2SymTensorVarCounts += stride;
    (*(R2SymVars[R2SymTensorVarCounts]))[elemNum] = epsel; R2SymTensorVarCounts += stride;
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
    fd0 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    omega = (*(R2SymVars[R2SymTensorVarCounts]))[elemNum]; R2SymTensorVarCounts += stride;
    epsilon = (*(R2SymVars[R2SymTensorVarCounts]))[elemNum]; R2SymTensorVarCounts += stride;
    epsid = (*(R2SymVars[R2SymTensorVarCounts]))[elemNum]; R2SymTensorVarCounts += stride;
    epsE = (*(R2SymVars[R2SymTensorVarCounts]))[elemNum]; R2SymTensorVarCounts += stride;
    epsel = (*(R2SymVars[R2SymTensorVarCounts]))[elemNum]; R2SymTensorVarCounts += stride;
  }
  inline DSIDDamageStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    fd0 *= factor;
    omega *= factor;
    epsilon *= factor;
    epsid *= factor;
    epsE *= factor;
    epsel *= factor;
    return *this;
  }

  inline DSIDDamageStateData&
  operator=(const DSIDDamageStateData& datum)
  {
    base::operator=(datum);
    fd0 = datum.fd0;
    omega = datum.omega;
    epsilon = datum.epsilon;
    epsid = datum.epsid;
    epsE = datum.epsE;
    epsel = datum.epsel;
    return *this;
  }

  inline DSIDDamageStateData&
  operator+=(const DSIDDamageStateData& datum)
  {
    base::operator+=(datum);
    fd0 += datum.fd0;
    omega += datum.omega;
    epsilon += datum.epsilon;
    epsid += datum.epsid;
    epsE += datum.epsE;
    epsel += datum.epsel;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, DSIDDamageStateData& p0, DSIDDamageStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, fd0, p0.fd0, p1.fd0);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, omega, p0.omega, p1.omega);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, epsilon, p0.epsilon, p1.epsilon);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, epsid, p0.epsid, p1.epsid);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, epsE, p0.epsE, p1.epsE);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, epsel, p0.epsel, p1.epsel);

  }

  void MapFromRegion(const DSIDDamageStateData& p0, const DSIDDamageStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.fd0, p1.fd0, fct0, fct1, fd0);
    GeometryUtilities::MapFromRegion(p0.omega, p1.omega, fct0, fct1, omega);
    GeometryUtilities::MapFromRegion(p0.epsilon, p1.epsilon, fct0, fct1, epsilon);
    GeometryUtilities::MapFromRegion(p0.epsid, p1.epsid, fct0, fct1, epsid);
    GeometryUtilities::MapFromRegion(p0.epsE, p1.epsE, fct0, fct1, epsE);
    GeometryUtilities::MapFromRegion(p0.epsel, p1.epsel, fct0, fct1, epsel);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class DSIDDamage: public MaterialBase
{
public:

  typedef DSIDDamageParameterData ParameterClass;
  typedef DSIDDamageStateData     StateClass;

  typedef Array1dT<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "DSIDDamage"; }

  DSIDDamage();
  virtual ~DSIDDamage();

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
  { SetVariableParametersFromDerived<DSIDDamage>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<DSIDDamage>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<DSIDDamage>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<DSIDDamage>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<DSIDDamage>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<DSIDDamage>( num ); }
 
  void GetVariableNames( sArray1d& intVars, sArray1d& realVars, sArray1d& R1TensorVars, sArray1d& R2TensorVars, sArray1d& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<DSIDDamage>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<DSIDDamage>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<DSIDDamage>(name, type ); }
  
  bool GetStateValues( const std::string& name, rArray1d& values ) const
  { return GetStateValuesFromDerived<DSIDDamage>(name, values); }

  bool GetParameterValues( const std::string& name, rArray1d& values ) const
  { return GetParameterValuesFromDerived<DSIDDamage>(name, values); }

  bool SetStateValues( const std::string& name, const rArray1d& values )
  { return SetStateValuesFromDerived<DSIDDamage>(name, values); }

  bool SetParameterValues( const std::string& name, const rArray1d& values )
  { return SetParameterValuesFromDerived<DSIDDamage>(name, values); }

  virtual void Serialize( Array1dT<iArray1d*>& intVars, Array1dT<rArray1d*>& realVars, Array1dT<Array1dT<R1Tensor>*>& R1Vars, Array1dT<Array1dT<R2Tensor>*>& R2Vars, Array1dT<Array1dT<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<DSIDDamage>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const Array1dT<iArray1d*>& intVars, const Array1dT<rArray1d*>& realVars, const Array1dT<Array1dT<R1Tensor>*>& R1Vars, const Array1dT<Array1dT<R2Tensor>*>& R2Vars, const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<DSIDDamage>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<DSIDDamage>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<DSIDDamage>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<DSIDDamage>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<DSIDDamage>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<DSIDDamage>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<DSIDDamage>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  
  void
  InitializeStates( const localIndex index );

  void
  StrainDrivenUpdateMember( const localIndex index0,
                              const localIndex index1,
                              const R2SymTensorT < 3 >& Ddt,
                              const R2TensorT < 3 >& L,
                              const R2Tensor& Rot,
                              const realT& volume_n,
                              const realT& volume_np1,
                              const realT dt );

  void
  StrainDrivenUpdateMember( const localIndex index0,
                                        const localIndex index1,
                                        const R2SymTensorT < 3 >& Ddt,
                                        const R2TensorT < 3 >& L,
                                        const R2Tensor& Rot,
                                        const realT dt );

  void
  EffectiveElasticStiffness(const realT Nu0, const realT E0,
                                        const realT a1, const realT a2,
                                        const realT a3,const realT a4,
                                        const R2SymTensorT < 3 > E,
                                        R2SymTensorT < 3 > omega,
                                        R4minSymTensorT < 3 >& De,
                                        R4minSymTensorT < 3 >& De0,
                                        R4minSymTensorT < 3 >& Se,
                                        R4minSymTensorT < 3 >& Se0);

  void
  DamageDrivingForce(const realT a1, const realT a2,
                                 const realT a3,const realT a4,
                                 const R2SymTensorT < 3 > E,
                                 R2SymTensorT < 3 > stresstotal,
                                 R2SymTensorT < 3 >& Yd);

  void
  ProjectionTensorP1(R2SymTensorT < 3 > stresstotal, R4minSymTensorT < 3 >& P1);

  inline realT
  DamageFunction(const realT a1, const realT a2,
                             const realT a3, const realT a4,
                             const realT c0, const realT c1, const realT alpha,
                             const R2SymTensorT < 3 > E,
                             R2SymTensorT < 3 > stresstotal,
                             R2SymTensorT < 3 > omega);

  void
  ProjectionTensorP2(R2SymTensorT < 3 > stresstotal, R4minSymTensorT < 3 >& P2);

  void
  Others(const realT Nu0, const realT E0, const realT a1, const realT a2, const realT a3, const realT a4,
                     const realT c0, const realT c1, const realT alpha, const R2SymTensorT < 3 > E,
                     R2SymTensorT < 3 > epsilon, R2SymTensorT < 3 > epsid,
                     R2SymTensorT < 3 > stress, R2SymTensorT < 3 > omega, realT& lambda, realT& fd,R2SymTensorT < 3 > depsilon,
                     R2SymTensorT < 3 >& depsid, R2SymTensorT < 3 >& omegaT, R2SymTensorT < 3 >& stressTT);



private:
  DSIDDamage(const DSIDDamage&);
  DSIDDamage& operator=(const DSIDDamage&);
  

};
#endif /* DSIDDAMAGE_H_ */
