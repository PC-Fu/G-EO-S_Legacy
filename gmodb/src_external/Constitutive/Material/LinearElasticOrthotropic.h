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

#ifndef LINEARELASTICORTHOTROPIC_H_
#define LINEARELASTICORTHOTROPIC_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/Material/MaterialBase.h"

/*
 * LinearElasticOrthotropic.h
 *
 *  Created on: Wed Jun 25 13:21:30 PDT 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticOrthotropicParameterData : public MaterialBaseParameterData
{

public:

  typedef MaterialBaseParameterData base;
  realT c11;
  realT c22;
  realT c33;
  realT c12;
  realT c13;
  realT c23;
  realT c44;
  realT c55;
  realT c66;


  LinearElasticOrthotropicParameterData():
    base(),
    c11(0),
    c22(0),
    c33(0),
    c12(0),
    c13(0),
    c23(0),
    c44(0),
    c55(0),
    c66(0)
  {}

  LinearElasticOrthotropicParameterData( const LinearElasticOrthotropicParameterData& source):
    base( source ),
    c11(source.c11),
    c22(source.c22),
    c33(source.c33),
    c12(source.c12),
    c13(source.c13),
    c23(source.c23),
    c44(source.c44),
    c55(source.c55),
    c66(source.c66)
  {}

  ~LinearElasticOrthotropicParameterData() {}
  friend class ConstitutiveBase;
  friend class LinearElasticOrthotropic;

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
    realVarCounts += 9;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("c11");
    realNames.push_back("c22");
    realNames.push_back("c33");
    realNames.push_back("c12");
    realNames.push_back("c13");
    realNames.push_back("c23");
    realNames.push_back("c44");
    realNames.push_back("c55");
    realNames.push_back("c66");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["c11"] = (char*)(&c11) - (char*)this;
    realOffsets["c22"] = (char*)(&c22) - (char*)this;
    realOffsets["c33"] = (char*)(&c33) - (char*)this;
    realOffsets["c12"] = (char*)(&c12) - (char*)this;
    realOffsets["c13"] = (char*)(&c13) - (char*)this;
    realOffsets["c23"] = (char*)(&c23) - (char*)this;
    realOffsets["c44"] = (char*)(&c44) - (char*)this;
    realOffsets["c55"] = (char*)(&c55) - (char*)this;
    realOffsets["c66"] = (char*)(&c66) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["c11"] = c11;
    realValues["c22"] = c22;
    realValues["c33"] = c33;
    realValues["c12"] = c12;
    realValues["c13"] = c13;
    realValues["c23"] = c23;
    realValues["c44"] = c44;
    realValues["c55"] = c55;
    realValues["c66"] = c66;
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
    (*(realVars[realVarCounts]))[index] = c11; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c22; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c33; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c12; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c13; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c23; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c44; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c55; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c66; realVarCounts++;
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
    c11 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c22 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c33 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c12 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c13 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c23 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c44 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c55 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c66 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline LinearElasticOrthotropicParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    c11 *= factor;
    c22 *= factor;
    c33 *= factor;
    c12 *= factor;
    c13 *= factor;
    c23 *= factor;
    c44 *= factor;
    c55 *= factor;
    c66 *= factor;
    return *this;
  }

  inline LinearElasticOrthotropicParameterData&
  operator=(const LinearElasticOrthotropicParameterData& datum)
  {
    base::operator=(datum);
    c11 = datum.c11;
    c22 = datum.c22;
    c33 = datum.c33;
    c12 = datum.c12;
    c13 = datum.c13;
    c23 = datum.c23;
    c44 = datum.c44;
    c55 = datum.c55;
    c66 = datum.c66;
    return *this;
  }

  inline LinearElasticOrthotropicParameterData&
  operator+=(const LinearElasticOrthotropicParameterData& datum)
  {
    base::operator+=(datum);
    c11 += datum.c11;
    c22 += datum.c22;
    c33 += datum.c33;
    c12 += datum.c12;
    c13 += datum.c13;
    c23 += datum.c23;
    c44 += datum.c44;
    c55 += datum.c55;
    c66 += datum.c66;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearElasticOrthotropicParameterData& p0, LinearElasticOrthotropicParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c11, p0.c11, p1.c11);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c22, p0.c22, p1.c22);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c33, p0.c33, p1.c33);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c12, p0.c12, p1.c12);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c13, p0.c13, p1.c13);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c23, p0.c23, p1.c23);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c44, p0.c44, p1.c44);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c55, p0.c55, p1.c55);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c66, p0.c66, p1.c66);

  }

  void MapFromRegion(const LinearElasticOrthotropicParameterData& p0, const LinearElasticOrthotropicParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.c11, p1.c11, fct0, fct1, c11);
    GeometryUtilities::MapFromRegion(p0.c22, p1.c22, fct0, fct1, c22);
    GeometryUtilities::MapFromRegion(p0.c33, p1.c33, fct0, fct1, c33);
    GeometryUtilities::MapFromRegion(p0.c12, p1.c12, fct0, fct1, c12);
    GeometryUtilities::MapFromRegion(p0.c13, p1.c13, fct0, fct1, c13);
    GeometryUtilities::MapFromRegion(p0.c23, p1.c23, fct0, fct1, c23);
    GeometryUtilities::MapFromRegion(p0.c44, p1.c44, fct0, fct1, c44);
    GeometryUtilities::MapFromRegion(p0.c55, p1.c55, fct0, fct1, c55);
    GeometryUtilities::MapFromRegion(p0.c66, p1.c66, fct0, fct1, c66);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    MaterialBaseParameterData::ReadXML( node );

    const realT E1 = node.GetAttributeOrDefault("E1", 0.0);
    const realT E2 = node.GetAttributeOrDefault("E2", 0.0);
    const realT E3 = node.GetAttributeOrDefault("E3", 0.0);

    const realT nu12 = node.GetAttributeOrDefault("nu12", 0.0);
    const realT nu23 = node.GetAttributeOrDefault("nu23", 0.0);
    const realT nu13 = node.GetAttributeOrDefault("nu13", 0.0);

    const realT nu21 = nu12 *( E2/E1 );
    const realT nu31 = nu13 *( E3/E1 );
    const realT nu32 = nu23 *( E3/E2 );

    const realT G12 = node.GetAttributeOrDefault("G12", 0.0);
    const realT G23 = node.GetAttributeOrDefault("G23", 0.0);
    const realT G13 = node.GetAttributeOrDefault("G13", 0.0);

    const realT gamma = 1.0 / ( 1.0 - nu21*nu21 - nu23*nu32 - nu13*nu31 - 2*nu21*nu32*nu13 );

    c11 = E1*( 1.0 - nu23*nu32)*gamma;
    c12 = E1*(nu21 + nu23*nu31)*gamma;
    c13 = E1*(nu31 + nu21*nu32)*gamma;

    c22 = E2*( 1.0 - nu13*nu31)*gamma;
    c23 = E2*(nu12*nu31 + nu32)*gamma;

    c33 = E3*( 1.0 - nu12*nu21)*gamma;

    c44 = 2*G23;
    c55 = 2*G13;
    c66 = 2*G12;

    // hack for incorrect timestep calculation
    init_shearModulus = std::max(std::max(G12,G13),G23);
    Lame = std::max(std::max(E1,E2),E3);
    E = (E1 + E2 + E3) / 3.0;	

    PostReadXML( node );

  }


};

//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticOrthotropicStateData : public MaterialBaseStateData
{

public:

  typedef MaterialBaseStateData base;


  LinearElasticOrthotropicStateData():
    base()
  {}

  LinearElasticOrthotropicStateData( const LinearElasticOrthotropicStateData& source):
    base( source )
  {}

  ~LinearElasticOrthotropicStateData() {}
  friend class ConstitutiveBase;
  friend class LinearElasticOrthotropic;

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

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
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
  }
  inline LinearElasticOrthotropicStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline LinearElasticOrthotropicStateData&
  operator=(const LinearElasticOrthotropicStateData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline LinearElasticOrthotropicStateData&
  operator+=(const LinearElasticOrthotropicStateData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearElasticOrthotropicStateData& p0, LinearElasticOrthotropicStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const LinearElasticOrthotropicStateData& p0, const LinearElasticOrthotropicStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class LinearElasticOrthotropic: public MaterialBase
{
public:

  typedef LinearElasticOrthotropicParameterData ParameterClass;
  typedef LinearElasticOrthotropicStateData     StateClass;

  typedef Array1dT<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "LinearElasticOrthotropic"; }

  LinearElasticOrthotropic();
  virtual ~LinearElasticOrthotropic();

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
  { SetVariableParametersFromDerived<LinearElasticOrthotropic>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<LinearElasticOrthotropic>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<LinearElasticOrthotropic>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<LinearElasticOrthotropic>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<LinearElasticOrthotropic>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<LinearElasticOrthotropic>( num ); }
 
  void GetVariableNames( sArray1d& intVars, sArray1d& realVars, sArray1d& R1TensorVars, sArray1d& R2TensorVars, sArray1d& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<LinearElasticOrthotropic>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<LinearElasticOrthotropic>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<LinearElasticOrthotropic>(name, type ); }
  
  bool GetStateValues( const std::string& name, rArray1d& values ) const
  { return GetStateValuesFromDerived<LinearElasticOrthotropic>(name, values); }

  bool GetParameterValues( const std::string& name, rArray1d& values ) const
  { return GetParameterValuesFromDerived<LinearElasticOrthotropic>(name, values); }

  bool SetStateValues( const std::string& name, const rArray1d& values )
  { return SetStateValuesFromDerived<LinearElasticOrthotropic>(name, values); }

  bool SetParameterValues( const std::string& name, const rArray1d& values )
  { return SetParameterValuesFromDerived<LinearElasticOrthotropic>(name, values); }

  virtual void Serialize( Array1dT<iArray1d*>& intVars, Array1dT<rArray1d*>& realVars, Array1dT<Array1dT<R1Tensor>*>& R1Vars, Array1dT<Array1dT<R2Tensor>*>& R2Vars, Array1dT<Array1dT<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<LinearElasticOrthotropic>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const Array1dT<iArray1d*>& intVars, const Array1dT<rArray1d*>& realVars, const Array1dT<Array1dT<R1Tensor>*>& R1Vars, const Array1dT<Array1dT<R2Tensor>*>& R2Vars, const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<LinearElasticOrthotropic>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticOrthotropic>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticOrthotropic>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticOrthotropic>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticOrthotropic>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<LinearElasticOrthotropic>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<LinearElasticOrthotropic>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  

  virtual void
  StrainDrivenUpdateMember(const localIndex index0,
                           const localIndex index1,
                           const R2SymTensorT<3>& Ddt,
                           const R2TensorT < 3 >& L,
                           const R2Tensor& Rot,
                           const realT dt);

  virtual void
  StrainDrivenUpdateMember(const localIndex index0,
                           const localIndex index1,
                           const R2SymTensorT<3>& Ddt,
                           const R2TensorT < 3 >& L,
                           const R2Tensor& Rot,
                           const realT& volume_n,
                           const realT& volume_np1,
                           const realT dt);

  void MeanPressureDevStress( const localIndex index, realT& pressure,
                              R2SymTensor& devStress) const;

private:
  LinearElasticOrthotropic(const LinearElasticOrthotropic&);
  LinearElasticOrthotropic& operator=(const LinearElasticOrthotropic&);
  

};
#endif /* LINEARELASTICORTHOTROPIC_H_ */
