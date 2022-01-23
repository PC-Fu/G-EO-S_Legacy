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

#ifndef LINEARELASTICGENERALANISOTROPIC_H_
#define LINEARELASTICGENERALANISOTROPIC_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/Material/MaterialBase.h"

/*
 * LinearElasticGeneralAnisotropic.h
 *
 *  Created on: Sun Jun 29 15:11:06 PDT 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticGeneralAnisotropicParameterData : public MaterialBaseParameterData
{

public:

  typedef MaterialBaseParameterData base;
  realT c11;
  realT c22;
  realT c33;
  realT c44;
  realT c55;
  realT c66;
  realT c12;
  realT c13;
  realT c14;
  realT c15;
  realT c16;
  realT c23;
  realT c24;
  realT c25;
  realT c26;
  realT c34;
  realT c35;
  realT c36;
  realT c45;
  realT c46;
  realT c56;


  LinearElasticGeneralAnisotropicParameterData():
    base(),
    c11(0),
    c22(0),
    c33(0),
    c44(0),
    c55(0),
    c66(0),
    c12(0),
    c13(0),
    c14(0),
    c15(0),
    c16(0),
    c23(0),
    c24(0),
    c25(0),
    c26(0),
    c34(0),
    c35(0),
    c36(0),
    c45(0),
    c46(0),
    c56(0)
  {}

  LinearElasticGeneralAnisotropicParameterData( const LinearElasticGeneralAnisotropicParameterData& source):
    base( source ),
    c11(source.c11),
    c22(source.c22),
    c33(source.c33),
    c44(source.c44),
    c55(source.c55),
    c66(source.c66),
    c12(source.c12),
    c13(source.c13),
    c14(source.c14),
    c15(source.c15),
    c16(source.c16),
    c23(source.c23),
    c24(source.c24),
    c25(source.c25),
    c26(source.c26),
    c34(source.c34),
    c35(source.c35),
    c36(source.c36),
    c45(source.c45),
    c46(source.c46),
    c56(source.c56)
  {}

  ~LinearElasticGeneralAnisotropicParameterData() {}
  friend class ConstitutiveBase;
  friend class LinearElasticGeneralAnisotropic;

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
    realVarCounts += 21;

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
    realNames.push_back("c44");
    realNames.push_back("c55");
    realNames.push_back("c66");
    realNames.push_back("c12");
    realNames.push_back("c13");
    realNames.push_back("c14");
    realNames.push_back("c15");
    realNames.push_back("c16");
    realNames.push_back("c23");
    realNames.push_back("c24");
    realNames.push_back("c25");
    realNames.push_back("c26");
    realNames.push_back("c34");
    realNames.push_back("c35");
    realNames.push_back("c36");
    realNames.push_back("c45");
    realNames.push_back("c46");
    realNames.push_back("c56");
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
    realOffsets["c44"] = (char*)(&c44) - (char*)this;
    realOffsets["c55"] = (char*)(&c55) - (char*)this;
    realOffsets["c66"] = (char*)(&c66) - (char*)this;
    realOffsets["c12"] = (char*)(&c12) - (char*)this;
    realOffsets["c13"] = (char*)(&c13) - (char*)this;
    realOffsets["c14"] = (char*)(&c14) - (char*)this;
    realOffsets["c15"] = (char*)(&c15) - (char*)this;
    realOffsets["c16"] = (char*)(&c16) - (char*)this;
    realOffsets["c23"] = (char*)(&c23) - (char*)this;
    realOffsets["c24"] = (char*)(&c24) - (char*)this;
    realOffsets["c25"] = (char*)(&c25) - (char*)this;
    realOffsets["c26"] = (char*)(&c26) - (char*)this;
    realOffsets["c34"] = (char*)(&c34) - (char*)this;
    realOffsets["c35"] = (char*)(&c35) - (char*)this;
    realOffsets["c36"] = (char*)(&c36) - (char*)this;
    realOffsets["c45"] = (char*)(&c45) - (char*)this;
    realOffsets["c46"] = (char*)(&c46) - (char*)this;
    realOffsets["c56"] = (char*)(&c56) - (char*)this;
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
    realValues["c44"] = c44;
    realValues["c55"] = c55;
    realValues["c66"] = c66;
    realValues["c12"] = c12;
    realValues["c13"] = c13;
    realValues["c14"] = c14;
    realValues["c15"] = c15;
    realValues["c16"] = c16;
    realValues["c23"] = c23;
    realValues["c24"] = c24;
    realValues["c25"] = c25;
    realValues["c26"] = c26;
    realValues["c34"] = c34;
    realValues["c35"] = c35;
    realValues["c36"] = c36;
    realValues["c45"] = c45;
    realValues["c46"] = c46;
    realValues["c56"] = c56;
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
    (*(realVars[realVarCounts]))[index] = c44; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c55; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c66; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c12; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c13; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c14; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c15; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c16; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c23; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c24; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c25; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c26; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c34; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c35; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c36; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c45; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c46; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = c56; realVarCounts++;
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
    c44 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c55 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c66 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c12 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c13 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c14 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c15 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c16 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c23 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c24 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c25 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c26 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c34 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c35 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c36 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c45 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c46 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    c56 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline LinearElasticGeneralAnisotropicParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    c11 *= factor;
    c22 *= factor;
    c33 *= factor;
    c44 *= factor;
    c55 *= factor;
    c66 *= factor;
    c12 *= factor;
    c13 *= factor;
    c14 *= factor;
    c15 *= factor;
    c16 *= factor;
    c23 *= factor;
    c24 *= factor;
    c25 *= factor;
    c26 *= factor;
    c34 *= factor;
    c35 *= factor;
    c36 *= factor;
    c45 *= factor;
    c46 *= factor;
    c56 *= factor;
    return *this;
  }

  inline LinearElasticGeneralAnisotropicParameterData&
  operator=(const LinearElasticGeneralAnisotropicParameterData& datum)
  {
    base::operator=(datum);
    c11 = datum.c11;
    c22 = datum.c22;
    c33 = datum.c33;
    c44 = datum.c44;
    c55 = datum.c55;
    c66 = datum.c66;
    c12 = datum.c12;
    c13 = datum.c13;
    c14 = datum.c14;
    c15 = datum.c15;
    c16 = datum.c16;
    c23 = datum.c23;
    c24 = datum.c24;
    c25 = datum.c25;
    c26 = datum.c26;
    c34 = datum.c34;
    c35 = datum.c35;
    c36 = datum.c36;
    c45 = datum.c45;
    c46 = datum.c46;
    c56 = datum.c56;
    return *this;
  }

  inline LinearElasticGeneralAnisotropicParameterData&
  operator+=(const LinearElasticGeneralAnisotropicParameterData& datum)
  {
    base::operator+=(datum);
    c11 += datum.c11;
    c22 += datum.c22;
    c33 += datum.c33;
    c44 += datum.c44;
    c55 += datum.c55;
    c66 += datum.c66;
    c12 += datum.c12;
    c13 += datum.c13;
    c14 += datum.c14;
    c15 += datum.c15;
    c16 += datum.c16;
    c23 += datum.c23;
    c24 += datum.c24;
    c25 += datum.c25;
    c26 += datum.c26;
    c34 += datum.c34;
    c35 += datum.c35;
    c36 += datum.c36;
    c45 += datum.c45;
    c46 += datum.c46;
    c56 += datum.c56;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearElasticGeneralAnisotropicParameterData& p0, LinearElasticGeneralAnisotropicParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c11, p0.c11, p1.c11);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c22, p0.c22, p1.c22);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c33, p0.c33, p1.c33);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c44, p0.c44, p1.c44);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c55, p0.c55, p1.c55);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c66, p0.c66, p1.c66);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c12, p0.c12, p1.c12);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c13, p0.c13, p1.c13);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c14, p0.c14, p1.c14);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c15, p0.c15, p1.c15);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c16, p0.c16, p1.c16);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c23, p0.c23, p1.c23);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c24, p0.c24, p1.c24);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c25, p0.c25, p1.c25);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c26, p0.c26, p1.c26);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c34, p0.c34, p1.c34);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c35, p0.c35, p1.c35);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c36, p0.c36, p1.c36);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c45, p0.c45, p1.c45);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c46, p0.c46, p1.c46);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, c56, p0.c56, p1.c56);

  }

  void MapFromRegion(const LinearElasticGeneralAnisotropicParameterData& p0, const LinearElasticGeneralAnisotropicParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.c11, p1.c11, fct0, fct1, c11);
    GeometryUtilities::MapFromRegion(p0.c22, p1.c22, fct0, fct1, c22);
    GeometryUtilities::MapFromRegion(p0.c33, p1.c33, fct0, fct1, c33);
    GeometryUtilities::MapFromRegion(p0.c44, p1.c44, fct0, fct1, c44);
    GeometryUtilities::MapFromRegion(p0.c55, p1.c55, fct0, fct1, c55);
    GeometryUtilities::MapFromRegion(p0.c66, p1.c66, fct0, fct1, c66);
    GeometryUtilities::MapFromRegion(p0.c12, p1.c12, fct0, fct1, c12);
    GeometryUtilities::MapFromRegion(p0.c13, p1.c13, fct0, fct1, c13);
    GeometryUtilities::MapFromRegion(p0.c14, p1.c14, fct0, fct1, c14);
    GeometryUtilities::MapFromRegion(p0.c15, p1.c15, fct0, fct1, c15);
    GeometryUtilities::MapFromRegion(p0.c16, p1.c16, fct0, fct1, c16);
    GeometryUtilities::MapFromRegion(p0.c23, p1.c23, fct0, fct1, c23);
    GeometryUtilities::MapFromRegion(p0.c24, p1.c24, fct0, fct1, c24);
    GeometryUtilities::MapFromRegion(p0.c25, p1.c25, fct0, fct1, c25);
    GeometryUtilities::MapFromRegion(p0.c26, p1.c26, fct0, fct1, c26);
    GeometryUtilities::MapFromRegion(p0.c34, p1.c34, fct0, fct1, c34);
    GeometryUtilities::MapFromRegion(p0.c35, p1.c35, fct0, fct1, c35);
    GeometryUtilities::MapFromRegion(p0.c36, p1.c36, fct0, fct1, c36);
    GeometryUtilities::MapFromRegion(p0.c45, p1.c45, fct0, fct1, c45);
    GeometryUtilities::MapFromRegion(p0.c46, p1.c46, fct0, fct1, c46);
    GeometryUtilities::MapFromRegion(p0.c56, p1.c56, fct0, fct1, c56);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    MaterialBaseParameterData::ReadXML( node );
    c11 = node.GetAttributeOrDefault("c11", 0.0);
    c22 = node.GetAttributeOrDefault("c22", 0.0);
    c33 = node.GetAttributeOrDefault("c33", 0.0);
    c44 = node.GetAttributeOrDefault("c44", 0.0);
    c55 = node.GetAttributeOrDefault("c55", 0.0);
    c66 = node.GetAttributeOrDefault("c66", 0.0);
    c12 = node.GetAttributeOrDefault("c12", 0.0);
    c13 = node.GetAttributeOrDefault("c13", 0.0);
    c14 = node.GetAttributeOrDefault("c14", 0.0);
    c15 = node.GetAttributeOrDefault("c15", 0.0);
    c16 = node.GetAttributeOrDefault("c16", 0.0);
    c23 = node.GetAttributeOrDefault("c23", 0.0);
    c24 = node.GetAttributeOrDefault("c24", 0.0);
    c25 = node.GetAttributeOrDefault("c25", 0.0);
    c26 = node.GetAttributeOrDefault("c26", 0.0);
    c34 = node.GetAttributeOrDefault("c34", 0.0);
    c35 = node.GetAttributeOrDefault("c35", 0.0);
    c36 = node.GetAttributeOrDefault("c36", 0.0);
    c45 = node.GetAttributeOrDefault("c45", 0.0);
    c46 = node.GetAttributeOrDefault("c46", 0.0);
    c56 = node.GetAttributeOrDefault("c56", 0.0);
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class LinearElasticGeneralAnisotropicStateData : public MaterialBaseStateData
{

public:

  typedef MaterialBaseStateData base;


  LinearElasticGeneralAnisotropicStateData():
    base()
  {}

  LinearElasticGeneralAnisotropicStateData( const LinearElasticGeneralAnisotropicStateData& source):
    base( source )
  {}

  ~LinearElasticGeneralAnisotropicStateData() {}
  friend class ConstitutiveBase;
  friend class LinearElasticGeneralAnisotropic;

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
  inline LinearElasticGeneralAnisotropicStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline LinearElasticGeneralAnisotropicStateData&
  operator=(const LinearElasticGeneralAnisotropicStateData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline LinearElasticGeneralAnisotropicStateData&
  operator+=(const LinearElasticGeneralAnisotropicStateData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, LinearElasticGeneralAnisotropicStateData& p0, LinearElasticGeneralAnisotropicStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const LinearElasticGeneralAnisotropicStateData& p0, const LinearElasticGeneralAnisotropicStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class LinearElasticGeneralAnisotropic: public MaterialBase
{
public:

  typedef LinearElasticGeneralAnisotropicParameterData ParameterClass;
  typedef LinearElasticGeneralAnisotropicStateData     StateClass;

  typedef Array1dT<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "LinearElasticGeneralAnisotropic"; }

  LinearElasticGeneralAnisotropic();
  virtual ~LinearElasticGeneralAnisotropic();

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
  { SetVariableParametersFromDerived<LinearElasticGeneralAnisotropic>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<LinearElasticGeneralAnisotropic>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<LinearElasticGeneralAnisotropic>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<LinearElasticGeneralAnisotropic>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<LinearElasticGeneralAnisotropic>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<LinearElasticGeneralAnisotropic>( num ); }
 
  void GetVariableNames( sArray1d& intVars, sArray1d& realVars, sArray1d& R1TensorVars, sArray1d& R2TensorVars, sArray1d& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<LinearElasticGeneralAnisotropic>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<LinearElasticGeneralAnisotropic>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<LinearElasticGeneralAnisotropic>(name, type ); }
  
  bool GetStateValues( const std::string& name, rArray1d& values ) const
  { return GetStateValuesFromDerived<LinearElasticGeneralAnisotropic>(name, values); }

  bool GetParameterValues( const std::string& name, rArray1d& values ) const
  { return GetParameterValuesFromDerived<LinearElasticGeneralAnisotropic>(name, values); }

  bool SetStateValues( const std::string& name, const rArray1d& values )
  { return SetStateValuesFromDerived<LinearElasticGeneralAnisotropic>(name, values); }

  bool SetParameterValues( const std::string& name, const rArray1d& values )
  { return SetParameterValuesFromDerived<LinearElasticGeneralAnisotropic>(name, values); }

  virtual void Serialize( Array1dT<iArray1d*>& intVars, Array1dT<rArray1d*>& realVars, Array1dT<Array1dT<R1Tensor>*>& R1Vars, Array1dT<Array1dT<R2Tensor>*>& R2Vars, Array1dT<Array1dT<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<LinearElasticGeneralAnisotropic>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const Array1dT<iArray1d*>& intVars, const Array1dT<rArray1d*>& realVars, const Array1dT<Array1dT<R1Tensor>*>& R1Vars, const Array1dT<Array1dT<R2Tensor>*>& R2Vars, const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<LinearElasticGeneralAnisotropic>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticGeneralAnisotropic>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticGeneralAnisotropic>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticGeneralAnisotropic>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<LinearElasticGeneralAnisotropic>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<LinearElasticGeneralAnisotropic>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<LinearElasticGeneralAnisotropic>( localIndices, buffer ); }

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
  LinearElasticGeneralAnisotropic(const LinearElasticGeneralAnisotropic&);
  LinearElasticGeneralAnisotropic& operator=(const LinearElasticGeneralAnisotropic&);
  

};
#endif /* LINEARELASTICGENERALANISOTROPIC_H_ */
