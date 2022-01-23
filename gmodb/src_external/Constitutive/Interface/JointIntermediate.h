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

#ifndef JOINTINTERMEDIATE_H_
#define JOINTINTERMEDIATE_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/Interface/InterfaceBase.h"

/*
 * JointIntermediate.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class JointIntermediateParameterData : public InterfaceBaseParameterData
{

public:

  typedef InterfaceBaseParameterData base;
  realT xsResidual;
  realT muResidual;
  realT cohesion;
  realT normalApproachYield;
  realT kplasticLoad;
  realT kplasticUnload;
  realT kdilation;
  realT dilationCoefficient0;
  realT xsDilationInit;
  realT xsDilationLimit;
  realT normalApproachDilationInit;
  realT stressDilationLimit;
  realT kCoefficientElastic;
  realT kshear;


  JointIntermediateParameterData():
    base(),
    xsResidual(0),
    muResidual(0),
    cohesion(0),
    normalApproachYield(0),
    kplasticLoad(0),
    kplasticUnload(0),
    kdilation(0),
    dilationCoefficient0(0),
    xsDilationInit(0),
    xsDilationLimit(0),
    normalApproachDilationInit(0),
    stressDilationLimit(0),
    kCoefficientElastic(0),
    kshear(0)
  {}

  JointIntermediateParameterData( const JointIntermediateParameterData& source):
    base( source ),
    xsResidual(source.xsResidual),
    muResidual(source.muResidual),
    cohesion(source.cohesion),
    normalApproachYield(source.normalApproachYield),
    kplasticLoad(source.kplasticLoad),
    kplasticUnload(source.kplasticUnload),
    kdilation(source.kdilation),
    dilationCoefficient0(source.dilationCoefficient0),
    xsDilationInit(source.xsDilationInit),
    xsDilationLimit(source.xsDilationLimit),
    normalApproachDilationInit(source.normalApproachDilationInit),
    stressDilationLimit(source.stressDilationLimit),
    kCoefficientElastic(source.kCoefficientElastic),
    kshear(source.kshear)
  {}

  ~JointIntermediateParameterData() {}
  friend class ConstitutiveBase;
  friend class JointIntermediate;

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
    realVarCounts += 14;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("shearAtResidual");
    realNames.push_back("coefficientOfFrictionAtResidual");
    realNames.push_back("cohesion");
    realNames.push_back("normalApproachAtYield");
    realNames.push_back("plasticStiffness");
    realNames.push_back("unloadStiffness");
    realNames.push_back("dilationStiffness");
    realNames.push_back("dilationCoefficient");
    realNames.push_back("shearAtInitiationOfDilation");
    realNames.push_back("shearAtLimitOfDilation");
    realNames.push_back("normalApproachAtInitiationOfDilation");
    realNames.push_back("stressAtLimitOfDilation");
    realNames.push_back("elasticStiffnessCoefficient");
    realNames.push_back("shearStiffness");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["shearAtResidual"] = (char*)(&xsResidual) - (char*)this;
    realOffsets["coefficientOfFrictionAtResidual"] = (char*)(&muResidual) - (char*)this;
    realOffsets["cohesion"] = (char*)(&cohesion) - (char*)this;
    realOffsets["normalApproachAtYield"] = (char*)(&normalApproachYield) - (char*)this;
    realOffsets["plasticStiffness"] = (char*)(&kplasticLoad) - (char*)this;
    realOffsets["unloadStiffness"] = (char*)(&kplasticUnload) - (char*)this;
    realOffsets["dilationStiffness"] = (char*)(&kdilation) - (char*)this;
    realOffsets["dilationCoefficient"] = (char*)(&dilationCoefficient0) - (char*)this;
    realOffsets["shearAtInitiationOfDilation"] = (char*)(&xsDilationInit) - (char*)this;
    realOffsets["shearAtLimitOfDilation"] = (char*)(&xsDilationLimit) - (char*)this;
    realOffsets["normalApproachAtInitiationOfDilation"] = (char*)(&normalApproachDilationInit) - (char*)this;
    realOffsets["stressAtLimitOfDilation"] = (char*)(&stressDilationLimit) - (char*)this;
    realOffsets["elasticStiffnessCoefficient"] = (char*)(&kCoefficientElastic) - (char*)this;
    realOffsets["shearStiffness"] = (char*)(&kshear) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["shearAtResidual"] = xsResidual;
    realValues["coefficientOfFrictionAtResidual"] = muResidual;
    realValues["cohesion"] = cohesion;
    realValues["normalApproachAtYield"] = normalApproachYield;
    realValues["plasticStiffness"] = kplasticLoad;
    realValues["unloadStiffness"] = kplasticUnload;
    realValues["dilationStiffness"] = kdilation;
    realValues["dilationCoefficient"] = dilationCoefficient0;
    realValues["shearAtInitiationOfDilation"] = xsDilationInit;
    realValues["shearAtLimitOfDilation"] = xsDilationLimit;
    realValues["normalApproachAtInitiationOfDilation"] = normalApproachDilationInit;
    realValues["stressAtLimitOfDilation"] = stressDilationLimit;
    realValues["elasticStiffnessCoefficient"] = kCoefficientElastic;
    realValues["shearStiffness"] = kshear;
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
    (*(realVars[realVarCounts]))[index] = xsResidual; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = muResidual; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = cohesion; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = normalApproachYield; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = kplasticLoad; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = kplasticUnload; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = kdilation; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = dilationCoefficient0; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = xsDilationInit; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = xsDilationLimit; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = normalApproachDilationInit; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = stressDilationLimit; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = kCoefficientElastic; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = kshear; realVarCounts++;
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
    xsResidual = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    muResidual = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    cohesion = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    normalApproachYield = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    kplasticLoad = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    kplasticUnload = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    kdilation = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    dilationCoefficient0 = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    xsDilationInit = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    xsDilationLimit = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    normalApproachDilationInit = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    stressDilationLimit = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    kCoefficientElastic = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    kshear = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline JointIntermediateParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    xsResidual *= factor;
    muResidual *= factor;
    cohesion *= factor;
    normalApproachYield *= factor;
    kplasticLoad *= factor;
    kplasticUnload *= factor;
    kdilation *= factor;
    dilationCoefficient0 *= factor;
    xsDilationInit *= factor;
    xsDilationLimit *= factor;
    normalApproachDilationInit *= factor;
    stressDilationLimit *= factor;
    kCoefficientElastic *= factor;
    kshear *= factor;
    return *this;
  }

  inline JointIntermediateParameterData&
  operator=(const JointIntermediateParameterData& datum)
  {
    base::operator=(datum);
    xsResidual = datum.xsResidual;
    muResidual = datum.muResidual;
    cohesion = datum.cohesion;
    normalApproachYield = datum.normalApproachYield;
    kplasticLoad = datum.kplasticLoad;
    kplasticUnload = datum.kplasticUnload;
    kdilation = datum.kdilation;
    dilationCoefficient0 = datum.dilationCoefficient0;
    xsDilationInit = datum.xsDilationInit;
    xsDilationLimit = datum.xsDilationLimit;
    normalApproachDilationInit = datum.normalApproachDilationInit;
    stressDilationLimit = datum.stressDilationLimit;
    kCoefficientElastic = datum.kCoefficientElastic;
    kshear = datum.kshear;
    return *this;
  }

  inline JointIntermediateParameterData&
  operator+=(const JointIntermediateParameterData& datum)
  {
    base::operator+=(datum);
    xsResidual += datum.xsResidual;
    muResidual += datum.muResidual;
    cohesion += datum.cohesion;
    normalApproachYield += datum.normalApproachYield;
    kplasticLoad += datum.kplasticLoad;
    kplasticUnload += datum.kplasticUnload;
    kdilation += datum.kdilation;
    dilationCoefficient0 += datum.dilationCoefficient0;
    xsDilationInit += datum.xsDilationInit;
    xsDilationLimit += datum.xsDilationLimit;
    normalApproachDilationInit += datum.normalApproachDilationInit;
    stressDilationLimit += datum.stressDilationLimit;
    kCoefficientElastic += datum.kCoefficientElastic;
    kshear += datum.kshear;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, JointIntermediateParameterData& p0, JointIntermediateParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, xsResidual, p0.xsResidual, p1.xsResidual);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, muResidual, p0.muResidual, p1.muResidual);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, cohesion, p0.cohesion, p1.cohesion);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, normalApproachYield, p0.normalApproachYield, p1.normalApproachYield);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kplasticLoad, p0.kplasticLoad, p1.kplasticLoad);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kplasticUnload, p0.kplasticUnload, p1.kplasticUnload);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kdilation, p0.kdilation, p1.kdilation);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dilationCoefficient0, p0.dilationCoefficient0, p1.dilationCoefficient0);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, xsDilationInit, p0.xsDilationInit, p1.xsDilationInit);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, xsDilationLimit, p0.xsDilationLimit, p1.xsDilationLimit);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, normalApproachDilationInit, p0.normalApproachDilationInit, p1.normalApproachDilationInit);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressDilationLimit, p0.stressDilationLimit, p1.stressDilationLimit);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kCoefficientElastic, p0.kCoefficientElastic, p1.kCoefficientElastic);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kshear, p0.kshear, p1.kshear);

  }

  void MapFromRegion(const JointIntermediateParameterData& p0, const JointIntermediateParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.xsResidual, p1.xsResidual, fct0, fct1, xsResidual);
    GeometryUtilities::MapFromRegion(p0.muResidual, p1.muResidual, fct0, fct1, muResidual);
    GeometryUtilities::MapFromRegion(p0.cohesion, p1.cohesion, fct0, fct1, cohesion);
    GeometryUtilities::MapFromRegion(p0.normalApproachYield, p1.normalApproachYield, fct0, fct1, normalApproachYield);
    GeometryUtilities::MapFromRegion(p0.kplasticLoad, p1.kplasticLoad, fct0, fct1, kplasticLoad);
    GeometryUtilities::MapFromRegion(p0.kplasticUnload, p1.kplasticUnload, fct0, fct1, kplasticUnload);
    GeometryUtilities::MapFromRegion(p0.kdilation, p1.kdilation, fct0, fct1, kdilation);
    GeometryUtilities::MapFromRegion(p0.dilationCoefficient0, p1.dilationCoefficient0, fct0, fct1, dilationCoefficient0);
    GeometryUtilities::MapFromRegion(p0.xsDilationInit, p1.xsDilationInit, fct0, fct1, xsDilationInit);
    GeometryUtilities::MapFromRegion(p0.xsDilationLimit, p1.xsDilationLimit, fct0, fct1, xsDilationLimit);
    GeometryUtilities::MapFromRegion(p0.normalApproachDilationInit, p1.normalApproachDilationInit, fct0, fct1, normalApproachDilationInit);
    GeometryUtilities::MapFromRegion(p0.stressDilationLimit, p1.stressDilationLimit, fct0, fct1, stressDilationLimit);
    GeometryUtilities::MapFromRegion(p0.kCoefficientElastic, p1.kCoefficientElastic, fct0, fct1, kCoefficientElastic);
    GeometryUtilities::MapFromRegion(p0.kshear, p1.kshear, fct0, fct1, kshear);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    InterfaceBaseParameterData::ReadXML( node );
    xsResidual = node.GetAttributeOrDefault("shearAtResidual", 0.0);
    muResidual = node.GetAttributeOrDefault("coefficientOfFrictionAtResidual", 0.0);
    cohesion = node.GetAttributeOrDefault("cohesion", 0.0);
    normalApproachYield = node.GetAttributeOrDefault("normalApproachAtYield", 0.0);
    kplasticLoad = node.GetAttributeOrDefault("plasticStiffness", 0.0);
    kplasticUnload = node.GetAttributeOrDefault("unloadStiffness", 0.0);
    kdilation = node.GetAttributeOrDefault("dilationStiffness", 0.0);
    dilationCoefficient0 = node.GetAttributeOrDefault("dilationCoefficient", 0.0);
    xsDilationInit = node.GetAttributeOrDefault("shearAtInitiationOfDilation", 0.0);
    xsDilationLimit = node.GetAttributeOrDefault("shearAtLimitOfDilation", 0.0);
    normalApproachDilationInit = node.GetAttributeOrDefault("normalApproachAtInitiationOfDilation", 0.0);
    stressDilationLimit = node.GetAttributeOrDefault("stressAtLimitOfDilation", 0.0);
    kCoefficientElastic = node.GetAttributeOrDefault("elasticStiffnessCoefficient", 0.0);
    kshear = node.GetAttributeOrDefault("shearStiffness", 0.0);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class JointIntermediateStateData : public InterfaceBaseStateData
{

public:

  typedef InterfaceBaseStateData base;
  int ifail;
  realT gap;
  realT mu;
  realT kcurrent;
  realT stressDilation;


  JointIntermediateStateData():
    base(),
    ifail(0),
    gap(0),
    mu(0),
    kcurrent(0),
    stressDilation(0)
  {}

  JointIntermediateStateData( const JointIntermediateStateData& source):
    base( source ),
    ifail(source.ifail),
    gap(source.gap),
    mu(source.mu),
    kcurrent(source.kcurrent),
    stressDilation(source.stressDilation)
  {}

  ~JointIntermediateStateData() {}
  friend class ConstitutiveBase;
  friend class JointIntermediate;

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
    intVarCounts += 1;
    realVarCounts += 4;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    intNames.push_back("ifail");
    realNames.push_back("gap");
    realNames.push_back("currentFrictionCoefficient");
    realNames.push_back("currentStiffness");
    realNames.push_back("dilationStress");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    intOffsets["ifail"] = (char*)(&ifail) - (char*)this;
    realOffsets["gap"] = (char*)(&gap) - (char*)this;
    realOffsets["currentFrictionCoefficient"] = (char*)(&mu) - (char*)this;
    realOffsets["currentStiffness"] = (char*)(&kcurrent) - (char*)this;
    realOffsets["dilationStress"] = (char*)(&stressDilation) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    intValues["ifail"] = ifail;
    realValues["gap"] = gap;
    realValues["currentFrictionCoefficient"] = mu;
    realValues["currentStiffness"] = kcurrent;
    realValues["dilationStress"] = stressDilation;
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
    (*(intVars[intVarCounts]))[elemNum] = ifail; intVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = gap; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mu; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = kcurrent; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = stressDilation; realVarCounts += stride;
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
    ifail = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    gap = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mu = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    kcurrent = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    stressDilation = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
  }
  inline JointIntermediateStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    ifail *= factor;
    gap *= factor;
    mu *= factor;
    kcurrent *= factor;
    stressDilation *= factor;
    return *this;
  }

  inline JointIntermediateStateData&
  operator=(const JointIntermediateStateData& datum)
  {
    base::operator=(datum);
    ifail = datum.ifail;
    gap = datum.gap;
    mu = datum.mu;
    kcurrent = datum.kcurrent;
    stressDilation = datum.stressDilation;
    return *this;
  }

  inline JointIntermediateStateData&
  operator+=(const JointIntermediateStateData& datum)
  {
    base::operator+=(datum);
    ifail += datum.ifail;
    gap += datum.gap;
    mu += datum.mu;
    kcurrent += datum.kcurrent;
    stressDilation += datum.stressDilation;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, JointIntermediateStateData& p0, JointIntermediateStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, ifail, p0.ifail, p1.ifail);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, gap, p0.gap, p1.gap);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu, p0.mu, p1.mu);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, kcurrent, p0.kcurrent, p1.kcurrent);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressDilation, p0.stressDilation, p1.stressDilation);

  }

  void MapFromRegion(const JointIntermediateStateData& p0, const JointIntermediateStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.ifail, p1.ifail, fct0, fct1, ifail);
    GeometryUtilities::MapFromRegion(p0.gap, p1.gap, fct0, fct1, gap);
    GeometryUtilities::MapFromRegion(p0.mu, p1.mu, fct0, fct1, mu);
    GeometryUtilities::MapFromRegion(p0.kcurrent, p1.kcurrent, fct0, fct1, kcurrent);
    GeometryUtilities::MapFromRegion(p0.stressDilation, p1.stressDilation, fct0, fct1, stressDilation);

  }



};


//**********************************************************************************************************************
//**********************************************************************************************************************


class JointIntermediate: public InterfaceBase
{
public:
  
  typedef JointIntermediateParameterData ParameterClass;
  typedef JointIntermediateStateData     StateClass;
  

  JointIntermediate( const int paramSize, const int stateSize );

  virtual ~JointIntermediate();
  
  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;
  
  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;
  
  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;
  
  virtual const JointIntermediateStateData* StateData( const localIndex index0,
                                                  const localIndex index1 ) const = 0;
  virtual       JointIntermediateStateData* StateData( const localIndex index0,
                                                  const localIndex index1 )  = 0;

  virtual const JointIntermediateParameterData* ParameterData( const localIndex index ) const = 0;
  virtual       JointIntermediateParameterData* ParameterData( const localIndex index ) = 0;
  
  inline void IncrementPtr( const JointIntermediateStateData* ptr ) const
  {
    ptr = reinterpret_cast<const JointIntermediateStateData*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const JointIntermediateParameterData* ptr ) const
  {
    ptr = reinterpret_cast<const JointIntermediateParameterData*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
  }
  
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
  
  virtual void ZeroStates() = 0;
    
  virtual localIndex NumStateIndex0() const = 0;
  virtual localIndex NumStateIndex1() const = 0;

  virtual localIndex NumParameterIndex0() const = 0;
  virtual localIndex NumParameterIndex1() const = 0;

  virtual void
  Initialize( const localIndex index,
                                 const realT stressNormal,
                                 const realT stressShear);

  virtual void
  StrainDrivenUpdate( const localIndex index );

  virtual realT
  StiffnessProjected(const localIndex index);


protected:
  virtual void
  UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                     InterfaceBaseStateData& matStateBase) const;

  virtual realT
  ShearStrength(const InterfaceBaseParameterData& ,
                                   InterfaceBaseStateData& matStateBase) const;

  virtual realT
  DilationalStressIncrement( const InterfaceBaseParameterData&,
                                                InterfaceBaseStateData& ,
                                                const realT ) const;

  virtual realT
  NormalApproachAtInitialization( const realT stressNormal,
                                                     const realT knCoefficientElastic,
                                                     const realT knPlasticLoading,
                                                     const realT knPlasticUnloading,
                                                     realT& kcurrent,
                                                     realT& normalGap,
                                                     const realT normalApproachNormalYield) const;

  realT
  NormalStiffnessRock( const realT knCoefficientElastic,
                                          const realT knPlasticLoading,
                                          const realT knPlasticUnloading,
                                          const realT normalApproach,
                                          realT& normalGap,
                                          const realT normalApproachNormalYield) const;

  realT
  ShearStrength(const realT stress,
                                   const realT normalStressAtDilationLimit,
                                   const realT tanFrictionCoefficientInitial,
                                   const realT tanFrictionCoefficientResidual,
                                   const realT cohesion,
                                   const int ifail) const;

  
private:
  JointIntermediate();
  JointIntermediate( const JointIntermediate& );
  JointIntermediate& operator=( const JointIntermediate& );
  
  
};
#endif /* JOINTINTERMEDIATE_H_ */
