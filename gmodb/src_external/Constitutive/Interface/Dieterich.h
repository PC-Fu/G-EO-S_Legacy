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

#ifndef DIETERICH_H_
#define DIETERICH_H_

#include "Utilities/GeometryUtilities.h"
#include "RateAndStateIntermediate.h"
#include "PhysicsSolvers/Seismicity/BoundaryElementDataManagerT.h"
#include "DataStructures/Tables/Table.h"
#if GPAC_MPI
#include <mpi.h>
#endif


/*
 * Dieterich.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class DieterichParameterData : public RateAndStateIntermediateParameterData
{

public:

  typedef RateAndStateIntermediateParameterData base;
  realT tFail;
  realT maxThetaPin;
  realT KShearSelf;
  realT KNormalSelf;
  realT dxsdtAB;
  realT stressShearFail;
  realT AreductionFactor;
  realT stressOvershootFactor;
  R1Tensor rake;


  DieterichParameterData():
    base(),
    tFail(0),
    maxThetaPin(0),
    KShearSelf(0),
    KNormalSelf(0),
    dxsdtAB(0),
    stressShearFail(0),
    AreductionFactor(0),
    stressOvershootFactor(0),
    rake(0)
  {}

  DieterichParameterData( const DieterichParameterData& source):
    base( source ),
    tFail(source.tFail),
    maxThetaPin(source.maxThetaPin),
    KShearSelf(source.KShearSelf),
    KNormalSelf(source.KNormalSelf),
    dxsdtAB(source.dxsdtAB),
    stressShearFail(source.stressShearFail),
    AreductionFactor(source.AreductionFactor),
    stressOvershootFactor(source.stressOvershootFactor),
    rake(source.rake)
  {}

  ~DieterichParameterData() {}
  friend class ConstitutiveBase;
  friend class Dieterich;

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
    realVarCounts += 8;
    R1TensorVarCounts += 1;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("tFail");
    realNames.push_back("maxStatePin");
    realNames.push_back("KShearSelf");
    realNames.push_back("KNormalSelf");
    realNames.push_back("shearSlipRateAB");
    realNames.push_back("shearStressAtFailure");
    realNames.push_back("AReductionFactor");
    realNames.push_back("stressOverShootFactor");
    R1TensorNames.push_back("unitRakeVector");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["tFail"] = (char*)(&tFail) - (char*)this;
    realOffsets["maxStatePin"] = (char*)(&maxThetaPin) - (char*)this;
    realOffsets["KShearSelf"] = (char*)(&KShearSelf) - (char*)this;
    realOffsets["KNormalSelf"] = (char*)(&KNormalSelf) - (char*)this;
    realOffsets["shearSlipRateAB"] = (char*)(&dxsdtAB) - (char*)this;
    realOffsets["shearStressAtFailure"] = (char*)(&stressShearFail) - (char*)this;
    realOffsets["AReductionFactor"] = (char*)(&AreductionFactor) - (char*)this;
    realOffsets["stressOverShootFactor"] = (char*)(&stressOvershootFactor) - (char*)this;
    R1TensorOffsets["unitRakeVector"] = (char*)(&rake) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["tFail"] = tFail;
    realValues["maxStatePin"] = maxThetaPin;
    realValues["KShearSelf"] = KShearSelf;
    realValues["KNormalSelf"] = KNormalSelf;
    realValues["shearSlipRateAB"] = dxsdtAB;
    realValues["shearStressAtFailure"] = stressShearFail;
    realValues["AReductionFactor"] = AreductionFactor;
    realValues["stressOverShootFactor"] = stressOvershootFactor;
    R1TensorValues["unitRakeVector"] = rake;
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
    (*(realVars[realVarCounts]))[index] = tFail; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = maxThetaPin; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = KShearSelf; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = KNormalSelf; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = dxsdtAB; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = stressShearFail; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = AreductionFactor; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = stressOvershootFactor; realVarCounts++;
    (*(R1Vars[R1TensorVarCounts]))[index] = rake; R1TensorVarCounts++;
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
    tFail = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    maxThetaPin = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    KShearSelf = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    KNormalSelf = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    dxsdtAB = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    stressShearFail = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    AreductionFactor = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    stressOvershootFactor = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    rake = (*(R1Vars[R1TensorVarCounts]))[index]; R1TensorVarCounts++;
  }
  inline DieterichParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    tFail *= factor;
    maxThetaPin *= factor;
    KShearSelf *= factor;
    KNormalSelf *= factor;
    dxsdtAB *= factor;
    stressShearFail *= factor;
    AreductionFactor *= factor;
    stressOvershootFactor *= factor;
    rake *= factor;
    return *this;
  }

  inline DieterichParameterData&
  operator=(const DieterichParameterData& datum)
  {
    base::operator=(datum);
    tFail = datum.tFail;
    maxThetaPin = datum.maxThetaPin;
    KShearSelf = datum.KShearSelf;
    KNormalSelf = datum.KNormalSelf;
    dxsdtAB = datum.dxsdtAB;
    stressShearFail = datum.stressShearFail;
    AreductionFactor = datum.AreductionFactor;
    stressOvershootFactor = datum.stressOvershootFactor;
    rake = datum.rake;
    return *this;
  }

  inline DieterichParameterData&
  operator+=(const DieterichParameterData& datum)
  {
    base::operator+=(datum);
    tFail += datum.tFail;
    maxThetaPin += datum.maxThetaPin;
    KShearSelf += datum.KShearSelf;
    KNormalSelf += datum.KNormalSelf;
    dxsdtAB += datum.dxsdtAB;
    stressShearFail += datum.stressShearFail;
    AreductionFactor += datum.AreductionFactor;
    stressOvershootFactor += datum.stressOvershootFactor;
    rake += datum.rake;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, DieterichParameterData& p0, DieterichParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, tFail, p0.tFail, p1.tFail);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, maxThetaPin, p0.maxThetaPin, p1.maxThetaPin);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, KShearSelf, p0.KShearSelf, p1.KShearSelf);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, KNormalSelf, p0.KNormalSelf, p1.KNormalSelf);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dxsdtAB, p0.dxsdtAB, p1.dxsdtAB);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressShearFail, p0.stressShearFail, p1.stressShearFail);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, AreductionFactor, p0.AreductionFactor, p1.AreductionFactor);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressOvershootFactor, p0.stressOvershootFactor, p1.stressOvershootFactor);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, rake, p0.rake, p1.rake);

  }

  void MapFromRegion(const DieterichParameterData& p0, const DieterichParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.tFail, p1.tFail, fct0, fct1, tFail);
    GeometryUtilities::MapFromRegion(p0.maxThetaPin, p1.maxThetaPin, fct0, fct1, maxThetaPin);
    GeometryUtilities::MapFromRegion(p0.KShearSelf, p1.KShearSelf, fct0, fct1, KShearSelf);
    GeometryUtilities::MapFromRegion(p0.KNormalSelf, p1.KNormalSelf, fct0, fct1, KNormalSelf);
    GeometryUtilities::MapFromRegion(p0.dxsdtAB, p1.dxsdtAB, fct0, fct1, dxsdtAB);
    GeometryUtilities::MapFromRegion(p0.stressShearFail, p1.stressShearFail, fct0, fct1, stressShearFail);
    GeometryUtilities::MapFromRegion(p0.AreductionFactor, p1.AreductionFactor, fct0, fct1, AreductionFactor);
    GeometryUtilities::MapFromRegion(p0.stressOvershootFactor, p1.stressOvershootFactor, fct0, fct1, stressOvershootFactor);
    GeometryUtilities::MapFromRegion(p0.rake, p1.rake, fct0, fct1, rake);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    RateAndStateIntermediateParameterData::ReadXML( node );
    tFail = node.GetAttributeOrDefault("tFail", 0.0);
    maxThetaPin = node.GetAttributeOrDefault("maxStatePin", 0.0);
    KShearSelf = node.GetAttributeOrDefault("KShearSelf", 0.0);
    KNormalSelf = node.GetAttributeOrDefault("KNormalSelf", 0.0);
    dxsdtAB = node.GetAttributeOrDefault("shearSlipRateAB", 0.0);
    stressShearFail = node.GetAttributeOrDefault("shearStressAtFailure", 0.0);
    AreductionFactor = node.GetAttributeOrDefault("AReductionFactor", 0.0);
    stressOvershootFactor = node.GetAttributeOrDefault("stressOverShootFactor", 0.0);
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class DieterichStateData : public RateAndStateIntermediateStateData
{

public:

  typedef RateAndStateIntermediateStateData base;
  int nextState;
  int apFail;
  int pinned;
  int slowSlip;
  int neighborInRuptureState;
  int currentState;
  int Areduced;
  realT stressReference;
  realT dstressdtDrive;
  realT dxsdtDrive;
  realT stressPin;
  realT stressShearReference;
  realT dstressShear;
  realT dstress;
  realT dstressShearDt;
  realT dstressShearDtDrive;
  realT dxs;
  realT pp;
  realT dppdt;
  realT H;
  realT mu2;
  realT mu2High;
  realT mu2aLow;
  realT mu2Low;
  realT mu3High;
  realT mu3Low;


  DieterichStateData():
    base(),
    nextState(0),
    apFail(0),
    pinned(0),
    slowSlip(0),
    neighborInRuptureState(0),
    currentState(0),
    Areduced(0),
    stressReference(0),
    dstressdtDrive(0),
    dxsdtDrive(0),
    stressPin(0),
    stressShearReference(0),
    dstressShear(0),
    dstress(0),
    dstressShearDt(0),
    dstressShearDtDrive(0),
    dxs(0),
    pp(0),
    dppdt(0),
    H(0),
    mu2(0),
    mu2High(0),
    mu2aLow(0),
    mu2Low(0),
    mu3High(0),
    mu3Low(0)
  {}

  DieterichStateData( const DieterichStateData& source):
    base( source ),
    nextState(source.nextState),
    apFail(source.apFail),
    pinned(source.pinned),
    slowSlip(source.slowSlip),
    neighborInRuptureState(source.neighborInRuptureState),
    currentState(source.currentState),
    Areduced(source.Areduced),
    stressReference(source.stressReference),
    dstressdtDrive(source.dstressdtDrive),
    dxsdtDrive(source.dxsdtDrive),
    stressPin(source.stressPin),
    stressShearReference(source.stressShearReference),
    dstressShear(source.dstressShear),
    dstress(source.dstress),
    dstressShearDt(source.dstressShearDt),
    dstressShearDtDrive(source.dstressShearDtDrive),
    dxs(source.dxs),
    pp(source.pp),
    dppdt(source.dppdt),
    H(source.H),
    mu2(source.mu2),
    mu2High(source.mu2High),
    mu2aLow(source.mu2aLow),
    mu2Low(source.mu2Low),
    mu3High(source.mu3High),
    mu3Low(source.mu3Low)
  {}

  ~DieterichStateData() {}
  friend class ConstitutiveBase;
  friend class Dieterich;

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
    intVarCounts += 7;
    realVarCounts += 19;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    intNames.push_back("nextState");
    intNames.push_back("apFail");
    intNames.push_back("pinned");
    intNames.push_back("slowSlip");
    intNames.push_back("neighborInRuptureState");
    intNames.push_back("currentState");
    intNames.push_back("Areduced");
    realNames.push_back("stressReference");
    realNames.push_back("stressRateDrive");
    realNames.push_back("shearRateDrive");
    realNames.push_back("stressPin");
    realNames.push_back("shearStressReference");
    realNames.push_back("dstressShear");
    realNames.push_back("dstress");
    realNames.push_back("shearStressRate");
    realNames.push_back("shearStressRateDrive");
    realNames.push_back("dshearOverTimestep");
    realNames.push_back("porePressure");
    realNames.push_back("dPorePressureDt");
    realNames.push_back("H");
    realNames.push_back("frictionCoefficient2");
    realNames.push_back("frictionCoefficient2High");
    realNames.push_back("frictionCoefficient2aLow");
    realNames.push_back("frictionCoefficient2Low");
    realNames.push_back("frictionCoefficient3High");
    realNames.push_back("frictionCoefficient3Low");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    intOffsets["nextState"] = (char*)(&nextState) - (char*)this;
    intOffsets["apFail"] = (char*)(&apFail) - (char*)this;
    intOffsets["pinned"] = (char*)(&pinned) - (char*)this;
    intOffsets["slowSlip"] = (char*)(&slowSlip) - (char*)this;
    intOffsets["neighborInRuptureState"] = (char*)(&neighborInRuptureState) - (char*)this;
    intOffsets["currentState"] = (char*)(&currentState) - (char*)this;
    intOffsets["Areduced"] = (char*)(&Areduced) - (char*)this;
    realOffsets["stressReference"] = (char*)(&stressReference) - (char*)this;
    realOffsets["stressRateDrive"] = (char*)(&dstressdtDrive) - (char*)this;
    realOffsets["shearRateDrive"] = (char*)(&dxsdtDrive) - (char*)this;
    realOffsets["stressPin"] = (char*)(&stressPin) - (char*)this;
    realOffsets["shearStressReference"] = (char*)(&stressShearReference) - (char*)this;
    realOffsets["dstressShear"] = (char*)(&dstressShear) - (char*)this;
    realOffsets["dstress"] = (char*)(&dstress) - (char*)this;
    realOffsets["shearStressRate"] = (char*)(&dstressShearDt) - (char*)this;
    realOffsets["shearStressRateDrive"] = (char*)(&dstressShearDtDrive) - (char*)this;
    realOffsets["dshearOverTimestep"] = (char*)(&dxs) - (char*)this;
    realOffsets["porePressure"] = (char*)(&pp) - (char*)this;
    realOffsets["dPorePressureDt"] = (char*)(&dppdt) - (char*)this;
    realOffsets["H"] = (char*)(&H) - (char*)this;
    realOffsets["frictionCoefficient2"] = (char*)(&mu2) - (char*)this;
    realOffsets["frictionCoefficient2High"] = (char*)(&mu2High) - (char*)this;
    realOffsets["frictionCoefficient2aLow"] = (char*)(&mu2aLow) - (char*)this;
    realOffsets["frictionCoefficient2Low"] = (char*)(&mu2Low) - (char*)this;
    realOffsets["frictionCoefficient3High"] = (char*)(&mu3High) - (char*)this;
    realOffsets["frictionCoefficient3Low"] = (char*)(&mu3Low) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    intValues["nextState"] = nextState;
    intValues["apFail"] = apFail;
    intValues["pinned"] = pinned;
    intValues["slowSlip"] = slowSlip;
    intValues["neighborInRuptureState"] = neighborInRuptureState;
    intValues["currentState"] = currentState;
    intValues["Areduced"] = Areduced;
    realValues["stressReference"] = stressReference;
    realValues["stressRateDrive"] = dstressdtDrive;
    realValues["shearRateDrive"] = dxsdtDrive;
    realValues["stressPin"] = stressPin;
    realValues["shearStressReference"] = stressShearReference;
    realValues["dstressShear"] = dstressShear;
    realValues["dstress"] = dstress;
    realValues["shearStressRate"] = dstressShearDt;
    realValues["shearStressRateDrive"] = dstressShearDtDrive;
    realValues["dshearOverTimestep"] = dxs;
    realValues["porePressure"] = pp;
    realValues["dPorePressureDt"] = dppdt;
    realValues["H"] = H;
    realValues["frictionCoefficient2"] = mu2;
    realValues["frictionCoefficient2High"] = mu2High;
    realValues["frictionCoefficient2aLow"] = mu2aLow;
    realValues["frictionCoefficient2Low"] = mu2Low;
    realValues["frictionCoefficient3High"] = mu3High;
    realValues["frictionCoefficient3Low"] = mu3Low;
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
    (*(intVars[intVarCounts]))[elemNum] = nextState; intVarCounts += stride;
    (*(intVars[intVarCounts]))[elemNum] = apFail; intVarCounts += stride;
    (*(intVars[intVarCounts]))[elemNum] = pinned; intVarCounts += stride;
    (*(intVars[intVarCounts]))[elemNum] = slowSlip; intVarCounts += stride;
    (*(intVars[intVarCounts]))[elemNum] = neighborInRuptureState; intVarCounts += stride;
    (*(intVars[intVarCounts]))[elemNum] = currentState; intVarCounts += stride;
    (*(intVars[intVarCounts]))[elemNum] = Areduced; intVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = stressReference; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dstressdtDrive; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dxsdtDrive; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = stressPin; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = stressShearReference; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dstressShear; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dstress; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dstressShearDt; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dstressShearDtDrive; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dxs; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = pp; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dppdt; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = H; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mu2; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mu2High; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mu2aLow; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mu2Low; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mu3High; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mu3Low; realVarCounts += stride;
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
    nextState = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    apFail = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    pinned = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    slowSlip = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    neighborInRuptureState = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    currentState = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    Areduced = (*(intVars[intVarCounts]))[elemNum]; intVarCounts += stride;
    stressReference = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dstressdtDrive = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dxsdtDrive = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    stressPin = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    stressShearReference = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dstressShear = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dstress = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dstressShearDt = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dstressShearDtDrive = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dxs = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    pp = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dppdt = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    H = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mu2 = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mu2High = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mu2aLow = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mu2Low = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mu3High = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mu3Low = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
  }
  inline DieterichStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    nextState *= factor;
    apFail *= factor;
    pinned *= factor;
    slowSlip *= factor;
    neighborInRuptureState *= factor;
    currentState *= factor;
    Areduced *= factor;
    stressReference *= factor;
    dstressdtDrive *= factor;
    dxsdtDrive *= factor;
    stressPin *= factor;
    stressShearReference *= factor;
    dstressShear *= factor;
    dstress *= factor;
    dstressShearDt *= factor;
    dstressShearDtDrive *= factor;
    dxs *= factor;
    pp *= factor;
    dppdt *= factor;
    H *= factor;
    mu2 *= factor;
    mu2High *= factor;
    mu2aLow *= factor;
    mu2Low *= factor;
    mu3High *= factor;
    mu3Low *= factor;
    return *this;
  }

  inline DieterichStateData&
  operator=(const DieterichStateData& datum)
  {
    base::operator=(datum);
    nextState = datum.nextState;
    apFail = datum.apFail;
    pinned = datum.pinned;
    slowSlip = datum.slowSlip;
    neighborInRuptureState = datum.neighborInRuptureState;
    currentState = datum.currentState;
    Areduced = datum.Areduced;
    stressReference = datum.stressReference;
    dstressdtDrive = datum.dstressdtDrive;
    dxsdtDrive = datum.dxsdtDrive;
    stressPin = datum.stressPin;
    stressShearReference = datum.stressShearReference;
    dstressShear = datum.dstressShear;
    dstress = datum.dstress;
    dstressShearDt = datum.dstressShearDt;
    dstressShearDtDrive = datum.dstressShearDtDrive;
    dxs = datum.dxs;
    pp = datum.pp;
    dppdt = datum.dppdt;
    H = datum.H;
    mu2 = datum.mu2;
    mu2High = datum.mu2High;
    mu2aLow = datum.mu2aLow;
    mu2Low = datum.mu2Low;
    mu3High = datum.mu3High;
    mu3Low = datum.mu3Low;
    return *this;
  }

  inline DieterichStateData&
  operator+=(const DieterichStateData& datum)
  {
    base::operator+=(datum);
    nextState += datum.nextState;
    apFail += datum.apFail;
    pinned += datum.pinned;
    slowSlip += datum.slowSlip;
    neighborInRuptureState += datum.neighborInRuptureState;
    currentState += datum.currentState;
    Areduced += datum.Areduced;
    stressReference += datum.stressReference;
    dstressdtDrive += datum.dstressdtDrive;
    dxsdtDrive += datum.dxsdtDrive;
    stressPin += datum.stressPin;
    stressShearReference += datum.stressShearReference;
    dstressShear += datum.dstressShear;
    dstress += datum.dstress;
    dstressShearDt += datum.dstressShearDt;
    dstressShearDtDrive += datum.dstressShearDtDrive;
    dxs += datum.dxs;
    pp += datum.pp;
    dppdt += datum.dppdt;
    H += datum.H;
    mu2 += datum.mu2;
    mu2High += datum.mu2High;
    mu2aLow += datum.mu2aLow;
    mu2Low += datum.mu2Low;
    mu3High += datum.mu3High;
    mu3Low += datum.mu3Low;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, DieterichStateData& p0, DieterichStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, nextState, p0.nextState, p1.nextState);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, apFail, p0.apFail, p1.apFail);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, pinned, p0.pinned, p1.pinned);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, slowSlip, p0.slowSlip, p1.slowSlip);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, neighborInRuptureState, p0.neighborInRuptureState, p1.neighborInRuptureState);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, currentState, p0.currentState, p1.currentState);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, Areduced, p0.Areduced, p1.Areduced);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressReference, p0.stressReference, p1.stressReference);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dstressdtDrive, p0.dstressdtDrive, p1.dstressdtDrive);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dxsdtDrive, p0.dxsdtDrive, p1.dxsdtDrive);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressPin, p0.stressPin, p1.stressPin);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, stressShearReference, p0.stressShearReference, p1.stressShearReference);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dstressShear, p0.dstressShear, p1.dstressShear);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dstress, p0.dstress, p1.dstress);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dstressShearDt, p0.dstressShearDt, p1.dstressShearDt);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dstressShearDtDrive, p0.dstressShearDtDrive, p1.dstressShearDtDrive);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dxs, p0.dxs, p1.dxs);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, pp, p0.pp, p1.pp);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dppdt, p0.dppdt, p1.dppdt);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, H, p0.H, p1.H);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu2, p0.mu2, p1.mu2);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu2High, p0.mu2High, p1.mu2High);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu2aLow, p0.mu2aLow, p1.mu2aLow);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu2Low, p0.mu2Low, p1.mu2Low);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu3High, p0.mu3High, p1.mu3High);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu3Low, p0.mu3Low, p1.mu3Low);

  }

  void MapFromRegion(const DieterichStateData& p0, const DieterichStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.nextState, p1.nextState, fct0, fct1, nextState);
    GeometryUtilities::MapFromRegion(p0.apFail, p1.apFail, fct0, fct1, apFail);
    GeometryUtilities::MapFromRegion(p0.pinned, p1.pinned, fct0, fct1, pinned);
    GeometryUtilities::MapFromRegion(p0.slowSlip, p1.slowSlip, fct0, fct1, slowSlip);
    GeometryUtilities::MapFromRegion(p0.neighborInRuptureState, p1.neighborInRuptureState, fct0, fct1, neighborInRuptureState);
    GeometryUtilities::MapFromRegion(p0.currentState, p1.currentState, fct0, fct1, currentState);
    GeometryUtilities::MapFromRegion(p0.Areduced, p1.Areduced, fct0, fct1, Areduced);
    GeometryUtilities::MapFromRegion(p0.stressReference, p1.stressReference, fct0, fct1, stressReference);
    GeometryUtilities::MapFromRegion(p0.dstressdtDrive, p1.dstressdtDrive, fct0, fct1, dstressdtDrive);
    GeometryUtilities::MapFromRegion(p0.dxsdtDrive, p1.dxsdtDrive, fct0, fct1, dxsdtDrive);
    GeometryUtilities::MapFromRegion(p0.stressPin, p1.stressPin, fct0, fct1, stressPin);
    GeometryUtilities::MapFromRegion(p0.stressShearReference, p1.stressShearReference, fct0, fct1, stressShearReference);
    GeometryUtilities::MapFromRegion(p0.dstressShear, p1.dstressShear, fct0, fct1, dstressShear);
    GeometryUtilities::MapFromRegion(p0.dstress, p1.dstress, fct0, fct1, dstress);
    GeometryUtilities::MapFromRegion(p0.dstressShearDt, p1.dstressShearDt, fct0, fct1, dstressShearDt);
    GeometryUtilities::MapFromRegion(p0.dstressShearDtDrive, p1.dstressShearDtDrive, fct0, fct1, dstressShearDtDrive);
    GeometryUtilities::MapFromRegion(p0.dxs, p1.dxs, fct0, fct1, dxs);
    GeometryUtilities::MapFromRegion(p0.pp, p1.pp, fct0, fct1, pp);
    GeometryUtilities::MapFromRegion(p0.dppdt, p1.dppdt, fct0, fct1, dppdt);
    GeometryUtilities::MapFromRegion(p0.H, p1.H, fct0, fct1, H);
    GeometryUtilities::MapFromRegion(p0.mu2, p1.mu2, fct0, fct1, mu2);
    GeometryUtilities::MapFromRegion(p0.mu2High, p1.mu2High, fct0, fct1, mu2High);
    GeometryUtilities::MapFromRegion(p0.mu2aLow, p1.mu2aLow, fct0, fct1, mu2aLow);
    GeometryUtilities::MapFromRegion(p0.mu2Low, p1.mu2Low, fct0, fct1, mu2Low);
    GeometryUtilities::MapFromRegion(p0.mu3High, p1.mu3High, fct0, fct1, mu3High);
    GeometryUtilities::MapFromRegion(p0.mu3Low, p1.mu3Low, fct0, fct1, mu3Low);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class Dieterich: public RateAndStateIntermediate
{
public:

  typedef DieterichParameterData ParameterClass;
  typedef DieterichStateData     StateClass;

  typedef Array1dT<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "Dieterich"; }

  Dieterich();
  virtual ~Dieterich();

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
  { SetVariableParametersFromDerived<Dieterich>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<Dieterich>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<Dieterich>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<Dieterich>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<Dieterich>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<Dieterich>( num ); }
 
  void GetVariableNames( sArray1d& intVars, sArray1d& realVars, sArray1d& R1TensorVars, sArray1d& R2TensorVars, sArray1d& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<Dieterich>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<Dieterich>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<Dieterich>(name, type ); }
  
  bool GetStateValues( const std::string& name, rArray1d& values ) const
  { return GetStateValuesFromDerived<Dieterich>(name, values); }

  bool GetParameterValues( const std::string& name, rArray1d& values ) const
  { return GetParameterValuesFromDerived<Dieterich>(name, values); }

  bool SetStateValues( const std::string& name, const rArray1d& values )
  { return SetStateValuesFromDerived<Dieterich>(name, values); }

  bool SetParameterValues( const std::string& name, const rArray1d& values )
  { return SetParameterValuesFromDerived<Dieterich>(name, values); }

  virtual void Serialize( Array1dT<iArray1d*>& intVars, Array1dT<rArray1d*>& realVars, Array1dT<Array1dT<R1Tensor>*>& R1Vars, Array1dT<Array1dT<R2Tensor>*>& R2Vars, Array1dT<Array1dT<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<Dieterich>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const Array1dT<iArray1d*>& intVars, const Array1dT<rArray1d*>& realVars, const Array1dT<Array1dT<R1Tensor>*>& R1Vars, const Array1dT<Array1dT<R2Tensor>*>& R2Vars, const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<Dieterich>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<Dieterich>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<Dieterich>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<Dieterich>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<Dieterich>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<Dieterich>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<Dieterich>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  
  realT
  Advance(EarthquakeSimulation::EarthquakeSimulationTimestep& timestep,
                     const realT dtGlobal,
                     const realT dxsdtEQ,
                     const Array1dT<R1Tensor>& centers,
                     const Table<4, realT>* porePressure);

  void
  SetStressRates(const realT time,
                            const EarthquakeSimulation::BoundaryElementDataManagerT& KGD,
                            const Table<4, realT>* porePressure,
                            const Array1dT<R1Tensor> centers,
                            const bool resetStress);

  void
  SetStress(const realT time, const R1Tensor& center, const Table<4, realT>* porePressure, DieterichStateData& matState);

  void
  ResetStates(const gArray1d& localToGlobal,
                         const EarthquakeSimulation::BoundaryElementDataManagerT& KGD,
                         const realT time,
                         const realT dxsdtEQ,
                         const globalIndex maxGlobalIndex,
                         const Array1dT<R1Tensor>& centers,
                         const Table<4, realT>* porePressure);

  realT
  NextTransitionTime( EarthquakeSimulation::EarthquakeSimulationTimestep& timestep, 
                                 const realT dxsdtEQ,
                                 const Array1dT<R1Tensor>& centers,
                                 const Table<4, realT>* porePressure);

  bool
  PostTimestepSynchronization( EarthquakeSimulation::EarthquakeSimulationTimestep& timestep,
                                          EarthquakeSimulation::TransitionState& current,
                                          realT& dxsdt);

  void
  TransitionLockedToNucleate(const localIndex iRupture, const realT dxsdtEQ);

  realT
  TransitionNucleateToRupture(const localIndex iRupture, const realT dxsdtEQ);

  void
  TransitionUpdateContributions(const localIndex iRupture,
                                           const globalIndex giRupture,
                                           const bool isQuiescent,
                                           const EarthquakeSimulation::BoundaryElementDataManagerT& KGD,
                                           const realT dxsdt,
                                           const int checkLevel);

  void
  TransitionNucleateToSlowSlip2A(const localIndex iRupture);

  realT
  TransitionRuptureToLocked(const localIndex iRupture, const realT dxsdtEQ);

  realT
  TransitionSlowSlip2AToLocked(const localIndex iRupture, const realT dxsdtEQ);

  realT
  TransitionSlowSlip2AToSlowSlip2B(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep);

  realT
  TransitionSlowSlip2BToSlowSlip2B(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep);

  realT
  TransitionSlowSlip2BToSlowSlip2A(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep);

  realT
  TransitionSlowSlip2BToSlowSlip2C(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep);

  realT
  TransitionSlowSlip2CToSlowSlip2B(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep);

  realT
  TransitionCreepToCreep(const localIndex iRupture, const realT dxsdtEQ, const realT dMuCreep);

  void
  TransitionLowSigma(const localIndex iRupture,
                                const globalIndex giRupture,
                                const EarthquakeSimulation::BoundaryElementDataManagerT& KGD);


protected:
  virtual void
  UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                             InterfaceBaseStateData& matStateBase) const;

  bool
  NextTransitionTime( const DieterichParameterData& matParams,
                                 DieterichStateData& matState,
                                 EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep,
                                 const realT dxsdtEQ,
                                 const R1Tensor& center,
                                 const localIndex i) const;

  virtual realT
  ShearStrength(const InterfaceBaseParameterData& ,
                           InterfaceBaseStateData& matStateBase) const;

  realT
  IncrementalShearDisplacement(const EarthquakeSimulation::TransitionState state,
                                          const realT dt,
                                          const realT dxsdtEQ,
                                          const DieterichParameterData& matParams,
                                          DieterichStateData& matState) const;

  static realT
  IncrementalShearDisplacement(const EarthquakeSimulation::TransitionState state,
                                          const realT dt,
                                          const realT sigma0,
                                          const realT sigmaDot0,
                                          const realT tau0,
                                          const realT tauDot0,
                                          const realT slipRateEQ,
                                          const realT Dc,
                                          const realT alpha,
                                          const realT KShearSelf,
                                          const realT A_,
                                          const realT B,
                                          realT& ddot,
                                          realT& H);

  static realT
  Theta(const EarthquakeSimulation::TransitionState state,
                   const DieterichParameterData& matParams,
                   DieterichStateData& matState);

  static realT
  Theta(const EarthquakeSimulation::TransitionState state,
                   const realT dt,
                   const realT sigmaDot,
                   const realT sigma0,
                   const realT alpha,
                   const realT B,
                   const realT Dc,
                   const realT dd,
                   const realT theta0);

  static realT
  ThetaLock(const realT dt,
                       const realT deffectiveStressDt,
                       const realT effectiveStress,
                       const realT alpha,
                       const realT B,
                       const realT theta0);

  static realT
  ThetaNucleate(const realT tau,
                           const realT sigma,
                           const realT ddot,
                           const realT mu0,
                           const realT ddotStar,
                           const realT A_,
                           const realT B,
                           const realT Dc);

  static realT
  ThetaNucleateAlternate( const realT theta0,
                                     const realT d,
                                     const realT Dc,
                                     const realT alpha,
                                     const realT B,
                                     const realT sigma0,
                                     const realT sigma);

  static void
  TransitionTimeAPriori(const realT time,
                                   const realT timeFail,
                                   const EarthquakeSimulation::FailureState apFail,
                                   const localIndex local,
                                   EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  TransitionTimeLocked(const realT dEffectiveStressDt,
                                  const realT effectiveStress,
                                  const realT dStressShearDt,
                                  const realT stressShear,
                                  const realT dShearSlipDtStar,
                                  const realT alpha, const realT A_, const realT B,
                                  const realT Dc, const realT frictionCoefficient0,
                                  const realT theta0, const localIndex local,
                                  EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  TransitionTimeNucleate(const realT dEffectiveStressDt,
                                    const realT effectiveStress,
                                    const realT dStressShearDt,
                                    const realT stressShear,
                                    const realT dShearSlipDt,
                                    const realT dShearSlipDtStar,
                                    const realT dShearSlipDtEQ,
                                    const realT alpha,
                                    const realT A_,
                                    const realT B,
                                    const realT Dc,
                                    const realT KShearSelf,
                                    const realT theta0,
                                    const int slowSlip,
                                    const int neighborInRuptureState,
                                    realT& H,
                                    const localIndex local,
                                    EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  TransitionTimeRupture(const realT dEffectiveStressDt,
                                   const realT effectiveStress,
                                   const realT dStressShearDt,
                                   const realT stressShear,
                                   const realT frictionCoefficient2,
                                   const localIndex local,
                                   EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  TransitionTimeSlowSlipA(const realT dEffectiveStressDt,
                                     const realT effectiveStress,
                                     const realT dStressShearDt,
                                     const realT stressShear,
                                     const realT frictionCoefficient2High,
                                     const realT frictionCoefficient2aLow,
                                     const localIndex local,
                                     EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  TransitionTimeSlowSlipB(const realT dEffectiveStressDt,
                                     const realT effectiveStress,
                                     const realT dStressShearDt,
                                     const realT stressShear,
                                     const realT frictionCoefficient0,
                                     const realT frictionCoefficient2High,
                                     const realT frictionCoefficient2Low,
                                     const realT dShearSlipDtAB,
                                     const realT dShearSlipDtStar,
                                     const realT dShearSlipDtEQ,
                                     const realT A_,
                                     const localIndex local,
                                     EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  TransitionTimeSlowSlipC(const realT dEffectiveStressDt,
                                     const realT effectiveStress,
                                     const realT dStressShearDt,
                                     const realT stressShear,
                                     const realT frictionCoefficient2High,
                                     const realT frictionCoefficient2Low,
                                     const localIndex local,
                                     EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  TransitionTimeCreep( const realT dEffectiveStressDt,
                                  const realT effectiveStress,
                                  const realT dStressShearDt,
                                  const realT stressShear,
                                  const realT frictionCoefficient3Low,
                                  const realT frictionCoefficient3High,
                                  const localIndex local,
                                  EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  LowStressCheck( const realT dEffectiveStressDt,
                             const realT effectiveStress,
                             const realT stressNormalPin,
                             const realT dStressShearDt,
                             const realT stressShear,
                             const localIndex local,
                             EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static void
  HighThetaCheck( const realT dt,
                             const realT theta,
                             const realT maxThetaPin,
                             const localIndex local,
                             EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static bool
  PorePressureChangeCheck( const Table<4, realT>* porePressure,
                                      EarthquakeSimulation::EarthquakeSimulationTimestep& ruptureTimestep);

  static realT
  f(const realT dt,
               const realT* params);

  static realT
  zbrent(const realT* params,
                    const realT lower,
                    const realT upper,
                    const realT tol);

  static void
  APrioriFail(const realT dxsdtEQ,
                         const DieterichParameterData& matParams,
                         DieterichStateData& matState);

  realT
  A(const localIndex iRupture) const;

  static realT
  A(const DieterichParameterData& matParams,
               const DieterichStateData& matState);

  void
  CheckA(const localIndex i);

  void
  CheckA(const DieterichParameterData& matParams,
                    DieterichStateData& matState) const;

  bool
  UpdateStressDriveRates(const globalIndex giRupture,
                                    const localIndex iCurrent,
                                    const realT ddotDrive,
                                    const EarthquakeSimulation::BoundaryElementDataManagerT& KGD);

  bool
  UpdateStressRates(const globalIndex giRupture,
                               const localIndex iCurrent,
                               const realT dxsdt,
                               const EarthquakeSimulation::BoundaryElementDataManagerT& KGD);


private:
  Dieterich(const Dieterich&);
  Dieterich& operator=(const Dieterich&);
  

};
#endif /* DIETERICH_H_ */
