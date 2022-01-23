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

#ifndef RATEANDSTATEINTERMEDIATE_H_
#define RATEANDSTATEINTERMEDIATE_H_

#include "Utilities/GeometryUtilities.h"
#include "Constitutive/Interface/PenaltyCoulombIntermediate.h"

/*
 * RateAndStateIntermediate.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class RateAndStateIntermediateParameterData : public PenaltyCoulombIntermediateParameterData
{

public:

  typedef PenaltyCoulombIntermediateParameterData base;
  realT A;
  realT B;
  realT vstar;
  realT thetastar;
  realT Dc;
  realT alpha;


  RateAndStateIntermediateParameterData():
    base(),
    A(0),
    B(0),
    vstar(0),
    thetastar(0),
    Dc(0),
    alpha(0)
  {}

  RateAndStateIntermediateParameterData( const RateAndStateIntermediateParameterData& source):
    base( source ),
    A(source.A),
    B(source.B),
    vstar(source.vstar),
    thetastar(source.thetastar),
    Dc(source.Dc),
    alpha(source.alpha)
  {}

  ~RateAndStateIntermediateParameterData() {}
  friend class ConstitutiveBase;
  friend class RateAndStateIntermediate;

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
    realVarCounts += 6;

  }

  static void GetVariableNames( sArray1d& intNames,
                                sArray1d& realNames,
                                sArray1d& R1TensorNames,
                                sArray1d& R2TensorNames,
                                sArray1d& R2SymTensorNames )
  {
    base::GetVariableNames( intNames, realNames, R1TensorNames, R2TensorNames, R2SymTensorNames);
    realNames.push_back("A");
    realNames.push_back("B");
    realNames.push_back("shearRateStar");
    realNames.push_back("stateStar");
    realNames.push_back("Dc");
    realNames.push_back("alpha");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["A"] = (char*)(&A) - (char*)this;
    realOffsets["B"] = (char*)(&B) - (char*)this;
    realOffsets["shearRateStar"] = (char*)(&vstar) - (char*)this;
    realOffsets["stateStar"] = (char*)(&thetastar) - (char*)this;
    realOffsets["Dc"] = (char*)(&Dc) - (char*)this;
    realOffsets["alpha"] = (char*)(&alpha) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["A"] = A;
    realValues["B"] = B;
    realValues["shearRateStar"] = vstar;
    realValues["stateStar"] = thetastar;
    realValues["Dc"] = Dc;
    realValues["alpha"] = alpha;
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
    (*(realVars[realVarCounts]))[index] = A; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = B; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = vstar; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = thetastar; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = Dc; realVarCounts++;
    (*(realVars[realVarCounts]))[index] = alpha; realVarCounts++;
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
    A = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    B = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    vstar = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    thetastar = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    Dc = (*(realVars[realVarCounts]))[index]; realVarCounts++;
    alpha = (*(realVars[realVarCounts]))[index]; realVarCounts++;
  }
  inline RateAndStateIntermediateParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    A *= factor;
    B *= factor;
    vstar *= factor;
    thetastar *= factor;
    Dc *= factor;
    alpha *= factor;
    return *this;
  }

  inline RateAndStateIntermediateParameterData&
  operator=(const RateAndStateIntermediateParameterData& datum)
  {
    base::operator=(datum);
    A = datum.A;
    B = datum.B;
    vstar = datum.vstar;
    thetastar = datum.thetastar;
    Dc = datum.Dc;
    alpha = datum.alpha;
    return *this;
  }

  inline RateAndStateIntermediateParameterData&
  operator+=(const RateAndStateIntermediateParameterData& datum)
  {
    base::operator+=(datum);
    A += datum.A;
    B += datum.B;
    vstar += datum.vstar;
    thetastar += datum.thetastar;
    Dc += datum.Dc;
    alpha += datum.alpha;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, RateAndStateIntermediateParameterData& p0, RateAndStateIntermediateParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, A, p0.A, p1.A);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, B, p0.B, p1.B);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, vstar, p0.vstar, p1.vstar);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, thetastar, p0.thetastar, p1.thetastar);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, Dc, p0.Dc, p1.Dc);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, alpha, p0.alpha, p1.alpha);

  }

  void MapFromRegion(const RateAndStateIntermediateParameterData& p0, const RateAndStateIntermediateParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.A, p1.A, fct0, fct1, A);
    GeometryUtilities::MapFromRegion(p0.B, p1.B, fct0, fct1, B);
    GeometryUtilities::MapFromRegion(p0.vstar, p1.vstar, fct0, fct1, vstar);
    GeometryUtilities::MapFromRegion(p0.thetastar, p1.thetastar, fct0, fct1, thetastar);
    GeometryUtilities::MapFromRegion(p0.Dc, p1.Dc, fct0, fct1, Dc);
    GeometryUtilities::MapFromRegion(p0.alpha, p1.alpha, fct0, fct1, alpha);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    PenaltyCoulombIntermediateParameterData::ReadXML( node );
    A = node.GetAttributeOrDefault("A", 0.0);
    B = node.GetAttributeOrDefault("B", 0.0);
    vstar = node.GetAttributeOrDefault("shearRateStar", 0.0);
    thetastar = node.GetAttributeOrDefault("stateStar", 0.0);
    Dc = node.GetAttributeOrDefault("Dc", 0.0);
    alpha = node.GetAttributeOrDefault("alpha", 0.0);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class RateAndStateIntermediateStateData : public PenaltyCoulombIntermediateStateData
{

public:

  typedef PenaltyCoulombIntermediateStateData base;
  realT theta;
  realT mu;
  realT dstressdt;


  RateAndStateIntermediateStateData():
    base(),
    theta(0),
    mu(0),
    dstressdt(0)
  {}

  RateAndStateIntermediateStateData( const RateAndStateIntermediateStateData& source):
    base( source ),
    theta(source.theta),
    mu(source.mu),
    dstressdt(source.dstressdt)
  {}

  ~RateAndStateIntermediateStateData() {}
  friend class ConstitutiveBase;
  friend class RateAndStateIntermediate;

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
    realNames.push_back("state");
    realNames.push_back("currentFrictionCoefficient");
    realNames.push_back("stressRate");
  }

  virtual void GetVariableOffsets( std::map<std::string, size_t>& intOffsets,
                                std::map<std::string, size_t>& realOffsets,
                                std::map<std::string, size_t>& R1TensorOffsets,
                                std::map<std::string, size_t>& R2TensorOffsets,
                                std::map<std::string, size_t>& R2SymTensorOffsets ) const
  {
    base::GetVariableOffsets( intOffsets, realOffsets, R1TensorOffsets, R2TensorOffsets, R2SymTensorOffsets);
    realOffsets["state"] = (char*)(&theta) - (char*)this;
    realOffsets["currentFrictionCoefficient"] = (char*)(&mu) - (char*)this;
    realOffsets["stressRate"] = (char*)(&dstressdt) - (char*)this;
  }

  virtual void GetVariableValues( std::map<std::string, int>& intValues,
                                std::map<std::string, realT>& realValues,
                                std::map<std::string, R1Tensor>& R1TensorValues,
                                std::map<std::string, R2Tensor>& R2TensorValues,
                                std::map<std::string, R2SymTensor>& R2SymTensorValues )
  {
    base::GetVariableValues( intValues, realValues, R1TensorValues, R2TensorValues, R2SymTensorValues);
    realValues["state"] = theta;
    realValues["currentFrictionCoefficient"] = mu;
    realValues["stressRate"] = dstressdt;
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
    (*(realVars[realVarCounts]))[elemNum] = theta; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = mu; realVarCounts += stride;
    (*(realVars[realVarCounts]))[elemNum] = dstressdt; realVarCounts += stride;
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
    theta = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    mu = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
    dstressdt = (*(realVars[realVarCounts]))[elemNum]; realVarCounts += stride;
  }
  inline RateAndStateIntermediateStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    theta *= factor;
    mu *= factor;
    dstressdt *= factor;
    return *this;
  }

  inline RateAndStateIntermediateStateData&
  operator=(const RateAndStateIntermediateStateData& datum)
  {
    base::operator=(datum);
    theta = datum.theta;
    mu = datum.mu;
    dstressdt = datum.dstressdt;
    return *this;
  }

  inline RateAndStateIntermediateStateData&
  operator+=(const RateAndStateIntermediateStateData& datum)
  {
    base::operator+=(datum);
    theta += datum.theta;
    mu += datum.mu;
    dstressdt += datum.dstressdt;
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, RateAndStateIntermediateStateData& p0, RateAndStateIntermediateStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, theta, p0.theta, p1.theta);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, mu, p0.mu, p1.mu);
    GeometryUtilities::MapToRegion(fctNormal, fct0, fct1, dstressdt, p0.dstressdt, p1.dstressdt);

  }

  void MapFromRegion(const RateAndStateIntermediateStateData& p0, const RateAndStateIntermediateStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);
    GeometryUtilities::MapFromRegion(p0.theta, p1.theta, fct0, fct1, theta);
    GeometryUtilities::MapFromRegion(p0.mu, p1.mu, fct0, fct1, mu);
    GeometryUtilities::MapFromRegion(p0.dstressdt, p1.dstressdt, fct0, fct1, dstressdt);

  }



};


//**********************************************************************************************************************
//**********************************************************************************************************************


class RateAndStateIntermediate: public PenaltyCoulombIntermediate
{
public:
  
  typedef RateAndStateIntermediateParameterData ParameterClass;
  typedef RateAndStateIntermediateStateData     StateClass;
  

  RateAndStateIntermediate( const int paramSize, const int stateSize );

  virtual ~RateAndStateIntermediate();
  
  virtual void ReadXML( TICPP::HierarchicalDataNode& node ) = 0;

  virtual void resize( const localIndex num ) = 0;
  
  virtual void resize( const localIndex num0,
                       const localIndex num1 ) = 0;
  
  virtual void insert( const localIndex num ) = 0;

  virtual void erase( const localIndex num ) = 0;
  
  virtual const RateAndStateIntermediateStateData* StateData( const localIndex index0,
                                                  const localIndex index1 ) const = 0;
  virtual       RateAndStateIntermediateStateData* StateData( const localIndex index0,
                                                  const localIndex index1 )  = 0;

  virtual const RateAndStateIntermediateParameterData* ParameterData( const localIndex index ) const = 0;
  virtual       RateAndStateIntermediateParameterData* ParameterData( const localIndex index ) = 0;
  
  inline void IncrementPtr( const RateAndStateIntermediateStateData* ptr ) const
  {
    ptr = reinterpret_cast<const RateAndStateIntermediateStateData*>( reinterpret_cast<const char*>(ptr) + m_stateSize );
  }

  inline void IncrementPtr( const RateAndStateIntermediateParameterData* ptr ) const
  {
    ptr = reinterpret_cast<const RateAndStateIntermediateParameterData*>( reinterpret_cast<const char*>(ptr) + m_paramSize );
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

  realT
  CriticalLength(const realT G,
                                           const InterfaceBaseParameterData& matParamsBase,
                                           InterfaceBaseStateData& matStateBase) const;


protected:
  virtual void
  UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                            InterfaceBaseStateData& matStateBase) const;

  virtual realT
  ShearStrength(const InterfaceBaseParameterData& ,
                                          InterfaceBaseStateData& matStateBase) const;

  
private:
  RateAndStateIntermediate();
  RateAndStateIntermediate( const RateAndStateIntermediate& );
  RateAndStateIntermediate& operator=( const RateAndStateIntermediate& );
  
  
};
#endif /* RATEANDSTATEINTERMEDIATE_H_ */
