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

#ifndef HEUZEJOINT_H_
#define HEUZEJOINT_H_

#include "Utilities/GeometryUtilities.h"
#include "JointIntermediate.h"

/*
 * HeuzeJoint.h
 *
 *  Created on: Tue Jan  7 22:46:45 PST 2014
 *      Author: johnson346, settgast
 */
 


//**********************************************************************************************************************
//**********************************************************************************************************************


class HeuzeJointParameterData : public JointIntermediateParameterData
{

public:

  typedef JointIntermediateParameterData base;


  HeuzeJointParameterData():
    base()
  {}

  HeuzeJointParameterData( const HeuzeJointParameterData& source):
    base( source )
  {}

  ~HeuzeJointParameterData() {}
  friend class ConstitutiveBase;
  friend class HeuzeJoint;

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
  }
  inline HeuzeJointParameterData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline HeuzeJointParameterData&
  operator=(const HeuzeJointParameterData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline HeuzeJointParameterData&
  operator+=(const HeuzeJointParameterData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, HeuzeJointParameterData& p0, HeuzeJointParameterData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const HeuzeJointParameterData& p0, const HeuzeJointParameterData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }
  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  {
    JointIntermediateParameterData::ReadXML( node );
    PostReadXML( node );

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************


class HeuzeJointStateData : public JointIntermediateStateData
{

public:

  typedef JointIntermediateStateData base;


  HeuzeJointStateData():
    base()
  {}

  HeuzeJointStateData( const HeuzeJointStateData& source):
    base( source )
  {}

  ~HeuzeJointStateData() {}
  friend class ConstitutiveBase;
  friend class HeuzeJoint;

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
  inline HeuzeJointStateData&
  operator*=(const realT factor)
  {
    base::operator*=(factor);
    return *this;
  }

  inline HeuzeJointStateData&
  operator=(const HeuzeJointStateData& datum)
  {
    base::operator=(datum);
    return *this;
  }

  inline HeuzeJointStateData&
  operator+=(const HeuzeJointStateData& datum)
  {
    base::operator+=(datum);
    return *this;
  }
  void MapToRegion(const realT fctNormal, const realT fct0,
                   const realT fct1, HeuzeJointStateData& p0, HeuzeJointStateData& p1)
  {
    base::MapToRegion(fctNormal, fct0, fct1, p0, p1);

  }

  void MapFromRegion(const HeuzeJointStateData& p0, const HeuzeJointStateData& p1, const realT fct0,
                     const realT fct1)
  {
    base::MapFromRegion(p0, p1, fct0, fct1);

  }



};

//**********************************************************************************************************************
//**********************************************************************************************************************
class HeuzeJoint: public JointIntermediate
{
public:

  typedef HeuzeJointParameterData ParameterClass;
  typedef HeuzeJointStateData     StateClass;

  typedef Array1dT<ParameterClass> ParameterArrayType;
  typedef Array2dT<StateClass>     StateArrayType;


  StateArrayType m_stateData;
  ParameterArrayType m_parameterData;

  localIndex NumStateIndex0() const { return m_stateData.Dimension(0); }
  localIndex NumStateIndex1() const { return m_stateData.Dimension(1); }

  localIndex NumParameterIndex0() const { return m_parameterData.size(); }
  localIndex NumParameterIndex1() const { return 1; }

  static std::string Name() { return "HeuzeJoint"; }

  HeuzeJoint();
  virtual ~HeuzeJoint();

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
  { SetVariableParametersFromDerived<HeuzeJoint>(varParams, newSize); }

  virtual void ReadXML( TICPP::HierarchicalDataNode& node )
  { ReadXMLFromDerived<HeuzeJoint>( node ); }

  virtual void resize( const localIndex num )
  { ResizeFromDerived<HeuzeJoint>( num ); }

  virtual void resize( const localIndex num0, const localIndex num1 )
  {
    m_stateData.resize2(num0, num1);
    ResizeFromDerived<HeuzeJoint>( num0 );
  }
  
  virtual void insert( const localIndex num )
  { InsertFromDerived<HeuzeJoint>( num ); }

  virtual void erase( const localIndex num )
  { EraseFromDerived<HeuzeJoint>( num ); }
 
  void GetVariableNames( sArray1d& intVars, sArray1d& realVars, sArray1d& R1TensorVars, sArray1d& R2TensorVars, sArray1d& R2SymTensorVars ) const
  { GetVariableNamesFromDerived<HeuzeJoint>(intVars, realVars, R1TensorVars, R2TensorVars, R2SymTensorVars ); }

  size_t GetStateOffset( const std::string& name, const int type ) const
  { return GetStateOffsetFromDerived<HeuzeJoint>(name, type); }
  
  size_t GetParameterOffset( const std::string& name, const int type ) const
  { return GetParameterOffsetFromDerived<HeuzeJoint>(name, type ); }
  
  bool GetStateValues( const std::string& name, rArray1d& values ) const
  { return GetStateValuesFromDerived<HeuzeJoint>(name, values); }

  bool GetParameterValues( const std::string& name, rArray1d& values ) const
  { return GetParameterValuesFromDerived<HeuzeJoint>(name, values); }

  bool SetStateValues( const std::string& name, const rArray1d& values )
  { return SetStateValuesFromDerived<HeuzeJoint>(name, values); }

  bool SetParameterValues( const std::string& name, const rArray1d& values )
  { return SetParameterValuesFromDerived<HeuzeJoint>(name, values); }

  virtual void Serialize( Array1dT<iArray1d*>& intVars, Array1dT<rArray1d*>& realVars, Array1dT<Array1dT<R1Tensor>*>& R1Vars, Array1dT<Array1dT<R2Tensor>*>& R2Vars, Array1dT<Array1dT<R2SymTensor>*>& R2SymVars ) const
  { SerializeFromDerived<HeuzeJoint>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual void Deserialize( const Array1dT<iArray1d*>& intVars, const Array1dT<rArray1d*>& realVars, const Array1dT<Array1dT<R1Tensor>*>& R1Vars, const Array1dT<Array1dT<R2Tensor>*>& R2Vars, const Array1dT<Array1dT<R2SymTensor>*>& R2SymVars  )
  { DeserializeFromDerived<HeuzeJoint>( intVars, realVars, R1Vars, R2Vars, R2SymVars ); }

  virtual unsigned int Pack( const lArray1d& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<HeuzeJoint>( localIndices, buffer, doBufferPacking ); }
  unsigned int Pack( const lSet& localIndices, bufvector& buffer, const bool doBufferPacking )
  { return PackFromDerived<HeuzeJoint>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Pack( const lArray1d& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<HeuzeJoint>( localIndices, buffer, doBufferPacking ); }
  virtual unsigned int Pack( const lSet& localIndices, char*& buffer, const bool doBufferPacking )
  { return PackFromDerived<HeuzeJoint>( localIndices, buffer, doBufferPacking ); }

  virtual unsigned int Unpack( const lArray1d& localIndices, const char*& buffer )
  { return UnpackFromDerived<HeuzeJoint>( localIndices, buffer ); }

  virtual unsigned int Unpack( const lSet& localIndices, const char*& buffer )
  { return UnpackFromDerived<HeuzeJoint>( localIndices, buffer ); }

  const StateClass* StateData( const localIndex index0, const localIndex index1 ) const
  { return &( m_stateData(index0,index1) );  }
  StateClass* StateData( const localIndex index0, const localIndex index1 )
  { return &( m_stateData(index0,index1) );  }

  const ParameterClass* ParameterData( const localIndex index ) const
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  ParameterClass* ParameterData( const localIndex index )
  { return &( m_parameterData(hasVariableParams ? index : 0) ); }
  
  virtual void
  Initialize( const localIndex index,
                          const realT stressNormal,
                          const realT stressShear);


protected:
  virtual void
  UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                     InterfaceBaseStateData& matStateBase) const;

  virtual realT
  DilationalStressIncrement( const InterfaceBaseParameterData& matParamsBase,
                                         InterfaceBaseStateData& matStateBase,
                                         const realT dxs) const;

  virtual realT
  NormalApproachAtInitialization( const realT stressNormal,
                                              const realT knCoefficientElastic,
                                              const realT knPlasticLoading,
                                              const realT knPlasticUnloading,
                                              realT& kcurrent,
                                              realT& normalGap,
                                              const realT normalApproachNormalYield) const;

  static realT
  ZerosFunctionElasticElastic(const realT x, const realT* params );


private:
  HeuzeJoint(const HeuzeJoint&);
  HeuzeJoint& operator=(const HeuzeJoint&);
  

};
#endif /* HEUZEJOINT_H_ */
