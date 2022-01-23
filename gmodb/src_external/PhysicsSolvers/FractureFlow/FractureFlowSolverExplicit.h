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
/**
 * @file FractureFlowSolverExplicit.h
 * @author hao1
 * @date Apr. 21, 2014
 */

#ifndef FRACTURE_FLOW_SOLVER_EXPLICIT_H_
#define FRACTURE_FLOW_SOLVER_EXPLICIT_H_

#include "PhysicsSolvers/SolverBase.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "ObjectManagers/ProblemManagerT.h"

#include "Common/Common.h"
#include "Utilities/TrilinosUtilities.h"
#include "Utilities/RCVSparse.h"

#include "PhysicsSolvers/PhysicsSolverStrings.h"

#include <set>

#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

using namespace std;

class FractureFlowSolverExplicit : public SolverBase
{

public:

  typedef Array1dT<rArray1d> rSArray2d; 
  typedef Array1dT<iArray1d*> iArrayPtrs;
  typedef Array1dT<rArray1d*> rArrayPtrs;

  FractureFlowSolverExplicit(const std::string &name, ProblemManagerT* const problemManager);

  virtual ~FractureFlowSolverExplicit();
  
  void ReadXML( TICPP::HierarchicalDataNode* hdn ) ;
  void RegisterFields( PhysicalDomainT& domain );
  void Initialize( PhysicalDomainT& domain, SpatialPartition& partition);

  void InitializeCommunications( PartitionBase& partition );

  double TimeStep( const realT& time, const realT& dt,
                 const int cycleNumber,
                 PhysicalDomainT& domain,
                 const sArray1d& namesOfSolverRegions, 
                 SpatialPartition& partition,
                 FractunatorBase* const fractunator );
                 
  static const char* SolverName(){return "FractureFlowSolverExplicit";};

private:

  realT Solve(const realT& dt, PhysicalDomainT& domain, SpatialPartition& partition);
  
  void SetDirichletNodeBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set,realT time);

  void DirichletBoundaryCondition(PhysicalDomainT& domain, ObjectDataStructureBaseT& object, BoundaryConditionBase* bcBase, const lSet& set, realT time);

  void SetDirichletBCs(PhysicalDomainT& domain);

  void UpdateDirichletBCs(PhysicalDomainT& domain, const realT &time);

  void UpdateFluidRockProps(PhysicalDomainT& domain);

  void MakeParameters(PhysicalDomainT& domain);


  // Flags
  bool m_doApertureUpdate;
  bool m_doMFCoupling;
  
  std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d> syncedFields;

  std::string m_TrilinosIndexStr;
  static int m_instances;

  realT m_dt;

  struct 
  {
    //liquid phase properties

    realT m_mu; 
    realT m_rho_o; 
    realT m_press_o;
    realT m_compress;



    //rock proprties

    realT m_rockCompress;

    realT  m_maxPorosity;
    realT  m_minPorosity;

    realT m_maxAperture;
    realT m_minAperture;


  } m_fluidRockProperty;

  struct {

    int ncomp;
    int nvar;

  } m_eqt;
  

  Array1dT<rSArray2d> m_femCoefs;
  rArray1d m_volume;

  rArray1d m_density;

  //trilinos index

  iArrayPtrs m_elem_is_ghost;

  //global variables

  sArray1d m_fieldName;

  iArray1d  m_elemRegionIndex;
  Array1dT<ElementRegionT *>  m_elemRegion;

  sArray1d m_fracture_regionName;

  int m_fracture_regionNumber; 
  int m_matrix_regionNumber; 

};

int FractureFlowSolverExplicit::m_instances = 0;

#endif /* FRACTUREFLOWSOLVEREXPLICIT_H_ */
