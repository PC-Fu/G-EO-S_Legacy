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
 * @file ProblemManagerT.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */

#ifndef PROBLEMMANAGERT_H_
#define PROBLEMMANAGERT_H_

#include "DataStructures/Tables/Table.h"
#include "Common/GPException.h"
#include "PhysicalDomainT.h"

#include "IO/FileManagerDataT.h"
#include "IO/FileManagerT.h"
#include "IO/silo/SiloFile.h"
#include "IO/ticpp/TinyXMLParser.h"
#include "PhysicsSolvers/SolverApplicationSet.h"
#include "PhysicsSolvers/SolverBase.h"
#include "MPI_Communications/SpatialPartition.h"
#include "SurfaceGeneration/Fractunator.h"

#if 1
// This is for writing out node-based graph.  Pengcheng is working on it.
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#endif

#include "Constitutive/ConstitutivePropertiesTable.h"

#include "MeshUtilities/MeshGenerator.h"

#include "Epetra_ConfigDefs.h"
#if GPAC_MPI
  #include "Epetra_MpiComm.h"
#else
  #include "Epetra_SerialComm.h"
#endif

using namespace TICPP;
class InitialConditionBase;

/**
 * @author Randolph Settgast
 */
class ProblemManagerT
{
public:
  ProblemManagerT();
  ~ProblemManagerT();

  void ParseCommandLineInput( const int& argc, char* const argv[]) throw(GPException);

  void ProblemSetup();

  void ReadGeometryInput( HierarchicalDataNode& hdn );
  
  void RegisterFields();
  
  void DisplayFields();

  void DisplayUnits();

  void SetInitialConditions();

  void SetInitialConstitutive();

  void CompleteObjectInitialization( );

  void OutputMeshInformation();

  void WriteSilo( const bool isRestart );
  
  void WriteFlowTxt (const int cycleNum, const realT time);

  void ReadSilo( const bool isRestart );

  void CalculateNextSiloTime( const std::string& intervalTableName,
                              realT& siloIntervals,
                              realT& nextSiloTime,
                              globalIndex& nextSiloIndex );



  void DisplayUsage( );
  
  void DisplayVersion( );
  
  void DisplaySplash( );

  void DisplaySolvers( );
  
  void VerifySolvers( );

  void UpdateOwnership( );

  void WriteNodeGraph( const localIndex nodeID,
                       PhysicalDomainT& domain);


  /// Run the problem
  void Run( double& outputTime, realT t_start );

  std::map<std::string,SolverBase*> m_solvers;
  std::vector<InitialConditionBase*> m_initialConditions;
  ConstitutivePropertiesTable m_initialConstitutive;

  GPAC_IO::FileManagerT m_FileManager;

  MeshGenerator m_MeshGenerator;

  int m_numDomains;

  realT m_problemTime;
  realT m_dt;
  realT m_maxWallTime;

  int m_size;
  int m_rank;
  int m_cycleNumber;
  int m_cycleReportFreq;

  SiloFile m_siloFile;

  bool m_writeFEM;
  bool m_writeFEMFaces;
  bool m_writeFEMEdges;
  bool m_writeFlowText;
  sArray1d m_nameFieldsToPlot;

  int m_fractureFlag;
  sArray1d m_preFractureSets;

  bool m_writePlot;
  std::string m_visitFileGroupFile; // ".visit" file with list of time series files (usually "geos.visit")
  realT m_plotIntervals;
  globalIndex m_nextPlotIndex;
  std::string m_plotIntervalTableName;
  realT m_nextPlotTime;

  bool m_writeRestart;
  realT m_restartIntervals;
  globalIndex m_nextRestartIndex;
  std::string m_restartIntervalTableName;
  realT m_nextRestartTime;


  //Fractunator3 m_surfaceSeparation;
  FractunatorBase* m_surfaceSeparation;

  std::vector<SolverApplicationSet> m_solverApplicationSets;
  
  SpatialPartition m_partition;
  
  std::map<std::string, std::string> m_simulationParameterMap;

  PhysicalDomainT m_Domains;

  // trilinos communicators
  #if GPAC_MPI
  Epetra_MpiComm m_epetraComm;
  #else
  Epetra_SerialComm m_epetraComm;
  #endif


private:
  bool m_beginFromRestart;
  std::string m_beginFromRestartFileName;

  bool m_doWriteXML;
  std::string m_xmlOutputFileName;

  bool m_hackInitialStress;


  int m_trackEnergy;
  EnergyT m_energy;
  realT m_initialEnergy;

  bool m_displayFields;
  bool m_displaySplash;

  bool m_echoParameters;

  ProblemManagerT(const ProblemManagerT&);
  ProblemManagerT& operator=(const ProblemManagerT&);
  void RegisterFilesIncludedFromCommandLine(HierarchicalDataNode* hdn, sArray1d& includedFiles);
  void ParseIncludedFiles(HierarchicalDataNode* hdn);
  void BuildSimulationParameterMap(HierarchicalDataNode* hdn);
  bool ReplaceMathematicalExpressions(std::string& valueStr);
  void ParseMetadata(HierarchicalDataNode* hdn,bool isRoot);
  bool EvaluateXMLIfThenElseStatement(HierarchicalDataNode* &hdn, HierarchicalDataNode* parentNode);
  void ReadXML(HierarchicalDataNode& hdn);
  void WriteXML(std::string& filename, HierarchicalDataNode& hdn);


  realT SetInitialTimeStep( std::vector<SolverApplicationSet>::iterator solverSet );

  //CompleteObjectInitialization helpers
  void SetDomainBoundaryObjects();
  void InitializeObjectManagers();
  void ResetGlobalToLocal();
  void SetExternal();
  void SetupNeighborListsPackUnpack(std::map<PhysicalDomainT::ObjectDataStructureKeys, sArray1d>& syncedFields);
  void InitializeSurfaceSeparation();

  void UpdateEnergy( const bool set_initial = false );

  static realT SetLocalMaxTimeStep( const realT& hard_dt, const realT& solver_dt, const realT& current_max_dt );
  static realT SetGlobalMaxTimeStep( const realT& local_dt );

  void MarkFieldsToWrite();

};


#endif /* PROBLEMMANAGERT_H_ */
