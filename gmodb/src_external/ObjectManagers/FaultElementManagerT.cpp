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
/*
 * FaultElementManagerT.cpp
 *
 *  Created on: Dec 2, 2012
 *      Author: johnson346
 */

#include "FaultElementManagerT.h"
#include "Utilities/GeometryUtilities.h"
#include "ObjectManagers/PhysicalDomainT.h"
#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "BoundaryConditions/BoundaryConditions.h"

#include "SurfaceGeneration/StatisticalDistributionBaseT.h"

/**
 * @brief Constructor to set internal pointers to the external node and face managers
 * @author Scott Johnson
 *
 * @param nm External node manager pointer
 * @param fm External face manager pointer
 */
FaultElementManagerT::FaultElementManagerT(NodeManagerT* nm = 0, FaceManagerT* fm = 0) :
      ObjectDataStructureBaseT(ObjectDataStructureBaseT::FaultElementManager),
      m_nodeManager(NULL),
      m_faceManager(NULL),
      m_maxGlobalNumber(std::numeric_limits<localIndex>::max()),
      m_eventIndex(1),
      m_magnitudeWriteMin(std::numeric_limits<realT>::max()),
      m_neighborList(m_VariableOneToManyMaps["neighborList"]),
      m_neighborListInverse(m_UnorderedVariableOneToManyMaps["neighborListInverse"]),
      m_contactBufferOffset(0.0),
      m_contactBufferFactor(0.0),
      m_sorted(false),
      m_sorter(NULL),
      m_useInfiniteSupport(true),
      m_isSerial(true),
      m_rupture(),
      m_contact(),
      m_boundaryElementMaterialParameters(),
      m_properties(),
      m_transformStressFrame(),
      m_dislocation(),
      m_dstressDz(),
      m_porePressure(0),
      m_initialConstitutiveFields(),
      m_initialConstitutiveTables()
      //m_initialConstitutive(0)
{
  m_faceManager = fm;
  m_nodeManager = nm;

  m_contact.SetVariableParameters(true);
  m_contact.resize(fm->DataLengths(),1);

  //orientation attributes
  m_transformStressFrame.PlusIdentity(1.0);

  //geometric attributes
  m_nodeManager->AddKeyedDataField<FieldInfo::displacement>();
  m_nodeManager->AddKeyedDataField<FieldInfo::referencePosition>();

  m_faceManager->AddKeylessDataField<realT> ( "area", true, true);
  m_faceManager->AddKeylessDataField<R1Tensor>( "shearSlip", true, true);

  m_faceManager->AddKeylessDataField<realT> ( "boundingRadiusLastSort" , true, false );
  m_faceManager->AddKeylessDataField<realT> ( "boundingRadius" , true, true );

  m_faceManager->AddKeylessDataField<R1Tensor> ( "centerLastSort" , true, false );
  m_faceManager->AddKeyedDataField<FieldInfo::referencePosition>();

  m_faceManager->AddKeylessDataField<R1Tensor> ( "normal" , true, true );

  //AddKeyedDataField<FieldInfo::velocity>();
  AddKeyedDataField<FieldInfo::referencePosition>();

  m_sorter = SpatialSorting::SpatialSorterFactory::NewSpatialSorter("N2");
}

FaultElementManagerT::~FaultElementManagerT()
{
  if(m_sorter)
    delete m_sorter;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ACCOUNTING
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


globalIndex
FaultElementManagerT::insert(const localIndex i, const bool assignGlobals )
{
  const globalIndex gi = ObjectDataStructureBaseT::insert(i);
  m_contact.insert(i);
  m_boundaryElementMaterialParameters.insert(i);
  return gi;
}
void
FaultElementManagerT::erase( const localIndex i )
{
  ObjectDataStructureBaseT::erase(i);
  m_contact.erase(i);
  m_boundaryElementMaterialParameters.erase(i);
}
globalIndex
FaultElementManagerT::resize( const localIndex size, const bool assignGlobals )
{
  const globalIndex gi = ObjectDataStructureBaseT::resize( size, assignGlobals );
  m_contact.resize(size);
  m_boundaryElementMaterialParameters.resize(size);
  return gi;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I/O
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void
FaultElementManagerT::WriteNonManagedDataMembersToSilo( SiloFile& siloFile,
                                                        const std::string& ,
                                                        const std::string& ,
                                                        const int ,
                                                        const int ,
                                                        const realT ,
                                                        const bool ,
                                                        const std::string& ,
                                                        const std::string& ,
                                                        const lArray1d& )
{
  m_rupture.WriteSilo(siloFile);
  siloFile.DBWriteWrapper("m_eventIndex", m_eventIndex );

  m_boundaryElementMaterialParameters.WriteSilo(siloFile);
}

void FaultElementManagerT::ReadNonManagedDataMembersFromSilo( const SiloFile& siloFile,
                                                              const std::string& ,
                                                              const std::string& ,
                                                              const int ,
                                                              const int ,
                                                              const realT ,
                                                              const bool ,
                                                              const std::string& ,
                                                              const lArray1d& )
{
  m_rupture.ReadSilo(siloFile);
  siloFile.DBReadWrapper("m_eventIndex", m_eventIndex );

  m_boundaryElementMaterialParameters.ReadSilo(siloFile);
}

void FaultElementManagerT::ReadXML(TableManager& tableManager, TICPP::HierarchicalDataNode*  hdn)
{
  //Read material properties
  {
    TICPP::HierarchicalDataNode* matNode = hdn->GetChild("Material");
    if(!matNode)
      throw GPException("FaultElementManagerT::ReadXML: You must define a material for the fault element manager");
    m_properties.ReadXML(*matNode);
  }

  //Read rupture properties
  m_rupture.ReadXML(hdn);

  //Get SILO event write parameter
  m_magnitudeWriteMin = hdn->GetAttributeOrDefault<realT>("magnitudeWriteMin",
                                                          std::numeric_limits<realT>::max());

  //Get orientation
  {
    m_transformStressFrame = 0.0;

    //Get maximum horizontal stress direction vector
    R1Tensor x(0); x(0) = 1;
    R1Tensor vmax = hdn->GetAttributeTensorOrDefault("maximumHorizontalStressDirection", x);
    vmax.Normalize();
    m_transformStressFrame.AddToColumn(0, vmax);

    //Get minimum horizontal stress direction vector
    R1Tensor z(0); z(2) = 1;
    R1Tensor vup = hdn->GetAttributeTensorOrDefault("upVector", z);
    vup.Normalize();
    m_transformStressFrame.AddToColumn(2, vup);

    //Get the orthogonal vector (right-handed coordinate system)
    z.Cross(vup, vmax);
    m_transformStressFrame.AddToColumn(1, z);

    //Get the rake direction
    m_dislocation = hdn->GetAttributeTensorOrDefault("rakeVector", x);
    m_dislocation.Normalize();
  }

  //Get change in stress with respect to vertical distance
  {
    m_dstressDz(0) = hdn->GetAttributeOrDefault<realT>("dHdu", 0.0);
    m_dstressDz(1) = hdn->GetAttributeOrDefault<realT>("dhdu", 0.0);
    m_dstressDz(2) = hdn->GetAttributeOrDefault<realT>("dvdu", 0.0);
  }

  //Get pore pressure table
  std::string str = hdn->GetAttributeString("porePressureTableName");
  if(str.length() != 0)
  {
    std::map<std::string,Table<4, realT> >::iterator it = tableManager.Tables<4>().find(str);
    if(it == tableManager.Tables<4>().end())
      throw GPException("Cannot find requested table in the table manager while attempting to set pore pressure");
    m_porePressure = &it->second;
  }
//  str = hdn->GetAttributeString("porePressureFileName");
//  if(str.length() != 0)
//  {
//    tableManager.Tables<4>()
//  }
}

void FaultElementManagerT::WriteSiloFaultElements( SiloFile& siloFile,
                                                   const std::string& siloDirName,
                                                   const std::string& meshname,
                                                   const int centering,
                                                   const int cycleNum,
                                                   const realT problemTime,
                                                   const bool isRestart,
                                                   const std::string& regionName,
                                                   const lArray1d& mask )
{
  std::string subDirectory = siloDirName;
  std::string rootDirectory = "/" + siloDirName;
  siloFile.MakeSubDirectory( subDirectory, rootDirectory );
  DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());

  //-----------------------------------------------------
  //HANDLE THE VARIABLES ASSOCIATED WITH THE JOINT STATE
  sArray1d intVarNames;
  sArray1d realVarNames;
  sArray1d R1TensorVarNames;
  sArray1d R2TensorVarNames;
  sArray1d R2SymTensorVarNames;

  Array1dT<iArray1d*> intVars;
  Array1dT<rArray1d*> realVars;
  Array1dT<Array1dT<R1Tensor>*> R1Vars;
  Array1dT<Array1dT<R2Tensor>*> R2Vars;
  Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

  m_contact.GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  AllocateTemporaryFields( intVarNames, intVars );
  AllocateTemporaryFields( realVarNames, realVars );
  AllocateTemporaryFields( R1TensorVarNames, R1Vars );
  AllocateTemporaryFields( R2TensorVarNames, R2Vars );
  AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

  m_contact.Serialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

  ObjectDataStructureBaseT::WriteSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DeallocateTemporaryFields<int>(intVarNames);
  DeallocateTemporaryFields<realT>(realVarNames);
  DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
  DeallocateTemporaryFields<R2Tensor>(R2TensorVarNames);
  DeallocateTemporaryFields<R2SymTensor>(R2SymTensorVarNames);

  WriteNonManagedDataMembersToSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DBSetDir(siloFile.m_dbFilePtr, "..");
}

void FaultElementManagerT::WriteSiloEvent(const int cycleNum)
{
  SiloFile siloFile;
  siloFile.m_fileRoot = "event";
  siloFile.Initialize(PMPIO_WRITE);
#if GPAC_MPI
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    siloFile.WaitForBaton(rank, cycleNum, false );
  }
#endif

  WriteSiloEvent(siloFile, "fault_event", "event_mesh", cycleNum, LastEvent().m_time);

#if GPAC_MPI
  {
    siloFile.HandOffBaton();
    siloFile.ClearEmptiesFromMultiObjects(cycleNum);
  }
#endif

  siloFile.Finish();
}

void FaultElementManagerT::WriteSiloEvent( SiloFile& siloFile,
                                           const std::string& siloDirName,
                                           const std::string& meshname,
                                           const int cycleNum,
                                           const realT problemTime,
                                           const std::string& regionName,
                                           const lArray1d& mask)
{
  const EarthquakeSimulation::FaultRupture& lastEvent = m_rupture.m_rupture;
  std::map< std::string, rArray1d > realData;

  //-----------------------------
  //Write facial mesh
  //-----------------------------
  {
    const OrderedVariableOneToManyRelation& elementToNodeMap = m_faceManager->m_toNodesRelation;
    const OrderedVariableOneToManyRelation::size_type numElements = lastEvent.m_events.size();
    const localIndex numTotalPoints = 4 * numElements;

    //set the nodal coordinate data structure
    realT* coords[3];
    dvector xcoords(numTotalPoints);
    dvector ycoords(numTotalPoints);
    dvector zcoords(numTotalPoints);
    {
      const Array1dT<R1Tensor>& npos = m_nodeManager->GetFieldData<FieldInfo::referencePosition>();
      localIndex nodeIndex = 0;
      for(Array1dT<EarthquakeSimulation::FaultRuptureData>::const_iterator it = lastEvent.m_events.begin();
          it != lastEvent.m_events.end(); ++it)
      {
        const EarthquakeSimulation::FaultRuptureData& d = *it;
        if(elementToNodeMap[d.m_elementIndex].size() != 4)
          throw GPException("Cannot handle non-quads\n");

        realData["SeismicBeginTime"].push_back(d.m_time);
        realData["SeismicRiseTime"].push_back(d.m_riseTime);
        realData["Slip"].push_back(d.m_slip);
        realData["Area"].push_back(d.m_area);
        realData["SeismicMagnitude"].push_back(d.m_magnitude);
        realData["SeismicMoment"].push_back(d.m_moment);

        const lArray1d& nodeIndices = elementToNodeMap[d.m_elementIndex];
        for(lArray1d::const_iterator itt = nodeIndices.begin(); itt != nodeIndices.end(); ++itt, ++nodeIndex)
        {
          xcoords[nodeIndex]   = npos[*itt](0);
          ycoords[nodeIndex]   = npos[*itt](1);
          zcoords[nodeIndex]   = npos[*itt](2);
        }
      }
      coords[0] = xcoords.data();
      coords[1] = ycoords.data();
      coords[2] = zcoords.data();
    }

    const int numFaceTypes = 1;
    Array1dT<localIndex*> meshConnectivity(numFaceTypes);
    ivector shapecnt(numFaceTypes);
    ivector shapetype(numFaceTypes);
    ivector shapesize(numFaceTypes);
    Array1dT<lArray2d> faceToNodeMap(numFaceTypes);
    {
      for (int faceType = 0; faceType < numFaceTypes; ++faceType)
      {
        faceToNodeMap[faceType].resize2(numElements, 4);
        localIndex count = 0;

        for (localIndex k = 0; k < numElements; ++k)
        {
          for (localIndex a = 0; a < 4; ++a)
          {
            faceToNodeMap[faceType][k][a] = count;
            ++count;
          }
        }

        meshConnectivity[faceType] = faceToNodeMap[faceType].data();

        shapecnt[faceType]  = numElements;
        shapetype[faceType] = DB_ZONETYPE_QUAD;
        shapesize[faceType] = 4;
      }
    }

    siloFile.WriteMeshObject(meshname, numTotalPoints, coords, NULL, numFaceTypes,
                             shapecnt.data(), meshConnectivity.data(), NULL, NULL,
                             shapetype.data(), shapesize.data(), cycleNum, problemTime);
  }

  //-----------------------------
  //Write facial data
  //-----------------------------
  {
    std::string subDirectory = siloDirName;
    std::string rootDirectory = "/" + siloDirName;
    siloFile.MakeSubDirectory( subDirectory, rootDirectory );
    DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());

    siloFile.WriteFieldMapToSilo<realT>( meshname, realData, DB_ZONECENT, cycleNum,
                                         problemTime, false, rootDirectory,
                                         regionName, mask );

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }
}


void FaultElementManagerT::ReadSiloFaultElements( const SiloFile& siloFile,
                                                  const std::string& siloDirName,
                                                  const std::string& meshname,
                                                  const int centering,
                                                  const int cycleNum,
                                                  const realT problemTime,
                                                  const bool isRestart,
                                                  const std::string& ,
                                                  const lArray1d& mask)

{
  if( DBSetDir(siloFile.m_dbFilePtr, siloDirName.c_str()) != -1 )
  {
    sArray1d intVarNames;
    sArray1d realVarNames;
    sArray1d R1TensorVarNames;
    sArray1d R2TensorVarNames;
    sArray1d R2SymTensorVarNames;

    Array1dT<iArray1d*> intVars;
    Array1dT<rArray1d*> realVars;
    Array1dT<Array1dT<R1Tensor>*> R1Vars;
    Array1dT<Array1dT<R2Tensor>*> R2Vars;
    Array1dT<Array1dT<R2SymTensor>*> R2SymVars;

    m_contact.GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

    AllocateTemporaryFields( intVarNames, intVars );
    AllocateTemporaryFields( realVarNames, realVars );
    AllocateTemporaryFields( R1TensorVarNames, R1Vars );
    AllocateTemporaryFields( R2TensorVarNames, R2Vars );
    AllocateTemporaryFields( R2SymTensorVarNames, R2SymVars );

    ObjectDataStructureBaseT::ReadSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, "none", mask);

    ReadNonManagedDataMembersFromSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, "none", mask );

    m_contact.Deserialize(intVars, realVars, R1Vars, R2Vars, R2SymVars);

    DeallocateTemporaryFields<int>(intVarNames);
    DeallocateTemporaryFields<realT>(realVarNames);
    DeallocateTemporaryFields<R1Tensor>(R1TensorVarNames);
    DeallocateTemporaryFields<R2Tensor>(R2TensorVarNames);
    DeallocateTemporaryFields<R2SymTensor>(R2SymTensorVarNames);

    DBSetDir(siloFile.m_dbFilePtr, "..");
  }
}

//void FaultElementManagerT::ReadPorePressureNUFT( const std::string& filename)
//{
//  m_porePressure = 0;
//  if(filename.empty())
//    return;
//
//  std::ifstream in;
//  in.open(filename.data());
//  if (!in.is_open())
//    return;
//  std::string inputline;
//
//  //--------------------------------------------------
//  // ASSUMED FILE FORMAT (REGULAR GRID IN 4D)
//  //
//  // nt nx ny nz
//  // origint originx originy originz
//  // dt_(1) dt_(2) ... dt_(nt)
//  // dx_(1) dx_(2) ... dx_(nx)
//  // dy_(1) dy_(2) ... dy_(ny)
//  // dz_(1) dz_(2) ... dz_(nz)
//  // index_t index_x index_y index_z pore_pressure
//  // ...
//  //--------------------------------------------------
//
//  //--------------------------------------------------
//  //get the number of entries along each dimension
//  //as well as the total number of spatial entries per
//  //timestep (n_xyz)
//  std::vector<std::vector<realT> > gridValues;
//  {
//    //get the total number of spatial entries as well as
//    //the components (stored in gridValues)
//    gridValues.resize(4);//t,x,y,z
//
//    //get number of entries
//    for (localIndex i = 0; i < 4; i++)
//    {
//      getline(in, inputline);
//      std::istringstream linestream(inputline);
//      localIndex n;
//      linestream >> n;
//      gridValues[i].resize(n+1);
//    }
//  }
//
//  //--------------------------------------------------
//  //allocate the porepressure table
//  //localIndex n_xyz = 0;
//  {
//    //get offsets
//    realT offset[4] =
//    { 0.0, 0.0, 0.0, 0.0 };
//    for (localIndex i = 0; i < 4; i++)
//    {
//      getline(in, inputline);
//      std::istringstream linestream(inputline);
//      linestream >> offset[i];
//    }
//
//    //get element sizes -> convert to values
//    for (localIndex i = 0; i < 4; i++)
//    {
//      getline(in, inputline);
//      std::istringstream linestream(inputline);
//      gridValues[i][0] = offset[i];
//      for (localIndex j = 1; j < gridValues[i].size(); ++j)
//      {
//        linestream >> gridValues[i][j];
//        gridValues[i][j] += gridValues[i][j-1];
//      }
//    }
//    m_porePressure->SetGrid(gridValues);
//    //n_xyz = gridValues[1].size()*gridValues[2].size()*gridValues[3].size();
//  }
//
//  //--------------------------------------------------
//  //fill the pore pressure table from the file
//  {
//    realT pp = 0.0;
//    iArray1d index(4);
//    while (!in.eof())
//    {
//      getline(in, inputline);
//      std::istringstream linestream(inputline);
//      linestream >> index[0] >> index[1] >> index[2] >> index[3] >> pp;
//      m_porePressure->set_value(index, pp);
//    }
//  }
//
//  in.close();
//}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// INTEGRATION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


realT
FaultElementManagerT::CalculateTimestep(const bool steadyState)
{
  if(steadyState)
    m_rupture.m_timestep.m_currentTime = 0.0;
  if(m_porePressure)
    m_porePressure->SetZeroGradient(steadyState);

  //Apply the boundary conditions
  ApplyBoundaryConditions();

  //Get new timestep and update rupture fields
  m_rupture.Reset();
  const realT dt = m_contact.NextTransitionTime( m_rupture.m_timestep,
                                m_rupture.dxsdtEQ,
                                GetFieldData<FieldInfo::referencePosition>(),
                                m_porePressure);

  //set the global index of the patch controlling the timestep
  m_rupture.m_global = m_faceManager->m_localToGlobalMap[m_rupture.m_timestep.m_local];

  //(2) Synchronize the rupture state across processes (including rupture.ddot)
  if(!m_useInfiniteSupport)
    throw GPException("cannot yet handle finite support");
  const bool onNeighbor = m_rupture.Synchronize(dt, m_faceManager->m_localToGlobalMap,
                                                m_faceManager->m_globalToLocalMap);

  if(!onNeighbor)
    m_contact.PostTimestepSynchronization(m_rupture.m_timestep,
                                          m_rupture.m_current,
                                          m_rupture.m_ddot);
  return m_rupture.m_timestep.m_dt;
}

bool FaultElementManagerT::TimeStep(const realT dt, const bool steadyState)
{
  //Synchronize the cached time with that of the problem manager
  m_rupture.m_timestep.m_currentTime = steadyState ? 0 : m_rupture.m_timestep.m_currentTime + dt;
  if (m_porePressure)
    m_porePressure->SetZeroGradient(steadyState);

  //(1) Update element states and rupture.ddot
  m_contact.Advance(m_rupture.m_timestep, dt, m_rupture.dxsdtEQ,
                    m_faceManager->GetFieldData<FieldInfo::referencePosition>(), m_porePressure);

  //  //(2) Synchronize the rupture state across processes (including rupture.ddot)
  //  if(!m_useInfiniteSupport)
  //    throw GPException("cannot yet handle finite support");
  //  m_rupture.Synchronize(dt, m_faceManager->m_localToGlobalMap,
  //                        m_faceManager->m_globalToLocalMap);
  //<-- NOT SURE WHY THIS WAS HERE:  MOVED IT TO CALCULATETIMESTEP

  //(3) Transition ruptured element
  m_rupture.Transition(m_contact, m_boundaryElementMaterialParameters);

  //(4) End earthquake if necessary
  const bool endOfEQ = m_rupture.EarthquakeJustFinished();
  if (endOfEQ)
  {
    //finalize the event
    {
      lArray1d index;
      m_rupture.CurrentRuptureFaultElementIndex(index);
      rArray1d currentShear(index.size(), 0.0);
      {
        rArray1d::iterator itr = currentShear.begin();
        for (lArray1d::const_iterator it = index.begin(); it != index.end(); ++it, ++itr)
          *itr = m_contact.StateData(*it, 0)->xs;
      }
      m_rupture.FinalizeJustFinishedEarthquake(
          m_properties.init_shearModulus, currentShear,
          m_faceManager->GetFieldData<FieldInfo::referencePosition>(),
          m_faceManager->GetFieldData<realT>("area"));
    }

    //output event details
    {
      const EarthquakeSimulation::FaultRuptureData& lastEvent = LastEvent();
      if (steadyState)
        std::cout << "ss_";
      std::cout << "eq_event: " << lastEvent.m_time << " " <<
      lastEvent.m_magnitude << " " <<
      lastEvent.m_moment << " " <<
      lastEvent.m_area << " " <<
      lastEvent.m_slip << " " <<
      lastEvent.m_hypocenter(0) << " " <<
      lastEvent.m_hypocenter(1) << " " <<
      lastEvent.m_hypocenter(2) << std::endl;

      //write silo if requested
      if(lastEvent.m_magnitude > m_magnitudeWriteMin)
      {
        m_rupture.FinalizeJustFinishedEarthquake2();
        WriteSiloEvent(m_eventIndex++);
      }
    }
  }
  return endOfEQ;
}

realT
FaultElementManagerT::LastMagnitude()
{
  //note: must be called after TimeStep and ONLY if EndAnyFinishedEarthquake returns true
  return m_rupture.m_rupture.m_main.m_magnitude;//UpdateTotalMagnitude(m_properties.init_shearModulus);
}

const EarthquakeSimulation::FaultRuptureData&
FaultElementManagerT::LastEvent() const
{
  return m_rupture.m_rupture.m_main;
}

void
FaultElementManagerT::ApplyBoundaryConditions()
{
  // iterate over all boundary conditions.
  for( Array1dT<BoundaryConditionBase*>::const_iterator bcItr=m_bcData.begin() ; bcItr!=m_bcData.end() ; ++ bcItr )
  {
    // check if the requested field has a wall boundary condition applied to it.
    SimpleBoundaryCondition* bc = dynamic_cast<SimpleBoundaryCondition*> (*bcItr);
    if( bc )
    {
      if(streq(bc->GetFieldName(0), "shearRateDrive"))
      {
        //const R1Tensor& n = bc->GetDirection(0);
        for(localIndex a = 0; a < m_DataLengths; a++)
        {
          //Array1dT<R1Tensor>& v = GetFieldData<FieldInfo::velocity>();
          //v[a] = n;
          //v[a] *= bc->m_scale;
          m_contact.StateData(a, 0)->dxsdtDrive = bc->m_scale;
        }
      }
      else
        throw GPException("Only shearRateDrive boundaries supported for EQ simulation");
    }
    else
      throw GPException("FaultRuptureBEMSolver::ApplyBoundaryConditions: Invalid BC for EQ simulation");
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// INITIALIZATION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void FaultElementManagerT::Initialize()
{
  if(m_faceManager->DataLengths() == 0)
    return;

  m_maxGlobalNumber = m_faceManager->m_maxGlobalNumber;
  if(!m_useInfiniteSupport)
  {
    throw GPException("cannot handle FMM yet\n");
  }

  //Calculate the stiffness matrix once
  CalculateStiffnessMatrix();
}

void FaultElementManagerT::DeserializeObjectField(const std::string& name, const rArray1d& field)
{
  if(m_DataLengths == 0)
    return;

  //states
  std::map<std::string, size_t> intOffsets;
  std::map<std::string, size_t> realOffsets;
  std::map<std::string, size_t> R1TensorOffsets;
  std::map<std::string, size_t> R2TensorOffsets;
  std::map<std::string, size_t> R2SymTensorOffsets;
  {
    DieterichStateData& matState = *m_contact.StateData(0,0);
    matState.GetVariableOffsets(intOffsets,
                                realOffsets,
                                R1TensorOffsets,
                                R2TensorOffsets,
                                R2SymTensorOffsets);

    const std::map<std::string, size_t>::const_iterator ir0 = realOffsets.find(name);
    if(ir0 != realOffsets.end())
    {
      const size_t sz = ir0->second;
      if(field.size() != m_contact.NumStateIndex0() || field.size() != m_DataLengths)
        throw GPException("DeserializeObjectField: error material state length is not the same as field length or datalength is not");
      for(localIndex i = 0; i < m_contact.NumStateIndex0(); i++)
      {
        for(localIndex j = 0; j < m_contact.NumStateIndex1(); j++)
        {
          memcpy(((char*)m_contact.StateData(i, j) + sz), &field[i], sizeof(realT));
        }
      }
      return;
    }
  }

  //parameters
  intOffsets.clear();
  realOffsets.clear();
  R1TensorOffsets.clear();
  R2TensorOffsets.clear();
  R2SymTensorOffsets.clear();
  {
    m_contact.ParameterData(0)->GetVariableOffsets(intOffsets,
                                                   realOffsets,
                                                   R1TensorOffsets,
                                                   R2TensorOffsets,
                                                   R2SymTensorOffsets);

    std::map<std::string, size_t>::const_iterator ir0 = realOffsets.find(name);
    if(ir0 != realOffsets.end())
    {
      const size_t sz = ir0->second;
      if(field.size() != m_contact.NumParameterIndex0() || field.size() != m_DataLengths)
        throw GPException("DeserializeObjectField: error material parameter length is not the same as field length or datalength is not");
      for(localIndex i = 0; i < m_contact.NumStateIndex0(); i++)
        memcpy(((char*)m_contact.ParameterData(i) + sz), &field[i], sizeof(realT));
      return;
    }
  }

  throw GPException("DeserializeObjectField: could not find field name");
}

void FaultElementManagerT::CalculateFaceNormalsAndCenters(rArray1d& centroidsnormals)
{
  const gArray1d& localToGlobal = m_faceManager->m_localToGlobalMap;
  const Array1dT<R1Tensor>& disp = m_nodeManager->GetFieldData<FieldInfo::displacement>();
  const Array1dT<R1Tensor>& refpos = m_nodeManager->GetFieldData<FieldInfo::referencePosition>();

  Array1dT<R1Tensor>& faceCenters = m_faceManager->GetFieldData<FieldInfo::referencePosition>();
  Array1dT<R1Tensor>& faceNormals = m_faceManager->GetFieldData<R1Tensor>("normal");
  rArray1d& faceAreas = m_faceManager->GetFieldData<realT>("area");

  const OrderedVariableOneToManyRelation& faceToNodes = m_faceManager->m_toNodesRelation;

  const R1Tensor up = Up();
  R1Tensor nnrm;
  nnrm.Cross(up, m_dislocation);

#if GPAC_MPI
  if(!m_isSerial)
  {
    rArray1d tmpcentroidsnormals(centroidsnormals.size(), 0.0);
    {
      localIndex a = 0;
      for (gArray1d::const_iterator iter = localToGlobal.begin();
          iter != localToGlobal.end(); ++iter, ++a)
      {
        R1Tensor center, normal;
        faceAreas[a] = GeometryUtilities::Centroid_3DPolygon(faceToNodes[a], refpos, disp, center, normal);
        if(Dot(normal, nnrm) < 0)
          normal *= -1.0;

        faceNormals[a] = normal;
        faceCenters[a] = center;

        const globalIndex gi = *iter;

        tmpcentroidsnormals[6 * gi] = center(0);
        tmpcentroidsnormals[6 * gi + 1] = center(1);
        tmpcentroidsnormals[6 * gi + 2] = center(2);
        tmpcentroidsnormals[6 * gi + 3] = normal(0);
        tmpcentroidsnormals[6 * gi + 4] = normal(1);
        tmpcentroidsnormals[6 * gi + 5] = normal(2);
      }

      MPI_Allreduce(tmpcentroidsnormals.data(), centroidsnormals.data(), centroidsnormals.size(),
                 MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
  }
  else
  {
#endif
  localIndex a = 0;
  for (gArray1d::const_iterator iter = localToGlobal.begin();
      iter != localToGlobal.end(); ++iter, ++a)
  {
    R1Tensor center, normal;
    faceAreas[a] = GeometryUtilities::Centroid_3DPolygon(faceToNodes[a], refpos, disp, center, normal);
    if(Dot(normal, nnrm) < 0)
      normal *= -1.0;

    faceNormals[a] = normal;
    faceCenters[a] = center;

    const globalIndex gi = *iter;

    centroidsnormals[6 * gi] = center(0);
    centroidsnormals[6 * gi + 1] = center(1);
    centroidsnormals[6 * gi + 2] = center(2);
    centroidsnormals[6 * gi + 3] = normal(0);
    centroidsnormals[6 * gi + 4] = normal(1);
    centroidsnormals[6 * gi + 5] = normal(2);
  }
#if GPAC_MPI
  }
#endif
}


void FaultElementManagerT::ResetStatesAndParameters( const bool resetStress)
{
  m_rupture.ClearState();

  const gArray1d& localToGlobal = m_faceManager->m_localToGlobalMap;
  const realT dd = m_rupture.dxsdtEQ < 1e-6 ? m_rupture.dxsdtEQ : 1e-6;
  const Array1dT<R1Tensor>& faceCenters = m_faceManager->GetFieldData<FieldInfo::referencePosition>();
  const Array1dT<R1Tensor>& faceNormals = m_faceManager->GetFieldData<R1Tensor>("normal");

  for(localIndex i = 0; i < DataLengths(); i++)
  {
    DieterichParameterData& matParams = *m_contact.ParameterData(i);
    DieterichStateData& matState = *m_contact.StateData(i,0);

    //cache the driving stress before clearing the states
    {
      const realT dssdt = matState.dstressShearDtDrive;
      const realT dnsdt = matState.dstressdtDrive;

      matState *= 0.0;
      matParams *= 0.0;

      matState.dstressShearDtDrive = dssdt;
      matState.dstressdtDrive = dnsdt;
    }

    matParams.mu0 = StatisticalDistributionBaseT::UniformSample(0.5, 0.8);
    matParams.tFail = -1.0;
    matParams.stressShearFail = -1.0;
    matParams.alpha = 0.25;
    matParams.maxThetaPin = 1.0e10;
    matParams.A = StatisticalDistributionBaseT::UniformSample(0.009, 0.011);
    matParams.B = StatisticalDistributionBaseT::UniformSample(0.015, 0.02);
    matParams.Dc = StatisticalDistributionBaseT::UniformSample(1e-5, 2e-5);
    matParams.vstar = dd;
    matParams.dxsdtAB = dd;
    matParams.thetastar = 0;
    matParams.stressOvershootFactor = 0.1;
    matParams.AreductionFactor = 0.1;
    matParams.rake = faceNormals[i];
    matParams.rake *= -Dot(m_dislocation, faceNormals[i]);
    matParams.rake += m_dislocation;
    matParams.rake.Normalize();

    matState.apFail = 0;
    matState.neighborInRuptureState = 0;
    matState.pinned = 0;
    matState.theta = 1.5e8;
    matState.currentState = EarthquakeSimulation::LOCK;
    matState.Areduced = 0;

    globalIndex gi = localToGlobal[i];
    R1Tensor& kgd = m_boundaryElementMaterialParameters(i,gi);
    matParams.KNormalSelf = kgd(0);
    matParams.KShearSelf = kgd(1);

    SigmaTau(faceCenters[i], faceNormals[i], matParams.rake, matState.stressReference, matState.stressShearReference);

    realT xyzt[] = {faceCenters[i](0), faceCenters[i](1), faceCenters[i](2), 0};
    matState.dstress = 0.0;
    matState.stress = matState.stressReference + matState.dstress -
        (m_porePressure ? m_porePressure->Lookup(xyzt) : 0);
    matState.stressPin = 0;//std::numeric_limits<realT>::max();//TODO: add setting via stressPinFactor
  }

  //Apply the initial conditions (if you have set any)
  ApplyInitialConstitutive();

  //Apply the boundary conditions (override the IC's if necessary)
  ApplyBoundaryConditions();

  m_contact.SetStressRates(m_rupture.m_timestep.m_currentTime,
                           m_boundaryElementMaterialParameters,
                           m_porePressure,
                           faceCenters,
                           resetStress);
}

void FaultElementManagerT::ResetStates()
{
  m_rupture.ClearState();

  const Array1dT<R1Tensor>& faceCenters = m_faceManager->GetFieldData<FieldInfo::referencePosition>();

  for(localIndex i = 0; i < DataLengths(); i++)
  {
    DieterichStateData& matState = *m_contact.StateData(i,0);

    matState.apFail = 0;
    matState.neighborInRuptureState = 0;
    matState.pinned = 0;
    matState.theta = 1.5e8;
    matState.currentState = EarthquakeSimulation::LOCK;
    matState.Areduced = 0;

    realT xyzt[] = {faceCenters[i](0), faceCenters[i](1), faceCenters[i](2), 0};
    matState.dstress = 0.0;
    matState.stress = matState.stressReference + matState.dstress -
        (m_porePressure ? m_porePressure->Lookup(xyzt) : 0);
    matState.stressPin = 0;//std::numeric_limits<realT>::max();//TODO: add setting via stressPinFactor
  }

  //Apply the initial conditions (if you have set any)
  ApplyInitialConstitutive();

  //Apply the boundary conditions (override the IC's if necessary)
  ApplyBoundaryConditions();

  m_contact.SetStressRates(m_rupture.m_timestep.m_currentTime,
                           m_boundaryElementMaterialParameters,
                           m_porePressure, faceCenters, false);
}


void FaultElementManagerT::CalculateStiffnessMatrix()
{
  if(!m_useInfiniteSupport)
    throw GPException("Currently, we only support infinite support :-)");

  //RECALCULATE THE FACE GEOMETRY AND CACHE EVERY FACE'S CENTROID (GLOBALLY)
  rArray1d centroidsnormals(6*(m_faceManager->m_maxGlobalNumber + 1), 0.0);
  CalculateFaceNormalsAndCenters(centroidsnormals);

  //GET DERIVED MECHANICAL PROPERTIES
  const realT alphaTmp = m_properties.Nu/(1.0 - 2.0*m_properties.Nu);
  const realT alpha = (alphaTmp + 0.5) / (alphaTmp + 1.0);
  const realT lambda = m_properties.Lame;//init_bulkModulus - (2.0/3.0)*m_properties.init_shearModulus;

  //dislocationStrike = A^T * dislocationSpace
  R1Tensor dislocationStrike, dislocation;
  R2Tensor transformPlaneFrame, transformStrikeFrame;
  R1Tensor up = Up();

  //GO THROUGH EACH ELEMENT AND GET THE CONTRIBUTION OF ALL OTHERS
  const OrderedVariableOneToManyRelation& faceToNodes = m_faceManager->m_toNodesRelation;
  const gArray1d& localToGlobal = m_faceManager->m_localToGlobalMap;
  const Array1dT<R1Tensor>& disp = m_nodeManager->GetFieldData<FieldInfo::displacement>();
  const Array1dT<R1Tensor>& refpos = m_nodeManager->GetFieldData<FieldInfo::referencePosition>();

  m_boundaryElementMaterialParameters.resize(faceToNodes.size(), m_maxGlobalNumber+1);

  for(localIndex i = 0; i < faceToNodes.size(); i++)
  {
    if(faceToNodes[i].size() != 4)
      throw GPException("Cannot process a non-quad face for Okada elements");

    //Get normal and centroid of current face in a more usable form
    R1Tensor xi, ni;
    const globalIndex global0 = localToGlobal[i];
    for(localIndex a = 0; a < 3; a++)
    {
      xi(a) = centroidsnormals[6*global0+a];
      ni(a) = centroidsnormals[6*global0+a+3];
    }
    GeometryUtilities::TransformPlaneFrame(ni, up, transformPlaneFrame);
    GeometryUtilities::TransformStrikeFrame(ni, up, transformStrikeFrame);
    R1Tensor dip;
    dip(0) = transformPlaneFrame(0,2);
    dip(1) = transformPlaneFrame(1,2);
    dip(2) = transformPlaneFrame(2,2);
    R1Tensor strike;
    strike(0) = transformPlaneFrame(0,0);
    strike(1) = transformPlaneFrame(1,0);
    strike(2) = transformPlaneFrame(2,0);

    //project the dislocation to the plane of the fault element
    {
      dislocation = ni;
      dislocation *= -Dot(m_dislocation, ni);
      dislocation += m_dislocation;
      dislocation.Normalize();
      dislocationStrike.AijBi(transformStrikeFrame, dislocation);
    }

    const realT dipRadians = GeometryUtilities::Dip(ni,up);

    //Get the depth, length, and width of the fault patch
    realT depth, length, width;
    {
      R1TensorT<2> lwMax(-std::numeric_limits<realT>::max()), lwMin(std::numeric_limits<realT>::max());
      {
        R1Tensor xx;
        for(localIndex b = 0; b < 4; b++)
        {
          xx = refpos[faceToNodes[i][b]];
          xx += disp[faceToNodes[i][b]];
          xx -= xi;
          const realT ll = Dot(strike, xx);
          const realT ww = Dot(dip, xx);
          lwMax(0) = lwMax(0) < ll ? ll : lwMax(0);
          lwMax(1) = lwMax(1) < ww ? ww : lwMax(1);
          lwMin(0) = lwMin(0) > ll ? ll : lwMin(0);
          lwMin(1) = lwMin(1) > ww ? ww : lwMin(1);
        }
      }
      length = lwMax(0) - lwMin(0);
      width = lwMax(1) - lwMin(1);
      depth = -Dot(xi, up);
    }

    //----------------------------
    // Use Okada, 1992, element formulation
    //----------------------------
    //realT ktmp = 0.0, gtmp = 0.0;
    for (globalIndex global1=0; global1 <= m_maxGlobalNumber; global1++) /* loop over receiver patches */
    {
      R1Tensor xneighbor, nneighbor;
      for(localIndex a = 0; a < nsdof; a++)
      {
        xneighbor(a) = centroidsnormals[6*global1+a];
        nneighbor(a) = centroidsnormals[6*global1+a+3];
      }

      if(isEqual(Dot(nneighbor, nneighbor),0.0))
        continue;

      //get the relative position of the neighbor in the directions
      //orthogonal to vertical; in the vertical direction, retain the
      //original coordinate
      R1Tensor dxneighbor(xneighbor);
      {
        const realT neighbordepth = -Dot(xneighbor,up);
        dxneighbor -= xi;
        R1Tensor ttmp(up);
        ttmp *= Dot(dxneighbor,up) + neighbordepth;
        dxneighbor -= ttmp;
      }

      //Calculate the Okada contributions
      R2Tensor sigmaSpace;
      {
        R1Tensor slipStrike, dxnStrike;
        dxnStrike.AijBi(transformStrikeFrame, dxneighbor);

        R2Tensor strainStrike;
        CalculateOkadaQuadStrains(dxnStrike,
                                  alpha,
                                  depth,
                                  length,
                                  width,
                                  dipRadians,
                                  dislocationStrike,
                                  slipStrike,
                                  strainStrike);

        R2SymTensor sigmaStrike;
        DeformationToStress(strainStrike, lambda, m_properties.init_shearModulus, sigmaStrike);

        //sigmaSpaceFrame = Ts * sigmaStrikeFrame * Ts^T
        R2Tensor temp0, temp1;
        for(localIndex ii = 0; ii < nsdof; ii++)
          for(localIndex jj = 0; jj < nsdof; jj++)
            temp0(ii,jj) = sigmaStrike(ii,jj);
        temp1.AijBjk( transformStrikeFrame, temp0 );
        sigmaSpace.AijBkj( temp1, transformStrikeFrame );
        sigmaSpace *= -1.0;
      }

      //Project the stress back to the fault ... this is the stress per unit strain (i.e., stiffness)
      R1Tensor& kgd = m_boundaryElementMaterialParameters(i,global1);
      SigmaTau(sigmaSpace, dislocation, ni, kgd(0), kgd(1));
      //ktmp += kgd(0);
      //gtmp += kgd(1);
    } /* fi->type == OKADA */
    //std::cout << "i: " << i << " k: " << ktmp << " g: " << gtmp << std::endl;

    //    if (faceToNodes[i].size() == 3)
    //    {
    //      //----------------------------
    //      //Use Meade elements
    //      //----------------------------
    //      xx.reserve(3);
    //      for(localIndex a = 0; a < 3; a++)
    //      {
    //        xx[a] = refpos[faceToNodes[i][a]];
    //        xx[a] += disp[faceToNodes[i][a]];
    //      }
    //
    //      for (localIndex j = 0; j < domain.m_faultPatchFaces.DataLengths(); j++)
    //      {
    //
    //        const localIndex j = *iter;
    //        R2SymTensor du;
    //        const realT nu = (3*K[j] - 2*G[j])/(2*(3*K[j]+G[j]));
    //        this->CalculateMeadeTriStrains(facectr[j], nu, slip[j], xx, du);
    //
    //        /* note that because of symmetry of c_{ijkl},
    //          sigma_{ij} = c_{ijkl} e_{kl} and also
    //          sigma_{ij} = c_{ijkl} d_{kl}, where d_{kl} = u_{k,l} = du_k/dx_l */
    //        /* for now am flipping the sign of the results from calcTriStrains()
    //           as it seems to be using a diff. sign convention, but maybe I should
    //           just be doing this for the off-diagonal comps? */
    //        du *= -1.0;
    ////        deformationToStress(du, lambda, mu, K);
    ////        projectStress(K, fj->nu, fj->u, &fj->this->K, &fj->this->G);
    ////
    ////        /* up until here, extension has been reckoned positive, but in the simulator
    ////           code, compression is reckoned positive, so need to switch sign of Ksigma */
    ////        fj->this->K = -fj->this->K;
    ////
    ////        if (allK != 0)
    ////        {
    ////          for (iK=0; iK<6; iK++) fj->K[i][iK] = K[iK];
    ////        }
    //      } /* for (j=startPatch; j<=stopPatch; j++) */
    //    } /* if fi->type == TRIANGULAR */
  } /* for (i=startPatch; i <= stopPatch; i++) */

  Initialize(centroidsnormals);
}

void FaultElementManagerT::Initialize(const rArray1d& centroidsnormals)
{
  const realT dd = m_rupture.dxsdtEQ < 1e-6 ? m_rupture.dxsdtEQ : 1e-6;
  const gArray1d& localToGlobal = m_faceManager->m_localToGlobalMap;

  Array1dT<R1Tensor>& faceCenters = m_faceManager->GetFieldData<FieldInfo::referencePosition>();
  Array1dT<R1Tensor>& faceNormals = m_faceManager->GetFieldData<R1Tensor>("normal");
  Array1dT<R1Tensor>& centers = GetFieldData<FieldInfo::referencePosition>();

  StatisticalDistributionBaseT::InitializeRandom();

  for(localIndex i = 0; i < DataLengths(); i++)
  {
    globalIndex gi = localToGlobal[i];
    for(localIndex j = 0; j < nsdof; j++)
    {
      const realT x = centroidsnormals[6*gi + j];
      centers[i](j) = x;
      faceCenters[i](j) = x;
      faceNormals[i](j) = centroidsnormals[6*gi + j + 3];
    }

    DieterichParameterData& matParams = *m_contact.ParameterData(i);
    DieterichStateData& matState = *m_contact.StateData(i,0);

    matState *= 0.0;
    matParams *= 0.0;

    matParams.mu0 = StatisticalDistributionBaseT::UniformSample(0.5, 0.8);
    matParams.tFail = -1.0;
    matParams.stressShearFail = -1.0;
    matParams.alpha = 0.25;
    matParams.maxThetaPin = 1.0e10;
    matParams.A = StatisticalDistributionBaseT::UniformSample(0.009, 0.011);
    matParams.B = StatisticalDistributionBaseT::UniformSample(0.015, 0.02);
    matParams.Dc = StatisticalDistributionBaseT::UniformSample(1e-5, 2e-5);
    matParams.vstar = dd;
    matParams.dxsdtAB = dd;
    matParams.thetastar = 0;
    matParams.stressOvershootFactor = 0.1;
    matParams.AreductionFactor = 0.1;
    matParams.rake = faceNormals[i];
    matParams.rake *= -Dot(m_dislocation, faceNormals[i]);
    matParams.rake += m_dislocation;
    matParams.rake.Normalize();

    matState.apFail = 0;
    matState.neighborInRuptureState = 0;
    matState.pinned = 0;
    matState.theta = 1.5e8;
    matState.currentState = EarthquakeSimulation::LOCK;
    matState.Areduced = 0;

    R1Tensor& kgd = m_boundaryElementMaterialParameters(i,gi);
    matParams.KNormalSelf = kgd(0);
    matParams.KShearSelf = kgd(1);

    SigmaTau(centers[i], faceNormals[i], matParams.rake, matState.stressReference, matState.stressShearReference);
    realT xyzt[] = {centers[i](0), centers[i](1), centers[i](2), 0};
    matState.dstress = 0.0;
    matState.stress = matState.stressReference + matState.dstress -
        (m_porePressure ? m_porePressure->Lookup(xyzt) : 0);
    matState.stressPin = 0;//std::numeric_limits<realT>::max();//TODO: add setting via stressPinFactor
  }

  //Apply the boundary conditions
  ApplyBoundaryConditions();

  m_contact.SetStressRates(m_rupture.m_timestep.m_currentTime,
                           m_boundaryElementMaterialParameters,
                           m_porePressure, faceCenters, true);
}

void
FaultElementManagerT::ApplyInitialConstitutive()
{
  const TableManager& tableManager = TableManager::Instance();

  const Array1dT<R1Tensor>& pos = GetFieldData<FieldInfo::referencePosition>();
  rArray1d field(pos.size(), 0.0);

  //foreach name and manager pair
  sArray1d::const_iterator itn = m_initialConstitutiveTables.begin();
  for(sArray1d::const_iterator ifn = m_initialConstitutiveFields.begin(); ifn != m_initialConstitutiveFields.end(); ++ifn, ++itn)
  {
    //fill the temporary array
    {
      //get the table to query
      const std::string& tableName = *itn;
      const std::map<std::string,Table<3, realT> >::const_iterator it3 = tableManager.Tables<3>().find(tableName);
      if(it3 == tableManager.Tables<3>().end())
        throw GPException("ConstitutivePropertiesTable::Apply : Cannot find requested table in the table manager: " + tableName);
      const Table<3, realT>& t3d = it3->second;

      //store the results in an array
      localIndex i = 0;
      for(Array1dT<R1Tensor>::const_iterator it = pos.begin(); it != pos.end(); ++it, ++i)
        field[i] = t3d.Lookup(*it,  TableInterpolation::linear);
    }

    //copy the field entries into the appropriate object members
    DeserializeObjectField(*ifn, field);
  }
}

void
FaultElementManagerT::SetInitialConstitutive(ConstitutivePropertiesTable& initialConstitutive)
{
//  m_initialConstitutive = &initialConstitutive;
  m_initialConstitutiveFields = initialConstitutive.m_fieldNames;
  m_initialConstitutiveTables = initialConstitutive.m_tableNames;
}

//  //DO THE SPATIAL DISTRIBUTION OF PROPERTIES
//  //TODO: at some point it would be nice to parallelize the fractal distribution algorithms
//  //since it is already spatially decomposed, it should not be too hard ... just ship the
//  //nearest kernels across process bounds
//  {
//    bool root = true;
//    bool serial = true;
//#if GPAC_MPI
//    {
//      int mpirank;
//      MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
//      root = mpirank == 0;
//      int mpisize;
//      MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
//      serial = mpisize < 2;
//      //std::cout << mpirank << " is setting fractal friction " << (root ? "as root" : "NOT root") << " and " << (serial ? " serial" : " parallel") << std::endl;
//    }
//    if(root)
//#endif
//    {
//      //I AM ROOT ... I HAVE TO DO ALL OF THE WORK
//
//      //GET PARAMETERS
//      realT* hurst = this->frictionCoefficientDistribution.GetParameter(StatisticalDistributionBaseT::HURST_EXPONENT);
//      if(!hurst)
//        throw GPException("Cannot initialize fractal distribution without Hurst exponent");
//      realT* mean = this->frictionCoefficientDistribution.GetParameter(StatisticalDistributionBaseT::MEAN);
//      if(!mean)
//        throw GPException("Cannot initialize fractal distribution without mean");
//      realT* stdev = this->frictionCoefficientDistribution.GetParameter(StatisticalDistributionBaseT::STANDARD_DEVIATION);
//      if(!stdev)
//        throw GPException("Cannot initialize fractal distribution without standard deviation");
//
//      //GET EXTREMA
//      R1Tensor e1, e2;
//      FractalSurface spatial;
//      {
//        R1Tensor minAll(std::numeric_limits<realT>::max());
//        R1Tensor maxAll(-std::numeric_limits<realT>::max());
//        for (globalIndex global1=0; global1 <= maxGlobalNumber; global1++) /* loop over receiver patches */
//        {
//          R1Tensor xneighbor, nneighbor;
//          for(localIndex a = 0; a < 3; a++)
//          {
//            xneighbor(a) = centroidsnormals[6*global1+a];
//            nneighbor(a) = centroidsnormals[6*global1+a+3];
//          }
//
//          if(isEqual(Dot(nneighbor, nneighbor),0.0))
//            continue;
//
//          minAll.SetMin(xneighbor);
//          maxAll.SetMax(xneighbor);
//        }
//
//        //GET DIP AND STRIKE VECTORS FROM EXTREMA
//        e1 = maxAll;
//        e1 -= minAll;
//        {
//          const realT ue1 = Dot(this->up, e1);
//          R1Tensor tmp(up);
//          tmp *= ue1;
//          e1 -= tmp;
//          e1.Normalize();
//          tmp.Cross(e1, up);
//          e2.Cross(tmp, e1);
//          e2.Normalize();
//        }
//
//        //INITIALIZE THE FRACTAL DISTRIBUTION
//        R1TensorT<2> min2;
//        min2(0) = Dot(e1,minAll);
//        min2(1) = Dot(e2,minAll);
//        R1TensorT<2> max2;
//        max2(0) = Dot(e1,maxAll);
//        max2(1) = Dot(e2,maxAll);
//
//        spatial.InitializeFractal(*mean, *stdev, min2, max2, *hurst, 5, 5, 5, 1.2);
//      }
//
//#if GPAC_MPI
//      if(serial)
//#endif
//      {
//        localIndex i = 0;
//        rArray1d& mu = faultFaces.GetFieldData<realT>("frictionCoefficient");
//        for(gArray1d::const_iterator iter = faultFaces.m_localToGlobalMap.begin();
//            iter != faultFaces.m_localToGlobalMap.end(); ++iter, ++i)
//        {
//          R1Tensor xneighbor, nneighbor;
//          for(localIndex a = 0; a < 3; a++)
//          {
//            xneighbor(a) = centroidsnormals[6*(*iter)+a];
//            nneighbor(a) = centroidsnormals[6*(*iter)+a+3];
//          }
//
//          if(isEqual(Dot(nneighbor, nneighbor),0.0))
//            continue;
//
//          R1TensorT<2> pos;
//          pos(0) = Dot(e1,xneighbor);
//          pos(1) = Dot(e2,xneighbor);
//
//          mu[i] = spatial.Aperture(pos);
//        }
//      }
//#if GPAC_MPI
//      else
//      {
//        rArray1d muG(this->maxGlobalNumber + 1, 0.0);
//        for (globalIndex global1=0; global1 <= maxGlobalNumber; global1++) /* loop over receiver patches */
//        {
//          R1Tensor xneighbor, nneighbor;
//          for(localIndex a = 0; a < 3; a++)
//          {
//            xneighbor(a) = centroidsnormals[6*global1+a];
//            nneighbor(a) = centroidsnormals[6*global1+a+3];
//          }
//
//          if(isEqual(Dot(nneighbor, nneighbor),0.0))
//            continue;
//
//          R1TensorT<2> pos;
//          pos(0) = Dot(e1,xneighbor);
//          pos(1) = Dot(e2,xneighbor);
//
//          muG[global1] = spatial.Aperture(pos);
//        }
//
//        MPI_Bcast(muG.data(), muG.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//        std::cout << "Fault rupture solver: root sent " << muG.size() << " bulk moduli sent\n";
//
//        rArray1d& mu = faultFaces.GetFieldData<realT>("frictionCoefficient");
//        localIndex i = 0;
//        for(gArray1d::const_iterator iter = faultFaces.m_localToGlobalMap.begin();
//            iter != faultFaces.m_localToGlobalMap.end(); ++iter, ++i)
//        {
//          mu[i] = muG[*iter];
//        }
//      }
//    }
//    else
//    {
//      //NOT ROOT ... WAIT FOR BROADCAST
//
//      rArray1d muG(this->maxGlobalNumber + 1, 0.0);
//      MPI_Bcast(muG.data(), muG.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//      //std::cout << "not root: " << muG.size() << " G's received from root\n";
//
//      rArray1d& mu = faultFaces.GetFieldData<realT>("frictionCoefficient");
//      localIndex i = 0;
//      for(gArray1d::const_iterator iter = faultFaces.m_localToGlobalMap.begin();
//          iter != faultFaces.m_localToGlobalMap.end(); ++iter, ++i)
//      {
//        mu[i] = muG[*iter];
//      }
//    }
//#endif
//  }






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// STRESS/STRAIN AND COORDINATE TRANSFORMS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void FaultElementManagerT::SigmaTau(const R1Tensor& x, const R1Tensor& n, const R1Tensor& d,
                                     realT& stressNormal, realT& stressShear) const
{
  R1Tensor vstrike, nstrike;
  R1Tensor up = Up();
  {
    GeometryUtilities::Strike(n, Up(), vstrike);
    vstrike.Normalize();
    nstrike.Cross(vstrike, up);
  }
  const R1Tensor sstrike = MaximumHorizontalStressDirection();
  const realT cs = Dot(vstrike, sstrike);
  const realT ss = Dot(nstrike, sstrike);
  const realT depth = Dot(x, up);

  R2Tensor A(0.0);
  A(0,0) = cs;
  A(0,1) = -ss;
  A(1,0) = ss;
  A(1,1) = cs;
  A(2,2) = 1.0;

  R2Tensor sigmaPrincipalStress(0.0);
  sigmaPrincipalStress(0,0) = m_dstressDz(0) * depth;
  sigmaPrincipalStress(1,1) = m_dstressDz(1) * depth;
  sigmaPrincipalStress(2,2) = m_dstressDz(2) * depth;

  R2Tensor sigmaSpace;
  {
    R2Tensor tmp;
    tmp.AijBkj(sigmaPrincipalStress, A);
    sigmaSpace.AijBjk(A, tmp);
  }

  SigmaTau(sigmaSpace, d, n, stressNormal, stressShear);
}


void FaultElementManagerT::SigmaTau(const R2Tensor& sigmaGlobal,
                                     const R1Tensor& dislocation,
                                     const R1Tensor& normal,
                                     realT& sigma,
                                     realT& tau)
{
  sigma = 0.0;
  tau = 0.0;
  for (localIndex i = 0; i < nsdof; i++)
  {
    for (localIndex j = 0; j < nsdof; j++)
    {
      const realT sigmaij_nj = sigmaGlobal(i,j) * normal(j);
      sigma += normal(i) * sigmaij_nj;
      tau += dislocation(i) * sigmaij_nj;
      //tau += sigmaij_nj * sigmaij_nj;
    }
  }
  //tau -= (sigma*sigma);
  //tau = sqrt(tau);
}

void FaultElementManagerT::DeformationToStress(const R2Tensor& strain, /* deformation tensor, see above */
                                               const realT lambda, const realT G, /* Lame's constants */
                                               R2SymTensor& sigma  /* stress tensor as a vector, see above */
                                               )
{
  const realT ltre = lambda * strain.Trace();

  for(localIndex i = 0; i < nsdof; i++)
  {
    for(localIndex j = 0; j < nsdof; j++)
      sigma(i,j) = G*(strain(i,j) + strain(j,i));
    sigma(i,i) += ltre;
  }
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//     MEADE'S TRI DISLOCATION ELEMENTS
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


//Translated to C++ by Scott Johnson, 2012
//From ...
//Translated to C by Keith Richards-Dinger, 2010
//From ...
//Brendan Meade's Matlab code, whose original
//   comments are reproduced below
//
//   note that I store the elements of the strain tensor
//   according to my convention in KRDOkada.h:
//   (0 1 2
//      3 4
//        5)
//
//function [S] = CalcTriStrains(sx, sy, sz, x, y, z, nu, ss, ts, ds)
//% CalcTriStrains.m
//%
//% Calculates strains due to slip on a triangular dislocation in an
//% elastic half space utilizing the symbolically differentiated
//% displacement gradient tensor derived from the expressions for
//% the displacements due to an angular dislocation in an elastic half
//% space (Comninou and Dunders, 1975).
//%
//% Arguments
//%  sx : x-coordinates of observation points
//%  sy : y-coordinates of observation points
//%  sz : z-coordinates of observation points
//%  x  : x-coordinates of triangle vertices.
//%  y  : y-coordinates of triangle vertices.
//%  z  : z-coordinates of triangle vertices.
//%  nu : Poisson's ratio
//%  ss : strike slip displacement
//%  ts : tensile slip displacement
//%  ds : dip slip displacement
//%
//% Returns
//%  S  : structure containing the strains (S.xx, S.yy, S.zz, S.xy, S.xz, S.yz)
//%
//% Implements algorithms described in the journal article:
//% Meade, B. J. Algorithms for calculating displacements,
//% strains, and stresses for triangular dislocation elements
//% in a uniform elastic half space
//% Computers and Geosciences, submitted, 2006.
//%
//% Use at your own risk and please let me know of any bugs/errors.
//%
//% Copyright (c) 2006 Brendan Meade
//%
//% Permission is hereby granted, free of charge, to any person obtaining a
//% copy of this software and associated documentation files (the
//% "Software"), to deal in the Software without restriction, including
//% without limitation the rights to use, copy, modify, merge, publish,
//% distribute, sublicense, and/or sell copies of the Software, and to permit
//% persons to whom the Software is furnished to do so, subject to the
//% following conditions:
//%
//% The above copyright notice and this permission notice shall be included
//% in all copies or substantial portions of the Software.
//%
//% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
//% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
//% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
//% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//% USE OR OTHER DEALINGS IN THE SOFTWARE.
//
void FaultElementManagerT::CalculateMeadeTriStrains(const R1Tensor& xobs,
                                                     const realT nu,
                                                     const R1Tensor& slipPlaneFrame,
                                                     Array1dT<R1Tensor>& xx,
                                                     R2SymTensor& S)
{
  const realT tol = 1e-6;

  if(xx.size()!=3)
    throw GPException("Attempting to evaluate Meade triangular elements with n!=3");
  const int n = 3;

  S = 0.0;

  //TODO: when you port the call make sure the arguments 1 & 2 are transposed!!!
  //i.e., in the original call (..ss,ts,ds..)
  //slipPlaneFrame[0] = ss;  slipPlaneFrame[1] = ds; slipPlaneFrame[2] = ts;

  //get the slip vector in global frame
  R1Tensor slip;
  {
    R2Tensor T;
    PlaneFrameTransform(xx, T);
    slip.AijBj(T, slipPlaneFrame);
  }

  //go through components
  for (int iTri=0; iTri<n; iTri++)
  {
    int iTriTgt = (iTri+1)%n;

    R1Tensor dx(xx[iTriTgt]);
    dx -= xx[iTri];

    /* Calculate strike and dip of current leg */
    const realT strike = (180/M_PI)*(atan2(dx(1), dx(0)));
    //const realT segMapLength = sqrt(dx(0)*dx(0) + dx(1)*dx(1));

    /* Get the strike rotation transform */
    R2Tensor R(0.0);
    {
      const realT alpha = -(M_PI/180)*strike;
      R(0,0) = cos(alpha);
      R(0,1) = -sin(alpha);
      R(1,0) = -R(0,1);
      R(1,1) = R(0,0);
      R(2,2) = 1.0;
    }

    /* Get the segment vector transformed to strike frame */
    R1Tensor dxs;
    dxs.AijBj(R, dx);

    const realT dip = (180/M_PI)*(atan2(dx(2), dxs(0)));
    const realT beta = dip >= 0 ? (((M_PI/180)*(90-dip)) > (0.5*M_PI) ? 0.5*M_PI-(M_PI/180)*(90-dip) : (M_PI/180)*(90-dip)) :
        ((-(M_PI/180)*(90+dip)) < (-0.5*M_PI) ? 0.5*M_PI-fabs(-(M_PI/180)*(90+dip)) : -(M_PI/180)*(90+dip));

    /* Get the slip vector in strike frame */
    R1Tensor slipg;
    slipg.AijBi(R, slip);

    if (fabs(beta) > tol && fabs(beta - M_PI) > tol)
    {

      /* First angular dislocation strain */
      R2SymTensor e1;
      {
        R1Tensor s1;
        {
          R1Tensor tmp(xobs);
          tmp -= xx[iTri];
          s1.AijBj(R, tmp);
        }
        TriStrainMeade(s1, xx[iTri](2), beta, nu, slipg, e1);
      }

      /* Second angular dislocation strain */
      R2SymTensor e2;
      {
        R1Tensor s2;
        {
          R1Tensor tmp(xobs);
          tmp -= xx[iTriTgt];
          s2.AijBj(R, tmp);
        }
        TriStrainMeade(s2, xx[iTriTgt](2), beta, nu, slipg, e2);
      }

      /* Rotate tensors to correct for strike */
      R2SymTensor en(0.0);
      {
        R2SymTensor e12 = e1;
        e12 -= e2;
        const realT strikeRadians = (M_PI/180)*strike;
        {
          const realT cg = cos(strikeRadians);
          const realT sg = sin(strikeRadians);
          en(0,0) = (cg*e12(0,0)-sg*e12(0,1))*cg-(cg*e12(0,1)-sg*e12(1,1))*sg;
          en(0,1) = (cg*e12(0,0)-sg*e12(0,1))*sg+(cg*e12(0,1)-sg*e12(1,1))*cg;
          en(0,2) = cg*e12(0,2)-sg*e12(1,2);
          en(1,1) = (sg*e12(0,0)+cg*e12(0,1))*sg+(sg*e12(0,1)+cg*e12(1,1))*cg;
          en(1,2) = sg*e12(0,2)+cg*e12(1,2);
          en(2,2) = e12(2,2);
        }
      }

      /* Add the strains from current leg */
      S += en;
    }
  } /* for (iTri=0; iTri<3; iTri++) */
}

void FaultElementManagerT::PlaneFrameTransform(Array1dT<R1Tensor>& xx,
                                                R2Tensor& T)
{
  /* Calculate the slip vector in XYZ coordinates */
  R1Tensor dx1(xx[1]);
  R1Tensor x0(xx[0]);
  dx1 -= x0;
  R1Tensor dx2(xx[2]);
  dx2 -= x0;
  R1Tensor normVec;
  normVec.Cross(dx1, dx2);
  normVec.Normalize();

  /* Enforce clockwise circulation */
  if (normVec(2) < 0)
  {
    R1Tensor tmp(xx[2]);
    xx[2] = xx[1];
    xx[1] = tmp;
    normVec *= -1;
  }

  /* Get the strike and slip vectors */
  R1Tensor strikeVec(-sin(atan2(normVec[1], normVec[0])), cos(atan2(normVec[1], normVec[0])), 0.0);

  R1Tensor dipVec;
  dipVec.Cross(normVec, strikeVec);

  T.AddToColumn(0, strikeVec);
  T.AddToColumn(1, dipVec);
  T.AddToColumn(2, normVec);
}

//void FaultElementManagerT::TriSlipVector(const R1Tensor& slipPlaneFrame,
//                                          Array1dT<R1Tensor>& xx,
//                                          R1Tensor& slip)
//{
//  //TODO: move this to a cached value on the surface ... no reason for this to be calculated EVERY time
//  R2Tensor T;
//  {
//    /* Calculate the slip vector in XYZ coordinates */
//    R1Tensor dx1(xx[1]);
//    R1Tensor x0(xx[0]);
//    dx1 -= x0;
//    R1Tensor dx2(xx[2]);
//    dx2 -= x0;
//    R1Tensor normVec;
//    normVec.Cross(dx1, dx2);
//    normVec.Normalize();
//
//    /* Enforce clockwise circulation */
//    if (normVec(2) < 0)
//    {
//      R1Tensor tmp(xx[2]);
//      xx[2] = xx[1];
//      xx[1] = tmp;
//      normVec *= -1;
//    }
//
//    /* Get the strike and slip vectors */
//    R1Tensor strikeVec (-sin(atan2(normVec[1],normVec[0])),
//                        cos(atan2(normVec[1],normVec[0])),
//                        0.0);
//
//    R1Tensor dipVec;
//    dipVec.Cross(normVec, strikeVec);
//
//    T.AddToColumn(0, strikeVec);
//    T.AddToColumn(1, dipVec);
//    T.AddToColumn(2, normVec);
//  }
//
//  slip.AijBj(T, slipPlaneFrame);
//}


/**
   These are the strains in a uniform elastic half space due to slip
   on an angular dislocation.  They were calculated by symbolically
   differentiating the expressions for the displacements (Comninou and
   Dunders, 1975, with typos noted by Thomas 1993) then combining the
   elements of the displacement gradient tensor to form the strain tensor.

   The six very long lines were converted by getting rid of the '.'s in '.*', './', etc.
   and then running them through my caret2pow utility.  Actually, there was a
   bit more: see ~/Papers/Others/Meade2007/README
*/
void FaultElementManagerT::TriStrainMeade(const R1Tensor& y,
                                           const realT a,
                                           const realT b,
                                           const realT nu,
                                           const R1Tensor& B,
                                           R2SymTensor& e)
{
  const realT ynorm = y.L2_Norm();
  R2Tensor yy;
  yy.dyadic_ab(y,y);
  const realT half = 0.5;
  const realT diff12nu = 1.0 - 2.0*nu;
  const realT diff22nu = 2.0 - 2.0*nu;
  const realT sb = sin( b );
  const realT cb = cos( b );
  const realT cotb = 1.0 / tan( b );
  const realT threehalves = 1.5;
  R1Tensor yplus2a = y;
  yplus2a += 2.0 * a;
  R1Tensor yplus2a2;
  yplus2a2.AiBi(yplus2a,yplus2a);

  e(0,0) = B(0)*(1/(realT)8*(diff22nu*(2*y(1)/(realT)yy(0,0)/(realT)(1+yy(1,1)/(realT)yy(0,0))-y(1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(y(1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)*y(0)-y(1)*ynorm*sb/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb-y(2)*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*y(0)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb+(yplus2a(2))*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))-y(1)*(1/(realT)ynorm/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))-y(0)*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(0)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(2),(realT)2)*y(0)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*y(0)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*y(0))-y(1)*cb*((1/(realT)ynorm*sb*y(0)-1)/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-(ynorm*sb-y(0))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)+(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb*y(0)-1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*((-2+2*nu)*diff12nu*(y(1)/(realT)yy(0,0)/(realT)(1+yy(1,1)/(realT)yy(0,0))-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*y(0)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb+(yplus2a(2))*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*Power((realT)cotb,(realT)2)-diff12nu*y(1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((1-2*nu-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb-y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)+diff12nu*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)*cotb-1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))-diff12nu*y(1)*cb*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-diff12nu*y(1)*cb*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-3*a*y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*y(0)-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((-1+2*nu)*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)-y(1)*(y(2)+a)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu)*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-2*a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2))-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)-a*(yplus2a(2))*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)-a*(yplus2a(2))*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(0)*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)*cotb+diff22nu*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb*y(0)-1)*cb)+2*a*(yplus2a(2))*cb*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(0)))/(realT)M_PI/(realT)(1-nu))+B(1)*(1/(realT)8*((-1+2*nu)*(1/(realT)ynorm*y(0)/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*((1/(realT)ynorm*y(0)-sb)/(realT)(ynorm-y(0)*sb-y(2)*cb)+(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))+2*y(0)*(1/(realT)ynorm/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))+yy(0,0)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(0)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(2),(realT)2)*y(0)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*y(0)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*y(0))+cb*(ynorm*sb-y(0))/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)+(y(0)*cb-y(2)*sb)*(1/(realT)ynorm*sb*y(0)-1)/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb-y(2)*sb)*(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-(y(0)*cb-y(2)*sb)*(ynorm*sb-y(0))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)+cb*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(y(0)*cb+(yplus2a(2))*sb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb*y(0)-1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*((diff22nu*Power((realT)cotb,(realT)2)+nu)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-(diff22nu*Power((realT)cotb,(realT)2)+1)*cb*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-diff12nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu)*y(0)*cotb+nu*(yplus2a(2))-a+a*y(0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)+diff12nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((-1+2*nu)*cotb+a*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-a*yy(0,0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+2*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-Power(y(0),(realT)3)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-Power(y(0),(realT)3)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))+diff12nu*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((y(0)*cb+(yplus2a(2))*sb)*cb-a*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-diff12nu*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(Power((realT)cb,(realT)2)-a*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb*y(0)-1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb+a*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*y(0))-a*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+3*a*yy(0,0)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*y(0)*cotb+a)-yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*y(0)*cotb+a)*y(0)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*diff12nu*cotb-2*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+Power(y(0),(realT)3)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+Power(y(0),(realT)3)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+Power(y(0),(realT)3)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a-2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)+3*a*Power(y(0),(realT)3)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2))-(y(2)+a)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(-cb*sb+a*y(0)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+(y(2)+a)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb-3*a*yy(0,0)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)/(realT)cb+(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb*y(0)-1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb))*y(0)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*y(0))))/(realT)M_PI/(realT)(1-nu))+B(2)*(1/(realT)8*y(1)*sb*((1/(realT)ynorm*sb*y(0)-1)/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-(ynorm*sb-y(0))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)+(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb*y(0)-1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*(-y(1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)+y(1)*cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+y(1)*cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0))+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))*y(0)-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(0)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0))-y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)-y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-2*a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(0)))/(realT)M_PI/(realT)(1-nu));
  e(1,1) = B(0)*(1/(realT)8*(diff12nu*(1/(realT)ynorm*y(1)/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*(1/(realT)ynorm*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))-2*y(1)*(1/(realT)ynorm/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*(1/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))-yy(1,1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(2),(realT)2)*y(1)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*y(1)-cb*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1))))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*((diff22nu*Power((realT)cotb,(realT)2)-nu)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-(diff22nu*Power((realT)cotb,(realT)2)+1-2*nu)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+diff12nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(y(0)*cotb*(1-2*nu-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+nu*(yplus2a(2))-a+yy(1,1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-diff12nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(a*y(0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)+2*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-Power(y(1),(realT)3)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-Power(y(1),(realT)3)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))+diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)+3*a*y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*y(0)-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(-2*nu+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*y(0)*cotb-a)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*y(0)*cotb-a)*y(1)+2*y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-Power(y(1),(realT)3)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-Power(y(1),(realT)3)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-Power(y(1),(realT)3)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a+2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-3*a*Power(y(1),(realT)3)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2))-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(Power((realT)cb,(realT)2)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb+a*cb)+a*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb+a*cb)*y(1)-3*a*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*y(1)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*y(1)+1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*y(1)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(1)*Power((realT)cb,(realT)2)+a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*y(1)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(1))))/(realT)M_PI/(realT)(1-nu))+B(1)*(1/(realT)8*(diff22nu*(-2/y(0)/(realT)(1+yy(1,1)/(realT)yy(0,0))+1/(realT)(y(0)*cb-y(2)*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)+yy(1,1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)-2*yy(1,1)*ynorm*sb/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))+1/(realT)(y(0)*cb+(yplus2a(2))*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)-2*yy(1,1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))+y(0)*(1/(realT)ynorm/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))+y(0)*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(2),(realT)2)*y(1)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*y(1))-(y(0)*cb-y(2)*sb)/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-y(1)*(-(y(0)*cb-y(2)*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-(y(0)*cb-y(2)*sb)/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)-(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*diff12nu*(-1/y(0)/(realT)(1+yy(1,1)/(realT)yy(0,0))+1/(realT)(y(0)*cb+(yplus2a(2))*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)-2*yy(1,1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*Power((realT)cotb,(realT)2)+diff12nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((-1+2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))-diff12nu*yy(1,1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+diff12nu*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)*cotb-y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0))-diff12nu*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+diff12nu*yy(1,1)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+diff12nu*yy(1,1)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb-a*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+3*a*yy(1,1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)+(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))))-yy(1,1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))))-yy(1,1)*(y(2)+a)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))))+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu*y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))*y(1)-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)))+(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((-2+2*nu)*cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb)-yy(1,1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((-2+2*nu)*cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb)-yy(1,1)*(y(2)+a)*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((-2+2*nu)*cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb)+y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*y(1)-2*a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)cb*y(1)))/(realT)M_PI/(realT)(1-nu))+B(2)*(1/(realT)8*(diff12nu*sb*(1/(realT)ynorm*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-2*y(1)*sb*(1/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-yy(1,1)*sb*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*(-sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+y(1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)+y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-(y(0)*cb+(yplus2a(2))*sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1))-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))*y(0)+y(0)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1))+(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(sb*(cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(sb*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(1)-2*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*a*(yplus2a(2))*y(1)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*y(1)+1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*y(1)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(1)*cb*sb+a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*y(1)-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(1))))/(realT)M_PI/(realT)(1-nu));
  e(2,2) = B(0)*(1/(realT)8*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)*y(2)+1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)-cb*((1/(realT)ynorm*cb*y(2)-1)/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-(ynorm*cb-y(2))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)-(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*(diff12nu*(-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(1/(realT)2*y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*(2*y(2)+4*a)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*cotb-y(1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-1/(realT)2*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+y(1)*cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+1/(realT)2*y(1)*cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a))+y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)2*y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(2*y(2)+4*a)+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-2*nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a))+y(1)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1-2*nu-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)2*y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1-2*nu-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(2*y(2)+4*a)-y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1-2*nu-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)-a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a)))/(realT)M_PI/(realT)(1-nu))+B(1)*(1/(realT)8*((-1+2*nu)*sb*((1/(realT)ynorm*y(2)-cb)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-y(0)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)*y(2)+1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a))-sb*(ynorm*cb-y(2))/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)+(y(0)*cb-y(2)*sb)*(1/(realT)ynorm*cb*y(2)-1)/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb-y(2)*sb)*(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-(y(0)*cb-y(2)*sb)*(ynorm*cb-y(2))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)-sb*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(y(0)*cb+(yplus2a(2))*sb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+1/(realT)2*(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*((-2+2*nu)*diff12nu*cotb*((1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+diff22nu*y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+1/(realT)2*diff22nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+diff22nu*sb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-diff22nu*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)2*diff22nu*(y(0)*cb+(yplus2a(2))*sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)2*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(2*y(2)+4*a)+(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*nu*y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a))-1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb*sb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(sb-(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))+(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb*sb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(sb-(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))*(2*y(2)+4*a)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb))-1/(realT)2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(sb-(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))*(2*y(2)+4*a)+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-(yplus2a(2))*sb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a)-sb*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(y(0)*cb+(yplus2a(2))*sb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+1/(realT)2*(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb))))/(realT)M_PI/(realT)(1-nu))+B(2)*(1/(realT)8*(diff22nu*(y(1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(y(1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)*y(2)+y(1)*ynorm*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))+y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))-(1/(realT)2*y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*(2*y(2)+4*a)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))+y(1)*sb*((1/(realT)ynorm*cb*y(2)-1)/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-(ynorm*cb-y(2))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)-(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*(-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(1/(realT)2*y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*(2*y(2)+4*a)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))-diff22nu*y(1)*sb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)2*diff22nu*y(1)*sb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+y(1)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)2*y(1)*(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(2*y(2)+4*a)-y(1)*(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+y(1)*(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a)))/(realT)M_PI/(realT)(1-nu));
  e(0,1) = 1/(realT)2*B(0)*(1/(realT)8*(diff22nu*(-2/y(0)/(realT)(1+yy(1,1)/(realT)yy(0,0))+1/(realT)(y(0)*cb-y(2)*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)+yy(1,1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)-2*yy(1,1)*ynorm*sb/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))+1/(realT)(y(0)*cb+(yplus2a(2))*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)-2*yy(1,1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))-y(0)*(1/(realT)ynorm/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))-y(0)*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(2),(realT)2)*y(1)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*y(1))-cb*((ynorm*sb-y(0))/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-y(1)*cb*(1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*sb*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-(ynorm*sb-y(0))/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)+1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*((-2+2*nu)*diff12nu*(-1/y(0)/(realT)(1+yy(1,1)/(realT)yy(0,0))+1/(realT)(y(0)*cb+(yplus2a(2))*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)-2*yy(1,1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*Power((realT)cotb,(realT)2)+diff12nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((1-2*nu-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb-y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))-diff12nu*yy(1,1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((1-2*nu-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb-y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+diff12nu*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)*cotb+y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0))+diff12nu*cb*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-diff12nu*yy(1,1)*cb*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-diff12nu*yy(1,1)*cb*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+a*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-3*a*yy(1,1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)+(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((-1+2*nu)*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((-1+2*nu)*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu)*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-2*a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1))+(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)-a*(yplus2a(2))*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)-a*(yplus2a(2))*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)-a*(yplus2a(2))*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(1)*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)*cotb+diff22nu/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb*y(1)*cb)+2*a*(yplus2a(2))*cb*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(1)*(1/(realT)8*((-1+2*nu)*(1/(realT)ynorm*y(1)/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*(1/(realT)ynorm*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))+yy(0,0)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(2),(realT)2)*y(1)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*y(1))+(y(0)*cb-y(2)*sb)/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*sb*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb-y(2)*sb)*(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-(y(0)*cb-y(2)*sb)*(ynorm*sb-y(0))/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)+(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*((diff22nu*Power((realT)cotb,(realT)2)+nu)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-(diff22nu*Power((realT)cotb,(realT)2)+1)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-diff12nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu)*y(0)*cotb+nu*(yplus2a(2))-a+a*y(0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+diff12nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-a*y(0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1))+diff12nu*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((y(0)*cb+(yplus2a(2))*sb)*cb-a*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-diff12nu*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*y(1)/(realT)cb+a*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*y(1))+3*a*y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*y(0)-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*y(0)*cotb+a)-yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*y(0)*cotb+a)*y(1)+yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*y(1)+yy(0,0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*y(1)+yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a*y(1)+3*a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*y(1))-(y(2)+a)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(-cb*sb+a*y(0)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+(y(2)+a)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-3*a*y(0)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)/(realT)cb*y(1)+1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*y(1)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb))*y(1)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*y(1))))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(2)*(1/(realT)8*sb*((ynorm*sb-y(0))/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))/(realT)M_PI/(realT)(1-nu)+1/(realT)8*y(1)*sb*(1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*sb*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-(ynorm*sb-y(0))/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)+1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*(1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-yy(1,1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-yy(1,1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+yy(1,1)*cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yy(1,1)*cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))-(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))+yy(1,1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1))+(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)*cb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))+y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-2*a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(0)*(1/(realT)8*(diff12nu*(1/(realT)ynorm*y(0)/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*((1/(realT)ynorm*y(0)-sb)/(realT)(ynorm-y(0)*sb-y(2)*cb)+(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))-yy(1,1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(0)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(2),(realT)2)*y(0)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*y(0)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*y(0)-cb*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-1/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb))))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*((diff22nu*Power((realT)cotb,(realT)2)-nu)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-(diff22nu*Power((realT)cotb,(realT)2)+1-2*nu)*cb*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+diff12nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(y(0)*cotb*(1-2*nu-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+nu*(yplus2a(2))-a+yy(1,1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-diff12nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((1-2*nu-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+a*yy(0,0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-yy(1,1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-yy(1,1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0))-diff12nu*cb*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-a*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+3*a*yy(0,0)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(-2*nu+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*y(0)*cotb-a)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*y(0)*cotb-a)*y(0)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*diff12nu*cotb-yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*y(0)-yy(1,1)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*y(0)-yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a*y(0)-3*a*yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*y(0))-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(Power((realT)cb,(realT)2)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb+a*cb)+a*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb+a*cb)*y(0)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*diff12nu*cb*cotb+a*(yplus2a(2))*cb*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-3*a*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*y(0)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*y(0)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-a*cb*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))+a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*y(0)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(0))))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(1)*(1/(realT)8*(diff22nu*(2*y(1)/(realT)yy(0,0)/(realT)(1+yy(1,1)/(realT)yy(0,0))-y(1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(y(1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)*y(0)-y(1)*ynorm*sb/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb-y(2)*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*y(0)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb+(yplus2a(2))*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))+y(1)*(1/(realT)ynorm/(realT)(ynorm-y(2))+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))+y(0)*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(0)-1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(2),(realT)2)*y(0)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*y(0)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*y(0))-y(1)*(cb/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb-y(2)*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-(y(0)*cb-y(2)*sb)/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)+cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*diff12nu*(y(1)/(realT)yy(0,0)/(realT)(1+yy(1,1)/(realT)yy(0,0))-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*y(0)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb+(yplus2a(2))*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*Power((realT)cotb,(realT)2)-diff12nu*y(1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)+diff12nu*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)*cotb+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))+diff12nu*y(1)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+diff12nu*y(1)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*y(0)+3*a*y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*y(0)-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))))*y(0)-y(1)*(y(2)+a)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))))*y(0)+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-2*nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))+2*nu*yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))+a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)))-y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((-2+2*nu)*cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb)*y(0)-y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((-2+2*nu)*cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*y(0)-2*a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)cb*y(0)))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(2)*(1/(realT)8*(diff12nu*sb*((1/(realT)ynorm*y(0)-sb)/(realT)(ynorm-y(0)*sb-y(2)*cb)+(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-yy(1,1)*sb*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-1/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*(-sb*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-(y(0)*cb+(yplus2a(2))*sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0))+(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))-yy(0,0)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))+y(0)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(0)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0))+(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(sb*(cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(sb*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)+cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)-2*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*a*(yplus2a(2))*y(0)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*y(0)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-a*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))+a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*y(0)-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(0))))/(realT)M_PI/(realT)(1-nu));
  e(0,2) = 1/(realT)2*B(0)*(1/(realT)8*(diff22nu*(y(1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(y(1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)*y(2)+y(1)*ynorm*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(1/(realT)2*y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*(2*y(2)+4*a)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))-y(0)*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(2)-1/(realT)ynorm/(realT)Power((realT)ynorm-y(2),(realT)2)*(1/(realT)ynorm*y(2)-1)-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*y(2)+4*a)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1))-y(1)*cb*(1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*sb*y(2)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-(ynorm*sb-y(0))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)+1/(realT)2/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*(2*y(2)+4*a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*((-2+2*nu)*diff12nu*(-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(1/(realT)2*y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*(2*y(2)+4*a)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*Power((realT)cotb,(realT)2)-diff12nu*y(1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((1-2*nu-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb-y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+diff12nu*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(1/(realT)2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)*cotb+y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+1/(realT)2*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a))-diff12nu*y(1)*cb*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)2*diff12nu*y(1)*cb*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)*cotb-3/(realT)2*a*y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*(2*y(2)+4*a)+y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((-1+2*nu)*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)2*y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*((-1+2*nu)*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(2*y(2)+4*a)-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu)*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-1/(realT)2*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a))+y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)-a*(yplus2a(2))*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)2*y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)-a*(yplus2a(2))*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(2*y(2)+4*a)-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)-a*(yplus2a(2))*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+diff22nu*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))*cb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)*(diff12nu*cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)*cotb+1/(realT)2*diff22nu/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb*(2*y(2)+4*a)*cb)-a*cb*cotb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+a*(yplus2a(2))*cb*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a)))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(1)*(1/(realT)8*((-1+2*nu)*((1/(realT)ynorm*y(2)-1)/(realT)(ynorm-y(2))+(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*((1/(realT)ynorm*y(2)-cb)/(realT)(ynorm-y(0)*sb-y(2)*cb)+(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))+yy(0,0)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(2)-1/(realT)ynorm/(realT)Power((realT)ynorm-y(2),(realT)2)*(1/(realT)ynorm*y(2)-1)-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*y(2)+4*a)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1))-sb*(ynorm*sb-y(0))/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)+(y(0)*cb-y(2)*sb)/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*sb*y(2)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb-y(2)*sb)*(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-(y(0)*cb-y(2)*sb)*(ynorm*sb-y(0))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)+sb*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+1/(realT)2*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*(2*y(2)+4*a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-1/(realT)2*(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*((diff22nu*Power((realT)cotb,(realT)2)+nu)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-(diff22nu*Power((realT)cotb,(realT)2)+1)*cb*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-diff12nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu)*y(0)*cotb+nu*(yplus2a(2))-a+a*y(0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+diff12nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu-1/(realT)2*a*y(0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)-yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-1/(realT)2*yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a))+diff12nu*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((y(0)*cb+(yplus2a(2))*sb)*cb-a*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-diff12nu*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb*sb-1/(realT)2*a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*(2*y(2)+4*a)/(realT)cb+1/(realT)2*a*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*(2*y(2)+4*a))-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)*cotb+3/(realT)2*a*y(0)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*(2*y(2)+4*a)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*y(0)*cotb+a)-yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*y(0)*cotb+a)-yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*y(0)*cotb+a)*(2*y(2)+4*a)+1/(realT)2*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(2*y(2)+4*a)+yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+1/(realT)2*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a*(2*y(2)+4*a)+3/(realT)2*a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*(2*y(2)+4*a))+cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-cb*sb+a*y(0)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)))-(y(2)+a)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(-cb*sb+a*y(0)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+(y(2)+a)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*y(0)-3/(realT)2*a*y(0)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)/(realT)cb*(2*y(2)+4*a)+1/(realT)2/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*(2*y(2)+4*a)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb))-1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb))*(2*y(2)+4*a)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*(2*y(2)+4*a))))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(2)*(1/(realT)8*y(1)*sb*(1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*sb*y(2)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-(ynorm*sb-y(0))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)+1/(realT)2/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*sb*(2*y(2)+4*a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb-y(0))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*(-y(1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-1/(realT)2*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+y(1)*cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+1/(realT)2*y(1)*cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a))-y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))+1/(realT)2*y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))*(2*y(2)+4*a)-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1))+y(1)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)2*y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(2*y(2)+4*a)-y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a)))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(0)*(1/(realT)8*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)*y(0)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-cb*(1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*cb*y(0)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-(ynorm*cb-y(2))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*(diff12nu*(y(1)/(realT)yy(0,0)/(realT)(1+yy(1,1)/(realT)yy(0,0))-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*y(0)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb+(yplus2a(2))*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*cotb-y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)+y(1)*cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+y(1)*cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0))-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-2*nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(0))-y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1-2*nu-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)-y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1-2*nu-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*yplus2a(0)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(0)))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(1)*(1/(realT)8*((-1+2*nu)*sb*((1/(realT)ynorm*y(0)-sb)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-1/(realT)ynorm+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)*y(0)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0))+cb*(ynorm*cb-y(2))/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)+(y(0)*cb-y(2)*sb)/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*cb*y(0)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb-y(2)*sb)*(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-(y(0)*cb-y(2)*sb)*(ynorm*cb-y(2))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)-cb*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*((-2+2*nu)*diff12nu*cotb*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-diff22nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+diff22nu*yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+diff22nu*yy(0,0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+diff22nu*cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-diff22nu*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-diff22nu*(y(0)*cb+(yplus2a(2))*sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)+(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-2*nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))+2*nu*yy(0,0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+2*a*yy(0,0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2))+(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb*sb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(sb-(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(0)*cotb*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))*y(0)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb))-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(sb-(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))*y(0)+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-(yplus2a(2))*cb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+2*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(0)-cb*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb))))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(2)*(1/(realT)8*(diff22nu*(-y(1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(y(1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)*y(0)-y(1)*ynorm*sb/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb-y(2)*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))+y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))-(y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*y(0)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb+(yplus2a(2))*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))+y(1)*sb*(1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*cb*y(0)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(0)-(ynorm*cb-y(2))/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(0)-sb)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(0)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*(y(1)/(realT)yy(0,0)/(realT)(1+yy(1,1)/(realT)yy(0,0))-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*cb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*y(0)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*(2*y(0)*cb+(yplus2a(2))*sb))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))-diff22nu*y(1)*sb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-diff22nu*y(1)*sb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-y(1)*(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(0)-y(1)*(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)+y(1)*(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(0)-sb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)-2*a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(0)))/(realT)M_PI/(realT)(1-nu));
  e(1,2) = 1/(realT)2*B(0)*(1/(realT)8*(diff12nu*((1/(realT)ynorm*y(2)-1)/(realT)(ynorm-y(2))+(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb*((1/(realT)ynorm*y(2)-cb)/(realT)(ynorm-y(0)*sb-y(2)*cb)+(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))-yy(1,1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(2)-1/(realT)ynorm/(realT)Power((realT)ynorm-y(2),(realT)2)*(1/(realT)ynorm*y(2)-1)-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*y(2)+4*a)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-cb*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-1/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb))))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*((diff22nu*Power((realT)cotb,(realT)2)-nu)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-(diff22nu*Power((realT)cotb,(realT)2)+1-2*nu)*cb*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+diff12nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(y(0)*cotb*(1-2*nu-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+nu*(yplus2a(2))-a+yy(1,1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-diff12nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(1/(realT)2*a*y(0)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+nu-yy(1,1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-1/(realT)2*yy(1,1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a))-diff12nu*sb*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+1/(realT)2*diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(0)*cotb+3/(realT)2*a*y(0)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*(2*y(2)+4*a)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-2*nu+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*y(0)*cotb-a)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(-2*nu+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*y(0)*cotb-a)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*y(0)*cotb-a)*(2*y(2)+4*a)-1/(realT)2*yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(2*y(2)+4*a)-yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-1/(realT)2*yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a*(2*y(2)+4*a)-3/(realT)2*a*yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*(2*y(2)+4*a))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(Power((realT)cb,(realT)2)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb+a*cb)+a*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))))-(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(Power((realT)cb,(realT)2)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb+a*cb)+a*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*(y(0)*cb+(yplus2a(2))*sb)*cotb+a*cb)*(2*y(2)+4*a)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*diff12nu*sb*cotb+a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+a*(yplus2a(2))*sb*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-3/(realT)2*a*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*(2*y(2)+4*a)+1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*(2*y(2)+4*a)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(yy(1,1)*Power((realT)cb,(realT)2)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-a*sb*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))+1/(realT)2*a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(2*y(2)+4*a)-a*(y(0)*cb+(yplus2a(2))*sb)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1))))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(1)*(1/(realT)8*(diff22nu*(y(1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(y(1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)*y(2)+y(1)*ynorm*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(1/(realT)2*y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*(2*y(2)+4*a)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))+y(0)*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(2))*y(2)-1/(realT)ynorm/(realT)Power((realT)ynorm-y(2),(realT)2)*(1/(realT)ynorm*y(2)-1)-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*y(2)+4*a)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1))-y(1)*(-sb/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb-y(2)*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-(y(0)*cb-y(2)*sb)/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)+sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)-1/(realT)2*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*diff12nu*(-y(1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2)*sb/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(1/(realT)2*y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)*(2*y(2)+4*a)-y(1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*y(0))/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*Power((realT)cotb,(realT)2)-diff12nu*y(1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*((-1+2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*cotb+y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+diff12nu*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(-1/(realT)2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)*cotb-y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)-1/(realT)2*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a))+diff12nu*y(1)*cotb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+1/(realT)2*diff12nu*y(1)*cotb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*(2*y(2)+4*a)-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)*cotb+3/(realT)2*a*y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)5/(realT)2)*(2*y(2)+4*a)+y(1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))))-1/(realT)2*y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))))*(2*y(2)+4*a)-y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu*y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+1/(realT)2*a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))*(2*y(2)+4*a)-a*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)))+y(1)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((-2+2*nu)*cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb)-1/(realT)2*y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((-2+2*nu)*cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb)*(2*y(2)+4*a)-y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*((-2+2*nu)*cb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)+y(1)*(y(2)+a)*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*((1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)cb)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)2*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)cb*(2*y(2)+4*a)+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)cb-a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)/(realT)cb*(2*y(2)+4*a)))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(2)*(1/(realT)8*(diff12nu*sb*((1/(realT)ynorm*y(2)-cb)/(realT)(ynorm-y(0)*sb-y(2)*cb)+(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-yy(1,1)*sb*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(2)-1/(realT)ynorm/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*(1/(realT)ynorm*y(2)-cb)-1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(2*y(2)+4*a)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff12nu*(-sb*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1)+1/(realT)2*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+sb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)2*(y(0)*cb+(yplus2a(2))*sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a))+y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))-1/(realT)2*y(0)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2)))*(2*y(2)+4*a)+y(0)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a)-1/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+1))-1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(sb*(cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))))+(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(sb*(cb-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)2*sb*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*y(2)+4*a)+sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-1/(realT)2*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(1+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*(2*y(2)+4*a)+(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*(2*y(2)+4*a))+1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*(2*y(2)+4*a)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(yy(1,1)*cb*sb-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2)))*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*y(2)+4*a)+cb)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-a*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))+1/(realT)2*a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*(2*y(2)+4*a)-a*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(1/(realT)2/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*(2*y(2)+4*a)+1))))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(0)*(1/(realT)8*(1/(realT)ynorm-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-cb*((ynorm*cb-y(2))/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))/(realT)M_PI/(realT)(1-nu)+1/(realT)8*y(1)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)*y(1)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-cb*(1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*cb*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-(ynorm*cb-y(2))/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*(diff12nu*(-1/y(0)/(realT)(1+yy(1,1)/(realT)yy(0,0))+1/(realT)(y(0)*cb+(yplus2a(2))*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)-2*yy(1,1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))*cotb+1/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-yy(1,1)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-yy(1,1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)-cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+yy(1,1)*cb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yy(1,1)*cb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves))+(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(2*nu/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))+a/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))+y(1)*(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-2*nu/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-2*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1))+(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1-2*nu-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1-2*nu-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)*cb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1-2*nu-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))+y(1)*(y(2)+a)*cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*yplus2a(1)*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(1)*(1/(realT)8*((-1+2*nu)*sb*(1/(realT)ynorm*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-y(0)*(-1/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)*y(1)+1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1))+(y(0)*cb-y(2)*sb)/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*cb*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(y(0)*cb-y(2)*sb)*(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-(y(0)*cb-y(2)*sb)*(ynorm*cb-y(2))/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)-(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*((-2+2*nu)*diff12nu*cotb*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+diff22nu*y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)*(2*nu+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)+diff22nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-diff22nu*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-diff22nu*(y(0)*cb+(yplus2a(2))*sb)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff12nu*cotb-2*nu*y(0)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2))-a*y(0)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))*y(1)+(y(2)+a)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*nu*y(0)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)+yplus2a(2),(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*yplus2a(1)*y(0)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1))+(y(2)+a)/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb*sb+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(sb-(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-(y(2)+a)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(1)*cotb*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(diff22nu*cb-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))*y(1)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))*cotb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(-cb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1))-a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*(sb-(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))-(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))*y(1)+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*(2*(yplus2a(2))*(y(0)*cb+(yplus2a(2))*sb)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1)-(y(0)*cb+(yplus2a(2))*sb)/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)+(y(0)*cb+(yplus2a(2))*sb)*(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1))))/(realT)M_PI/(realT)(1-nu))+1/(realT)2*B(2)*(1/(realT)8*(diff22nu*(1/(realT)(y(0)*cb-y(2)*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb-y(2)*sb,(realT)2))+(ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)+yy(1,1)/(realT)ynorm*sb/(realT)(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb)-2*yy(1,1)*ynorm*sb/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yy(2,2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb-y(2)*sb)+yy(1,1)*cb,(realT)2))-1/(realT)(y(0)*cb+(yplus2a(2))*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)-2*yy(1,1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))+sb*((ynorm*cb-y(2))/(realT)ynorm/(realT)(ynorm-y(0)*sb-y(2)*cb)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb))+y(1)*sb*(1/(realT)(yy(0,0)+yy(1,1)+yy(2,2))*cb*y(1)/(realT)(ynorm-y(0)*sb-y(2)*cb)-(ynorm*cb-y(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yy(2,2),threehalves)/(realT)(ynorm-y(0)*sb-y(2)*cb)*y(1)-(ynorm*cb-y(2))/(realT)(yy(0,0)+yy(1,1)+yy(2,2))/(realT)Power((realT)ynorm-y(0)*sb-y(2)*cb,(realT)2)*y(1)-1/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*y(1)+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu)+1/(realT)4*(diff22nu*(-1/y(0)/(realT)(1+yy(1,1)/(realT)yy(0,0))+1/(realT)(y(0)*cb+(yplus2a(2))*sb)/(realT)(1+yy(1,1)/(realT)Power(y(0)*cb+(yplus2a(2))*sb,(realT)2))+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)+yy(1,1)/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb)-2*yy(1,1)*Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*sb/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)*cb)/(realT)(1+yy(1,1)*(yy(0,0)+yy(1,1)+yplus2a2(2))*Power((realT)sb,(realT)2)/(realT)Power(y(0)*(y(0)*cb+(yplus2a(2))*sb)+yy(1,1)*cb,(realT)2)))+diff22nu*sb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-diff22nu*yy(1,1)*sb/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-diff22nu*yy(1,1)*sb/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)+(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))-yy(1,1)*(y(2)+a)*sb/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(1+(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))+a*(yplus2a(2))/(realT)(yy(0,0)+yy(1,1)+yplus2a2(2)))+y(1)*(y(2)+a)*sb/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(1/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb*y(1)/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)Power((realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb,(realT)2)*(cb+a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*y(1)-(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)*cb+yplus2a(2))/(realT)(Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),half)-y(0)*sb+(yplus2a(2))*cb)*a/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),threehalves)*y(1)-2*a*(yplus2a(2))/(realT)Power((realT)yy(0,0)+yy(1,1)+yplus2a2(2),(realT)2)*y(1)))/(realT)M_PI/(realT)(1-nu));
  return;
}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//     OKADA'S QUAD DISLOCATION ELEMENTS
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

///Calculate the strains from the Okada formulation
/**
 @param[in] xStrike Observation point in strike frame
 @param[in] alpha Alpha mechanical property
 @param[in] depth Depth that is positive in the direction of -up
 @param[in] length Dimension of fault patch in the strike direction
 @param[in] width Dimension of fault patch in the dip direction
 @param[in] dip Radian angle of dip
 @param[in] dislocationStrike Rake vector in the strike frame
 @param[out] slipStrike Slip vector in the strike frame
 @param[out] SStrike Strain tensor in the strike frame
 */
void FaultElementManagerT::CalculateOkadaQuadStrains(const R1Tensor& xStrike,
                                                      const realT alpha,
                                                      const realT depth,
                                                      const realT length,
                                                      const realT width,
                                                      const realT dip,
                                                      const R1Tensor& dislocationStrike,
                                                      R1Tensor& slipStrike,
                                                      R2Tensor& SStrike)
{
  const realT dipDegrees = 180 * dip / M_PI;
  const localIndex i = 0, j = 2, k = 1;
  const realT lhalf = 0.5 * length;
  const realT whalf = 0.5 * width;

  /* Calculate the strain tensor */
  dc3d(alpha, xStrike(0), xStrike(1), xStrike(2),
       depth, dipDegrees, -lhalf, lhalf, -whalf, whalf,
       dislocationStrike(i), dislocationStrike(j), dislocationStrike(k),//strike, dip, normal (the strike frame is strike, normal, dip)
       &slipStrike(0), &slipStrike(1), &slipStrike(2),
       &SStrike(0,0), &SStrike(1,0), &SStrike(2,0),
       &SStrike(0,1), &SStrike(1,1), &SStrike(2,1),
       &SStrike(0,2), &SStrike(1,2), &SStrike(2,2));
}


///*
// C********************************************************************   04680005
// C*****                                                          *****   04690005
// C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   04700005
// C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   04710005
// C*****              CODED BY  Y.OKADA ... SEP.1991              *****   04720005
// C*****              REVISED ... NOV.1991, APR.1992, MAY.1993,   *****   04730005
// C*****                          JUL.1993                        *****   04740005
// C********************************************************************   04750005
// C                                                                       04760005
// C***** INPUT                                                            04770005
// C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           04780005
// C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          04790005
// C*****   DEPTH : DEPTH OF REFERENCE POINT                               04800005
// C*****   DIP   : DIP-ANGLE (DEGREE)                                     04810005
// C*****   AL1,AL2   : FAULT LENGTH RANGE                                 04820005
// C*****   AW1,AW2   : FAULT WIDTH RANGE                                  04830005
// C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              04840005
// C                                                                       04850005
// C***** OUTPUT                                                           04860005
// C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF DISL)               04870005
// C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT=(UNIT OF DISL) /             04880005
// C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH,AL,AW) )04890005
// C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     04900005
// C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  04910005
// C
// */
//int FaultElementManagerT::dc3d(const realT alpha,
//                                const R1Tensor& xobs,
//                                const realT depth, const realT dip,
//                                const realT al1, const realT al2,
//                                const realT aw1, const realT aw2,
//                                const realT disl1, const realT disl2, const realT disl3,
//                                R1Tensor& ux,
//                                R2Tensor& uxx)
//{
//  const realT eps = 1e-6;
//
//  //DECLARE AND SET THE FOLLOWING
//  realT alp1, alp2, alp3, alp4, alp5;
//  realT sd, cd, sdsd, cdcd, sdcd, s2d, c2d;
//  dccon0(alpha, dip,
//         alp1, alp2, alp3, alp4, alp5,
//         sd, cd, sdsd, cdcd, sdcd, s2d, c2d);
//
//  const realT zz = xobs(2);
//
//  realT xi[2];
//  {
//    xi[0] = xobs(0) - al1;
//    xi[1] = xobs(0) - al2;
//    if (isEqual(xi[0],0.0,eps))
//      xi[0] = 0;
//    if (isEqual(xi[1],0.0,eps))
//      xi[1] = 0;
//  }
//
//  /*
//   C======================================                                 05170005
//   C=====  REAL-SOURCE CONTRIBUTION  =====                                 05180005
//   C======================================                                 05190005
//   */
//  realT d, p, q;
//  realT et[2];
//  {
//    d = depth + xobs(2);
//    p = xobs(1) * cd + d * sd;
//    q = xobs(1) * sd - d * cd;
//
//    et[0] = p - aw1;
//    et[1] = p - aw2;
//
//    if (isEqual(q, 0.0, eps))
//      q = 0;
//    if (isEqual(et[0], 0.0, eps))
//      et[0] = 0;
//    if (isEqual(et[1], 0.0, eps))
//      et[1] = 0;
//  }
//
//  /*
//   C--------------------------------                                       05280005
//   C----- REJECT SINGULAR CASE -----                                       05290005
//   C--------------------------------                                       05300005
//   */
//  /*
//   C----- ON FAULT EDGE                                                    05310014
//   */
//  if (isEqual(q,0.0) && ((xi[0] * xi[1] < 0 && isEqual(et[0] * et[1],0)) || (et[0] * et[1] < 0 && isEqual(xi[0] * xi[1],0))))
//    return 1;
//
//  /*
//   C----- ON NEGATIVE EXTENSION OF FAULT EDGE                              05360014
//   */
//  realT r12, r21, r22;
//  realT kxi[2];
//  realT ket[2];
//  {
//    r12 = sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q);
//    r21 = sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q);
//    r22 = sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q);
//
//    kxi[0] = 0;
//    kxi[1] = 0;
//
//    ket[0] = 0;
//    ket[1] = 0;
//    if (xi[0] < 0 && r21 + xi[1] < eps)
//      kxi[0] = 1;
//    if (xi[0] < 0 && r22 + xi[1] < eps)
//      kxi[1] = 1;
//    if (et[0] < 0 && r12 + et[1] < eps)
//      ket[0] = 1;
//    if (et[0] < 0 && r22 + et[1] < eps)
//      ket[1] = 1;
//  }
//
//  //temporary variables
//  int i, j, k;
//  realT xi2_c2, et2_c2, q2_c2, r_c2, r2_c2, r3_c2, r5_c2, y_c2, d_c2, tt_c2;
//  realT alx_c2, x11_c2, x32_c2, ale_c2, y11_c2, y32_c2;
//  realT ey_c2, ez_c2, fy_c2, fz_c2, gy_c2, gz_c2, hy_c2, hz_c2;
//  realT dua[12], du[12], u[12], dub[12], duc[12];
//
//  for (k = 0; k < 2; k++)
//  {
//    for (j = 0; j < 2; j++)
//    {
//      //dccon2(xi[j], et[k], q, sd, cd, kxi[k], ket[j]);
//      dccon2(xi[j], et[k], q, sd, cd, kxi[k], ket[j], xi2_c2, et2_c2, q2_c2, r_c2, r2_c2, r3_c2,
//             r5_c2, y_c2, d_c2, tt_c2, alx_c2, x11_c2, x32_c2, ale_c2, y11_c2, y32_c2, ey_c2, ez_c2,
//             fy_c2, fz_c2, gy_c2, gz_c2, hy_c2, hz_c2);
//
//      //ua(xi[j], et[k], q, disl1, disl2, disl3, dua);
//      ua(xi[j], et[k], q, disl1, disl2, disl3, dua,
//          alp1, alp2,
//          sd, cd,
//          xi2_c2,
//          q2_c2, r_c2,
//          r3_c2, y_c2,
//          d_c2, tt_c2, alx_c2,
//          x11_c2, ale_c2,
//          y11_c2, y32_c2, ey_c2,
//          ez_c2, fy_c2, fz_c2,
//          gy_c2, gz_c2, hy_c2, hz_c2);
//
//      for (i = 0; i <= 9; i += 3)
//      {
//        du[i] = -dua[i];
//        du[i + 1] = -dua[i + 1] * cd + dua[i + 2] * sd;
//        du[i + 2] = -dua[i + 1] * sd - dua[i + 2] * cd;
//        if (i == 9)
//        {
//          du[i] = -du[i];
//          du[i + 1] = -du[i + 1];
//          du[i + 2] = -du[i + 2];
//        }
//      }
//
//      for (i = 0; i < 12; i++)
//      {
//        if (j + k != 1)
//          u[i] += du[i];
//        else
//          u[i] -= du[i];
//      }
//    }
//  }
//
//  /*
//   C=======================================                                05700005
//   C=====  IMAGE-SOURCE CONTRIBUTION  =====                                05710005
//   C=======================================                                05720005
//   */
//
//  d = depth - xobs(2);
//  p = xobs(1) * cd + d * sd;
//  q = xobs(1) * sd - d * cd;
//  et[0] = p - aw1;
//  et[1] = p - aw2;
//  if (fabs(q) < eps)
//    q = 0;
//  if (fabs(et[0]) < eps)
//    et[0] = 0;
//  if (fabs(et[1]) < eps)
//    et[1] = 0;
//
//  /*
//   C--------------------------------                                       05810005
//   C----- REJECT SINGULAR CASE -----                                       05820005
//   C--------------------------------                                       05830005
//   */
//  /*
//   C----- ON FAULT EDGE                                                    05840015
//   */
//  if (isEqual(q, 0.0) && ((xi[0] * xi[1] < 0 && isEqual(et[0] * et[1], 0)) || (et[0] * et[1] < 0 && isEqual(xi[0] * xi[1], 0))))
//    return (1);
//  /*
//   C----- ON NEGATIVE EXTENSION OF FAULT EDGE                              05890015
//   */
//  kxi[0] = 0;
//  kxi[1] = 0;
//  ket[0] = 0;
//  ket[1] = 0;
//  r12 = sqrt(xi[0] * xi[0] + et[1] * et[1] + q * q);
//  r21 = sqrt(xi[1] * xi[1] + et[0] * et[0] + q * q);
//  r22 = sqrt(xi[1] * xi[1] + et[1] * et[1] + q * q);
//  if (xi[0] < 0 && r21 + xi[1] < eps)
//    kxi[0] = 1;
//  if (xi[0] < 0 && r22 + xi[1] < eps)
//    kxi[1] = 1;
//  if (et[0] < 0 && r12 + et[1] < eps)
//    ket[0] = 1;
//  if (et[0] < 0 && r22 + et[1] < eps)
//    ket[1] = 1;
//
//  for (k = 0; k < 2; k++)
//  {
//    for (j = 0; j < 2; j++)
//    {
//      //dccon2(xi[j], et[k], q, sd, cd, kxi[k], ket[j]);
//      dccon2(xi[j], et[k], q, sd, cd, kxi[k], ket[j],
//             xi2_c2, et2_c2, q2_c2,
//             r_c2, r2_c2, r3_c2, r5_c2,
//             y_c2, d_c2, tt_c2,
//             alx_c2, x11_c2, x32_c2,
//             ale_c2, y11_c2, y32_c2,
//             ey_c2, ez_c2, fy_c2, fz_c2,
//             gy_c2, gz_c2, hy_c2, hz_c2);
//
//      //ua(xi[j], et[k], q, disl1, disl2, disl3, dua);
//      ua(xi[j], et[k], q, disl1, disl2, disl3, dua,
//          alp1, alp2,
//          sd, cd,
//          xi2_c2,
//          q2_c2, r_c2,
//          r3_c2, y_c2,
//          d_c2, tt_c2, alx_c2,
//          x11_c2, ale_c2,
//          y11_c2, y32_c2, ey_c2,
//          ez_c2, fy_c2, fz_c2,
//          gy_c2, gz_c2, hy_c2, hz_c2);
//
//      //ub(xi[j], et[k], q, disl1, disl2, disl3, dub);
//      ub(xi[j], et[k], q, disl1, disl2, disl3, dub,
//         alp3,
//         sd, cd, sdsd, cdcd, sdcd,
//         xi2_c2, q2_c2,
//         r_c2, r3_c2,
//         y_c2, d_c2, tt_c2,
//         x11_c2,
//         ale_c2, y11_c2, y32_c2,
//         ey_c2, ez_c2, fy_c2, fz_c2,
//         gy_c2, gz_c2, hy_c2, hz_c2);
//
//      //uc(xi[j], et[k], q, zz, disl1, disl2, disl3, duc);
//      uc(xi[j], et[k], q, zz, disl1, disl2, disl3, duc,
//         alp4, alp5,
//         sd, cd, sdsd, cdcd, sdcd,
//         xi2_c2, et2_c2, q2_c2,
//         r_c2, r2_c2, r3_c2, r5_c2,
//         y_c2, d_c2, x11_c2, x32_c2,
//         y11_c2, y32_c2);
//
//      for (i = 0; i <= 9; i += 3)
//      {
//        du[i] = dua[i] + dub[i] + xobs(2) * duc[i];
//        du[i + 1] = (dua[i + 1] + dub[i + 1] + xobs(2) * duc[i + 1]) * cd - (dua[i + 2] + dub[i + 2] + xobs(2) * duc[i + 2]) * sd;
//        du[i + 2] = (dua[i + 1] + dub[i + 1] - xobs(2) * duc[i + 1]) * sd + (dua[i + 2] + dub[i + 2] - xobs(2) * duc[i + 2]) * cd;
//
//        if (i == 9)
//        {
//          du[9] = du[9] + duc[0];
//          du[10] = du[10] + duc[1] * cd - duc[2] * sd;
//          du[11] = du[11] - duc[1] * sd - duc[2] * cd;
//        }
//      }
//
//      for (i = 0; i < 12; i++)
//      {
//        if (j + k != 1)
//          u[i] += du[i];
//        else
//          u[i] -= du[i];
//      }
//
//    }
//  }
//
//  ux(0) = u[0];
//  ux(1) = u[1];
//  ux(2) = u[2];
//
//  uxx(0,0) = u[3];
//  uxx(1,0) = u[4];
//  uxx(2,0) = u[5];
//  uxx(0,1) = u[6];
//  uxx(1,1) = u[7];
//  uxx(2,1) = u[8];
//  uxx(0,2) = u[9];
//  uxx(1,2) = u[10];
//  uxx(2,2) = u[11];
//
//  return (0);
//}//OK 5/14/12
//
///*
// C********************************************************************   06640005
// C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   06650005
// C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   06660005
// C********************************************************************   06670005
// C                                                                       06680005
// C***** INPUT                                                            06690005
// C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  06700005
// C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              06710005
// C***** OUTPUT                                                           06720005
// C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     06730005
// */
//
//void FaultElementManagerT::ua(realT xi, realT et, realT q, realT disl1, realT disl2,
//                                   realT disl3, realT *u,
//                                   const realT alp1, const realT alp2,
//                                   const realT sd, const realT cd,
//                                   const realT xi2_c2,
//                                   const realT q2_c2, const realT r_c2,
//                                   const realT r3_c2, const realT y_c2,
//                                   const realT d_c2, const realT tt_c2, const realT alx_c2,
//                                   const realT x11_c2, const realT ale_c2,
//                                   const realT y11_c2, const realT y32_c2, const realT ey_c2,
//                                   const realT ez_c2, const realT fy_c2, const realT fz_c2,
//                                   const realT gy_c2, const realT gz_c2, const realT hy_c2,
//                                   const realT hz_c2)
//{
//  const realT pi2 = 6.283185307179586;
//  realT du[12];
//  realT xy, qx, qy;
//  int i;
//
//  for (i = 0; i < 12; i++)
//    u[i] = 0;
//
//  xy = xi * y11_c2;
//  qx = q * x11_c2;
//  qy = q * y11_c2;
//
//  /*
//   C======================================                                 06850005
//   C=====  STRIKE-SLIP CONTRIBUTION  =====                                 06860005
//   C======================================                                 06870005
//   */
//  if (!isEqual(disl1,0))
//  {
//    du[0] = tt_c2 / 2 + alp2 * xi * qy;
//    du[1] = alp2 * q / r_c2;
//    du[2] = alp1 * ale_c2 - alp2 * q * qy;
//    du[3] = -alp1 * qy - alp2 * xi2_c2 * q * y32_c2;
//    du[4] = -alp2 * xi * q / r3_c2;
//    du[5] = alp1 * xy + alp2 * xi * q2_c2 * y32_c2;
//    du[6] = alp1 * xy * sd + alp2 * xi * fy_c2 + d_c2 / 2 * x11_c2;
//    du[7] = alp2 * ey_c2;
//    du[8] = alp1 * (cd / r_c2 + qy * sd) - alp2 * q * fy_c2;
//    du[9] = alp1 * xy * cd + alp2 * xi * fz_c2 + y_c2 / 2 * x11_c2;
//    du[10] = alp2 * ez_c2;
//    du[11] = -alp1 * (sd / r_c2 - qy * cd) - alp2 * q * fz_c2;
//
//    for (i = 0; i < 12; i++)
//      u[i] += disl1 / pi2 * du[i];
//  }
//
//  /*
//   C======================================                                 07040005
//   C=====    DIP-SLIP CONTRIBUTION   =====                                 07050005
//   C======================================                                 07060005
//   */
//  if (!isEqual(disl2, 0))
//  {
//    du[0] = alp2 * q / r_c2;
//    du[1] = tt_c2 / 2 + alp2 * et * qx;
//    du[2] = alp1 * alx_c2 - alp2 * q * qx;
//    du[3] = -alp2 * xi * q / r3_c2;
//    du[4] = -qy / 2 - alp2 * et * q / r3_c2;
//    du[5] = alp1 / r_c2 + alp2 * q2_c2 / r3_c2;
//    du[6] = alp2 * ey_c2;
//    du[7] = alp1 * d_c2 * x11_c2 + xy / 2 * sd + alp2 * et * gy_c2;
//    du[8] = alp1 * y_c2 * x11_c2 - alp2 * q * gy_c2;
//    du[9] = alp2 * ez_c2;
//    du[10] = alp1 * y_c2 * x11_c2 + xy / 2 * cd + alp2 * et * gz_c2;
//    du[11] = -alp1 * d_c2 * x11_c2 - alp2 * q * gz_c2;
//
//    for (i = 0; i < 12; i++)
//      u[i] += disl2 / pi2 * du[i];
//  }
//
//  /*
//   C========================================                               07230005
//   C=====  TENSILE-FAULT CONTRIBUTION  =====                               07240005
//   C========================================                               07250005
//   */
//  if (!isEqual(disl3, 0))
//  {
//    du[0] = -alp1 * ale_c2 - alp2 * q * qy;
//    du[1] = -alp1 * alx_c2 - alp2 * q * qx;
//    du[2] = tt_c2 / 2 - alp2 * (et * qx + xi * qy);
//    du[3] = -alp1 * xy + alp2 * xi * q2_c2 * y32_c2;
//    du[4] = -alp1 / r_c2 + alp2 * q2_c2 / r3_c2;
//    du[5] = -alp1 * qy - alp2 * q * q2_c2 * y32_c2;
//    du[6] = -alp1 * (cd / r_c2 + qy * sd) - alp2 * q * fy_c2;
//    du[7] = -alp1 * y_c2 * x11_c2 - alp2 * q * gy_c2;
//    du[8] = alp1 * (d_c2 * x11_c2 + xy * sd) + alp2 * q * hy_c2;
//    du[9] = alp1 * (sd / r_c2 - qy * cd) - alp2 * q * fz_c2;
//    du[10] = alp1 * d_c2 * x11_c2 - alp2 * q * gz_c2;
//    du[11] = alp1 * (y_c2 * x11_c2 + xy * cd) + alp2 * q * hz_c2;
//    for (i = 0; i < 12; i++)
//      u[i] += disl3 / pi2 * du[i];
//  }
//
//  return;
//} //OK 5/14/12
//
///*
// C********************************************************************   07480005
// C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   07490005
// C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   07500005
// C********************************************************************   07510005
// C                                                                       07520005
// C***** INPUT                                                            07530005
// C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  07540005
// C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              07550005
// C***** OUTPUT                                                           07560005
// C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     07570005
// */
//
//void FaultElementManagerT::ub(realT xi, realT et, realT q, realT disl1, realT disl2,
//                                   realT disl3, realT *u,
//                                   const realT alp3,
//                                   const realT sd, const realT cd, const realT sdsd,
//                                   const realT cdcd, const realT sdcd,
//                                   const realT xi2_c2,
//                                   const realT q2_c2, const realT r_c2,
//                                   const realT r3_c2, const realT y_c2,
//                                   const realT d_c2, const realT tt_c2,
//                                   const realT x11_c2, const realT ale_c2,
//                                   const realT y11_c2, const realT y32_c2, const realT ey_c2,
//                                   const realT ez_c2, const realT fy_c2, const realT fz_c2,
//                                   const realT gy_c2, const realT gz_c2, const realT hy_c2,
//                                   const realT hz_c2)
//{
//  realT du[12], rd, d11, aj1, aj2, aj3, aj4, aj5, aj6, ai1, ai2, ai3, ai4;
//  realT ak1, ak2, ak3, ak4, x;
//  realT rd2, xy, qx, qy;
//  const realT pi2 = 6.283185307179586;
//  int i;
//
//  rd = r_c2 + d_c2;
//  d11 = 1 / (r_c2 * rd);
//  aj2 = xi * y_c2 / rd * d11;
//  aj5 = -(d_c2 + y_c2 * y_c2 / rd) * d11;
//  if (!isEqual(cd,0))
//  {
//    if (isEqual(xi,0))
//    {
//      ai4 = 0;
//    }
//    else
//    {
//      x = sqrt(xi2_c2 + q2_c2);
//      ai4 = 1 / cdcd * (xi / rd * sdcd + 2 * atan(
//          (et * (x + q * cd) + x * (r_c2 + x) * sd) / (xi * (r_c2 + x) * cd)));
//    }
//    ai3 = (y_c2 * cd / rd - ale_c2 + sd * log(rd)) / cdcd;
//    ak1 = xi * (d11 - y11_c2 * sd) / cd;
//    ak3 = (q * y11_c2 - y_c2 * d11) / cd;
//    aj3 = (ak1 - aj2 * sd) / cd;
//    aj6 = (ak3 - aj5 * sd) / cd;
//  }
//  else
//  {
//    rd2 = rd * rd;
//    ai3 = (et / rd + y_c2 * q / rd2 - ale_c2) / 2;
//    ai4 = xi * y_c2 / rd2 / 2;
//    ak1 = xi * q / rd * d11;
//    ak3 = sd / rd * (xi2_c2 * d11 - 1);
//    aj3 = -xi / rd2 * (q2_c2 * d11 - 0.5);
//    aj6 = -y_c2 / rd2 * (xi2_c2 * d11 - 0.5);
//  }
//
//  xy = xi * y11_c2;
//  ai1 = -xi / rd * cd - ai4 * sd;
//  ai2 = log(rd) + ai3 * sd;
//  ak2 = 1 / r_c2 + ak3 * sd;
//  ak4 = xy * cd - ak1 * sd;
//  aj1 = aj5 * cd - aj6 * sd;
//  aj4 = -xy - aj2 * cd + aj3 * sd;
//
//  for (i = 0; i < 12; i++)
//    u[i] = 0;
//  qx = q * x11_c2;
//  qy = q * y11_c2;
//
//  /*
//   C======================================                                 08030005
//   C=====  STRIKE-SLIP CONTRIBUTION  =====                                 08040005
//   C======================================                                 08050005
//   */
//  if (!isEqual(disl1,0))
//  {
//    du[0] = -xi * qy - tt_c2 - alp3 * ai1 * sd;
//    du[1] = -q / r_c2 + alp3 * y_c2 / rd * sd;
//    du[2] = q * qy - alp3 * ai2 * sd;
//    du[3] = xi2_c2 * q * y32_c2 - alp3 * aj1 * sd;
//    du[4] = xi * q / r3_c2 - alp3 * aj2 * sd;
//    du[5] = -xi * q2_c2 * y32_c2 - alp3 * aj3 * sd;
//    du[6] = -xi * fy_c2 - d_c2 * x11_c2 + alp3 * (xy + aj4) * sd;
//    du[7] = -ey_c2 + alp3 * (1 / r_c2 + aj5) * sd;
//    du[8] = q * fy_c2 - alp3 * (qy - aj6) * sd;
//    du[9] = -xi * fz_c2 - y_c2 * x11_c2 + alp3 * ak1 * sd;
//    du[10] = -ez_c2 + alp3 * y_c2 * d11 * sd;
//    du[11] = q * fz_c2 + alp3 * ak2 * sd;
//    for (i = 0; i < 12; i++)
//      u[i] += disl1 / pi2 * du[i];
//  }
//
//  /*
//   C======================================                                 08220005
//   C=====    DIP-SLIP CONTRIBUTION   =====                                 08230005
//   C======================================                                 08240005
//   */
//  if (!isEqual(disl2,0))
//  {
//    du[0] = -q / r_c2 + alp3 * ai3 * sdcd;
//    du[1] = -et * qx - tt_c2 - alp3 * xi / rd * sdcd;
//    du[2] = q * qx + alp3 * ai4 * sdcd;
//    du[3] = xi * q / r3_c2 + alp3 * aj4 * sdcd;
//    du[4] = et * q / r3_c2 + qy + alp3 * aj5 * sdcd;
//    du[5] = -q2_c2 / r3_c2 + alp3 * aj6 * sdcd;
//    du[6] = -ey_c2 + alp3 * aj1 * sdcd;
//    du[7] = -et * gy_c2 - xy * sd + alp3 * aj2 * sdcd;
//    du[8] = q * gy_c2 + alp3 * aj3 * sdcd;
//    du[9] = -ez_c2 - alp3 * ak3 * sdcd;
//    du[10] = -et * gz_c2 - xy * cd - alp3 * xi * d11 * sdcd;
//    du[11] = q * gz_c2 - alp3 * ak4 * sdcd;
//    for (i = 0; i < 12; i++)
//      u[i] += disl2 / pi2 * du[i];
//  }
//
//  /*
//   C========================================                               08410005
//   C=====  TENSILE-FAULT CONTRIBUTION  =====                               08420005
//   C========================================                               08430005
//   */
//  if (!isEqual(disl3,0))
//  {
//    du[0] = q * qy - alp3 * ai3 * sdsd;
//    du[1] = q * qx + alp3 * xi / rd * sdsd;
//    du[2] = et * qx + xi * qy - tt_c2 - alp3 * ai4 * sdsd;
//    du[3] = -xi * q2_c2 * y32_c2 - alp3 * aj4 * sdsd;
//    du[4] = -q2_c2 / r3_c2 - alp3 * aj5 * sdsd;
//    du[5] = q * q2_c2 * y32_c2 - alp3 * aj6 * sdsd;
//    du[6] = q * fy_c2 - alp3 * aj1 * sdsd;
//    du[7] = q * gy_c2 - alp3 * aj2 * sdsd;
//    du[8] = -q * hy_c2 - alp3 * aj3 * sdsd;
//    du[9] = q * fz_c2 + alp3 * ak3 * sdsd;
//    du[10] = q * gz_c2 + alp3 * xi * d11 * sdsd;
//    du[11] = -q * hz_c2 + alp3 * ak4 * sdsd;
//    for (i = 0; i < 12; i++)
//      u[i] += disl3 / pi2 * du[i];
//  }
//
//  return;
//} //OK 5/14/12
//
///*
// C********************************************************************   08660005
// C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-C)             *****   08670005
// C*****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****   08680005
// C********************************************************************   08690005
// C                                                                       08700005
// C***** INPUT                                                            08710005
// C*****   XI,ET,Q,Z   : STATION COORDINATES IN FAULT SYSTEM              08720005
// C*****   DISL1-DISL3 : STRIKE-, DIP-, TENSILE-DISLOCATIONS              08730005
// C***** OUTPUT                                                           08740005
// C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     08750005
// */
//
//void FaultElementManagerT::uc(const realT xi, const realT et, const realT q, const realT z,
//                               const realT disl1, const realT disl2,
//                               const realT disl3, realT *u,
//                               const realT alp4, const realT alp5,
//                               const realT sd, const realT cd, const realT sdsd,
//                               const realT cdcd, const realT sdcd,
//                               const realT xi2_c2, const realT et2_c2,
//                               const realT q2_c2, const realT r_c2, const realT r2_c2,
//                               const realT r3_c2, const realT r5_c2, const realT y_c2,
//                               const realT d_c2, const realT x11_c2, const realT x32_c2,
//                               const realT y11_c2, const realT y32_c2)
//
//{
//  realT du[12], c, x53, y53, z53, h, z32, ppy, ppz, y0, z0;
//  realT qq, qy, qr, qqy, qqz, xy;//qx
//  realT cdr, yy0;//cqx
//  const realT pi2 = 6.283185307179586;
//  int i;
//
//  c = d_c2 + z;
//  x53 = (8 * r2_c2 + 9 * r_c2 * xi + 3 * xi2_c2) * x11_c2 * x11_c2 * x11_c2 / r2_c2;
//  y53 = (8 * r2_c2 + 9 * r_c2 * et + 3 * et2_c2) * y11_c2 * y11_c2 * y11_c2 / r2_c2;
//  h = q * cd - z;
//  z32 = sd / r3_c2 - h * y32_c2;
//  z53 = 3 * sd / r5_c2 - h * y53;
//  y0 = y11_c2 - xi2_c2 * y32_c2;
//  z0 = z32 - xi2_c2 * z53;
//  ppy = cd / r3_c2 + q * y32_c2 * sd;
//  ppz = sd / r3_c2 - q * y32_c2 * cd;
//  qq = z * y32_c2 + z32 + z0;
//  qqy = 3 * c * d_c2 / r5_c2 - qq * sd;
//  qqz = 3 * c * y_c2 / r5_c2 - qq * cd + q * y32_c2;
//  xy = xi * y11_c2;
//  //qx = q * x11_c2;
//  qy = q * y11_c2;
//  qr = 3 * q / r5_c2;
//  //cqx = c * q * x53;
//  cdr = (c + d_c2) / r3_c2;
//  yy0 = y_c2 / r3_c2 - y0 * cd;
//
//  for (i = 0; i < 12; i++)
//    u[i] = 0;
//
//  /*
//   C======================================                                 09050005
//   C=====  STRIKE-SLIP CONTRIBUTION  =====                                 09060005
//   C======================================                                 09070005
//   */
//  if (!isEqual(disl1, 0))
//  {
//    du[0] = alp4 * xy * cd - alp5 * xi * q * z32;
//    du[1] = alp4 * (cd / r_c2 + 2 * qy * sd) - alp5 * c * q / r3_c2;
//    du[2] = alp4 * qy * cd - alp5 * (c * et / r3_c2 - z * y11_c2 + xi2_c2 * z32);
//    du[3] = alp4 * y0 * cd - alp5 * q * z0;
//    du[4] = -alp4 * xi * (cd / r3_c2 + 2 * q * y32_c2 * sd) + alp5 * c * xi * qr;
//    du[5] = -alp4 * xi * q * y32_c2 * cd + alp5 * xi * (3 * c * et / r5_c2 - qq);
//    du[6] = -alp4 * xi * ppy * cd - alp5 * xi * qqy;
//    du[7] = alp4 * 2 * (d_c2 / r3_c2 - y0 * sd) * sd - y_c2 / r3_c2 * cd - alp5 * (cdr * sd - et / r3_c2 - c * y_c2 * qr);
//    du[8] = -alp4 * q / r3_c2 + yy0 * sd + alp5 * (cdr * cd + c * d_c2 * qr - (y0 * cd + q * z0) * sd);
//    du[9] = alp4 * xi * ppz * cd - alp5 * xi * qqz;
//    du[10] = alp4 * 2 * (y_c2 / r3_c2 - y0 * cd) * sd + d_c2 / r3_c2 * cd - alp5 * (cdr * cd + c * d_c2 * qr);
//    du[11] = yy0 * cd - alp5 * (cdr * sd - c * y_c2 * qr - y0 * sdsd + q * z0 * cd);
//    for (i = 0; i < 12; i++)
//      u[i] += disl1 / pi2 * du[i];
//  }
//
//  /*
//   C======================================                                 09250005
//   C=====    DIP-SLIP CONTRIBUTION   =====                                 09260005
//   C======================================                                 09270005
//   */
//  if (!isEqual(disl2, 0))
//  {
//    du[0] = alp4 * cd / r_c2 - qy * sd - alp5 * c * q / r3_c2;
//    du[1] = alp4 * y_c2 * x11_c2 - alp5 * c * et * q * x32_c2;
//    du[2] = -d_c2 * x11_c2 - xy * sd - alp5 * c * (x11_c2 - q2_c2 * x32_c2);
//    du[3] = -alp4 * xi / r3_c2 * cd + alp5 * c * xi * qr + xi * q * y32_c2 * sd;
//    du[4] = -alp4 * y_c2 / r3_c2 + alp5 * c * et * qr;
//    du[5] = d_c2 / r3_c2 - y0 * sd + alp5 * c / r3_c2 * (1 - 3 * q2_c2 / r2_c2);
//    du[6] = -alp4 * et / r3_c2 + y0 * sdsd - alp5 * (cdr * sd - c * y_c2 * qr);
//    du[7] = alp4 * (x11_c2 - y_c2 * y_c2 * x32_c2) - alp5 * c * ((d_c2 + 2 * q * cd) * x32_c2 - y_c2 * et * q * x53);
//    du[8] = xi * ppy * sd + y_c2 * d_c2 * x32_c2 + alp5 * c * ((y_c2 + 2 * q * sd) * x32_c2 - y_c2 * q2_c2 * x53);
//    du[9] = -q / r3_c2 + y0 * sdcd - alp5 * (cdr * cd + c * d_c2 * qr);
//    du[10] = alp4 * y_c2 * d_c2 * x32_c2 - alp5 * c * ((y_c2 - 2 * q * sd) * x32_c2 + d_c2 * et * q * x53);
//    du[11] = -xi * ppz * sd + x11_c2 - d_c2 * d_c2 * x32_c2 - alp5 * c * ((d_c2 - 2 * q * cd) * x32_c2 - d_c2 * q2_c2 * x53);
//    for (i = 0; i < 12; i++)
//      u[i] += disl2 / pi2 * du[i];
//  }
//
//  /*
//   C========================================                               09440005
//   C=====  TENSILE-FAULT CONTRIBUTION  =====                               09450005
//   C========================================                               09460005
//   */
//  if (!isEqual(disl3, 0))
//  {
//    du[0] = -alp4 * (sd / r_c2 + qy * cd) - alp5 * (z * y11_c2 - q2_c2 * z32);
//    du[1] = alp4 * 2 * xy * sd + d_c2 * x11_c2 - alp5 * c * (x11_c2 - q2_c2 * x32_c2);
//    du[2] = alp4 * (y_c2 * x11_c2 + xy * cd) + alp5 * q * (c * et * x32_c2 + xi * z32);
//    du[3] = alp4 * xi / r3_c2 * sd + xi * q * y32_c2 * cd + alp5 * xi * (3 * c * et / r5_c2 - 2 * z32 - z0);
//    du[4] = alp4 * 2 * y0 * sd - d_c2 / r3_c2 + alp5 * c / r3_c2 * (1 - 3 * q2_c2 / r2_c2);
//    du[5] = -alp4 * yy0 - alp5 * (c * et * qr - q * z0);
//    du[6] = alp4 * (q / r3_c2 + y0 * sdcd) + alp5 * (z / r3_c2 * cd + c * d_c2 * qr - q * z0 * sd);
//    du[7] = -alp4 * 2 * xi * ppy * sd - y_c2 * d_c2 * x32_c2 + alp5 * c * ((y_c2 + 2 * q * sd) * x32_c2 - y_c2 * q2_c2 * x53);
//    du[8] = -alp4 * (xi * ppy * cd - x11_c2 + y_c2 * y_c2 * x32_c2) + alp5 * (c * ((d_c2 + 2 * q * cd) * x32_c2 - y_c2 * et * q * x53) + xi * qqy);
//    du[9] = -et / r3_c2 + y0 * cdcd - alp5 * (z / r3_c2 * sd - c * y_c2 * qr - y0 * sdsd + q * z0 * cd);
//    du[10] = alp4 * 2 * xi * ppz * sd - x11_c2 + d_c2 * d_c2 * x32_c2 - alp5 * c * ((d_c2 - 2 * q * cd) * x32_c2 - d_c2 * q2_c2 * x53);
//    du[11] = alp4 * (xi * ppz * cd + y_c2 * d_c2 * x32_c2) + alp5 * (c * ((y_c2 - 2 * q * sd) * x32_c2 + d_c2 * et * q * x53) + xi * qqz);
//    for (i = 0; i < 12; i++)
//      u[i] += disl3 / pi2 * du[i];
//  }
//
//  return;
//} //OK 5/14/12
//
///*
// C*******************************************************************    09720005
// C*****   CALCULATE MEDIUM CONSTANTS AND FAULT-DIP CONSTANTS    *****    09730005
// C*******************************************************************    09740005
// C                                                                       09750005
// C***** INPUT                                                            09760005
// C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           09770005
// C*****   DIP   : DIP-ANGLE (DEGREE)                                     09780005
// C### CAUTION ### IF COS(DIP) IS SUFFICIENTLY SMALL, IT IS SET TO ZERO   09790005
// */
//
//void FaultElementManagerT::dccon0(const realT alpha, const realT dip,
//                                   realT& alp1, realT& alp2, realT& alp3,
//                                   realT& alp4, realT& alp5,
//                                   realT& sd, realT& cd,
//                                   realT& sdsd, realT& cdcd, realT& sdcd,
//                                   realT& s2d, realT& c2d)
//{
//  const realT eps = 1e-6, pi2 = 6.283185307179586;
//
//  alp1 = (1 - alpha) / 2;
//  alp2 = alpha / 2;
//  alp3 = (1 - alpha) / alpha;
//  alp4 = 1 - alpha;
//  alp5 = alpha;
//
//  const realT p18 = pi2 / 360;
//  sd = sin(dip * p18);
//  cd = cos(dip * p18);
//  if (fabs(cd) < eps)
//  {
//    cd = 0;
//    if (sd > 0)
//      sd = 1;
//    if (sd < 0)
//      sd = -1;
//  }
//  sdsd = sd * sd;
//  cdcd = cd * cd;
//  sdcd = sd * cd;
//  s2d = 2 * sdcd;
//  c2d = cdcd - sdsd;
//
//  return;
//} //OK 5/14/12
//
///*
// C********************************************************************** 10090005
// C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR POINT SOURCE    ***** 10100005
// C********************************************************************** 10110005
// C                                                                       10120005
// C***** INPUT                                                            10130005
// C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    10140005
// C### CAUTION ### IF X,Y,D ARE SUFFICIENTLY SMALL, THEY ARE SET TO ZERO  10150005
// */
//
//void FaultElementManagerT::dccon1(const realT sd, const realT cd, realT x, realT y, realT d,
//                                       realT& p_c1, realT& q_c1, realT& s_c1, realT& t_c1,
//                                       realT& xy_c1, realT& x2_c1, realT& y2_c1, realT& d2_c1,
//                                       realT& r_c1, realT& r2_c1, realT& r3_c1, realT& r5_c1,
//                                       realT& a3_c1, realT& a5_c1, realT& b3_c1, realT& c3_c1,
//                                       realT& qr_c1, realT& qrx_c1, realT& uy_c1, realT& uz_c1,
//                                       realT& vy_c1, realT& vz_c1, realT& wy_c1, realT& wz_c1)
//{
//  const realT eps = 1e-6;
//
//  if (fabs(x) < eps)
//    x = 0;
//  if (fabs(y) < eps)
//    y = 0;
//  if (fabs(d) < eps)
//    d = 0;
//  p_c1 = y * cd + d * sd;
//  q_c1 = y * sd - d * cd;
//  s_c1 = p_c1 * sd + q_c1 * cd;
//  t_c1 = p_c1 * cd - q_c1 * sd;
//  xy_c1 = x * y;
//  x2_c1 = x * x;
//  y2_c1 = y * y;
//  d2_c1 = d * d;
//  r2_c1 = x2_c1 + y2_c1 + d2_c1;
//  r_c1 = sqrt(r2_c1);
//  if (isEqual(r_c1,0))
//    return;
//  r3_c1 = r_c1 * r2_c1;
//  r5_c1 = r3_c1 * r2_c1;
//  //const realT r7 = r5_c1 * r2_c1;
//
//  a3_c1 = 1 - 3 * x2_c1 / r2_c1;
//  a5_c1 = 1 - 5 * x2_c1 / r2_c1;
//  b3_c1 = 1 - 3 * y2_c1 / r2_c1;
//  c3_c1 = 1 - 3 * d2_c1 / r2_c1;
//
//  qr_c1 = 3 * q_c1 / r5_c1;
//  qrx_c1 = 5 * qr_c1 * x / r2_c1;
//
//  uy_c1 = sd - 5 * y * q_c1 / r2_c1;
//  uz_c1 = cd + 5 * d * q_c1 / r2_c1;
//  vy_c1 = s_c1 - 5 * y * p_c1 * q_c1 / r2_c1;
//  vz_c1 = t_c1 + 5 * d * p_c1 * q_c1 / r2_c1;
//  wy_c1 = uy_c1 + sd;
//  wz_c1 = uz_c1 + cd;
//
//  return;
//} //OK 5/14/12
//
///*
// C********************************************************************** 10590005
// C*****   CALCULATE STATION GEOMETRY CONSTANTS FOR FINITE SOURCE   ***** 10600005
// C********************************************************************** 10610005
// C                                                                       10620005
// C***** INPUT                                                            10630005
// C*****   XI,ET,Q : STATION COORDINATES IN FAULT SYSTEM                  10640005
// C*****   SD,CD   : SIN, COS OF DIP-ANGLE                                10650005
// C*****   KXI,KET : KXI=1, KET=1 MEANS R+XI<EPS, R+ET<EPS, RESPECTIVELY  10660005
// C                                                                       10670005
// C### CAUTION ### IF XI,ET,Q ARE SUFFICIENTLY SMALL, THey_c2 ARE SET TO ZER010680005
// */
//
//void FaultElementManagerT::dccon2(realT xi, realT et, realT q, const realT sd, const realT cd,
//                                   const realT kxi, const realT ket,
//                                   realT& xi2_c2, realT& et2_c2, realT& q2_c2,
//                                   realT& r_c2, realT& r2_c2, realT& r3_c2, realT& r5_c2,
//                                   realT& y_c2, realT& d_c2, realT& tt_c2, realT& alx_c2,
//                                   realT& x11_c2, realT& x32_c2, realT& ale_c2, realT& y11_c2,
//                                   realT& y32_c2,
//                                   realT& ey_c2, realT& ez_c2, realT& fy_c2,
//                                   realT& fz_c2, realT& gy_c2, realT& gz_c2, realT& hy_c2,
//                                   realT& hz_c2)
//{
//  const realT eps = 1e-6;
//
//  if (isEqual(xi, 0.0, eps))
//    xi = 0;
//  if (isEqual(et, 0.0, eps))
//    et = 0;
//  if (isEqual(q, 0.0, eps))
//    q = 0;
//
//  xi2_c2 = xi * xi;
//  et2_c2 = et * et;
//  q2_c2 = q * q;
//  r2_c2 = xi2_c2 + et2_c2 + q2_c2;
//  if (isEqual(r2_c2,0))
//  {
//    r_c2 = 0;
//    return;
//  }
//
//  r_c2 = sqrt(r2_c2);
//  r3_c2 = r_c2 * r2_c2;
//  r5_c2 = r3_c2 * r2_c2;
//  y_c2 = et * cd + q * sd;
//  d_c2 = et * sd - q * cd;
//
//  if (isEqual(q,0))
//    tt_c2 = 0;
//  else
//    tt_c2 = atan(xi * et / (q * r_c2));
//
//  if (isEqual(kxi,1))
//  {
//    alx_c2 = -log(r_c2 - xi);
//    x11_c2 = 0;
//    x32_c2 = 0;
//  }
//  else
//  {
//    const realT rxi = r_c2 + xi;
//    alx_c2 = log(rxi);
//    x11_c2 = 1 / (r_c2 * rxi);
//    x32_c2 = (r_c2 + rxi) * x11_c2 * x11_c2 / r_c2;
//  }
//
//  if (isEqual(ket,1))
//  {
//    ale_c2 = -log(r_c2 - et);
//    y11_c2 = 0;
//    y32_c2 = 0;
//  }
//  else
//  {
//    const realT ret = r_c2 + et;
//    ale_c2 = log(ret);
//    y11_c2 = 1 / (r_c2 * ret);
//    y32_c2 = (r_c2 + ret) * y11_c2 * y11_c2 / r_c2;
//  }
//
//  ey_c2 = sd / r_c2 - y_c2 * q / r3_c2;
//  ez_c2 = cd / r_c2 + d_c2 * q / r3_c2;
//
//  fy_c2 = d_c2 / r3_c2 + xi2_c2 * y32_c2 * sd;
//  fz_c2 = y_c2 / r3_c2 + xi2_c2 * y32_c2 * cd;
//
//  gy_c2 = 2 * x11_c2 * sd - y_c2 * q * x32_c2;
//  gz_c2 = 2 * x11_c2 * cd + d_c2 * q * x32_c2;
//
//  hy_c2 = d_c2 * q * x32_c2 + xi * q * y32_c2 * sd;
//  hz_c2 = y_c2 * q * x32_c2 + xi * q * y32_c2 * cd;
//
//  return;
//} //OK 5/14/12
//
///*
// C************************************************************           11310005
// C*****   CHECK SINGULARITIES RELATED TO R+XI AND R+ET   *****           11320005
// C************************************************************           11330005
// C                                                                       11340005
// C***** INPUT                                                            11350005
// C*****   X,P,Q : STATION COORDINATE                                     11360005
// C*****   AL1,AL2,AW1,AW2 : FAULT DIMENSIONS                             11370005
// C***** OUTPUT                                                           11380005
// C*****   KXI(2) : =1 IF R+XI<EPS                                        11390005
// C*****   KET(2) : =1 IF R+ET<EPS                                        11400005
// */
//
//void FaultElementManagerT::dccon3(realT x, realT p, realT q, realT al1, realT al2, realT aw1,
//                                       realT aw2, realT *kxi, realT *ket)
//{
//  const realT eps = 1e-6;
//  realT et, xi, r, rxi, ret;
//  int j, k;
//
//  kxi[0] = 0;
//  kxi[1] = 0;
//  ket[0] = 0;
//  ket[1] = 0;
//  if (x > al1 && p > aw2)
//    return;
//
//  for (k = 0; k < 2; k++)
//  {
//    if (k == 0)
//      et = p - aw1;
//    if (k == 1)
//      et = p - aw2;
//    for (j = 0; j < 2; j++)
//    {
//      if (j == 0)
//        xi = x - al1;
//      if (j == 1)
//        xi = x - al2;
//      r = sqrt(xi * xi + et * et + q * q);
//      rxi = r + xi;
//      ret = r + et;
//      if (x < al1 && rxi < eps)
//        kxi[k] = 1;
//      if (p < aw1 && ret < eps)
//        ket[j] = 1;
//    }
//  }
//  return;
//} //OK 5/14/12
