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
 * File: FileManagerT.h
 * Class provides file IO
 * created : RRS (10/11/2001)
 */
 
#ifndef _FILE_MANAGER_T_H_
#define _FILE_MANAGER_T_H_

// ***** Included Headers ******************************************************
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>
#include "ticpp/TinyXMLParser.h"
#include "Common/Common.h"
#include "AbaqusFileManagerDataT.h"
#include "EllipsoidFileManagerDataT.h"
#include <metis.h>
#include <vector>

class SpatialPartition;
class PartitionBase;

// ***** Forward Declarations **************************************************
class PhysicalDomainT;
class NodeManagerT;
class ElementManagerT;
class DiscreteElementManagerT;

// ***** Class Declaration *****************************************************
namespace GPAC_IO
{
  class FileManagerT
  {
    public:
    /**
    * Constructor sets refences
    * @param pman reference to problem manager
    */
    FileManagerT();

    /** Destructor */
    ~FileManagerT(void);

    /**
    * Set fileroot std::string
    * @param root file root
    */
    void SetRoot(const char* const root);

    /**
    * Set input filename std::string
    * @param filename input filename
    */
    void SetInputFilename(const char* const filename);

    std::string GetInputFilename(void) const{return outputfilename;};

    /**
    * Set geometry filename std::string
    * @param filename geometry filename
    */
    void SetGeometryFilename(const char* const filename);

    /** Set geometry units - later used to convert mesh units into problem units **/
    void SetGeometryUnits(const realT units){geometryUnits=units; };

    /** Set mpi mesh reading message size (size of each message block) **/
    void SetMessageSize(const globalIndex messageSize){geometryMessageSize=messageSize; };

    /**
    * Set discrete element geometry filename std::string
    * @param filename geometry filename
    */
    void SetDiscreteElementGeometryFilename(const char* const filename);

    /**
    * Set ellipsoidal discrete element geometry filename std::string
    * @param filename geometry filename
    */
    void SetEllipsoidalDiscreteElementGeometryFilename(const char* const filename);

#ifdef SRC_EXTERNAL
    /**
    * Set fault patch element geometry filename std::string
    * @param filename geometry filename
    */
    void SetFaultPatchElementGeometryFilename(const char* const filename);
#endif
  
    std::string GetGeometryFilename(void) const{return geometryfilename;};

    /** Read Input data from an Ascii file */
    void ReadAsciiInput();
  
    /** Read Input data from an XML file */
    bool ReadXMLInput(TICPP::HierarchicalDataNode& hdn) const;
    bool ReadXMLInput(const char* const filename, TICPP::HierarchicalDataNode& hdn) const;
  
    /** Determine mesh data from an XML node */
    bool ReadMeshXML(TICPP::HierarchicalDataNode* meshNode);
  
    /** Read Geometry data from an Ascii file */
    void ReadAsciiGeometry( PhysicalDomainT& domain );
  
    bool ExtractElementRegionFromAbaqus( PhysicalDomainT& domain, const PartitionBase& partition );
  
   private:
  #if 1
    static unsigned int ElementTypeToNumberOfNodes(const std::string& elementType);
  //  static unsigned int ElementTypeToNumberOfNodesPerFace(const std::string& elementType);
  //  static unsigned int ElementTypeToNumberOfFaces(const std::string& elementType);
   // static unsigned int ElementTypeToNumberOfDimensions(const std::string& elementType);
  #endif
    void ReadDiscreteAbaqusMeshA( AbaqusFileManagerDataT& fd);
    virtual void ReadAbaqusMeshA( AbaqusFileManagerDataT& fd);
    virtual void ReadAbaqusMeshA( AbaqusFileManagerDataT& fd, const globalIndex mpiNodeLimit, globalIndex &countGlobalNodeNum, std::map<globalIndex,R1Tensor> &tempNodalPositionsMap);
    virtual void ReadAbaqusMeshB( SpatialPartition& partition, AbaqusFileManagerDataT& fd);
    //virtual void ReadAbaqusMeshB( SpatialPartition& partition, AbaqusFileManagerDataT& fd, const globalIndex mpiElemLimit);
    virtual void ReadAbaqusNodeRegion(AbaqusFileManagerDataT& femData,
                                      std::map<globalIndex,R1Tensor> &tempNodalPositionsMap,
                                      std::map<globalIndex,globalIndex> &resortNodeMap);
    virtual void ReadAbaqusElement( PhysicalDomainT& domain , SpatialPartition& partition,
                                    AbaqusFileManagerDataT& femData, std::map<globalIndex,globalIndex> &resortNodeMap,
                                    std::map<globalIndex,R1Tensor> &tempNodalPositionsMap);
    virtual void ReadAbaqusNodeSet( PhysicalDomainT& domain , SpatialPartition& partition,
                                    AbaqusFileManagerDataT& femData, std::map<globalIndex,globalIndex> &resortNodeMap);
    void SetNumNodes(FileManagerDataT& fd);
    void SetNumNodesDiscrete(FileManagerDataT& fd);


    void ReadAbaqusMeshDEM_PartitionDEM( SpatialPartition& partition, AbaqusFileManagerDataT& fd, std::map<std::string, gSet>& dems);
    void ReadDiscreteElementAbaqusMesh(PhysicalDomainT& domain, SpatialPartition& partition, AbaqusFileManagerDataT& fd);
    void ReadFiniteElementAbaqusMesh( PhysicalDomainT& domain , SpatialPartition& partition, AbaqusFileManagerDataT& fd );
#ifdef SRC_EXTERNAL
    void ReadFaultPatchElementAbaqusMesh( PhysicalDomainT& domain , SpatialPartition& partition, AbaqusFileManagerDataT& fd );
#endif
    void ReadEllipsoidFileA( EllipsoidFileManagerDataT& fd);
    void ReadEllipsoidFileB( SpatialPartition& partition, EllipsoidFileManagerDataT& fd);
    void ReadEllipsoidFile( PhysicalDomainT& domain, SpatialPartition& partition, EllipsoidFileManagerDataT& fd);
  
   public:
    void ReadMesh(
            PhysicalDomainT& domain,
            SpatialPartition& partition );

    void ReadMeshforMetis(PhysicalDomainT& domain, SpatialPartition& partition );


    void WriteOutputFile();
  
    void WriteEnsight( const int cycle, const realT time, const PhysicalDomainT& domain );

    void WriteRestartFile( const realT time , const realT ensight_output_time ,
                           const realT th_output_time ,const int n , const int iset ,
                           const realT TotalEnergy,
                           const realT ElasticStrainEnergy,
                           const realT PlasticDissipationEnergy,
                           const realT KineticEnergy,
                           const realT SeparationEnergy,
                           const realT NodalSeparationEnergy );
  
    void ReadRestartFile( realT& time , realT& ensight_output_time ,
                          realT& th_output_time , int& n , int& iset ,
                          realT& TotalEnergy,
                          realT& ElasticStrainEnergy,
                          realT& PlasticDissipationEnergy,
                          realT& KineticEnergy,
                          realT& SeparationEnergy ,
                          realT& NodalSeparationEnergy);

    void SetTemp(const int filenum);
    
  
    //***** Data Member Declarations **********************************************
    protected:
    R2SymTensorT<nsdof> tempMat;
    //std::string fileroot;        /** root for all input/output files*/
    std::string inputfilename;    /** std::string containing name of the input file */

    std::string ensightoutputfilename;
    std::string restartfilename;

    std::string outputfilename;    /** std::string containing name of the output file */
    std::string geometryfilename;  /** std::string containing name of the geometry file */
    std::string degeometryfilename;  /** std::string containing name of the discrete element geometry file */
    std::string edegeometryfilename;  /** std::string containing name of the ellipsoidal discrete element geometry file */
    std::string fpgeometryfilename;  /** std::string containing name of the fault patch element geometry file */

    std::ifstream input;        /** input file stream */
    std::ofstream output;      /** output file stream */


    int ensight_sequence_num;

    std::string temp;

    int output_precision;
    int restart_flag;

    realT geometryUnits; /** scale of the mesh geometry **/
    globalIndex geometryMessageSize; /** size of message for mpi mesh reading**/

    static int s_var;
  
  public:
    int Restart_Flag(void)
    {  return restart_flag; }
  };
}
#endif
