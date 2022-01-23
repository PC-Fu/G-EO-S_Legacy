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
 * File: FileManagerDataT.h
 * Class provides file IO data structure
 * created : SJ (11/16/2013)
 * ported from work by RRS (10/2001)
 */
 
#ifndef _FILEMANAGERDATAT_H_
#define _FILEMANAGERDATAT_H_

#include <fstream>
#include <string>
#include <stdio.h>
#include "Common/Common.h"

namespace GPAC_IO
{
  /**
   * @author Scott Johnson
   * @brief Class to manager the data associated with the Abaqus geometry file reader
   */
  class FileManagerDataT
  {

  public:

    FileManagerDataT();
    ~FileManagerDataT(){}

    //file operations
    bool OpenFile();
    bool OK() {return !geometry.eof();}
    void CloseFile() {geometry.close();}
    void Write();

    //manage state
    virtual void Reset();

    virtual bool AdvanceLine(std::string& inputline);
    virtual void AddNodalPositionLine(std::string&, const realT){}
    virtual void AddElementRegionLine(std::string&){}

    static std::string ExtractValue(const std::string& key, const std::string& inputline);

    bool exist; // Does the manager exist - set ReadMesh
    localIndex numNodes; // Number of nodes - init (A) set (B)
    globalIndex maxGlobalNodeID; // maximum nodal ID - set (A)
    globalIndex maxGlobalElemID; // maximum element ID
    localIndex numElementRegions;
    localIndex numNodeSets; // Number of nodes - init (A) set (B)
    std::map<std::string, localIndex> numElementsInRegion; // Number of elements in a region
    std::map<std::string, std::string> elementRegionNameTypes; // Map of element region names to their types
    sArray1d elementRegionNames;
    sArray1d elementRegionTypes;
    gArray2d elemToNodes;

    Array1dT<R1Tensor> nodalPositions;
    std::map<globalIndex,R1Tensor> nodalPositionsMap;
    std::map<globalIndex,R1Tensor> potentialNodalPositionsMap;
    std::map<globalIndex,gArray1d> elemToNodesMap;


    R1Tensor spatialMin;
    R1Tensor spatialMax;
    // these arrays indicate if the objects status in the computational domain. The values are:
    //  0 = Not in the domain
    //  1 = In the domain
    iArray1d isElemInDomain;
    iArray1d isNodeInDomain;
    std::map<std::string, gArray1d> elemsInRegion;
    std::map<std::string, globalIndex> maxNumElementsInRegion;
    std::streampos beginElems;
    std::streampos beginNodes;
    std::streampos beginNodesets;
    std::streampos current;
    lvector GlobalToLocalNodeMap;
    std::ifstream geometry;      /** geometry file stream */
    std::string filename;

    globalIndex totalNumElem; //total number of elements
    gSet neighborList;

  };
}
#endif //_FILEMANAGERDATAT_H_
