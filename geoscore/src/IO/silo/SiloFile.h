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
 * @file SiloFile.h
 * @author settgast1
 * @date Jan 24, 2011
 */

#ifndef SILOFILE_H_
#define SILOFILE_H_

#include "Common/Common.h"
#include "silo.h"
#include <vector>

#if GPAC_MPI
#include <mpi.h>
#endif

#include "pmpio.h"


/////////////////////////////////////////////////


#include "DataStructures/InterObjectRelation.h"

/// typedef my callback function for use with WriteMultiXXXX. This function is essentially a placeholder
/// for DB,PutMultimesh, DBPutMultivar, etc.
typedef struct _PMPIO_baton_t PMPIO_baton_t;
//typedef int (*DBPutMultimeshType)(DBfile *, char DB_CONSTARR1, int, char DB_CONSTARR2, int DB_CONSTARR1, DBoptlist const *);
//typedef int (*DBPutMultivarType)(DBfile *, const char *, int, char **, int *, DBoptlist *);
class ElementManagerT;

// *********************************************************************************************************************
// *********************************************************************************************************************
/**
 * @author settgast
 * This class serves as a wrapper to isolate the code from the specifics of SILO output/input. Its members contain all
 * the necessary information for reading/writing a group of SILO files.
 */
class SiloFile
{

public:

  /// Default Constructor
  SiloFile();

  /// Destructor
  virtual ~SiloFile();

  /// Initializes the silo library
  void Initialize( const PMPIO_iomode_t readwrite );

  /// finishes up the silo library usage
  void Finish();

  /// Wait for the Baton when doing PMPIO
  void WaitForBaton( const int domainNumber, const int cycleNum, const bool isRestart );

  void WaitForBaton( const int domainNumber, const std::string& restartFileName );

  /// Hand off the Baton when doing PMPIO
  void HandOffBaton();


  void MakeSubDirectory( const std::string& subdir, const std::string& rootdir )
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
      DBMkDir(m_dbFilePtr, rootdir.c_str());
    }
    DBMkDir(m_dbFilePtr, subdir.c_str());
  }

  /// Write out a single mesh object
  void WriteMeshObject(const std::string& meshName,
                       const localIndex nnodes,
                       realT* coords[3],
                       const globalIndex*,
                       const int numRegions,
                       const int* shapecnt,
                       const localIndex* const * const meshConnectivity,
                       const globalIndex* const * const globalElementNum,
                       const int* const * const,
                       const int* const shapetype,
                       const int* const shapesize,
                       const int cycleNumber,
                       const realT problemTime);

  /**
   * @brief Write out a single discrete element mesh object
   * @author Scott Johnson
   */
  void WriteDiscreteElementMeshObject(const std::string& meshName,
				      const localIndex nnodes,
				      realT* coords[3],
				      const globalIndex* globalNodeNum,
				      const int nDiscreteElements,
				      const int nfaces,
				      int* nodecnts,
				      const int sumnodecnts,
				      int* nodelist,
				      int* facecnts,
				      const int sumfacecnts,
				      int* facelist,
				      const globalIndex* const * const globalElementNum,
				      const int* ghostFlag,
				      const int cycleNumber,
				      const realT problemTime);
  
  /**
   * @brief Write out a single ellipsoidal discrete element mesh object
   * @author Scott Johnson
   */
  void WriteDiscreteElementCSGObject(const std::string& meshName,
                                     const Array1dT<R1Tensor>& x,
                                     const Array1dT<R1Tensor>& r,
                                     const Array1dT<R2Tensor>& R,
                                     const int cycleNumber,
                                     const realT problemTime);

  void WritePointMesh( const std::string& meshName,
                       const localIndex numPoints,
                       realT* coords[3],
                       const int cycleNumber,
                       const realT problemTime );

  void WriteBeamMesh(const std::string& meshName,
                     const localIndex nnodes,
                     realT* coords[3],
                     const lArray1d& node1,
                     const lArray1d& node2,
                     const int cycleNumber,
                     const realT problemTime);

  void WriteBeamMesh(const std::string& meshName,
                     const localIndex nnodes,
                     realT* coords[3],
                     const std::map<int, int>& connectivity,
                     const int cycleNumber,
                     const realT problemTime);

  void WriteBeamMesh(const std::string& meshName,
                     const localIndex nnodes,
                     realT* coords[3],
                     iArray1d& nodelist,
                     const int cycleNumber,
                     const realT problemTime);

  /**
   * @brief Write out a single "quad" (i.e. hex) mesh object, and record the mesh dimensions for subsequent write data field calls.
   * @author walsh24
   */
  int  WriteQuadMeshObject(const std::string& meshName,
                            const localIndex nX,
                            const localIndex nY,
                            const localIndex nZ,
                            realT* coords[3],
                            const int cycleNumber,
                            const realT problemTime);

//  void TestWriteDiscreteElementMeshObject();

  void WriteRegionSpecifications(const ElementManagerT& elementManager,
                                 const std::string& meshName,
                                 const int cycleNumber,
                                 const realT problemTime);

  /// writes out fields in a data member map
  template< typename OUTPUTTYPE, typename TYPE >
  void WriteFieldMapToSilo( const std::string& meshname,
                            const std::map< std::string, TYPE>& member,
                            const int centering,
                            const int cycleNum,
                            const realT problemTime,
                            const bool isRestart,
                            const std::string& multiRoot,
                            const std::string& regionName,
                            const lArray1d& mask );

  template< typename INPUTTYPE, typename TYPE >
  void ReadFieldMapFromSilo( std::map< std::string, Array1dT<TYPE> >& member,
                             const std::string& meshname,
                             const int centering,
                             const int cycleNum,
                             const realT problemTime,
                             const bool isRestart,
                             const std::string& regionName,
                             const lArray1d& mask ) const;

  /// Write out a data field
  template<typename OUTTYPE, typename TYPE>
  void WriteDataField( const std::string& meshName,
                       const std::string& fieldName,
                       const Array1dT<TYPE>& field,
                       const int centering,
                       const int cycleNumber,
                       const realT problemTime,
                       const std::string& multiRoot,
                       const std::string& regionName );

  /// Write out a data field
  template< typename INPUTTYPE, typename TYPE>
  void ReadDataField( Array1dT<TYPE>& field,
                      const std::string& meshName,
                      const std::string& fieldName,
                      const int centering,
                      const int cycleNumber,
                      const realT problemTime,
                      const std::string& regionName ) const;

  int GetMeshType( const std::string& meshName ) const
  {
    int meshType = -1;
    {
      // in order to get mesh type, we might have to go up a few levels in the silo directory structure
      // before we can find the mesh.
      char pwd[256];
      DBGetDir(m_dbFilePtr, pwd);
      for( int i=0 ; i<3 ; ++i )
      {
        meshType = DBInqMeshtype(m_dbFilePtr,meshName.c_str());
        if( meshType != -1 && meshType != 610 )
          break;
        else
          DBSetDir(m_dbFilePtr,"..");
      }
      DBSetDir(m_dbFilePtr,pwd);
    }
    return meshType;
  }

  template <typename TYPE>
  void** GetDataVar( const std::string& fieldName,
                     const std::string& meshName,
                     const typename Array1dT<TYPE>::size_type nels,
                     const int centering,
                     const int cycleNumber,
                     const realT problemTime,
                     const std::string& ) const;

  /// Write out a multi-mesh object
  template< typename CBF >
  void WriteMultiXXXX(const DBObjectType type, CBF DBPutMultiCB,
                      const int centering, const std::string name, const int cycleNumber,
                      const std::string& multiRoot, const DBoptlist* optlist = NULL);


  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const TYPE& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const Array1dT<TYPE>& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const std::set<TYPE>& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const Array2dT<TYPE>& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const Array1dT<Array1dT<TYPE> >& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const Array1dT<Array2dT<TYPE> >& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const Array1dT<Array1dT<Array1dT<TYPE> > >& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const Array1dT<std::set<TYPE> >& data );

  template< typename T1, typename T2 >
  void DBWriteWrapper( const std::string& name, const std::map< T1, T2 >& datamap );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& subdir, const std::map< std::string, TYPE>& member );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& subdir, const std::map< std::string, InterObjectRelation<TYPE> >& member );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const TYPE* const data, const int size );





  template<typename TYPE>
  void DBReadWrapper( const std::string& name, TYPE& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, Array1dT<TYPE>& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, std::set<TYPE>& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, Array2dT<TYPE>& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, Array1dT<Array1dT<TYPE> >& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, Array1dT<Array2dT<TYPE> >& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, Array1dT<Array1dT<Array1dT<TYPE> > >& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, Array1dT<std::set<TYPE> >& data ) const;

  template< typename T1, typename T2 >
  void DBReadWrapper( const std::string& name, std::map< T1, T2 >& datamap ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& subdir, std::map< std::string, TYPE>& member ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& subdir, std::map< std::string, InterObjectRelation<TYPE> >& member ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, TYPE* const data, const int size ) const;








  /// dummy function to get rid of compiler warnings about unused functions from PMPIO's use of static function
  /// definitions in the header. This should go away as we fully utilize PMPIO.
  void StopSiloCompilerWarnings();

  void ClearEmptiesFromMultiObjects(const int cycleNum);

  //private:

  /// pointer to the DBfile that this class is working on
  DBfile* m_dbFilePtr;

  /// total number of "parallel" files to write out
  int m_numGroups;

  /// the pmpio baton. A processor needs this to write to the file.
  PMPIO_baton_t *m_baton;

  /// which database to use. DB_PDB or DB_HDF5
  int m_driver;

  /// root of the filename that we will be reading/writing
  std::string m_fileRoot;

  std::string m_restartFileRoot;

  std::string m_slaveDirectory;

  std::string m_fileName;

  std::string m_baseFileName;

  bool m_markGhosts;

  sArray1d m_emptyMeshes;
  sArray1d m_emptyVariables;


private:
  int m_quadMeshDims[3];
  int m_quadMeshNDims;

};

// *********************************************************************************************************************
// *********************************************************************************************************************
/**
 * Namespace to hold some utilities needed by the functions in SiloFile.
 */
namespace SiloFileUtilities
{
  /**
   * @author settgast
   * @tparam OUTTYPE the type of data to write out (int,float,realT)
   * @return the integer identifier associated with the enum DBdatatype defined in
   * the silo.h file.
   *
   * This templated function is a "specialization only" definition. There is no general
   * definition, only specializations for predetermined data types.
   */
  template<typename OUTTYPE>
  int DB_TYPE();

  /**
   * @author settgast
   * @tparam TYPE the data type in question
   * @return the number of "variables" in a data field. For instance, 1 for a int, float or realT.
   * 3 for a R1Tensor...etc.
   */
  template<typename TYPE>
  int GetNumberOfVariablesInField();

  template<typename TYPE>
  int GetTensorRank();

  template<typename OUTTYPE, typename TYPE>
    OUTTYPE CastField(const TYPE& field, const int i = 0); // avoids compiler warning

  template<typename OUTTYPE, typename TYPE>
  OUTTYPE CastField(const TYPE& field, const int i)
  {
    return static_cast<OUTTYPE>(field.Data()[i]);
  }

  template<> inline int CastField<int, int> (const int& field, const int )
  {
    return field;
  }


  template<> inline int CastField<int, localIndex> (const localIndex& field, const int )
  {
    return static_cast<int>(field);
  }

  template<> inline localIndex CastField<localIndex, localIndex> (const localIndex& field, const int )
  {
    return field;
  }


  template<> inline globalIndex CastField<globalIndex, globalIndex> (const globalIndex& field, const int )
  {
    return field;
  }

  template<> inline int CastField<int, long long unsigned int > (const long long unsigned int& field, const int )
  {
    return static_cast<int>(field);
  }


  template<> inline realT CastField<realT, realT> (const realT& field, const int )
  {
    return field;
  }
  template<> inline float CastField<float, realT> (const realT& field, const int )
  {
    return static_cast<float> (field);
  }


  template< typename TYPE, typename PTR_TYPE> inline PTR_TYPE* DataPtr( TYPE& data)
  {
    return (PTR_TYPE*)&data;
  }

  template<> inline realT* DataPtr( R1Tensor& data)        {    return data.begin();  }
  template<> inline realT* DataPtr( R2Tensor& data)        {    return data.begin();  }
  template<> inline realT* DataPtr( R2SymTensor& data)  {    return data.begin();  }



  template<typename TYPE>
  void SetVariableNames(const std::string& fieldName, sArray1d& varnamestring, char* varnames[]);

  template<typename OBJECT_TYPE>
  int FieldCentering();

  void SetCenteringSubdir(const int centering, std::string& subdir);

}

// *********************************************************************************************************************
// *********************************************************************************************************************

/**
 * @author Randy Settgast, Scott Johnson
 * @tparam OUTPUTTYPE type of data write to the output file
 * @tparam T native type of data being output
 * @param[in] siloFile file handle
 * @param[in] member member that we are writing
 * @param[in] centering the centering location of the data (i.e. Node, Element, Face)
 * @param[in] cycleNum the the cycle number
 * @param[in] problemTime the problem time
 * @param[in] mask list of indices of entries to plot
 *
 * This function writes all fields in a member map to a silo file.
 */
template< typename OUTPUTTYPE, typename T >
void SiloFile::WriteFieldMapToSilo( const std::string& meshname,
                                    const std::map< std::string, T>& member,
                                    const int centering,
                                    const int cycleNum,
                                    const realT problemTime,
                                    const bool isRestart,
                                    const std::string& multiRoot,
                                    const std::string& regionName,
                                    const lArray1d& mask )
{

  // iterate over all entries in the member map
  for( typename std::map< std::string, T >::const_iterator iter = member.begin() ; iter!=member.end() ; ++iter )
  {
    // the field name is the key to the map
    const std::string fieldName = iter->first;

    // check to see if the field should be written
    if( FieldInfo::AttributesByName.find(fieldName) != FieldInfo::AttributesByName.end() )
    {
      if( (  isRestart && FieldInfo::AttributesByName[fieldName]->m_WriteToRestart) ||
          ( !isRestart && FieldInfo::AttributesByName[fieldName]->m_WriteToPlot ) )
      {
        // the field data is mapped value
        const T& fieldData = iter->second;

        if( !(mask.empty()) && !isRestart )
        {
          T dataToWrite(mask.size());
          for( lArray1d::size_type i = 0; i < mask.size(); ++i)
            dataToWrite[i] = fieldData[mask[i]];

          // write the data field
          WriteDataField<OUTPUTTYPE>(meshname.c_str(), fieldName, dataToWrite, centering, cycleNum, problemTime, multiRoot, regionName );
        }
        else
        {
          WriteDataField<OUTPUTTYPE>(meshname.c_str(), fieldName, fieldData, centering, cycleNum, problemTime, multiRoot, regionName );
        }
      }
    }
  }
}


/**
 *
 * @param meshName
 * @param fieldName
 * @param field
 * @param centering
 * @param cycleNumber
 * @param problemTime
 */
template<typename OUTTYPE, typename TYPE>
void SiloFile::WriteDataField( const std::string& meshName,
                               const std::string& fieldName,
                               const Array1dT<TYPE>& field,
                               const int centering,
                               const int cycleNumber,
                               const realT problemTime,
                               const std::string& multiRoot,
                               const std::string& regionName )
{
  const int nvars = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();
  const int nels = field.size();

  const int meshType = GetMeshType( meshName );


  DBoptlist *optlist = DBMakeOptlist(5);
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

  char *regionpnames[2] =
  { NULL, NULL };




  if (regionName != "none")
  {
    regionpnames[0] = const_cast<char*> (regionName.c_str());
    regionpnames[1] = NULL;
    DBAddOption(optlist, DBOPT_REGION_PNAMES, &regionpnames);
  }

  // if the number of elements is zero, then record the path to the var. This will be used later to delete the entry
  // from the multivar.
  if (nels == 0)
  {
    char pwd[256];
    DBGetDir(m_dbFilePtr, pwd);
    std::string emptyObject = pwd;
    emptyObject += "/" + fieldName;
    m_emptyVariables.push_back(emptyObject);
  }
  else
  {
    Array1dT<char*> varnames(nvars);
    Array1dT<void*> vars(nvars);


    sArray1d varnamestring(nvars);
    std::vector<std::vector<OUTTYPE> > castedField(nvars);

    SiloFileUtilities::SetVariableNames<TYPE>(fieldName, varnamestring, varnames.data() );

    for (int i = 0; i < nvars; ++i)
    {
      castedField[i].resize(nels);
      vars[i] = static_cast<void*> (&(castedField[i][0]));
      for (int k = 0; k < nels; ++k)
      {
        castedField[i][k] = SiloFileUtilities::CastField<OUTTYPE>(field[k], i);
      }
    }

    int err = -2;
    if( meshType == DB_UCDMESH )
    {
      err = DBPutUcdvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, varnames.data(), (float**) vars.data(),
                         nels, NULL, 0, SiloFileUtilities::DB_TYPE<OUTTYPE>(), centering, optlist);
    }
    else if( meshType == DB_POINTMESH )
    {
      err = DBPutPointvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, (float**) vars.data(),
                           nels, SiloFileUtilities::DB_TYPE<OUTTYPE>(), optlist);
    }
    else if( meshType == DB_QUADCURV )
    {
      err = DBPutQuadvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, varnames.data(), (float**) vars.data(),
                          m_quadMeshDims,m_quadMeshNDims,NULL, 0,  SiloFileUtilities::DB_TYPE<OUTTYPE>() ,centering, optlist);
    }
    if(err < 0)
    {
      if(err < -1)
        throw GPException("unhandled case in SiloFile::WriteDataField A\n");
      else
        throw GPException("unhandled failure in adding variable during SiloFile::WriteDataField\n");
    }

  }

  // write multimesh object
  int rank = 0;
#if GPAC_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0)
  {
    int tensorRank = SiloFileUtilities::GetTensorRank<TYPE>();
    DBAddOption(optlist, DBOPT_TENSOR_RANK, const_cast<int*> (&tensorRank));
    DBAddOption(optlist, DBOPT_MMESH_NAME, const_cast<char*> (meshName.c_str()));

    DBObjectType vartype;

    if( meshType == DB_UCDMESH )
    {
      vartype = DB_UCDVAR;
    }
    else if( meshType == DB_POINTMESH )
    {
      vartype = DB_POINTVAR;
    }
    else if( meshType == DB_QUADCURV )
    {
      vartype = DB_QUADVAR;
    }
    else
    {
      throw GPException("unhandled case in SiloFile::WriteDataField B\n");
    }


    WriteMultiXXXX(vartype, DBPutMultivar, centering, fieldName.c_str(), cycleNumber, multiRoot,
                   optlist);
  }

  DBFreeOptlist(optlist);

}

template< typename INPUTTYPE, typename TYPE >
void SiloFile::ReadFieldMapFromSilo( std::map< std::string, Array1dT<TYPE> >& member,
                                     const std::string& meshname,
                                     const int centering,
                                     const int cycleNum,
                                     const realT problemTime,
                                     const bool isRestart,
                                     const std::string& regionName,
                                     const lArray1d& mask ) const
{
  // iterate over all entries in the member map
  for( typename std::map< std::string, Array1dT<TYPE> >::iterator iter = member.begin() ; iter!=member.end() ; ++iter )
  {
    // the field name is the key to the map
    const std::string fieldName = iter->first;

    // check to see if the field should have been written
    if( FieldInfo::AttributesByName.find(fieldName) != FieldInfo::AttributesByName.end() )
    {
      if( (  isRestart && FieldInfo::AttributesByName[fieldName]->m_WriteToRestart) ||
          ( !isRestart && FieldInfo::AttributesByName[fieldName]->m_WriteToPlot ) )
      {
        // the field data is mapped value
        Array1dT<TYPE> & fieldData = iter->second;

        if( !(mask.empty()) && !isRestart )
        {
          Array1dT<TYPE> dataToRead( mask.size() );
          // write the data field
          ReadDataField<INPUTTYPE>( dataToRead, meshname.c_str(), fieldName, centering, cycleNum, problemTime, regionName );

          for( lArray1d::size_type i = 0; i < mask.size(); ++i)
          {
            fieldData[mask[i]] = dataToRead[i];
          }

        }
        else
        {
          ReadDataField<INPUTTYPE>( fieldData, meshname.c_str(), fieldName, centering, cycleNum, problemTime, regionName );
        }
      }
    }
  }

}

template< typename INPUTTYPE, typename TYPE>
void SiloFile::ReadDataField( Array1dT<TYPE>& field,
                              const std::string& meshName,
                              const std::string& fieldName,
                              const int centering,
                              const int cycleNumber,
                              const realT problemTime,
                              const std::string& regionName ) const
{

  INPUTTYPE** var = (INPUTTYPE**) GetDataVar<TYPE>( fieldName, meshName, field.size(), centering, cycleNumber, problemTime, regionName );

  for( typename Array1dT<TYPE>::size_type a=0 ; a<field.size() ; ++a )
  {
    TYPE temp;
    INPUTTYPE* ptemp = SiloFileUtilities::DataPtr<TYPE,INPUTTYPE>( temp );

    for( int i=0 ; i<SiloFileUtilities::GetNumberOfVariablesInField<TYPE>() ; ++i  )
    {
      ptemp[i] = var[i][a];
    }
    field[a] = temp;
  }

}

/**
 *
 * @param type
 * @param DBPutMultiCB
 * @param centering
 * @param name
 * @param cycleNumber
 * @param optlist
 */
template< typename CBF >
void SiloFile::WriteMultiXXXX( const DBObjectType type,
                               CBF DBPutMultiCB,
                               const int centering ,
                               const std::string name,
                               const int,
                               const std::string& multiRoot,
                               const DBoptlist* optlist)
{
  int size = 1;
#if GPAC_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  sArray1d vBlockNames(size);
  std::vector<char*> BlockNames(size);
  ivector blockTypes(size);
  char tempBuffer[1024];
  char currentDirectory[256];

  DBGetDir(m_dbFilePtr, currentDirectory);
  DBSetDir(m_dbFilePtr, multiRoot.c_str());


  std::string multiRootString(multiRoot);
  if( !(multiRootString.compare("/")) )
  {
    multiRootString.clear();
  }

  for (int i = 0; i < size; ++i)
  {
    int groupRank = PMPIO_GroupRank(m_baton, i);

    if (groupRank == 0)
    {

      sprintf(tempBuffer, "/domain_%04d%s/%s", i, multiRootString.c_str(), name.c_str());
      //      sprintf(tempBuffer, "%s_%04d:/domain_%04d%s/%s", m_fileRoot.c_str(),
      //              cycleNumber, i, multiRoot.c_str(), name.c_str() );

    }
    else
    {
      if (m_slaveDirectory.empty())
      {
        sprintf(tempBuffer, "%s.%03d:/domain_%04d%s/%s", m_baseFileName.c_str(),
                groupRank, i, multiRootString.c_str(), name.c_str());
      }
      else
      {
        sprintf(tempBuffer, "%s%s%s.%03d:/domain_%04d%s/%s", m_slaveDirectory.c_str(), "/", m_baseFileName.c_str(),
                groupRank, i, multiRootString.c_str(), name.c_str());

      }

    }
    vBlockNames[i] = tempBuffer;
    BlockNames[i] = (char*) vBlockNames[i].c_str();
    blockTypes[i] = type;
  }

  std::string multiName = name;
  DBPutMultiCB(m_dbFilePtr, multiName.c_str(), size, BlockNames.data(), blockTypes.data(),
               const_cast<DBoptlist*> (optlist));

  DBSetDir(m_dbFilePtr, currentDirectory);

}

#endif /* SILOFILE_H_ */
