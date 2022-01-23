/*
 * BoundaryElementDataManagerT.cpp
 *
 *  Created on: Apr 25, 2013
 *      Author: johnson346
 */

#include "BoundaryElementDataManagerT.h"

namespace EarthquakeSimulation
{

  BoundaryElementDataManagerT::BoundaryElementDataManagerT()
  {
    // TODO Auto-generated constructor stub
  }

  BoundaryElementDataManagerT::~BoundaryElementDataManagerT()
  {
    // TODO Auto-generated destructor stub
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // CrUD
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  void
  BoundaryElementDataManagerT::insert( const localIndex i )
  {
    R1Tensor t;
    m_KGD.Insert( i, t );
  }

  void
  BoundaryElementDataManagerT::erase( const localIndex i )
  {
    m_KGD.erase(m_KGD.begin() + i);
  }

  void
  BoundaryElementDataManagerT::resize(const localIndex dim1)
  {
    m_KGD.resize(dim1);
  }

  void
  BoundaryElementDataManagerT::resize(const localIndex dim1, const globalIndex dim2)
  {
    m_KGD.resize2(dim1, dim2);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // I/O
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  void
  BoundaryElementDataManagerT::WriteSilo(SiloFile& siloFile)
  {
    if(m_KGD.Dimension(0) == 0)
      return;

    std::string subDirectory =   "FaultElementMatrix";
    if ( DBMkDir(siloFile.m_dbFilePtr, subDirectory.c_str()) != 0)
      throw GPException("Cannot make directory FaultElementMatrix");

    DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str());
    //size the m_KGD array
    {
      Array2dT<R1Tensor>::size_type n0 = m_KGD.Dimension(0), n1 = m_KGD.Dimension(1);
      localIndex nentries = n0;
      siloFile.DBWriteWrapper("n0_KGD", nentries);
      nentries = n1;
      siloFile.DBWriteWrapper("n1_KGD", nentries);
    }
    siloFile.DBWriteWrapper("m_KGD", m_KGD);
    DBSetDir(siloFile.m_dbFilePtr, "..");
  }

  void
  BoundaryElementDataManagerT::ReadSilo(const SiloFile& siloFile)
  {
    std::string subDirectory =   "FaultElementMatrix";
    if ( DBSetDir(siloFile.m_dbFilePtr, subDirectory.c_str()) == 0)
    {
      //size the m_KGD array
      {
        Array2dT<R1Tensor>::size_type n0, n1;
        localIndex nentries;
        siloFile.DBReadWrapper("n0_KGD", nentries);
        n0 = nentries;
        siloFile.DBReadWrapper("n1_KGD", nentries);
        n1 = nentries;
        m_KGD.resize2(n0, n1);
      }

      siloFile.DBReadWrapper("m_KGD", m_KGD);
      DBSetDir(siloFile.m_dbFilePtr, "..");
    }
  }
} /* namespace EarthquakeSimulation */
