/*
 * BoundaryElementDataManagerT.h
 *
 *  Created on: Apr 25, 2013
 *      Author: johnson346
 */

#ifndef BOUNDARYELEMENTDATAMANAGERT_H_
#define BOUNDARYELEMENTDATAMANAGERT_H_

#include "IO/silo/SiloFile.h"
#include "EarthquakeSimulationCommon.h"

namespace EarthquakeSimulation
{

  class BoundaryElementDataManagerT
  {
  public:
    BoundaryElementDataManagerT();
    virtual ~BoundaryElementDataManagerT();

    inline size_t size() const { return m_KGD.size(); }

    inline Array2dT<R1Tensor>::size_type Dimension( const int dimnum ) const {return m_KGD.Dimension(dimnum);}

    inline const R1Tensor* operator[]( const Array2dT<R1Tensor>::size_type index ) const { return m_KGD[index]; }

    void insert(const localIndex i);

    void erase(const localIndex i);

    void resize(const localIndex dim1);

    void resize(const localIndex dim1, const globalIndex dim2);

    void WriteSilo(SiloFile& siloFile);

    void ReadSilo(const SiloFile& siloFile);

    inline R1Tensor& operator()(const Array2dT<R1Tensor>::size_type dim1 , const Array2dT<R1Tensor>::size_type dim2)
    {
      return m_KGD(dim1, dim2);
    }

    inline const R1Tensor& operator() (const Array2dT<R1Tensor>::size_type dim1 , const Array2dT<R1Tensor>::size_type dim2) const
    {
      return m_KGD(dim1, dim2);
    }

  private:
    Array2dT<R1Tensor> m_KGD;
    size_t m_maxGlobal;
  };

} /* namespace EarthquakeSimulation */
#endif /* BOUNDARYELEMENTDATAMANAGERT_H_ */
