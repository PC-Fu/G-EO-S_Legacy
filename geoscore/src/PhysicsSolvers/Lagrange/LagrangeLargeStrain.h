/*
 * LagrangeLargeStrain.h
 *
 *  Created on: Nov 8, 2012
 *      Author: settgast1
 */

#ifndef LAGRANGELARGESTRAIN_H_
#define LAGRANGELARGESTRAIN_H_

#include "LagrangeSolverBase.h"

class LagrangeLargeStrain: public LagrangeSolverBase
{
public:
  LagrangeLargeStrain(  const std::string& name,
                        ProblemManagerT* const pm );
  virtual ~LagrangeLargeStrain();


  static const char* SolverName()
  {
    return "LagrangeLargeStrain";
  };

  void ApplyGravity( NodeManagerT& nodeManager,
                     ElementRegionT& elemRegion,
                     const realT dt );

private:
  LagrangeLargeStrain();


  void ProcessElementRegion( NodeManagerT& nodeManager,
                             ElementRegionT& elemRegion,
                             const realT dt );

  virtual realT CalculateElementResidualAndDerivative( const MaterialBaseParameterData& matParams,
                                                const FiniteElementBase& fe,
                                                const Array2dT<R1Tensor>& dNdX,
                                                const realT* const detJ,
                                                const Epetra_SerialDenseVector& dof_np1,
                                                Epetra_SerialDenseMatrix& dRdU,
                                                Epetra_SerialDenseVector& R ){ return 0;}






};

#endif /* LAGRANGELARGESTRAIN_H_ */
