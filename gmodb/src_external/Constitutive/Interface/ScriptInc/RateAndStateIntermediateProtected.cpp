//FUNCTION_BEGIN_PARSE
virtual_void
RateAndStateIntermediate::UpdateFriction( const InterfaceBaseParameterData& matParamsBase,
                                          InterfaceBaseStateData& matStateBase) const
{
  RateAndStateIntermediateStateData& matState = static_cast<RateAndStateIntermediateStateData&> (matStateBase);
  const RateAndStateIntermediateParameterData& matParams = static_cast<const RateAndStateIntermediateParameterData&> (matParamsBase);

  matState.mu = matParams.mu0 + matParams.A * log( 1.0 + matState.dxsdt / matParams.vstar) +
      matParams.B * log(1.0 + matState.theta / matParams.thetastar);
  const realT dthetadt = 1.0 - matState.theta * matState.dxsdt / matParams.Dc -
      (matParams.alpha * matState.theta * matState.dstressdt) / (matParams.B * matState.stress);
  matState.theta += dthetadt * matState.dt;
}

//FUNCTION_BEGIN_PARSE
virtual_realT
RateAndStateIntermediate::ShearStrength(const InterfaceBaseParameterData& ,
                                        InterfaceBaseStateData& matStateBase) const
{
  RateAndStateIntermediateStateData& matState = static_cast<RateAndStateIntermediateStateData&> (matStateBase);
  return matState.stress * matState.mu;
}
