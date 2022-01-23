//FUNCTION_BEGIN_PARSE
realT
RateAndStateIntermediate::CriticalLength(const realT G,
                                         const InterfaceBaseParameterData& matParamsBase,
                                         InterfaceBaseStateData& matStateBase) const
{
  RateAndStateIntermediateStateData& matState = static_cast<RateAndStateIntermediateStateData&> (matStateBase);
  const RateAndStateIntermediateParameterData& matParams = static_cast<const RateAndStateIntermediateParameterData&> (matParamsBase);
  return G * (4.0 * atan(1)) * matParams.Dc / (4 * matState.stress * (matParams.B - matParams.A));
}
