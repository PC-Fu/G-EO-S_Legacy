#ifndef EARTHQUAKESIMULATIONCOMMON_H_
#define EARTHQUAKESIMULATIONCOMMON_H_

namespace EarthquakeSimulation
{
  enum TransitionState
  {
    LOCK = 0,
    NUCLEATE = 1,
    RUPTURE = 2,
    CREEP = 3,
    LOW_SIGMA = 4,
    LOW_TAU = 5,
    HIGH_THETA = 6,
    SLOWSLIP_2A = 7,
    SLOWSLIP_2B = 8,
    SLOWSLIP_2C = 9,
    POREPRESSURERATECHANGE = 10
  };

  enum TransitionStateConsequence
  {
    LOCK_TO_NUCLEATE = 0,
    NUCLEATE_TO_LOCK = 1,
    NUCLEATE_TO_RUPTURE = 2,
    NUCLEATE_TO_SLOWSLIP_2A = 3,
    RUPTURE_TO_LOCK = 4,
    CREEP_TO_CREEP = 5,
    SLOWSLIP_2A_TO_LOCK = 6,
    SLOWSLIP_2A_TO_SLOWSLIP_2B = 7,
    SLOWSLIP_2B_TO_SLOWSLIP_2A = 8,
    SLOWSLIP_2B_TO_SLOWSLIP_2B = 9,
    SLOWSLIP_2B_TO_SLOWSLIP_2C = 10,
    SLOWSLIP_2C_TO_SLOWSLIP_2B = 11,
    ANY_TO_POREPRESSURERATECHANGE = 12
  };

  enum FailureState
  {
    NOTYET = 0,
    NOW = 1,
    ALREADY = 2
  };

  static TransitionState TransitionStateFromIndex(const int index)
  {
    switch (index)
    {
      case LOCK:
        return LOCK;
      case NUCLEATE:
        return NUCLEATE;
      case RUPTURE:
        return RUPTURE;
      case CREEP:
        return CREEP;
      case LOW_SIGMA:
        return LOW_SIGMA;
      case LOW_TAU:
        return LOW_TAU;
      case HIGH_THETA:
        return HIGH_THETA;
      case SLOWSLIP_2A:
        return SLOWSLIP_2A;
      case SLOWSLIP_2B:
        return SLOWSLIP_2B;
      case SLOWSLIP_2C:
        return SLOWSLIP_2C;
      default:
        throw GPException("Invalid TransitionState");
    }
  }

  static std::string TransitionStateToString(const TransitionState state)
  {
    switch (state)
    {
      case LOCK:
        return "LOCK";
      case NUCLEATE:
        return "NUCLEATE";
      case RUPTURE:
        return "RUPTURE";
      case CREEP:
        return "CREEP";
      case LOW_SIGMA:
        return "LOW_SIGMA";
      case LOW_TAU:
        return "LOW_TAU";
      case HIGH_THETA:
        return "HIGH_THETA";
      case SLOWSLIP_2A:
        return "SLOWSLIP_2A";
      case SLOWSLIP_2B:
        return "SLOWSLIP_2B";
      case SLOWSLIP_2C:
        return "SLOWSLIP_2C";
      case POREPRESSURERATECHANGE:
      default:
        throw GPException("Invalid TransitionState");
    }
  }

  static FailureState FailureStateFromIndex(const int index)
  {
    switch (index)
    {
      case NOTYET:
        return NOTYET;
      case NOW:
        return NOW;
      case ALREADY:
        return ALREADY;
      default:
        throw GPException("Invalid FailureState");
    }
  }


  struct EarthquakeSimulationTimestep
  {
    localIndex m_local;
    realT m_dt, m_currentTime;
    TransitionState m_next;
    bool m_isPorePressureRateChange;

    bool CheckForTimestepControl(const realT dt, const localIndex local, const TransitionState next)
    {
      if(dt < 0)
        throw GPException("EarthquakeSimulationTimestep::CheckForTimestepControl timestep is negative - exiting");
      if (dt >= m_dt)
        return false;

      m_dt = dt;
      m_local = local;
      m_isPorePressureRateChange = next == POREPRESSURERATECHANGE;
      if(!m_isPorePressureRateChange)
        m_next = next;
      return true;
    }

    void Reset()
    {
      m_local = 0;
      m_dt = std::numeric_limits<realT>::max();
      m_next = LOCK;
    }

    void Initialize()
    {
      Reset();
      m_currentTime = 0.0;
    }
  };
}
#endif
