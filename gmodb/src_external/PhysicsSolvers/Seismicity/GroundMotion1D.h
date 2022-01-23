/*
 * GroundMotion1D.h
 *
 *  Created on: Mar 15, 2013
 *      Author: johnson346
 *
 *      Fortran conversion from grmdlib.f aided by FABLE (LBL)
 *      port of gmsiglib.f and grmdlib.f
 */

#ifndef GROUNDMOTION1D_H_
#define GROUNDMOTION1D_H_

#include "IO/ticpp/HierarchicalDataNode.h"
#include "Utilities/StringUtilities.h"

namespace EarthquakeSimulation
{
  enum GroundMotionModel
  {
    BOORE_ATKINSON_RV = 0,
    TORRO_MCGUIRE_RV = 1,
    EXPERT2_RV = 2,
    EXPERT3_RV = 3,
    TRIFUNAC_ANDERSON = 4,
    MODEL_F = 5,
    EXPERT2_SPECTRAL_RV = 6,
    EXPERT3_RV_5SV = 7,
    EXPERT2_RV_5SVA = 8,
    ABRAHAMSON_SILVA = 11,
    CAMPBELL_1997 = 12,
    CAMPBELL_BOZORGNIA = 13,
    PGE_LTSP = 14,
    IDRISS_1987 = 15,
    XJB_82 = 16,
    CAMPBELL_1989 = 17,
    WCC_NPR = 18,
    ROMANIAN_1997 = 19,
    SABETTA_PUGLIESE = 20,
    SADIGH_YOUNGS_ROCK = 21,
    SADIGH_YOUNGS_SOIL = 22,
    ABRAHAMSON_TIP = 23,
    TORE_1997 = 24,
    YOUNGS_1997 = 25
  };

  class GroundMotion1D
  {
  public:
    GroundMotion1D();
    virtual ~GroundMotion1D();

    void ReadXML(TICPP::HierarchicalDataNode* hdn);

    inline void SetModel(const GroundMotionModel model) { m_model = model; }

    inline void SetModel(const std::string name) { SetModel(StringToModel(name)); }

    static inline GroundMotionModel StringToModel(const std::string name)
    {
      GroundMotionModel model;
      const char* cstr = name.c_str();
      if(streq(cstr, "BooreAtkinson"))
        model = BOORE_ATKINSON_RV;
      else if(streq(cstr, "TorroMcGuire"))
        model = TORRO_MCGUIRE_RV;
      else if(streq(cstr, "Expert2"))
        model = EXPERT2_RV;
      else if(streq(cstr, "Expert3"))
        model = EXPERT3_RV;
      else if(streq(cstr, "TrifunacAnderson"))
        model = TRIFUNAC_ANDERSON;
      else if(streq(cstr, "ModelF"))
        model = MODEL_F;
      else if(streq(cstr, "Expert2Spectral"))
        model = EXPERT2_SPECTRAL_RV;
      else if(streq(cstr, "Expert3Spectral"))
        model = EXPERT3_RV_5SV;
      else if(streq(cstr, "Expert2_5SVa"))
        model = EXPERT2_RV_5SVA;
      else if(streq(cstr, "AbrahamsonSilva"))
        model = ABRAHAMSON_SILVA;
      else if(streq(cstr, "Campbell1997"))
        model = CAMPBELL_1997;
      else if(streq(cstr, "CampbellBozorgnia"))
        model = CAMPBELL_BOZORGNIA;
      else if(streq(cstr, "PGE_LTSP"))
        model = PGE_LTSP;
      else if(streq(cstr, "Idriss1987"))
        model = IDRISS_1987;
      else if(streq(cstr, "XJB_82"))
        model = XJB_82;
      else if(streq(cstr, "Campbell1989"))
        model = CAMPBELL_1989;
      else if(streq(cstr, "WCC_NPR"))
        model = WCC_NPR;
      else if(streq(cstr, "Romanian1997"))
        model = ROMANIAN_1997;
      else if(streq(cstr, "SabettaPugliese"))
        model = SABETTA_PUGLIESE;
      else if(streq(cstr, "SadighYoungsRock"))
        model = SADIGH_YOUNGS_ROCK;
      else if(streq(cstr, "SadighYoungsSoil"))
        model = SADIGH_YOUNGS_SOIL;
      else if(streq(cstr, "AbrahamsonTIP"))
        model = ABRAHAMSON_TIP;
      else if(streq(cstr, "Tore1997"))
        model = TORE_1997;
      else if(streq(cstr, "Youngs1997"))
        model = YOUNGS_1997;
      else
        throw GPException("GroundModel1D: model name not recognized!");
      return model;
    }

    inline void SetCoefficient(const localIndex icoefficient, const realT value)
    {
      if(m_coefficients.size() <= icoefficient)
        m_coefficients.resize(icoefficient + 1);
      m_coefficients[icoefficient] = value;
    }

    void
    IncrementHazard(const realT distance,
                    const realT magnitude,
                    const rArray1d& logPeakAccelerationBins,
                    rArray1d& logPeakAccelerationHazards);

//    void GroundMotionHazard(const realT distance,
//                            const realT magnitude,
//                            const rArray1d& logPeakAccelerationBins,
//                            rArray1d& logPeakAccelerationHazards,
//                            iArray1d& logPeakAccelerationHazardCounts,
//                            const int nbins = 10,
//                            const int nsamples = 5);

  private:

    static void
    IncrementHazard2(const realT logPeakAccelerationMean,
                     const realT sigma,
                     const rArray1d& logPeakAccelerationBins,
                     rArray1d& logPeakAccelerationHazards);

//    static void
//    EpistemicHazard(const realT logPeakAccelerationMean,
//                    const realT sigma,
//                    const int nsigma,
//                    const rArray1d& logPeakAccelerationBins,
//                    rArray1d& logPeakAccelerationHazards,
//                    iArray1d& logPeakAccelerationHazardCounts,
//                    const int nbins,
//                    const int nsamples);


    // magnitude  = magnitude value
    // distance   = actual arithmetic distance value (km)
    // zzk = the output Neperian (base-e) logarithm of the
    //       calculated ground motion (log(cm/s/s))
    inline realT LogPeakAcceleration(const realT distance,
                                     const realT magnitude,
                                     realT& stdev) const
    {
      stdev = m_standardDeviation;
      return LogPeakAcceleration(distance, magnitude, m_coefficients, m_model, stdev);
    }

    static realT
    LogPeakAcceleration(const realT distance,
                        const realT magnitude,
                        const rArray1d& coefficient,
                        const GroundMotionModel model,
                        realT& stdev);

  private:
    GroundMotionModel m_model;
    rArray1d m_coefficients;
    realT m_standardDeviation;
    int m_numberStandardDeviations;

    //
    //******************************************************************************
    //                Library of Ground Motion Models
    //******************************************************************************


    ///boore-atkinson rv-model
    /**    model a (model index=1)
     * if model c(19)=1.55 new 11/4/92 Boore's model fit
     */
    static realT
    BooreAtkinsonRVModel(const realT xment,
                         const realT d,
                         const rArray1d& c);


    //
    //---------------------------------------
    //
    //    model b  (model index=2)
    // torro-mcguire rv-model
    //
    static realT
    TorroMcGuireRVModel(const realT xm,
                        const realT d,
                        const rArray1d& c);

    //
    //---------------------------------------
    //
    //    model c (index=3)
    //
    //    expert 2's rv-accel model
    //
    static realT
    SecondExpertRVModel(const realT xm,
                        const realT d,
                        const rArray1d& c);

    //
    //---------------------------------------
    //
    //    d model  (index=4)
    //
    //    expert 3's rv-5a accel model
    //
    static realT
    ThirdExpertRVModel(const realT xm,
                       const realT d,
                       const rArray1d& c);

    //
    //---------------------------------------
    //
    //      model e  (model index = 5)
    //      *******
    //    trifunac-anderson  accel model
    //
    static realT
    TrifunacAndersonModel(const realT xm,
                          const realT d,
                          const rArray1d& c,
                          const realT xi,
                          const int icat);

    //
    //---------------------------------------
    //
    //      model f (model index = 6)
    //    model f (model index = 6)
    //    *******
    //
    static realT
    ModelF(const realT xm,
           const realT d,
           const rArray1d& c);


    //
    //---------------------------------------
    //
    //    model g  (model index=7)
    //
    //    expert 2's rv spectral model
    //
    static realT
    SecondExpertSpectralRVModel(const realT xm,
                                const realT d,
                                const rArray1d& c);


    //
    //---------------------------------------
    //
    //     model  h   (index = 8 )
    //
    //    expert 3's rv-5sv model for spectra
    //
    static realT
    ThirdExpertRV5SVModel(const realT xm,
                          const realT d,
                          const rArray1d& c);


    //
    //---------------------------------------
    //
    //    model xi  (index = 9 )
    //
    //    expert 2's n-h spectral model using his rv-5a, v models
    //
    static realT
    SecondExpertRV5aModel(const realT xm,
                          const realT d,
                          const rArray1d& c,
                          const int l);


//
//    //
//    //==========================================
//    //
//    void
//    numark2(
//      const int maxattncoeff,
//      const realT xm,
//      const realT /* xi */,
//      const realT r,
//      const realT /* alr */,
//      rArray1d& attn,
//      realT& z,
//      int& iflagav)
//    {
//      attn(dimension(maxattncoeff));
//      arr_1d<2, int> mod(fem::fill0);
//      int m2 = 0;
//      arr_1d<2, realT> x(fem::fill0);
//      int i = 0;
//      int l = 0;
//      int m = 0;
//      //
//      //==========================================
//      //
//      //    this routine is used only when dealing with newmark spectra.
//      //    for a (any) frequency, the sv value is calculated in two ways.
//      //         1. the spectrum considered to be anchored on a velocity
//      //            attenuaion curve
//      //         2. the spectrum considered to be anchored on an acceleration
//      //            attenuation curve.
//      //    the value of sv used is taken as the smaller of the two above
//      //    values.
//      //    the velocity part of the eqn is given by the first 8 coefficients
//      //    the acceleration part is given by coeficients 11 to 18 inclusive.
//      //
//      //    the variable iflagav is a flag which indicates if we are in
//      //    the velocity/acceleration regime or in the acceleration regime.
//      //    iflagav = 1 vel/acc regime , we calculate both
//      //              2 acc regime only. calculate only acceleration.
//      //
//      mod(1) = attn(19);
//      mod(2) = attn(20);
//      m2 = 3 - iflagav;
//      //
//      //    loop over the two parts of the equation
//      x(2) = 0.0;
//      FEM_DO_SAFE(i, 1, m2) {
//        l = 20 - (i * 10);
//        m = mod(i);
//        //
//        switch (m) {
//          case 1: goto statement_100;
//          case 2: goto statement_100;
//          case 3: goto statement_100;
//          case 4: goto statement_100;
//          case 5: goto statement_100;
//          case 6: goto statement_200;
//          case 7: goto statement_100;
//          case 8: goto statement_100;
//          case 9: goto statement_300;
//          default: break;
//        }
//        statement_100:
//        goto statement_1000;
//        statement_200:
//        x(i) = fmodel(maxattncoeff, xm, r, attn, l);
//        goto statement_1000;
//        statement_300:
//        x(i) = ximod(maxattncoeff, xm, r, attn, l);
//        statement_1000:;
//      }
//      //
//      z = x(1);
//      if (iflagav == 2) {
//        goto statement_2000;
//      }
//      if (z < x(2)) {
//        goto statement_1500;
//      }
//      z = x(2);
//      return;
//      statement_1500:
//      iflagav = 2;
//      statement_2000:;
//      //
//    }
//
//    //
//    //==========================================
//    //
//    void
//    trifun(
//      const int maxattncoeff,
//      const realT xi,
//      const realT r,
//      const realT alr,
//      const realT vl1,
//      rArray1d& attn,
//      const int icat,
//      realT& p)
//    {
//      attn(dimension(maxattncoeff));
//      //
//      //==========================================
//      //
//      //      this routine returns the probability of the spectral velocity
//      //      being smaller than vl for a site located at distance r from the
//      //      source of intensity xi.
//      //      the site correction factor is given by its natural log., corsite.
//      //      the equations are taken from trifunac and anderson report ce 77-03
//      //      and trif. & lee  ce 85-04
//      //      of usc.
//      //
//      //  convert  trifunac & lee to cm/sec
//      realT attn3 = attn(3) + .405;
//      //
//      realT xsit = 0.0;
//      if (icat == 1) {
//        xsit = 2.0;
//      }
//      if (icat == 2 || icat == 3) {
//        xsit = 1.0;
//      }
//      if (icat == 6 || icat == 7) {
//        xsit = 1.0;
//      }
//      //
//      //      calculate the attenuated intensity via modified gupta-nuttli eq.
//      realT xis = attn(4) + xi + attn(5) * r + attn(6) * alr;
//      //
//      //      calculate the pl
//      realT pl = ((vl1 / 2.30259) - attn(2) * xis - attn3 - attn(14) *
//        xsit) / attn(1);
//      //
//      //      calculate the pa
//      realT z = exp(attn(7) * pl + attn(8));
//      p = 1.0 - exp(-z);
//      if (attn(11) > 1.0) {
//        p = fem::pow(p, attn(11));
//      }
//    }
//
    //
    //---------------------------------------
    //
    //    Model  xabsilva    index  =  20
    //
    //    added 4/24/2001
    //       Abrahamson and Silva model: Seismological Research Letters
    //       January/February 1997
    //
    static realT
    AbrahamsonSilvaModel(const realT xmm,
                         const realT d,
                         const rArray1d& c);

    static realT AbrahamsonSilvaStdev(const realT xm, const rArray1d& c)
    {
      //comment indicates the following is for "Abramson SLR97" ... probably for Abrahamson
      if(c.size() < 24)
        throw GPException("c array too short");
      const realT stdev = xm >= 7.0 ? c(22) - 2.0*c(23) :
          (xm > 5.0 ? c(22) - c(23)*(xm-5.0) :
              c(22));
      return stdev;
    }

    //
    //---------------------------------------
    //
    // model  xcampbell97    index  =  23
    //
    //   added 5/2/2001
    //  Campbell, SRL, Jan/Feb. 1997 model
    //
    static realT
    Campbell1997Model(const realT xm,
                     const realT d,
                     const rArray1d& c);

    static realT
    Campbell1997Stdev(const rArray1d& c)
    {
      if(c.size() < 29)
        throw GPException("c array too short");
      //comment indicates the following is for "Campbell SLR97"
      //const realT ggk = zzk - 6.889;
      //const realT zz = exp (ggk);
      //stdev = zz < coefficient(22) ? coefficient(23) : (zz <= coefficient(24) ?
      //    coefficient(25) - coefficient(26)*ggk : coefficient(27));
      realT stdev = 0.47;//<--this value was hard-coded in and overrides the calculation commented above
      stdev = sqrt (stdev*stdev + c(28)*c(28));
      return stdev;
    }


    //
    //---------------------------------------
    //    Campbell and Bozorgnia (2007-NGA) model
    //
    //    THIS IS MODEL NUMBER : 25
    //
    static realT
    CampbellBozorgniaModel(const realT xm,
                           const realT d,
                           const rArray1d& c);

    static realT
    CampbellBozorgniaStdev(const rArray1d& c)
    {
      if(c.size() < 25)
        throw GPException("c array too short");
      //comment indicates the following is for "Campbell and Bozorgnia NGA 2007"
      //For ANPP assume VS30 = 620 m/s for rock
      const realT stdev = c(14) > 620. ? sqrt (c(22)*c(22) + c(23)*c(23)) : c(24);
      return stdev;
    }

    //
    //---------------------------------------
    //
    // model xdabc   index=14
    //
    // model in PG&E LTSP rpt   more generally sadigh's model
    //
    //   c(8)   plays the role of fault type
    //
    //  c(20) plays the role of a dummy depth to compare to J&B model
    //
    static realT
    PGE_LTSP_Model(const realT xm,
                   const realT d,
                   const rArray1d& c,
                   const int ll);


    //
    //---------------------------------------
    //
    //   model xidrs   index= 13
    //
    // added 5/21/90
    //idriss's 1987 model  fron j&b 1988 summary paper
    //
    static realT
    IdrissModel(const realT xm,
                const realT dd,
                const rArray1d& c);


    //
    //---------------------------------------
    //
    // model  xjb82    index  =  11
    //
    //   added 5/21/90
    //  joyner and boore's 1982 model
    //
    //  c(5)  plays the role of the h parameter in their model
    //
    static realT
    XJB82Model(const realT xmm,
               const realT d,
               const rArray1d& c);



    //
    //---------------------------------------
    //
    //    added 5/21 /90
    //
    static realT
    Campbell1989Model(const realT xm,
                      const realT d,
                      const rArray1d& c);



    //
    //---------------------------------------
    //
    // model  xmod15    index=15
    //  generic model of the form
    // ln(a)=c1+c2*m+c3*m**2 +c4*ln(r) +c5*r
    //
    //  in particular form used for my fit to WCC model for the NPR at INEL
    //
    static realT
    WCC_NPR_Model(const realT xm,
                  const realT d,
                  const rArray1d& c);

    //
    //---------------------------------------
    //
    static realT
    Romanian1997Model(const realT xm,
                      const realT d,
                      const rArray1d& c);


    //
    //---------------------------------------
    //
    static realT
    SabettaPuglieseModel(const realT xm,
                         const realT d,
                         const rArray1d& c);

    //
    //---------------------------------------
    //
    //    Model  xSadighRock    index  =  21
    //
    //    added 5/1/2001
    //       Sadigh, Youngs et al. Rock model: Seismological Research Letters
    //       January/February 1997
    //
    static realT
    SadighYoungsRockModel(const realT xm,
                          const realT d,
                          const rArray1d& c);

    static realT
    SadighStdev(const realT xm,
                const rArray1d& c)
    {
      //comment indicates the following is for "Saddigh SLR97" .... not clear which model this applies to
      const realT stdev = xm < c(23) ? c(24) - c(25)*xm : c(22);
      return stdev;
    }


    //
    //---------------------------------------
    //
    //    Model  xSadighSoil    index  =  22
    //
    //    added 5/1/2001
    //       Sadigh, Youngs et al. Soil model: Seismological Research Letters
    //       January/February 1997
    //
    static realT
    SadighYoungsSoilModel(const realT xm,
                          const realT d,
                          const rArray1d& c);



    //
    //
    //
    //      model xtip98  index = 16
    //
    //      added july 28,1998, by JBS
    //
    //      Model developed for the TIP project, by N. Abrahamson
    //       c(19) is the model type (=16)
    //       c(10) is the magnitude threshold (6.25)
    //       c(11), c(12) and c(13) are for the Mblg to Mw
    //             conversion
    //       Values of PGA are in cm/s/s
    //       Values of PGV are in cm/s
    //       Log are Natural Logs
    //
    static realT
    AbrahamsonTIPModel(const realT xmm,
                       const realT d,
                       const rArray1d& c);

    //
    //---------------------------------------
    static realT
    Tore1997Model(const realT xmm,
                  const realT d,
                  const rArray1d& c);
    //
    //---------------------------------------
    static realT
    Youngs1997Model(const realT xmm,
                    const realT d,
                    const rArray1d& c);
  };
}
#endif
