/*
 * GroundMotion1D.cpp
 *
 *  Created on: Mar 15, 2013
 *      Author: johnson346
 */

#include "GroundMotion1D.h"
#include "SurfaceGeneration/StatisticalDistributionBaseT.h"

namespace EarthquakeSimulation
{

  GroundMotion1D::GroundMotion1D() : m_model(BOORE_ATKINSON_RV)
  {
    // TODO Auto-generated constructor stub
  }

  GroundMotion1D::~GroundMotion1D()
  {
    // TODO Auto-generated destructor stub
  }

  void
  GroundMotion1D::ReadXML(TICPP::HierarchicalDataNode* hdn)
  {
    //set model type
    const std::string modelName = hdn->GetAttributeString("type");
    SetModel(modelName);

    //get standard deviation info
    m_standardDeviation = hdn->GetAttributeOrDefault<realT>("defaultStandardDeviation", 0.0);
    m_numberStandardDeviations = hdn->GetAttributeOrDefault<int>("numberOfStandardDeviations", 1);

    //read coefficients
    for(TICPP::HierarchicalDataNode* child = hdn->Next(true); child; child = hdn->Next())
    {
      const int index = child->GetAttributeOrDefault<int>("index", -1);
      if(index < 0)
        throw GPException("GroundMotion1D: Cannot specify a negative coefficient index");
      if(m_coefficients.size() < ((unsigned int)index + 1))
        m_coefficients.resize(index + 1, 0.0);
      m_coefficients[index] = child->GetAttributeOrDefault("value", 0.0);
    }
  }

  void
  GroundMotion1D::IncrementHazard(const realT distance,
                                  const realT magnitude,
                                  const rArray1d& logPeakAccelerationBins,
                                  rArray1d& logPeakAccelerationHazards)
  {
    realT stdev = m_standardDeviation;

    //NOTE: input to LogPeakAcceleration should be km and output is cm/s/s
    //assuming input is SI: DO CONVERSION
    const realT distance_km = distance * 0.001;
    const realT zzk = LogPeakAcceleration(distance_km, magnitude, stdev);

    std::cout << "GM at distance=" << distance_km << " km for M0=" << magnitude << " is " << (exp(zzk)) << " cm/s/s (log=" << zzk << ")" << std::endl;

    IncrementHazard2(zzk, stdev, logPeakAccelerationBins, logPeakAccelerationHazards);
  }

  void
  GroundMotion1D::IncrementHazard2(const realT logPeakAccelerationMean,
                                   const realT sigma,
                                   const rArray1d& logPeakAccelerationBins,
                                   rArray1d& logPeakAccelerationHazards)
  {
    rArray1d::iterator ith = logPeakAccelerationHazards.begin();
    for(rArray1d::const_iterator it = logPeakAccelerationBins.begin();
        it != logPeakAccelerationBins.end(); ++it, ++ith)
    {
      //find the probability that the actual value of the measurement
      //with normally distributed error >= the bin value
      const realT z = (*it - logPeakAccelerationMean) / sigma;
      //where z is the normalized, signed distance from the mean

      //if d is positive, then the probability that a measurement, under the influence of
      //normally distributed errors with standard deviation sigma, has a distance less
      //than d from the mean value is erf(d/(sigma*sqrt(2)))
      //so the probability that a measurement has a distance > d from the mean is erfc(d/(sigma*sqrt(2)))
      //therefore, if d is negative, the (one-tailed) probability that a measurement is < (mean - fabs(d))
      //is 0.5*erfc(d/(sigma*sqrt(2))), therefore the complement is the probability of exceedence
      //conversely, if d is positive, the (one-tailed) probability that a measurement is > (mean + z)
      //is 0.5*erfc(d/(sigma*sqrt(2)))
      //using z = d/sigma = (value - mean)/sigma

      *ith += z < 0 ? 1.0 - StatisticalDistributionBaseT::PhiNormalDifferenceInf(fabs(z)) :
          StatisticalDistributionBaseT::PhiNormalDifferenceInf(z);
    }
  }

//  void
//  GroundMotion1D::GroundMotionHazard(const realT distance,
//                                     const realT magnitude,
//                                     const rArray1d& logPeakAccelerationBins,
//                                     rArray1d& logPeakAccelerationHazards,
//                                     iArray1d& logPeakAccelerationHazardCounts,
//                                     const int nbins,
//                                     const int nsamples)
//  {
//    realT stdev = m_standardDeviation;
//
//    //NOTE: be careful ... these are always given in cm/s/s!!!
//    const realT zzk = LogPeakAcceleration(distance, magnitude, stdev);
//
//    std::cout << "GM at distance=" << distance << " for M0=" << magnitude << " is " << (exp(zzk)) << " cm/s (log=" << zzk << ")" << std::endl;
//
//    EpistemicHazard(zzk, stdev, m_numberStandardDeviations,
//                    logPeakAccelerationBins,
//                    logPeakAccelerationHazards,
//                    logPeakAccelerationHazardCounts,
//                    nbins,
//                    nsamples);
//  }
//
//  void
//  GroundMotion1D::EpistemicHazard(const realT logPeakAccelerationMean,
//                                  const realT sigma,
//                                  const int nsigma,
//                                  const rArray1d& logPeakAccelerationBins,
//                                  rArray1d& logPeakAccelerationHazards,
//                                  iArray1d& logPeakAccelerationHazardCounts,
//                                  const int nbins,
//                                  const int nsamples)
//  {
//    //Latin-HyperCube sampling the distribution of GM
//    //-----------------------------------------------
//    //with "nsamples" in each bin, and using "nbins" bins
//    //Limit the excursion to nsigma (Read from input)
//    // nsigma = Max number of sigmas for sampling
//
//    //Bounds of the distribution
//    const realT zzkmin = zzk - nsigma * sigma;
//    const realT zzkmax = zzk + nsigma * sigma;
//    const realT dzzk = (zzkmax - zzkmin) / nbins;
//
//    //Normalization factor
//    const realT xnormal = 2.*StatisticalDistributionBaseT::ProbabilityOfExceedence(dzzk) - 1.0;
//
//    //Loop on number of LHC-sampled bins
//    realT zzknorm = StatisticalDistributionBaseT::ProbabilityOfExceedence( (zzkmin - zzk) / sigma );
//
//    const realT binweight = 1.0 / nbins;
//    for(int ibin = 0; ibin < nbins; ++ibin)
//    {
//      const realT binmin = zzkmin + dzzk * ibin;
//      const realT binmax = binmin + dzzk;
//      //Weight of the LH sampling bin
//      const realT binnorm = StatisticalDistributionBaseT::ProbabilityOfExceedence( (binmax - zzk ) / sigma );
//      const realT sampleweight = ( binnorm - zzknorm ) / (xnormal * nsamples);
//      zzknorm = binnorm;
//      for(int ii = 0; ii < nsamples; ++ii)
//      {
//        //now we have the bounds and weight of the bin
//        //Do the Uniform sampling within a bin
//        const realT binval = StatisticalDistributionBaseT::UniformSample(binmin, binmax);
//        iArray1d::iterator ithc = logPeakAccelerationHazardCounts.begin();
//        rArray1d::iterator ith = logPeakAccelerationHazards.begin();
//        for(rArray1d::const_iterator itb = logPeakAccelerationBins.begin(); itb != logPeakAccelerationBins.end(); ++itb, ++ith, ++ithc)
//        {
//          if(*itb > binval)
//          {
//            *ith += sampleweight * binweight;
//            *ithc += 1;
//            break;
//          }
//        }
//      }
//    }
//  }

  realT
  GroundMotion1D::LogPeakAcceleration(const realT distance,
                                      const realT magnitude,
                                      const rArray1d& coefficient,
                                      const GroundMotionModel model,
                                      realT& stdev)
  {
    //DISTANCES ARE KM, RETURNED VALUE IS CM/S/S!!!
    const int l = 0;
    if(coefficient(22) > 0)
      stdev = coefficient(22);
    switch(model)
    {
      case BOORE_ATKINSON_RV:
        return BooreAtkinsonRVModel(magnitude, distance, coefficient);
      case TORRO_MCGUIRE_RV:
        return TorroMcGuireRVModel(magnitude, distance, coefficient);
      case EXPERT2_RV:
        return SecondExpertRVModel(magnitude, distance, coefficient);
      case EXPERT3_RV:
        return ThirdExpertRVModel(magnitude, distance, coefficient);
      case TRIFUNAC_ANDERSON:
        //NOTE: magnitude is used as xi and 0 is hard-coded as icat in original
        return TrifunacAndersonModel(magnitude, distance, coefficient, magnitude, 0);
      case MODEL_F:
        return ModelF(magnitude, distance, coefficient);
      case EXPERT2_SPECTRAL_RV:
        return SecondExpertSpectralRVModel(magnitude, distance, coefficient);
      case EXPERT3_RV_5SV:
        return ThirdExpertRV5SVModel(magnitude, distance, coefficient);
      case EXPERT2_RV_5SVA:
        //implicit declaration of
        return SecondExpertRV5aModel(magnitude, distance, coefficient, l);
      case ABRAHAMSON_SILVA:
        stdev = AbrahamsonSilvaStdev(magnitude, coefficient);
        return AbrahamsonSilvaModel(magnitude, distance, coefficient);
      case CAMPBELL_1997:
        stdev = Campbell1997Stdev(coefficient);
        return Campbell1997Model(magnitude, distance, coefficient);
      case CAMPBELL_BOZORGNIA:
        stdev = CampbellBozorgniaStdev(coefficient);
        return CampbellBozorgniaModel(magnitude, distance, coefficient);
      case PGE_LTSP:
        return PGE_LTSP_Model(magnitude, distance, coefficient, l);
      case IDRISS_1987:
        return IdrissModel(magnitude, distance, coefficient);
      case XJB_82:
        return XJB82Model(magnitude, distance, coefficient);
      case CAMPBELL_1989:
        return Campbell1989Model(magnitude, distance, coefficient);
      case WCC_NPR:
        return WCC_NPR_Model(magnitude, distance, coefficient);
      case ROMANIAN_1997:
        return Romanian1997Model(magnitude, distance, coefficient);
      case SABETTA_PUGLIESE:
        return SabettaPuglieseModel(magnitude, distance, coefficient);
      case SADIGH_YOUNGS_ROCK:
        stdev = SadighStdev(magnitude, coefficient);
        return SadighYoungsRockModel(magnitude, distance, coefficient);
      case SADIGH_YOUNGS_SOIL:
        return SadighYoungsSoilModel(magnitude, distance, coefficient);
      case ABRAHAMSON_TIP:
        return AbrahamsonTIPModel(magnitude, distance, coefficient);
      case TORE_1997:
        return Tore1997Model(magnitude, distance, coefficient);
      case YOUNGS_1997:
        return Youngs1997Model(magnitude, distance, coefficient);
      default:
        throw GPException("unrecognized ground motion model!!");
        break;
    }
    return 0.0;
  }




  ///boore-atkinson rv-model
  realT
  GroundMotion1D::BooreAtkinsonRVModel(const realT xment,
                                       const realT d,
                                       const rArray1d& c)
  {
    if(c.size() < 20)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT r = 0.0;
    realT xlr = 0.0;
    realT xm2 = 0.0;
    realT xm3 = 0.0;
    int l = 0;
    realT v = 0.0;
    int jc10 = 0;
    realT yc10 = 0.0;
    int ic10 = 0;
    realT xc10 = 0.0;
    realT xm = xment;
    const realT pi = 0.25 * atan(1);
    //truncation for new Boore model c(19)=1.55
    if (c(19) > 1.54999 && c(19) < 1.55001) {
      if (xm > 7.5) {
        xm = 7.5;
      }
    }
    //truncation for new herrmann 11/11/92 model c(19)=1.45
    if (c(19) > 1.44999 && c(19) < 1.45001) {
      if (xm > 7.5) {
        xm = 7.5;
      }
    }
    r = sqrt(d*d + c(20)*c(20));
    xlr = log10(r);
    xm2 = xm * xm;
    xm3 = xm * xm2;
    l = 0;
    if (r > 100.0) {
      l = 10;
    }

    bool goto60 = false;
    if (c(19) > 1.9 && xm > 6.4 && r <= 100) {
      //goto statement_45
      v = -3.48 + 1.349 * xm - .15636 * xm2 + .010277 * xm3 - .0383 *
        xm * xlr + .00907 * xm2 * xlr - .000537 * xm3 * xlr - .001463 *
        r - xlr;
      goto60 = true;
    }
    else
    {
      //statement_30:
      v = c(1 + l) + c(2 + l) * xm + c(3 + l) * r + c(4 + l) * xm2 + c(5 +
        l) * xm3 + c(6 + l) * xm * xlr + c(7 + l) * xm2 * xlr + c(8 + l) *
        xm3 * xlr;
      if (r > 100.0) {
        //goto statement_50;
        v += c(9) * xlr;
        if (c(19) < 1.4) {
          v += c(10) * cos(2.0 * pi * (xlr - 2.375));
        }
      }
      else
      {
        v = v - xlr;
        if (c(19) > 1.5) {
          v += c(10) * cos(pi * (xm - 3.75));
        }
        // change added for fit to Boore's 11/4/92 model
        if (c(19) == 1.55) {
          jc10 = c(20);
          yc10 = c(20) - jc10;
          v += yc10 * cos(2.0 * pi * (xlr - 2.375));
        }
        // end of change
        goto60 = true;
      }
    }

    //
    //   next card corrects error in logic  added  11/16/91
    //
    //  most cases ok because c(20) was a number like 15.000 km
    //   but if c(20) = XX.X then a problem exists  in that xc10 .ne. 0.
    //  original coding overlooked this point
    //
    //if (c(19) < 1.5) {
    //  goto statement_60;
    //}

    //
    //   since coef of cos term is always small inclue in c(20) when
    // needed        c19=1.1 std model if c10 not zero then with
    //  c10cos(2pi(logr-2.375)) term in ff
    //  for cases when both ff and nf branches are involved  the coef for the
    //  ff branch is carried as the dec part of the depth term c20
    //   term
    //            c19=1.6  with cos(pi(m-3.75)) in near field branch
    //             c19=1.8 cos(pi(m-3.75)) in both nf & ff branches
    //            c19=1.85 cos(pi(m-3.75)) in nf branch & cos(2pi(logr-2.375)) in ff
    //  c19=1.91  aki's hardwired model for charleston large mag fit nf
    //
    //if(c(19) > 1.7) {
    //  if (c(19) > 1.9) {
    //    goto statement_60;
    //  }
    //}
    if (!goto60 && c(19) >= 1.5 && c(19) <= 1.7) //use this instead of goto statements
    {
      ic10 = c(20);
      xc10 = c(20) - ic10;
      if (c(19) < 1.84) {
        v += xc10 * cos(pi * (xm - 3.75));
      }
      else if (c(19) > 1.84) {
        v += xc10 * cos(2.0 * pi * (xlr - 2.375));
      }
    }
    //statement_60:
    return_value = v * 2.3026;
    return return_value;
  }


  /**
   *---------------------------------------
   *
   *    model b  (model index=2)
   * torro-mcguire rv-model
   *
   */
  realT
  GroundMotion1D::TorroMcGuireRVModel(const realT xm,
                                      const realT d,
                                      const rArray1d& c)
  {
    if(c.size() < 21)
      throw GPException("c array too short");
    realT return_value = 0.0;
    realT r = 0.0;
    realT sr = 0.0;
    realT xlr = 0.0;
    realT xm2 = 0.0;
    realT xm3 = 0.0;
    realT v = 0.0;
    //
    r = sqrt(d * d + c(20) * c(20));
    sr = sqrt(r);
    xlr = log10(r);
    xm2 = (xm * xm);
    xm3 = xm2 * xm;
    if (r > 100.0) {
      //goto statement_100;
      //statement_100:
      v = c(11) + c(12) * xm + c(13) * r + c(14) * xlr + c(15) * xm2 + c(
        16) * xm3 + c(17) * xm * r + c(18) * xm * xlr + (c(8) + c(9) *
        xm2 + c(10) * xm3) / sr;
    }
    else
    {
      v = c(1) + c(2) * xm + c(3) * r - xlr + c(5) * xm2 + c(6) * xm3 + c(
        7) * xm * r + c(4) * (xm - 4.25) * (xm - 8.0) * sin(
        3.1415926535898 * (xm - 4.5));
      //goto statement_200;
    }
    //statement_200:
    return_value = v * 2.3026;
    return return_value;
  }



  /**
   *---------------------------------------
   *
   *    model c (index=3)
   *
   *    expert 2's rv-accel model
   */
  realT
  GroundMotion1D::SecondExpertRVModel(const realT xm,
                                      const realT d,
                                      const rArray1d& c)
  {
    if(c.size() < 19)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT xm2 = 0.0;
    realT xm3 = 0.0;
    int l = 0;
    realT h = 0.0;
    realT r = 0.0;
    realT xlr = 0.0;
    realT v = 0.0;
    //
    xm2 = (xm * xm);
    xm3 = xm2 * xm;
    l = 0;
    h = 2.5 * (xm - 1.0);
    if (xm < 5.0) {
      h = 5.0 * (xm - 3.0);
    }
    r = sqrt((d * d) + (h * h));
    xlr = log10(r);
    if (r > 100.0) {
      //goto statement_100;
      //statement_100:
      v = 2.772 + .248 * xm - .00119 * r - 3.432 * xlr + .000135 *
        xm * r + .501 * xm * xlr - .0288 * xm2 * xlr + .00208 *
        cos(8.378 * (xlr - 2.0));
    }
    else
    {
      if (xm <= 4.5) {
        l = 10;
      }
      v = c(1 + l) + c(2 + l) * xm + c(3 + l) * r - xlr + c(4 + l) *
        xm2 + c(5 + l) * xm * r + c(6 + l) * xm2 * r + c(7 + l) * xm *
        xlr + c(8 + l) * xm3 * xlr;
      //goto statement_200;
    }
    //statement_200:
    return_value = v * 2.3026;
    return return_value;
  }



  /**
   *---------------------------------------
   *
   *    d model  (index=4)
   *
   *    expert 3's rv-5a accel model
   */
  realT
  GroundMotion1D::ThirdExpertRVModel(const realT xm,
                                     const realT d,
                                     const rArray1d& c)
  {
    realT return_value = 0.0;
    realT xm2 = 0.0;
    realT xm3 = 0.0;
    realT r = 0.0;
    realT xlr = 0.0;
    realT v = 0.0;
    //
    xm2 = (xm * xm);
    xm3 = xm2 * xm;
    r = sqrt(d * d + 64.0);
    xlr = log10(r);
    if (r > 100.0)
    {
      //goto statement_100;
      //statement_100:
      v = c(11) + c(12) * xm + c(13) * r + c(14) * xlr + c(15) * xm2 + c(16) * xm3 + c(17) * xm * r + c(
          18) * xm3 * xlr - .0197 * cos(8.378 * (xlr - 2.0));
    }
    else
    {
      v = c(1) + c(2) * xm + c(3) * r - xlr + c(4) * xm2 + c(5) * xm2 * r + c(6) * xm * xlr + c(7) * xm2 * xlr + c(
          8) * xm3 * xlr;
      //goto statement_200;
    }
    //statement_200:
    return_value = v * 2.3026;
    return return_value;
  }




  /**
   *---------------------------------------
   *
   *      model e  (model index = 5)
   *      *******
   *    trifunac-anderson  accel model
   */
  realT
  GroundMotion1D::TrifunacAndersonModel(const realT xm,
                                        const realT d,
                                        const rArray1d& c,
                                        const realT xi,
                                        const int icat)
  {
    if (c.size() < 20)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT r = 0.0;
    realT alr = 0.0;
    realT xsit = 0.0;
    realT xi0 = 0.0;
    realT xis = 0.0;
    //
    r = d;
    if (r < c(9))
    {
      r = c(9);
    }
    if (c(19) < 5.7)
    {
      alr = log(r);
    }
    else
    {
      alr = log10(r);
    }
    xsit = 0.0;
    if (icat == 1)
    {
      xsit = 2.0;
    }
    if (icat == 2 || icat == 3)
    {
      xsit = 1.0;
    }
    if (icat == 6 || icat == 7)
    {
      xsit = 1.0;
    }
    if (c(19) > 5.2)
    {
      //goto statement_20;
      //statement_20:
      xi0 = c(11) + c(12) * xm + c(13) * xm * xm;
      if (xi0 >= 12)
      {
        xi0 = 12;
      }
      xis = c(1) + c(2) * xi0 + c(3) * r + c(4) * alr;
      return_value = c(6) + c(7) * xis + c(5) * xsit;
      if (c(19) > 5.7)
      {
        return_value = return_value * 2.30259;
      }
    }
    else
    {
      return_value = c(1) + c(2) * xi + c(3) * r + c(4) * alr + c(5) * xsit;
      //goto statement_30;
    }
    //statement_30:
    return return_value;
  }


  /**
   *---------------------------------------
   *
   *      model f (model index = 6)
   *    model f (model index = 6)
   *    *******
   */
  realT
  GroundMotion1D::ModelF(const realT xm,
                         const realT d,
                         const rArray1d& c)
  {
    if (c.size() < 11)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT h = 0.0;
    realT r = 0.0;
    int l = 0;
    //
    //*************************** added 7/30/91 ************
    if (c(10) > 3.0)
    {
      if (xm > c(10))
      {
        l = 10;
      }
    }
    const unsigned int uitmp = l + 5;
    if (c.size() < uitmp)
      throw GPException("c array too short");

    //***************************************************************
    //
    if (c(5) > 1.0)
    {
      //goto statement_100;
      //statement_100:
      h = c(5);
    }
    else
    {
      if (xm < 5.0)
      {
        h = 5.0 * (xm - 3.0);
      }
      if (xm >= 5.0)
      {
        h = 2.5 * (xm - 1.0);
      }
      //goto statement_101;
    }

    //statement_101:
    r = sqrt((d * d) + (h * h));
    return_value = c(1 + l) + c(2 + l) * xm + c(3 + l) * log(r) + c(4 + l) * r;
    return return_value;
  }



 /**
   *---------------------------------------
   *
   *    model g  (model index=7)
   *
   *    expert 2's rv spectral model
   */
  realT
  GroundMotion1D::SecondExpertSpectralRVModel(const realT xm,
                                              const realT d,
                                              const rArray1d& c)
  {
    if (c.size() < 21)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT xm2 = 0.0;
    realT xm3 = 0.0;
    realT h = 0.0;
    realT r = 0.0;
    realT xlr = 0.0;
    realT v = 0.0;
    //
    xm2 = (xm * xm);
    xm3 = xm * xm2;
    h = 2.5 * (xm - 1.0);
    r = sqrt(d * d + h * h);
    xlr = log10(r);
    if (r > 100.0)
    {
      //goto statement_100;
      //statement_100:
      v = c(11) + c(12) * xm + c(13) * r + c(14) * xlr + c(15) * xm2 + c(16) * xm * r + c(17) * xm * xlr + c(
          18) * xm2 * xlr + c(9) * cos(8.378 * (xlr - 2.0)) + c(10) * cos(1.57 * (xm - 5.0));
    }
    else
    {
      v = c(1) + c(2) * xm + c(3) * r - xlr + c(4) * xm2 + c(5) * xm * r + c(6) * xm2 * r + c(7) * xm * xlr + c(
          8) * xm3 * xlr + c(20) * sin(1.57 * (xm - 6.5));
      //goto statement_200;
    }
    //statement_200:
    return_value = v * 2.3026;
    return return_value;
  }



  /**
   *---------------------------------------
   *
   *     model  h   (index = 8 )
   *
   *    expert 3's rv-5sv model for spectra
   */
  realT
  GroundMotion1D::ThirdExpertRV5SVModel(const realT xm,
                                        const realT d,
                                        const rArray1d& c)
  {
    if (c.size() < 21)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT xm2 = 0.0;
    realT xm3 = 0.0;
    realT r = 0.0;
    realT xlr = 0.0;
    realT cs = 0.0;
    realT cc = 0.0;
    realT v = 0.0;
    realT sr = 0.0;
    //
    xm2 = (xm * xm);
    xm3 = xm2 * xm;
    r = sqrt(d * d + c(20) * c(20));
    xlr = log10(r);
    if (r > 100.0)
    {
      //goto statement_200;
      //statement_200:
      sr = 1.0 / sqrt(r);
      v = c(11) + c(12) * r + c(13) * xlr + c(14) * xm2 + c(15) * xm3 + c(16) * sr;
      if (c(19) == 8.1)
      {
        v += .000202 * xm * r - .0225 * cos(8.378 * (xlr - 2.14));
      }
      if (c(19) == 8.2 || c(19) == 8.3)
      {
        v += c(17) * xm3 * r + c(18) * xm * xlr + c(10) * xm3 * xlr;
      }
      //if (c(19) < 8.4) {
      //  goto statement_300;
      //}
      if (c(19) >= 8.4)
      { //added in lieu of goto
        cs = 0.0;
        if (c(19) == 8.5)
        {
          cs = 1.0;
        }
        v += (c(17) * xm + c(18) * xm2 + c(10) * xm3) * xlr - cs * .0337 * cos(
            3.1415926535898 * (xm - 5.0));
      }
    }
    else
    {
      //todo: check this logic ... used if(c(19)-8.4) 70,70,71 ... ???
      if(c(19) <= 8.4)
      {
        cs = 0.0;
        cc = 1.0;
      }
      else
      {
        cs = 1.0;
        cc = 0.0;
      }
      //statement_75:
      v = c(1) + c(2) * r - xlr + c(3) * xm2 + c(4) * xm3 + c(5) * xm * r + (c(6) * xm + c(7) * xm2 + c(
          8) * xm3) * xlr + c(9) * (cs * sin(2.098 * (xm - 5.0)) + cc * cos(1.047 * (xm - 5.0)));
      //goto statement_300;
    }
    //statement_300:
    return_value = v * 2.303;
    //
    return return_value;
  }



  /**
   *---------------------------------------
   *
   *    model xi  (index = 9 )
   *
   *    expert 2's n-h spectral model using his rv-5a, v models
   */
  realT
  GroundMotion1D::SecondExpertRV5aModel(const realT xm,
                                        const realT d,
                                        const rArray1d& c,
                                        const int l)
  {
    if (c.size() < 10)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT h = 0.0;
    realT xm2 = 0.0;
//    realT xm3 = 0.0;
    realT r = 0.0;
    realT xlr = 0.0;
    realT v = 0.0;
    //
    if (l == 10)
    {
      //goto statement_100;
      //statement_100:

      //  accel part
      v = SecondExpertRVModel(xm, d, c);

      //  note v has been converted to log base e in cmodel
      return_value = v - c(9);
      return return_value;
    }
    h = 2.5 * (xm - 1.0);
    xm2 = (xm * xm);
//    xm3 = xm2 * xm;
    //
    // vel part
    r = sqrt((d * d) + (h * h));
    xlr = log10(r);
    if (r > 100.0)
    {
      //goto statement_50;
      //statement_50:
      v = -5.398 + 1.998 * xm - .00026 * r - 1.663 * xlr - .129 * xm2 + .0198 * xm2 * xlr + .0144 * cos(
          1.57 * (xm - 5.0));
      //
      // note the .22 is the log of the amp fact for n-h velto sv
    }
    else
    {
      v = -3.169 + 1.024 * xm - .0245 * r - xlr - .0206 * xm2 + .00587 * xm * r - .000348 * xm2 * r - .00087 * xm2 * xlr - .0162 * sin(
          1.57 * (xm - 6.5));
      //goto statement_75;
    }
    //statement_75:
    return_value = (v + .22) * 2.303;
    return return_value;
  }


  /**
   *---------------------------------------
   *
   *    Model  xabsilva    index  =  20
   *
   *    added 4/24/2001
   *       Abrahamson and Silva model: Seismological Research Letters
   *       January/February 1997
   */
  realT
  GroundMotion1D::AbrahamsonSilvaModel(const realT xmm,
                                       const realT d,
                                       const rArray1d& c)
  {
    if (c.size() < 17)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT r = sqrt(d * d + (c(8) * c(8)));
    realT xc = 0.0;
    if (xmm <= c(7))
    {
      xc = c(2);
    }
    else
    {
      xc = c(4);
    }
    //
    //       Fault type factor
    //       c(15) = F  = 1.0 for Reverse
    //                                0.5 for reverse/oblique
    //                                0.0 for others
    realT xf = c(5);
    if ((xmm > 5.8) && (xmm <= c(7)))
    {
      xf = c(5) + ((c(6) - c(5)) / (c(7) - 5.8)) * (xmm - 5.8);
    }
    //
    if (xmm > c(7))
    {
      xf = c(6);
    }
    //
    realT y = c(1) + xc * (xmm - c(7)) + c(12) * (8.5 - xmm) * (8.5 - xmm) + (c(3) + c(13) * (xmm - c(
        7))) * log(r) + xf * c(15);
    //
    //       Site properties factor
    //       c(16) = S  =  1 for Rock and  Shallow soil
    //                                 0 for Deep Soil
    //       In the SRL equation, the unit of y is g's. Since we have
    //       transformed everything in cm/s/s, need to go back to g's
    //       in following statement.
    realT xs = exp(y) / 981.0;
    xs = c(16) * (c(10) + c(11) * (log(xs + c(14))));
    //
    return_value = y + xs;
    //       xabsilva = y + xs + 0.16187 !with correction for HW at 7.5km and M6
    //
    return return_value;
  }



  /**
   *---------------------------------------
   *
   * model  xcampbell97    index  =  23
   *
   *   added 5/2/2001
   *  Campbell, SRL, Jan/Feb. 1997 model
   */
  realT
  GroundMotion1D::Campbell1997Model(const realT xm,
                                    const realT dd,
                                    const rArray1d& c)
  {
    if (c.size() < 16)
      throw GPException("c array too short");

    realT d = dd;
    realT return_value = 0.0;
    //
    //Style of faulting =0 for SS, =1 for reverse)
    realT f = c(11);
    //Local site conditions (ssr=shr=0) for alluv/firm soil
    realT ssr = c(12);
    //(ssr=1,shr=0)/soft rock and (ssr=0,shr=1)/hard rock
    realT shr = c(13);
    //Depth to basement (=0 for Rock site)
    realT xd = c(14);
    //Depth to the seismogenic zone. By definition, the
    realT xs = c(15);
    //distace REIS has to be greater than th edistance to t
    //seismigenic zone (d >+ xs)
    if (d < xs)
    {
      d = xs;
    }
    xd = 5.0;
    //
    //       Calculate the natural log of Ah expressed in g's
    //       (per SLR 97, page 164)
    realT y = -3.512 + 0.904 * xm - 0.664 * log(d * d + (.0222 * exp(1.294 * xm))) + (1.125 - 0.112 * log(
        d) - 0.0957 * xm) * f + (0.440 - 0.171 * log(d)) * ssr + (0.405 - 0.222 * log(d)) * shr;
    //
    //       Calculate the Spectral ordinates, in g's
    //       (per SLR 97, page 164)
    realT yy = 0.0;
    realT fsa = 0.0;
    if (fabs((c(1) + c(2) + c(3) + c(4) + c(5) + c(6))) > 1.e-6)
    {
      fsa = 0.0;
      if (xd < 1.0)
      {
        fsa = c(6) * (1.0 - shr) * (1.0 - xd) + 0.5 * c(6) * (1.0 - xd) * ssr;
      }
      yy = c(1) + c(2) * tanh(c(3) * (xm - 4.7)) + (c(4) + c(5) * xm) * d + 0.5 * c(6) * ssr + c(
          6) * shr + c(7) * tanh(c(8) * xd) * (1.0 - shr) + fsa;
    }
    y += yy;
    //
    //     convert from units of g's to cm/s/s
    return_value = y + 6.889;
    //
    return return_value;
  }



  /**
   *---------------------------------------
   *    Campbell and Bozorgnia (2007-NGA) model
   *
   *    THIS IS MODEL NUMBER : 25
   */
  realT
  GroundMotion1D::CampbellBozorgniaModel(const realT xm,
                                         const realT d,
                                         const rArray1d& c)
  {
    if (c.size() < 11)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    //      Setting some parameters for the specific case of the ANPP
    //Considering that the most important faults are nearby (SaF and NWF)
    //and that we do not know the rake for NWF and Saf is SS, then
    realT frv = 0.0;
    //Same considerations for Fnm
    realT fnm = 0.0;
    //The most damaging EQ from SaF or NWF would be for M>=6, for whi
    //the rupture plane would probably cut to the surface. then ZTOR
    realT ztor = 0.0;
    realT ffltz = 0.0;
    if (ztor >= 1.0)
    {
      ffltz = 1.0;
    }
    //Delta i sht edip of the rupture plane, in degrees.
    realT delta = 65.0;
    //Condidering that the faults are mostly SS vertical, Rjb = Rrup
    realT rjb = d;
    //
    //Fmag------------------------------------------------------
    realT fmag = c(1) + c(2) * xm;
    if (xm > 5.5)
    {
      fmag += c(3) * (xm - 5.5);
      if (xm > 6.5)
      {
        fmag += c(4) * (xm - 6.5);
      }
    }
    //
    //Fdis------------------------------------------------------
    realT fdis = (c(5) + c(6) * xm) * 0.5 * log(d * d + c(7) * c(7));
    //
    //Fflt------------------------------------------------------
    realT fflt = c(8) * frv * ffltz + c(9) * fnm;
    //
    //Fhng------------------------------------------------------
    realT fhngr = 1.0;
    if (rjb > 0.0)
    {
      const realT tval = sqrt(rjb * rjb + 1.0);
      fhngr = d > tval - rjb ? d : tval - rjb;
      fhngr /= d > tval ? d : tval;
    }
    else
    {
      fhngr = (d - rjb) / d;
    }
    //
    realT fhngm = .0;
    if (xm > 6.0)
    {
      fhngm = 2.0 * (xm - 6.0);
    }
    if (xm > 6.5)
    {
      fhngm = 1.0;
    }
    //
    realT fhngz = 0.0;
    if (ztor < 20.0)
    {
      fhngz = (20.0 - ztor) / 20.0;
    }
    //
    realT fhngd = 1.0;
    if (delta > 70.0)
    {
      fhngd = (90.0 - delta) / 20.0;
    }
    //
    realT fhng = c(10) * fhngr * fhngm * fhngz * fhngd;
    //
    //Fsite and Fsed are not considered since site specific studies will be done
    //for ANPP.
    //
    //Natural Log of the ground motion
    return_value = fmag + fdis + fflt + fhng;
    //Convert from G's to cm/s/s
    return_value += 6.889;
    //
    return return_value;
    //
  }

  /**
   *---------------------------------------
   *
   * model xdabc   index=14
   *
   * model in PG&E LTSP rpt   more generally sadigh's model
   *
   *   c(8)   plays the role of fault type
   *
   *  c(20) plays the role of a dummy depth to compare to J&B model
   */
  realT
  GroundMotion1D::PGE_LTSP_Model(const realT xm,
                                 const realT d,
                                 const rArray1d& c,
                                 const int ll)
  {
    int l = ll;
    if (c.size() < 21)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    //****************************************************************
    //
    //     added next 2 lines   7/30 /91
    if (xm > c(10))
    {
      l = 10;
    }
    realT r = sqrt(d * d + c(20) * c(20));
    //
    //******************* removed next 2 lines 7/30/91 **************
    //     r=d
    //     if(d.le.15.)r=dsqrt(d*d+25.)
    //
    //****************************************************************
    return_value = c(1 + l) + c(2 + l) * xm + c(3 + l) * pow((8.5 - xm), c(4 + l)) + c(5 + l) * log(
        r + c(6 + l) * exp(c(7 + l) * xm)) + c(8);
    return return_value;
  }



  /**
   *---------------------------------------
   *
   *   model xidrs   index= 13
   *
   * added 5/21/90
   *idriss's 1987 model  fron j&b 1988 summary paper
   */
  realT
  GroundMotion1D::IdrissModel(const realT xm,
                              const realT dd,
                              const rArray1d& c)
  {
    realT d = dd;
    realT return_value = 0.0;
    realT r = d;
    //   note next line added to account for area zone around the site
    //   more generally may want to delete this line
    if (d <= 15.0)
    {
      r = sqrt(d * d + 25.0);
    }
    realT xm2 = xm * xm;
    realT xm3 = xm2 * xm;
    realT xm4 = xm3 * xm;
    //
    realT xlna = c(1) + c(2) * xm + c(3) * xm2 + c(4) * xm3 + c(5) * xm4;
    d = c(11) + c(12) * xm + c(13) * xm2 + c(14) * xm3 + c(15) * xm4;
    return_value = xlna + d * log(r + 20.0);
    return return_value;
  }

  /**
   *---------------------------------------
   *
   * model  xjb82    index  =  11
   *
   *   added 5/21/90
   *  joyner and boore's 1982 model
   *
   *  c(5)  plays the role of the h parameter in their model
   */
  realT
  GroundMotion1D::XJB82Model(const realT xmm,
                             const realT d,
                             const rArray1d& c)
  {
    if (c.size() < 14)
      throw GPException("c array too short");

    realT return_value = 0.0;
    realT r = sqrt(d * d + (c(5) * c(5)));
    //
    // if c(11) .ne. 0. then using atkinson's model  so must convert xm to
    // moment magintude
    realT xm = xmm;
    if (c(11) > 0.0)
    {
      xm = c(11) + c(12) * xm + c(13) * xm * xm;
    }
    //
    return_value = c(1) + c(2) * (xm - 6.0) + c(3) * ((xm - 6.0) * (xm - 6.0)) + c(4) * log(r) + c(6) * r;
    //
    //    Case for very small Earthquakes scaling for M<6
    if (xmm < 4.5)
    {
      return_value = c(1) - c(2) * 1.5 + c(4) * log(r);
      //        Magnitude scaling
      return_value = return_value - 1.7 * (4.5 - xmm);
    }
    //
    //    At this point the model is for use parameters with NATURAL LOG
    //        convert to base e
    //     xjb82=xjb82*2.3026
    return return_value;
  }



  /**
   *---------------------------------------
   *
   *    added 5/21 /90
   */
  realT
  GroundMotion1D::Campbell1989Model(const realT xm,
                                    const realT d,
                                    const rArray1d& c)
  {
    if (c.size() < 19)
      throw GPException("c array too short");

    realT return_value = 0.0;

    //ken campbell's 1989  models
    //
    //  c(8) plays the role of fault type
    // c(16)  is depth to bedrock
    //  c(20) plays the role of a dummy depth to compare to the J&B model
    //
    // c(17) plays the role of rupture direction  in some of campbell's
    //  earlier models    c(18) = coef of rupture direction
    //
    //     r=dsqrt(d*d+c(20)*c(20))
    realT r = d;
    //
    //   note next line added to account for area zone around the site
    //   more generally may want to delete this line
    if (d <= 15.0)
    {
      r = sqrt(d * d + 25.0);
    }
    return_value = c(1) + c(2) * xm + c(3) * log(r + c(4) * exp(c(5) * xm)) + c(6) * r + c(7) * c(8) + c(
        17) * c(18);
    if (c(11) > 0)
    {
      return_value += c(11) * tanh(c(12) * (xm + c(13))) + c(14) * tanh(c(15) * c(16));
    }
    return return_value;
  }


  /**
   *---------------------------------------
   *
   * model  xmod15    index=15
   *  generic model of the form
   * ln(a)=c1+c2*m+c3*m**2 +c4*ln(r) +c5*r
   *
   *  in particular form used for my fit to WCC model for the NPR at INEL
   */
  realT
  GroundMotion1D::WCC_NPR_Model(const realT xm,
                                const realT d,
                                const rArray1d& c)
  {
    if (c.size() < 6)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    realT r = d;
    //
    //    this line added to take care of adding a depth term for the two
    //     area zones  in which the INEL NPR site sits
    //     for more general use remove
    //
    if (r <= 15.0)
    {
      r = sqrt(r * r + 25.0);
    }
    return_value = c(1) + c(2) * xm + c(3) * (xm * xm) + c(4) * log(r) + c(5) * r;
    return return_value;
  }

  /**
    *---------------------------------------
    */
  realT
  GroundMotion1D::Romanian1997Model(const realT xm,
                                    const realT d,
                                    const rArray1d& c)
  {
    if(c.size() < 5)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    //---------------------------------------
    //
    //       Model type:  Romanian study, 1997
    //       Index = 19
    //       For use in dealing with the Vrancea zone.
    //
    //       Ln(y) = c1 + c2.M + c3.Ln(Rh+c4)
    //       Natural logs
    //       cm/s/s, Local/Ms magnitudes, Hypocentral distance km
    //
    realT rh = sqrt((d * d) + (c(5) * c(5)));
    return_value = c(1) + c(2) * xm + c(3) * log(rh + c(4));
    return return_value;
  }


  /**
    *---------------------------------------
    */
  realT
  GroundMotion1D::SabettaPuglieseModel(const realT xm,
                                       const realT d,
                                       const rArray1d& c)
  {
    if (c.size() < 4)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    //       Model type: Sabetta and Pugliese, 1987
    //       Index = 18
    //       For use in generic environment. From Italian data
    //
    realT r = (d * d) + (c(3) * c(3));
    r = log(r) / 2.0;
    return_value = c(1) + c(2) * xm - r;
    return return_value;
  }


  /**
   *---------------------------------------
   *
   *    Model  xSadighRock    index  =  21
   *
   *    added 5/1/2001
   *       Sadigh, Youngs et al. Rock model: Seismological Research Letters
   *       January/February 1997
   */
  realT
  GroundMotion1D::SadighYoungsRockModel(const realT xm,
                                        const realT d,
                                        const rArray1d& c)
  {
    if (c.size() < 18)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    //Shift in attn indeces for m>6.5
    int k = 0;
    if (xm > 6.5)
    {
      k = 10;
    }
    //
    realT xkm = 8.5 - xm;
    if (xm > 8.5)
    {
      xkm = 0.0;
    }
    realT y = c(1 + k) + (c(2 + k) * xm) + c(3 + k) * xkm * xkm * sqrt(xkm) + c(4 + k) * log(
        d + exp(c(5 + k) + c(6 + k) * xm)) + c(7 + k) * log(d + 2.0);
    //
    //       Transform g's into cm/s/s
    return_value = y + 6.889;
    return return_value;
  }



  /**
   *---------------------------------------
   *
   *    Model  xSadighSoil    index  =  22
   *
   *    added 5/1/2001
   *       Sadigh, Youngs et al. Soil model: Seismological Research Letters
   *       January/February 1997
   */
  realT
  GroundMotion1D::SadighYoungsSoilModel(const realT xm,
                                        const realT d,
                                        const rArray1d& c)
  {
    if(c.size() < 8)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    realT x4 = 2.1863;
    realT x5 = 0.32;
    if (xm > 6.5) {
      x4 = 0.3825;
      x5 = 0.5882;
    }
    //
    realT xkm = 8.5 - xm;
    if (xkm < 0.0) {
      xkm = 0.0;
    }
    realT y = c(1) + c(2) * xm - c(3) * log(d + (x4 * exp(
      x5 * xm))) + c(6) + c(7) * xkm * xkm * sqrt(xkm);
    //
    //       Transform g's into cm/s/s
    return_value = y + 6.889;
    return return_value;
  }


  /**
   *
   *
   *      model xtip98  index = 16
   *
   *      added july 28,1998, by JBS
   *
   *      Model developed for the TIP project, by N. Abrahamson
   *       c(19) is the model type (=16)
   *       c(10) is the magnitude threshold (6.25)
   *       c(11), c(12) and c(13) are for the Mblg to Mw
   *             conversion
   *       Values of PGA are in cm/s/s
   *       Values of PGV are in cm/s
   *       Log are Natural Logs
   */
  realT
  GroundMotion1D::AbrahamsonTIPModel(const realT xmm,
                                     const realT d,
                                     const rArray1d& c)
  {
    if (c.size() < 16)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    realT r = sqrt(d * d + c(15) * c(15));
    realT xm = xmm;
    //       convert Mblg to Mw if c(11) non-zero
    if (c(11) > 0.0)
    {
      xm = c(11) + xm * c(12) + c(13) * xm * xm;
    }
    realT ccm = c(2);
    //       ccm is the coeff of (M-m1). Change for xm above
    //       magnitude threshold
    if (xm > c(10))
    {
      ccm = c(14);
    }
    //
    realT xmc = xm - c(10);
    return_value = c(1) + ccm * xmc + c(3) * ((8.5 - xm) * (8.5 - xm)) + (c(4) + c(5) * xmc) * log(
        r) + c(6);
    //
    return return_value;
  }




  //
  //---------------------------------------
  realT
  GroundMotion1D::Tore1997Model(const realT xm,
                                const realT d,
                                const rArray1d& c)
  {
    if (c.size() < 8)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    //       Mode; type Tore 97 for Midcontinents, Moment Magnitude
    //       Index = 17
    //
    const realT xm6 = xm - 6.0;
    const realT rm = sqrt(d * d + (c(7) * c(7)));
    const realT rx = log(rm);
    realT ry = log(rm / 100.0);
    if (ry <= 0.0)
      ry = 0.0;
    return_value = c(1) + (c(2) + c(3) * xm6) * xm6 - c(4) * rx - (c(4) - c(5)) * ry - c(6) * rm;
    return return_value;
  }

  //
  //---------------------------------------
  realT
  GroundMotion1D::Youngs1997Model(const realT xm,
                                  const realT d,
                                  const rArray1d& c)
  {
    if (c.size() < 16)
      throw GPException("c array too short");

    realT return_value = 0.0;
    //
    //       Model type Youngs SRL 1997 for Subduction zones
    //               xm = M = Moment Magnitude
    //               d  = distance to rupture (km)
    //               c  = coefficients
    //               c(6) = H = dominant depth
    //               c(7) = Zt = 0 for Interface
    //                                       1 for Intraslab
    //
    //       Model Index = 24
    //
    const realT xh = c(6);
    const realT zt = c(7);
    realT p3 = (10.0 - xm);
    p3 *= p3 * p3;
    return_value = c(8) + c(11) * xm + c(1) + c(2) * p3;
    return_value += c(3) * log(d + c(12) * exp(c(13) * xm));
    return_value += c(14) * xh + c(15) * zt;
    //
    return return_value;
  }

} /* namespace EarthquakeSimulation */
