//
// Created by Kyowon Jeong on 4/22/20.
//


#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>

namespace OpenMS
{
  class DeconvolutedSpectrum;
  class OPENMS_DLLAPI PeakGroup
  {
  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;
    typedef FLASHDeconvHelperStructs::PrecalcularedAveragine PrecalcularedAveragine;

    DeconvolutedSpectrum *deconvSpec;
    std::vector<LogMzPeak> peaks;
    double monoisotopicMass = .0;
    double avgMass = .0;
    double intensity = .0;
    Size massBinIndex = 0;

    float isotopeCosineScore = .0;
    float chargeCosineScore = .0;

    int massIndex;
    int maxCharge, minCharge;
    int maxSNRcharge = 0;
    std::unordered_map<int, float> perChargeSNR;
    float maxSNR = 0;
    float totalSNR = 0;
    double maxSNRmaxMz, maxSNRminMz;

    ~PeakGroup();

    void push_back(LogMzPeak &p);

    void reserve(Size n);

    void emptyChargeSNRMap();

    bool operator<(const PeakGroup &a) const;

    bool operator>(const PeakGroup &a) const;

    bool operator==(const PeakGroup &a) const;

    void updateMassesAndIntensity(PrecalcularedAveragine &averagines, int offset = 0, int maxIsoIndex = 0);
  };
}