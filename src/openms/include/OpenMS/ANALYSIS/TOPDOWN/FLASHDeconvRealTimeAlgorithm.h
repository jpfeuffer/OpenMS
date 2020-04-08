// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <queue>

namespace OpenMS
{
  class OPENMS_DLLAPI FLASHDeconvRealTimeAlgorithm
  {/*
    struct OPENMS_DLLAPI MzRange
    {
      double mz1, mz2;
      int nominalMass;
      double rt;
      MzRange();
      ~MzRange();
      explicit MzRange(double mz1, double mz2, double rt, int mass);
    };

    struct OPENMS_DLLAPI ShortMassTrace
    {
      double lastRt;
      std::queue<double> intensities; // keep the length to 3 // last rt -rt >30 -> reset
      char color; // 'r' ''b' 'g' 'x'

      ShortMassTrace();
      ~ShortMassTrace();
      explicit ShortMassTrace(double rt);

      void setColor(char c);
      void updateColor();
      void addIntensity(double intensity);
      void trim();
    };

  public:
    typedef FLASHDeconvHelperStructs::Parameter Parameter;
    typedef FLASHDeconvHelperStructs::PeakGroup PeakGroup;
    typedef FLASHDeconvHelperStructs::PrecalcularedAveragine PrecalcularedAveragine;
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    typedef FLASHDeconvRealTimeAlgorithm::MzRange mzRange;
    typedef FLASHDeconvRealTimeAlgorithm::ShortMassTrace ShortMassTrace;

    FLASHDeconvRealTimeAlgorithm(MSSpectrum &s, std::vector<int> greenMasses, Parameter &p);
    ~FLASHDeconvRealTimeAlgorithm();

    std::vector<MzRange>& Deconvolution(FLASHDeconvHelperStructs::PrecalcularedAveragine &avg);

  protected:
    MSSpectrum &spec;
    Parameter &param;

  private:
    //std::queue<mzRange> mzQueue;
    bool compareIntensity(const FLASHDeconvHelperStructs::PeakGroup &pg1,
                          const FLASHDeconvHelperStructs::PeakGroup &pg2);

    static std::map<int, ShortMassTrace> massMemory; // static TODO
    double deltaRtForMassMemory;
    double deltaRtForQueue;
*/
  };
}