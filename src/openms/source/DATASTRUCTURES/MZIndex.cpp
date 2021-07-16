// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Utils/MZIndex.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak2D.h>

namespace OpenMS
{
  MzIndex::MzIndex(const MSExperiment& exp, Size nr_bins):
      exp_(std::make_shared<MSExperiment>(exp))
  {
    if (!exp.isSorted())
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION, "MSExperiment not sorted. Aborting.");
    }
    index_.resize(exp.size());
    for (auto& spec_bins : index_)
    {
      spec_bins.resize(nr_bins);
    }
    double maxmz = exp.getMaxMZ();
    double minmz = exp.getMinMZ();
    // add to maxmz, so the last peak falls into the last bin and not after
    bin_width_ = (maxmz - minmz) / double(nr_bins);
    Size s = 0;
    // TODO use OpenMP since everything is preallocated
    for (const auto& spec : exp)
    {
      std::vector<Size>& bin_limit_idcs = index_[s];
      bin_limit_idcs[0] = 0; //TODO probably possible without this bin
      Size b = 1;
      Size n = 0; //start at one, so bin end points to one-after last included index
      double currmz = minmz + bin_width_;
      for (const auto& peak : spec)
      {
        //we can assume that every peak mz is unique, so we can open a bin after the first equals encounter
        if (peak.getMZ() > currmz+0.00001)
        {
          Size bemptybins = b+Size((peak.getMZ()-(currmz+0.00001)) / bin_width_);
          for (; b < bemptybins; b++)
          {
            bin_limit_idcs[b] = n;
            currmz += bin_width_;
          }
          bin_limit_idcs[b] = n;
          currmz += bin_width_;
          ++b;
        }
        ++n;
      }

      // fill rest with last peak idx
      for (Size m = b; m < nr_bins; ++m)
      {
        bin_limit_idcs[m] = n;
      }
      ++s;
    }
  }

  Size MzIndex::findNearest(Size spec_idx, double mz)
  {
    double diff = mz - exp_->getMinMZ();
    if (diff < 0) return 0;
    Size bin = std::trunc(diff / bin_width_);

    if (bin >= index_[spec_idx].size()) return exp_->getSpectrum(spec_idx).size() - 1;
    //TODO size should be allowed as endidx. begin + size = (one-past)end. This should be allowed for lower_bound
    Size start = index_[spec_idx][bin];
    if (start > 1)
    {
      start--;
    }
    //std::min(exp_->getSpectrum(spec_idx).size(), index_[spec_idx][bin+1]+3)
    return exp_->getSpectrum(spec_idx).findNearestInRange(mz, start, index_[spec_idx][bin+1]);
  }
}
