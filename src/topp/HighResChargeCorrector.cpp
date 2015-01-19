// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/MISC/LinearInterpolation.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <boost/lambda/lambda.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <math.h>
#include <unsupported/Eigen/FFT>

using namespace OpenMS;
using namespace std;

/**
 @page UTILS_HighResChargeCorrector HighResChargeCorrector
 
 @brief Corrects the charge of high resolution data.
*/

class TOPPHighResChargeCorrector :
public TOPPBase
{
public:
  TOPPHighResChargeCorrector() :
  TOPPBase("HighResChargeCorrector", "Corrects the precursor charge determined by the instrument software.")
  {
  }
  
protected:
  void registerOptionsAndFlags_()
  {
    // input files
    registerInputFile_("in", "<file>", "", "input file (centroided data)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file");
    setValidFormats_("out", ListUtils::create<String>("mzML"));
    registerOutputFile_("out_csv", "<file>", "", "Optional csv output file containing scores for different precursor charges\n", false);
    setValidFormats_("out_csv", ListUtils::create<String>("csv"));
  }
  
  void getPrecursors_(const PeakMap & exp, vector<Precursor> & precursors, vector<double> & precursors_rt, vector<OpenMS::String> & precursors_id)
  {
    for (Size i = 0; i != exp.size(); ++i)
    {
      vector<Precursor> pcs = exp[i].getPrecursors();
      if (pcs.empty())
      {
        continue;
      }
      vector<double> pcs_rt(pcs.size(), exp[i].getRT());
      vector<OpenMS::String> pcs_id(pcs.size(), exp[i].getNativeID());
      
      //for (Precursor c : pcs)
      //  std::cout << c.getMZ() << ' ';
      //cout << endl;
      
      copy(pcs.begin(), pcs.end(), back_inserter(precursors));
      copy(pcs_rt.begin(), pcs_rt.end(), back_inserter(precursors_rt));
      copy(pcs_id.begin(), pcs_id.end(), back_inserter(precursors_id));
    }
  }
  
  String printRanks(const vector<double> & ranks){
    String result = "";
    for (vector<double>::const_iterator i = ranks.begin(); i != ranks.end(); ++i){
      result += *i;
      result += ", ";
    }
    
    return result.prefix(result.length()-2);
  }
  
  void writeHist(String out_csv, const vector<int> & ids, const vector<int> & oldchgs, const vector<int> & newchgs, const vector<vector<double> > & fourranks, const vector<vector<double> > & pattranks, const vector<vector<double> > & combranks, const vector<double> & rts, const vector<double> & itys, const vector<double> & mzs)
  {
    //cout << "writing data" << endl;
    ofstream csv_file(out_csv.c_str());
    csv_file << setprecision(9);
    
    // header
    csv_file << "ScanID\tuncorrectedChg\tcorrChg\tFourRanks\tPattRanks\tCombRanks\tRT\tIntensity\tMZ" << endl;
    
    // entries
    for (vector<int>::const_iterator it = ids.begin(); it != ids.end(); ++it)
    {
      UInt index = it - ids.begin();
      String f = printRanks(fourranks.at(index));
      String p = printRanks(pattranks.at(index));
      String c = printRanks(combranks.at(index));
      csv_file << *it << "\t" << oldchgs[index] << "\t" << newchgs[index] << "\t" << f << "\t" << p << "\t" << c << "\t" << rts[index] << "\t" << itys[index] << "\t" << mzs[index] << endl;
    }
    csv_file.close();
  }
  
protected:
  void correct(PeakMap & exp, vector<int> & ids, vector<int> & oldcharge, vector<int> & newcharge, vector<vector<double> > & fourrank, vector<vector<double> > & pattrank, vector<vector<double> > & combrank, vector<double> & rts, vector<double> & maxint, vector<double> & mzs)
  {
    // load experiment and extract precursors
    vector<Precursor> precursors;  // precursor
    vector<double> precursors_rt;  // RT of precursor MS2 spectrum
    vector<OpenMS::String> precursors_id;  // Scan Number of precursor MS2 spectrum
    getPrecursors_(exp, precursors, precursors_rt, precursors_id);
    
    for (Size i = 0; i != precursors_rt.size(); ++i)
    {
      // get precursor rt
      double rt = precursors_rt[i];
      
      // get scan number of MS2
      vector<OpenMS::String> splits;
      precursors_id[i].split(String("scan="), splits);
      int id = splits[1].toInt();
      
      // get precursor MZ
      double mz = precursors[i].getMZ();
      
      // get precursor charge (assigned by instrument)
      double chg = precursors[i].getCharge();
      
      double inty = precursors[i].getIntensity();

      
      //cout << rt << " " << mz << " " << chg << endl;
      //cout << id << " " << chg;
      
      ids.push_back(id);
      oldcharge.push_back(chg);
      rts.push_back(rt);
      maxint.push_back(inty);
      mzs.push_back(mz);
      
      
      
      
      // get precursor spectrum
      MSExperiment<Peak1D>::ConstIterator rt_it = exp.RTBegin(rt);
      
      // store index of MS2 spectrum
      UInt precursor_spectrum_idx = rt_it - exp.begin();
      
      // get parent (MS1) of precursor spectrum
      rt_it = exp.getPrecursorSpectrum(rt_it);
      
      
      /*/########## TEMP TEST EXPERIMENT############
      MSExperiment<> exper;
      std::vector<MSSpectrum<Peak1D> > spectra;
      MSSpectrum<> spectrum;
      Peak1D peak;
      double tempint = 1000.0;
      for (double mz = 600.0; mz <= 608.0; mz += 1.0)
      {
        peak.setMZ(mz);
        peak.setIntensity(tempint/(mz-599.0));
        spectrum.push_back(peak);
      }
      spectrum.setRT(1.0);
      spectra.push_back(spectrum);
      exper.setSpectra(spectra);
      mz = 600.0;
      rt_it=exper.RTBegin(1.0);
      //########################################*/
      
      if (rt_it->getMSLevel() != 1)
      {
        cout << "Error: no MS1 spectrum for this precursor" << endl;
      }
      

      //get window around precursor
      double startmz = rt_it->front().getMZ();
      double endmz = rt_it->back().getMZ();
      
      double rel_int = 0.001; //make it a specifiable value
      double window_startmz = mz - mz * rel_int;
      double window_endmz = mz + mz * rel_int;
      
      // Create a LIP object
      OpenMS::Math::LinearInterpolation<double,double> lip;
      
      double spacing = 1.0/100.0; //make it a specifiable value
      lip.getData().resize((window_endmz-window_startmz)/spacing);
      lip.setScale(spacing);
      lip.setOffset(window_startmz);
      
      
      MSSpectrum<Peak1D>::ConstIterator it;
      
      //get max peak of window (for normalization)
      double max(0);
      for (it = rt_it->MZBegin(window_startmz); it != rt_it->MZEnd(window_endmz); ++it)
      {
        if (max < it->getIntensity())
        {
          max = it->getIntensity();
        }
      }

      //Fill LIP Object
      for (it = rt_it->MZBegin(window_startmz); it != rt_it->MZEnd(window_endmz);++it)
      {
        lip.addValue(it->getMZ(), it->getIntensity()*100/max);
      }
      
      //Get interpolated values with constant spacing
      int maxstep = (window_endmz - window_startmz)/spacing;
      vector<float> timevec;
      
      
      double currmz;
      for (int nstep = 0; nstep <= maxstep; ++nstep)
      {
        currmz = (nstep * spacing) + startmz;
        timevec.push_back(lip.value(currmz));
      }

      //Up to which charge should the calculations go (reduces time only for PattTransform)
      int maxchg = 5;
      vector<double> patt_charge_map = pattersonTransform(timevec, maxchg, spacing);
      vector<double> four_charge_map = fourierTransform(timevec, maxchg, spacing);
      vector<double> comb_charge_map(patt_charge_map.size());
      
      std::transform(patt_charge_map.begin(), patt_charge_map.end(), four_charge_map.begin(), comb_charge_map.begin(), std::multiplies<float>());
      
      int pred_chg = round((std::distance(comb_charge_map.begin(), std::max_element(comb_charge_map.begin(), comb_charge_map.end()))+1));

      
      newcharge.push_back(pred_chg);
      
      
      //fourrank.push_back(getRankFor(four_charge_map,chg));
      //pattrank.push_back(getRankFor(patt_charge_map,chg));
      //combrank.push_back(getRankFor(comb_charge_map,chg));
      
      fourrank.push_back(four_charge_map);
      pattrank.push_back(patt_charge_map);
      combrank.push_back(comb_charge_map);
      //cout << rt_it->getRT() << " " << rt_it->size() << endl;

    }
  }
  
//  vector<double> pattersonTransform(MSExperiment<Peak1D>::ConstIterator rt_it, double mz, double rel_int, double thresh)
//  {
//    MSSpectrum<Peak1D>::ConstIterator it;
//    double startmz = mz - mz * rel_int;
//    double endmz = mz + mz * rel_int;
//    
//    double lastmz;
//    double nextmz;
//    double lastint;
//    double patt_val = 0;
//    
//    for (double chg = 2/3; chg <= 8 + 1/3; chg += 1/3)
//    {
//      it = rt_it->MZBegin(startmz);
//      lastmz = it->getMZ();
//      lastint = it->getIntensity();
//      ++it;
//      
//      for (it = it; it != rt_it->MZEnd(endmz);)
//      {
//        nextmz = lastmz + 1/chg;
//        cout << it->getMZ()<< " ";
//        if(std::abs(it->getMZ() - nextmz) < nextmz * thresh/1000000)
//        {
//          patt_val += lastint * it->getIntensity();
//          lastmz = it->getMZ();
//        }
//        ++it;
//      }
//      
//      cout << endl;
//    }
//  }
  
  vector<double> pattersonTransform(vector<float> &timevec, int maxchg, double spacing)
  {
    double patt_val = 0.0;
    vector<double> chargemap; //result
    unsigned int step;
    
    for (int chg = 1; chg <= maxchg; ++chg)
    {
      step = static_cast<unsigned int>((1/spacing)/chg); // works well for the usual charge states
      
      //cout << step << endl;
      for (unsigned i = 0; i < timevec.size()-step; ++i)
      {
        patt_val = patt_val + (timevec[i] * timevec[i+step]);
      }
      //cout << patt_val << endl;
      chargemap.push_back(patt_val);
      patt_val = 0.0;
      
    }
    return chargemap;
    
  }
  
  
//  vector<double> pattersonTransform(MSExperiment<Peak1D>::ConstIterator rt_it, double mz, double rel_int, OpenMS::Math::LinearInterpolation<double,double> lip)
//  {
//    
//    double startmz = mz - mz * rel_int;
//    double endmz = mz + mz * rel_int;
//    double patt_val = 0.0;
//    int maxchg = 5;
//    
//    vector<double> chargemap;
//    
//    double step;
//    
//    MSSpectrum<Peak1D>::ConstIterator it;
//    
//    for (int chg = 1; chg <= maxchg; ++chg)
//    {
//      step = chg;
//      //cout << step << endl;
//      for (it = rt_it->MZBegin(startmz); it != rt_it->MZEnd(endmz);++it)
//      //for (int currmz = startmz; currmz <= endmz; ++currmz)
//      {
//        patt_val = patt_val + lip.value(/*currmz*/it->getMZ() /*- 1/step/2*/) * lip.value(/*currmz*/it->getMZ() + 1/step);
//      }
//      //cout << patt_val << endl;
//      chargemap.push_back(patt_val);
//      patt_val = 0.0;
//    }
//    return chargemap;
//    
//  }

  
  vector<double> fourierTransform(vector<float> &timevec, int maxchg, double spacing)
  {
    
    //zero pad to next power of 2 for dft speed
    unsigned length = timevec.size();
    double y = floor(log2(length));
    
    unsigned nextpow2 = (unsigned int)pow(2, y + 1);
    vector<float> zeropadtimevec = timevec;
    zeropadtimevec.resize(nextpow2, 0.0);
    
    Eigen::FFT<float> fft;
    vector<complex<float> > freqvec;
    fft.fwd(freqvec, zeropadtimevec);
    
    vector<double> chargemap;
    int idx;
    double chg;
    for (int chgstp = 1; chgstp <= maxchg; ++chgstp)
    {
      chg = chgstp;
      idx = chg * nextpow2/*length*/ * spacing;
      chargemap.push_back(abs(freqvec[idx].real()));
    }
    
    return chargemap;
    
  }

//  vector<double> fourierTransform(MSExperiment<Peak1D>::ConstIterator rt_it, double mz, double rel_int, OpenMS::Math::LinearInterpolation<double,double> lip)
//  {
//    
//    double startmz = mz - mz * rel_int;
//    double endmz = mz + mz * rel_int;
//    int maxchg = 5;
//    //double step = 1.0/(maxchg*2);
//    double step = 1.0/128;
//    double currmz = 0.0;
//    
//    
//    vector<float> timevec;
//    
//    int maxstep = (endmz - startmz)/step;
//    
//    for (int nstep = 0; nstep <= maxstep; ++nstep)
//    {
//      currmz = (nstep * step) + startmz;
//      timevec.push_back(lip.value(currmz));
//    }
//    
//    
//    Eigen::FFT<float> fft;
//    vector<complex<float> > freqvec;
//    fft.fwd(freqvec, timevec);
//    
//    vector<double> chargemap;
//    int idx;
//    double chg;
//    for (int chgstp = 1; chgstp <= maxchg; ++chgstp)
//    {
//      chg = chgstp;
//      idx = chg * freqvec.size() * step;
//      chargemap.push_back(abs(freqvec[idx].real()));
//    }
//
//    return chargemap;
//    
//  }

  
  vector<int> getChargeRanks(std::vector<double> intVector)
  {
  
  //#################   Temporary Ranks   ##########################
  ;
  std::vector<int> rank;
  
  
  for(Size i = 1; i < intVector.size(); ++i )
  {
    rank.push_back( i );
  }
  
  using namespace boost::lambda;
  std::sort(
            rank.begin(), rank.end(),
            var( intVector )[ _1 ] < var( intVector )[ _2 ]
            );
  

  return rank;
  
  //#################################################
  }
  
  int getRankFor(std::vector<double> intVector, int chg)
  {
    
    //#################   Temporary Ranks   ##########################
    ;
    std::vector<int> rank;
    
    
    for( Size i = 1; i < intVector.size(); ++i )
    {
      rank.push_back( i );
    }
    
    using namespace boost::lambda;
    std::sort(
              rank.begin(), rank.end(),
              var( intVector )[ _1 ] > var( intVector )[ _2 ]
              );
    
    int index = 0;             //Index
    std::vector<int>::iterator it;
    it = find(rank.begin(), rank.end(), chg-1);
    
    if(it != rank.end())
      index = it-rank.begin();
    
    return index+1;
    
    //#################################################
  }
  
  ExitCodes main_(int, const char **)
  {
    const string in_mzml(getStringOption_("in"));
    const string out_mzml(getStringOption_("out"));
    const string out_csv = getStringOption_("out_csv");
    
    PeakMap exp;
    MzMLFile().load(in_mzml, exp);
    
    cout << setprecision(12);
    
    // determine accuracy
    vector<int> ids;
    vector<int> oldchgs;
    vector<int> newchgs;
    vector<vector<double> > fourranks;
    vector<vector<double> > pattranks;
    vector<vector<double> > combranks;
    vector<double> rts;
    vector<double> itys;
    vector<double> mzs;
    
    correct(exp, ids, oldchgs, newchgs, fourranks, pattranks, combranks, rts, itys, mzs);
    
    MzMLFile().store(out_mzml, exp);
    
    if (out_csv != "")
    {
      writeHist(out_csv, ids, oldchgs, newchgs, fourranks, pattranks, combranks, rts, itys, mzs  );
    }
    
    return EXECUTION_OK;
  }
  
};

int main(int argc, const char ** argv)
{
  TOPPHighResChargeCorrector tool;
  return tool.main(argc, argv);
}

/// @endcond
