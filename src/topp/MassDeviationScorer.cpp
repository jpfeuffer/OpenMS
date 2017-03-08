// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <boost/math/distributions/normal.hpp>
#include <OpenMS/CONCEPT/Constants.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_MassDeviationScorer MassDeviationScorer

 @brief Uses deviation from theoretical precursor masses to calculate another probabilistic score for PSMs.
<CENTER>
 <table>
  <tr>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tool </td>
   <td VALIGN="middle" ROWSPAN=5> \f$ \longrightarrow \f$ MassDeviationScorer \f$ \longrightarrow \f$</td>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MapRTTransformer (with trafoXML from InternalCalibration) </td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideRescorer (or TOPPPerc) </td>
  </tr>
 </table>
</CENTER>

 Use InternalCalibration before to correct for the influence of retention time.

 @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_MassDeviationScorer.cli
 <B>INI file documentation of this tool:</B>
 @htmlinclude TOPP_MassDeviationScorer.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPMassDeviationScorer :
    public TOPPBase
{
public:
  TOPPMassDeviationScorer() :
      TOPPBase("MassDeviationScorer", "Uses deviation from theoretical precursor masses to calculate another probabilistic score for PSMs.")
  {

  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerFlag_("use_all_psms", "If enabled, uses not only the best matching PSM per spectra but all.");
    registerFlag_("combine_isotope_errors", "If enabled, the tool does not use different distributions for"
        " deviations of PSMs per isotope error (if enabled in the search).");
  }

  ExitCodes main_(int, const char**)
  {
    // Reading
    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
    IdXMLFile().load(inputfile_name, proteins, peptides);

    //TODO generalize for multiple ProteinID objects
    double tolerance = proteins[0].getSearchParameters().precursor_mass_tolerance;

    //E.g. MSGF+ searches for isotope errors and that needs to be corrected for here
    int isotopeErrorRange = 1;
    if(proteins[0].getSearchParameters().metaValueExists("MaxIsotopeError"))
    {
      isotopeErrorRange = String(proteins[0].getSearchParameters().getMetaValue("MaxIsotopeError")).toInt() -
                          String(proteins[0].getSearchParameters().getMetaValue("MinIsotopeError")).toInt() + 1;
    }

    //We treat deviations different isotopes seperately
    vector<RunningStatistic> stats(isotopeErrorRange);

    //TODO Alternatives: Fitting, Gauss around mean 0.0 ppm


    // Adding information to datastructures
    bool all_psms = getFlag_("use_all_psms");
    bool combine_isotope_errors = getFlag_("combine_isotope_errors");
    std::cout << "considering all psms?" << all_psms << std::endl;

    double obs_mz = 0.0;
    for (vector<OpenMS::PeptideIdentification>::iterator pep_it = peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      obs_mz = pep_it->getMZ();
      vector<PeptideHit>& hits = pep_it->getHits();
      processHit_(hits[0], stats, obs_mz, combine_isotope_errors);

      if (all_psms)
      {
        for (Size i = 1; i < hits.size(); ++i)
        {
          processHit_(hits[i], stats, obs_mz, combine_isotope_errors);
        }
      }
    }

    // Statistics for isotope error = 0, or all if combined
    std::cout << "Mean error: " << stats[0].getMean() << std::endl;
    std::cout << "SD error: " << stats[0].getSD() << std::endl;

    // Fitting/Estimating
    vector<boost::math::normal_distribution<double> > error_dists(isotopeErrorRange);

    for (int j = 0; j < stats.size(); ++j)
    {
     error_dists[j] = boost::math::normal_distribution<double>( stats[j].getMean(), stats[j].getSD() );
    }

    // Writing
    std::ofstream tmp_outfile ("devs_all_relppm_transformed.txt", std::ofstream::out | std::ofstream::trunc);
    for (vector<OpenMS::PeptideIdentification>::iterator pep_it = peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      //std::cout << "Peptide: " << std::endl;
      obs_mz = pep_it->getMZ();
      //std::cout << obs_mz << std::endl;
      vector<PeptideHit>& hits = pep_it->getHits();
      //std::cout << hits[0].getSequence().getMonoWeight(Residue::Full, hits[0].getCharge()) / hits[0].getCharge() << std::endl;
      double curr_dev = hits[0].getMetaValue("ppmDeviation");
      int isotope_error = 0;
      if (!combine_isotope_errors && hits[0].metaValueExists("IsotopeError"))
      {
        isotope_error = String(hits[0].getMetaValue("IsotopeError")).toInt();
      }
      tmp_outfile << std::setprecision(16) << curr_dev << "\t" << hits[0].getScore() << "\t" << obs_mz << "\t" << hits[0].getSequence().getMonoWeight(Residue::Full, hits[0].getCharge()) / (double) hits[0].getCharge() << "\t" << hits[0].getMetaValue("target_decoy") << "\t" << hits[0].getCharge() << "\t" << pep_it->getRT() << "\t" << hits[0].getMetaValue("calcMZ") << std::endl;
      //std::cout << boost::math::pdf(error_dist, curr_dev) << std::endl;
      hits[0].setMetaValue("mz_dev_prob_correct",
                           boost::math::pdf(error_dists[isotope_error],curr_dev)/(boost::math::pdf(error_dists[isotope_error],curr_dev) + 1/(2*tolerance)));
      if (all_psms)
      {
        for (Size i = 1; i < hits.size(); ++i)
        {
          curr_dev = hits[i].getMetaValue("ppmDeviation");
          tmp_outfile << std::setprecision(16) << curr_dev << "\t" << hits[i].getScore() << "\t" << obs_mz << "\t" << hits[i].getSequence().getMonoWeight(Residue::Full, hits[i].getCharge()) / (double) hits[i].getCharge() << "\t" << hits[i].getMetaValue("target_decoy") << "\t" << hits[i].getCharge() << "\t" << pep_it->getRT() << "\t" << hits[i].getMetaValue("calcMZ") << std::endl;
          isotope_error = 0;
          if (!combine_isotope_errors && hits[i].metaValueExists("IsotopeError"))
          {
            isotope_error = String(hits[i].getMetaValue("IsotopeError")).toInt();
          }
          hits[i].setMetaValue("mz_dev_prob_correct",
                               boost::math::pdf(error_dists[isotope_error],curr_dev)/(boost::math::pdf(error_dists[isotope_error],curr_dev) + 1/(2*tolerance)));
        }
      }
    }
    tmp_outfile.close();

    IdXMLFile().store(outputfile_name, proteins, peptides);

    return EXECUTION_OK;
  }

private:

  // TODO maybe put to Statistics file and templatize
  struct RunningStatistic{
    // Running sum of totals t_i to the i.
    int t_0;
    double t_1;
    double t_2;

    void addValue(double val)
    {
      ++t_0;
      t_1 += val;
      t_2 += val * val;
    }
    double getMean() { return t_1/t_0; }
    double getSD() { return (1.0/t_0) * std::sqrt((t_0 * t_2) - (t_1 * t_1)); }
  };


  void inline processHit_(PeptideHit& hit, vector<RunningStatistic>& stats, double obs_mz, bool combine_isotope_errors)
  {
    int isotope_error = 0;

    //TODO if fitting, we have to store every value. Totals are not enough
    if (!combine_isotope_errors && hit.metaValueExists("IsotopeError"))
      isotope_error = String(hit.getMetaValue("IsotopeError")).toInt();

    double dev = calcPPMDevForHit_(hit, obs_mz, isotope_error);
    hit.setMetaValue("ppmDeviation", dev);
    stats[isotope_error].addValue(dev);
  }

  double inline calcPPMDevForHit_(PeptideHit const & hit, double obs_mz, int isotope_error = 0)
  {
    // Usually we should not use the calculated MZ from e.g. MSGFPlus (metavalue "calcMZ")
    // since it is only single precision (float cast to double)
    if (hit.metaValueExists("ppmDeviation"))
    {
      return hit.getMetaValue("ppmDeviation");
    }
    else
    {
      // Start with an offset of error*neutron_mass
      double calc_mz = isotope_error * Constants::NEUTRON_MASS_U;
      // Add uncharged mass + charge*proton_mass (yes, getMonoWeight works like this).
      calc_mz += hit.getSequence().getMonoWeight(Residue::Full, hit.getCharge());
      // Divide by charge
      calc_mz /= static_cast<double>(hit.getCharge());
      // Return deviation [ppm]
      return ((obs_mz * 1e6) - (calc_mz * 1e6))/obs_mz;
    }

  }

};


int main(int argc, const char** argv)
{
  TOPPMassDeviationScorer tool;

  return tool.main(argc, argv);
}

/// @endcond
