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

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
 @page TOPP_PSMPrescorer PSMPrescorer

 @brief Combines different kinds of scores to calculate a probability for a PSM being correct.
<CENTER>
 <table>
  <tr>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tool </td>
   <td VALIGN="middle" ROWSPAN=5> \f$ \longrightarrow \f$ PSMPrescorer \f$ \longrightarrow \f$</td>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MassDeviationScorer</td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FidoAdapter</td>
  </tr>
 </table>
</CENTER>

 Use MassDeviationScorer and bootstrapped RTPredict to aggregate more types of scores first.

 @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_PSMPrescorer.cli
 <B>INI file documentation of this tool:</B>
 @htmlinclude TOPP_PSMPrescorer.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPPSMPrescorer :
    public TOPPBase
{
public:
  TOPPPSMPrescorer() :
      TOPPBase("PSMPrescorer", "Combines different kinds of scores to calculate a probability for a PSM being correct.")
  {

  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("idXML"));
    registerDoubleOption_("pe1_given_RT1","",0.7,"Peptide present given its RT deviation comes from the positive distribution.",
                          false);
    registerDoubleOption_("pe1_given_RT0","",0.0,"Peptide present given its RT deviation comes from the negative distribution.",
                          false);
    registerDoubleOption_("pe1_given_mz1","",0.9,"Peptide present given its m/z deviation comes from the positive distribution.",
                          false);
    registerDoubleOption_("pe1_given_mz0","",0.0,"Peptide present given its m/z deviation comes from the negative distribution.",
                          false);
    registerDoubleOption_("pe1_given_score1","",0.9,"Peptide present given its m/z deviation comes from the positive distribution.",
                          false);
    registerDoubleOption_("pe1_given_score0","",0.0,"Peptide present given its m/z deviation comes from the negative distribution.",
                          false);
  }

  ExitCodes main_(int, const char**)
  {
    // Reading
    String inputfile_name = getStringOption_("in");
    String outputfile_name = getStringOption_("out");

    vector<ProteinIdentification> proteins;
    vector<PeptideIdentification> peptides;
    IdXMLFile().load(inputfile_name, proteins, peptides);

    pe1RT1 = getDoubleOption_("pe1_given_RT1");
    pe1RT0 = getDoubleOption_("pe1_given_RT0");
    pe1mz1 = getDoubleOption_("pe1_given_mz1");
    pe1mz0 = getDoubleOption_("pe1_given_mz0");
    pe1score1 = getDoubleOption_("pe1_given_score1");
    pe1score0 = getDoubleOption_("pe1_given_score0");


    for (vector<OpenMS::PeptideIdentification>::iterator pep_it = peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      vector<PeptideHit>& hits = pep_it->getHits();

      if (pep_it->getScoreType() == "Posterior Probability" /*&& pep_it->isHigherScoreBetter()*/)
      {
        for (vector<OpenMS::PeptideHit>::iterator hit_it = hits.begin(); hit_it != hits.end(); ++hit_it)
        {
          processHit_(*hit_it);
        }
      }
      else
      {
        // TODO Exception? Wrong score type in at least one PepID. Right predecessor tool?
      }
    }

    IdXMLFile().store(outputfile_name, proteins, peptides);

    return EXECUTION_OK;
  }

private:

  double pe1RT1, pe1RT0, pe1mz1, pe1mz0, pe1score1, pe1score0;

  /// Rescores a PeptideHit if the two metavalues
  /// mz_dev_prob_correct and
  /// rt_dev_prob_correct exist.
  /// Implicitly assumes that the score is a probability coming from the distribution of correct scores.
  /// E.g. from IDPep with prob_correct option
  void inline processHit_(PeptideHit& hit)
  {
    // Assumes score of hit is a probability of being from the correct distribution
    if (hit.metaValueExists("mz_dev_prob_correct") &&
        hit.metaValueExists("predicted_RT_prob_correct_alt"))
    {
      double new_score_pos = 1.0;
      double new_score_neg = 1.0;
      double evidence = 0.0;
      // multiply in marginal probabilities for the peptide with different types of evidence given
      evidence = hit.getMetaValue("mz_dev_prob_correct");
      new_score_pos *= multiplyWithCPTRow_(pe1mz1, evidence);
      new_score_neg *= multiplyWithCPTRow_(pe1mz0, evidence);
      evidence = hit.getMetaValue("predicted_RT_prob_correct_alt");
      new_score_pos *= multiplyWithCPTRow_(pe1RT1, evidence);
      new_score_neg *= multiplyWithCPTRow_(pe1RT0, evidence);
      evidence = hit.getScore();
      new_score_pos *= multiplyWithCPTRow_(pe1score1, evidence);
      new_score_neg *= multiplyWithCPTRow_(pe1score0, evidence);

      //Normalize
      new_score_pos = new_score_pos / (new_score_pos + new_score_neg);
      //new_score_neg = 1 - new_score_pos;

      //Save old score as metavalue
      hit.setMetaValue("score_prob_correct", evidence);
      hit.setScore(new_score_pos);
      //score type is already assumed to be probability with higher=better
    }
    else
    {
      std::cout << "Warning: PeptideHit without the mz_dev/predicted_RT_prob_correct metavalues. Skipping."<< std::endl;
    }
  }

  /// Simulates marginalizing out evidence (defined by param p) for Y in a
  /// conditional probability table P(X|Y) partly specified by the parameter a
  /// @param  a = P(X=1|Y=y)
  /// @param  p = P(X=1)
  /// @returns unnormalized P(Y=y)
  double inline multiplyWithCPTRow_(double a, double p)
  {
    return a * p + (1 - a) * (1 - p);
  }

};


int main(int argc, const char** argv)
{
  TOPPPSMPrescorer tool;

  return tool.main(argc, argv);
}

/// @endcond
