// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#ifndef OPENMS_ANALYSIS_ID_BAYESIANPROTEININFERENCE_H
#define OPENMS_ANALYSIS_ID_BAYESIANPROTEININFERENCE_H

//#define INFERENCE_BENCH

#include <OpenMS/ANALYSIS/ID/MessagePasserFactory.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/MATH/MISC/GridSearch.h>
#include <vector>
#include <functional>

namespace OpenMS
{
  class ConsensusMap;
  class IDBoostGraph;
  class PeptideIdentification;
  class ProteinIdentification;

  class OPENMS_DLLAPI BayesianProteinInferenceAlgorithm :
      public DefaultParamHandler,
      public ProgressLogger
  {
  public:
    /// Constructor
    explicit BayesianProteinInferenceAlgorithm(unsigned int debug_lvl = 0);

    /// Destructor
    ~BayesianProteinInferenceAlgorithm() override = default;

    /// A function object to pass into the IDBoostGraph class to perform algorithms on
    /// connected components
    class GraphInferenceFunctor;

    /// A function object to pass into the IDBoostGraph class to perform algorithms on
    /// connected components
    class ExtendedGraphInferenceFunctor;

    /// A function object to pass into the IDBoostGraph class to annotate and add
    /// indistinguishable groups to the underlying ID objects based on the graph.
    class AnnotateIndistGroupsFunctor;

    /// A function object to pass into the GridSearch;
    struct GridSearchEvaluator;

    /// Perform inference. Writes its results into protein and (optionally) peptide hits (overwrites score).
    /// Optionally adds indistinguishable protein groups with seperate scores, too.
    /// Currently only takes first proteinID run.
    /// TODO loop over all runs
    void inferPosteriorProbabilities(std::vector<ProteinIdentification>& proteinIDs, std::vector<PeptideIdentification>& peptideIDs);

    void inferPosteriorProbabilities(ConsensusMap& cmap);

  private:

    void inferPosteriorProbabilities_(IDBoostGraph& ibg);

    GridSearch<double,double,double> initGridSearchFromParams_(
        std::vector<double>& alpha_search,
        std::vector<double>& beta_search,
        std::vector<double>& gamma_search
        );

    std::function<void(PeptideIdentification&)> checkConvertAndFilterPepHits;

    /// The grid search object initialized with a default grid
    GridSearch<double,double,double> grid{{0.008,0.032,0.128},{0.001},{0.5}};

    unsigned int debug_lvl_;

    #ifdef INFERENCE_BENCH
    std::vector<std::pair<double,Size>> debug_times_;
    #endif

  };
}
#endif // OPENMS_ANALYSIS_ID_BAYESIANPROTEININFERENCE_H
