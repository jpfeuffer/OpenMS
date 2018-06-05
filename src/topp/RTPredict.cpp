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
// $Maintainer: Mathias Walzer $
// $Authors: Nico Pfeifer, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/SVM/SVMWrapper.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/LibSVMEncoder.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

#include <map>
#include <cmath>
#include <boost/math/distributions/normal.hpp>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_RTPredict RTPredict

    @brief This application is used to predict retention times
                 for peptides or peptide separation.

    This methods and applications of this model are described
    in several publications:

    Nico Pfeifer, Andreas Leinenbach, Christian G. Huber and Oliver Kohlbacher
    Statistical learning of peptide retention behavior in chromatographic separations: A new kernel-based approach for computational proteomics.
    BMC Bioinformatics 2007, 8:468

    Nico Pfeifer, Andreas Leinenbach, Christian G. Huber and Oliver Kohlbacher
    Improving Peptide Identification in Proteome Analysis by a Two-Dimensional Retention Time Filtering Approach
    J. Proteome Res. 2009, 8(8):4109-15


    The input of this application
    is an svm model and a file with peptide identifications (idXML or text).
    The svm model file is specified
    by the <b>svm_model</b> parameter in the command line or the INI file.
    This file should have been produced by the @ref TOPP_RTModel application.
    <br>
    For retention time prediction the peptide sequences are extracted
    from the idXML/text inputfile
    and passed to the svm. The svm then predicts retention times
    according to the trained model. The predicted retention times
    are stored as @code <userParam name="predicted_retention_time" value="<predicted retention time>" />
    @endcode inside the peptide entities in the idXML output file.

    For separation prediction you have to specify two output file names.
    'out_id:positive' is the filename of the peptides which are predicted
    to be collected by the column and 'out_id:negative' is the file
    of the predicted flowthrough peptides.

    Retention time prediction and separation prediction cannot be combined!

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_RTPredict.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_RTPredict.html

    @todo This needs serious clean up! Combining certain input and output options will
          result in strange behaviour, especially when using text output/input.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPRTPredict :
  public TOPPBase
{
public:
  TOPPRTPredict() :
    TOPPBase("RTPredict", "Predicts retention times for peptides using a model trained by RTModel.")
  {
  }

protected:
  struct OligoModel
  {
    bool first_dim_rt;
    UInt border_length;
    UInt k_mer_length;
    double sigma;
    double sigma_0;
    double sigma_max;
  };

  void registerOptionsAndFlags_()
  {
    // input
    registerInputFile_("in_id", "<file>", "", "Peptides with precursor information", false);
    setValidFormats_("in_id", ListUtils::create<String>("idXML"));
    registerInputFile_("in_text", "<file>", "", "Peptides as text-based file", false);
    setValidFormats_("in_text", ListUtils::create<String>("txt"));


    registerInputFileList_("svm_model", "<files>", StringList(), "SVM model(s) in libsvm format (can be produced by RTModel). Multiple models can be used to infer bootstrap statistics.");
    setValidFormats_("svm_model", ListUtils::create<String>("txt"));
    registerInputFileList_("in_oligo_params", "<files>", StringList(), "Input file(s) with additional model parameters when using the OLIGO kernel. Same number and order as in svm_model expected.", false);
    setValidFormats_("in_oligo_params", ListUtils::create<String>("paramXML"));
    registerInputFileList_("in_oligo_trainset", "<files>", StringList(), "input file(s) with the used training dataset when using the OLIGO kernel. Same number and order as in svm_model expected.", false);
    setValidFormats_("in_oligo_trainset", ListUtils::create<String>("txt"));

    registerDoubleOption_("total_gradient_time", "<time>", 1.0, "The time (in seconds) of the gradient (peptide RT prediction)", false);
    setMinFloat_("total_gradient_time", 0.00001);
    registerIntOption_("max_number_of_peptides", "<int>", 100000, "The maximum number of peptides considered at once (bigger number will lead to faster results but needs more memory).", false, true);

    // output
    registerTOPPSubsection_("out_id", "Output files in idXML format");
    registerOutputFile_("out_id:file", "<file>", "", "Output file with peptide RT prediction", false);
    setValidFormats_("out_id:file", ListUtils::create<String>("idXML"));
    registerOutputFile_("out_id:positive", "<file>", "", "Output file in idXML format containing positive predictions for peptide separation prediction - requires 'out_id:negative' to be present as well.", false);
    setValidFormats_("out_id:positive", ListUtils::create<String>("idXML"));
    registerOutputFile_("out_id:negative", "<file>", "", "Output file in idXML format containing negative predictions for peptide separation prediction - requires 'out_id:positive' to be present as well.", false);
    setValidFormats_("out_id:negative", ListUtils::create<String>("idXML"));
    registerFlag_("out_id:rewrite_peptideidentification_rtmz", "Rewrites each peptideidentification's rt and mz from prediction and calculation (according to the best hit)", true);

    registerTOPPSubsection_("out_text", "Output files in text format");
    registerOutputFile_("out_text:file", "<file>", "", "Output file with predicted RT values", false);
    setValidFormats_("out_text:file", ListUtils::create<String>("csv"));

    registerTOPPSubsection_("bootstrap", "Options for bootstrapping when given multiple files");
    registerFlag_("bootstrap:use_mean", "Use the mean predicted RT instead of the initial prediction", true);

  }

  static inline double computeSquare (double x) { return x*x; }

  void loadStrings_(String filename, std::vector<String>& sequences)
  {
    TextFile text_file(filename.c_str(), true);
    TextFile::ConstIterator it;

    sequences.clear();

    it = text_file.begin();
    while (it != text_file.end())
    {
      sequences.push_back(*it);
      sequences.back().trim();
      ++it;
    }
  }

  void writeStringLabelLines_(String filename, map<String, double> predicted_data)
  {
    ofstream os;
    map<String, double>::iterator it;

    os.open(filename.c_str(), ofstream::out);

    for (it = predicted_data.begin(); it != predicted_data.end(); ++it)
    {
      os << it->first << " " << it->second << "\n";
    }
    os.flush();
    os.close();
  }


  // TODO maybe put whole code for Oligo kernel in SVM Wrapper or a subclass of it)
  ExitCodes readOligoModel_(String const & params_filename, String const & trainset_filename, bool separation_prediction, SVMWrapper & svm, OligoModel & oligomodel)
  {
    Param additional_parameters;
    ParamXMLFile paramFile;
    paramFile.load(params_filename, additional_parameters);
    if (additional_parameters.exists("first_dim_rt")
        && additional_parameters.getValue("first_dim_rt") != DataValue::EMPTY) {
      oligomodel.first_dim_rt = additional_parameters.getValue("first_dim_rt").toBool();
    }
    if (additional_parameters.getValue("kernel_type") != DataValue::EMPTY) {
      svm.setParameter(SVMWrapper::KERNEL_TYPE, ((String) additional_parameters.getValue("kernel_type")).toInt());
    }

    if (additional_parameters.getValue("border_length") == DataValue::EMPTY
        && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
      writeLog_("No border length saved in additional parameters file. Aborting!");
      cout << "No border length saved in additional parameters file. Aborting!" << endl;
      return ILLEGAL_PARAMETERS;
    }
    oligomodel.border_length = ((String) additional_parameters.getValue("border_length")).toInt();
    if (additional_parameters.getValue("k_mer_length") == DataValue::EMPTY
        && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
      writeLog_("No k-mer length saved in additional parameters file. Aborting!");
      cout << "No k-mer length saved in additional parameters file. Aborting!" << endl;
      return ILLEGAL_PARAMETERS;
    }
    oligomodel.k_mer_length = ((String) additional_parameters.getValue("k_mer_length")).toInt();
    if (additional_parameters.getValue("sigma") == DataValue::EMPTY
        && svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
      writeLog_("No sigma saved in additional parameters file. Aborting!");
      cout << "No sigma saved in additional parameters file. Aborting!" << endl;
      return ILLEGAL_PARAMETERS;
    }
    oligomodel.sigma = ((String) additional_parameters.getValue("sigma")).toDouble();
    if (!separation_prediction && additional_parameters.getValue("sigma_0") == DataValue::EMPTY) {
      writeLog_("No sigma_0 saved in additional parameters file. Aborting!");
      cout << "No sigma_0 length saved in additional parameters file. Aborting!" << endl;
      return ILLEGAL_PARAMETERS;
    }
    if (!separation_prediction && additional_parameters.getValue("sigma_0") != DataValue::EMPTY) {
      oligomodel.sigma_0 = additional_parameters.getValue("sigma_0");
    }
    if (!separation_prediction && additional_parameters.getValue("sigma_max") == DataValue::EMPTY) {
      writeLog_("No sigma_max saved in additional parameters file. Aborting!");
      cout << "No sigma_max length saved in additional parameters file. Aborting!" << endl;
      return ILLEGAL_PARAMETERS;
    }
    if (!separation_prediction && additional_parameters.getValue("sigma_max") != DataValue::EMPTY) {
      oligomodel.sigma_max = additional_parameters.getValue("sigma_max");
    }
    SVMData training_samples;
    training_samples.load(trainset_filename);
    svm.setTrainingSample(training_samples);

    svm.setParameter(SVMWrapper::BORDER_LENGTH, (Int) oligomodel.border_length);
    svm.setParameter(SVMWrapper::SIGMA, oligomodel.sigma);
    return EXECUTION_OK;
  }

  ExitCodes main_(int, const char**) override
  {
    //CONSTANTS
    static const UInt maximum_length = 50;

    // For IdXML in and output
    IdXMLFile idXML_file;
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    vector<PeptideHit> temp_peptide_hits;
    vector<PeptideIdentification> identifications_positive;
    vector<PeptideIdentification> identifications_negative;
    // Input files
    StringList svmfile_names = StringList();
    StringList oligo_model_file_names = StringList();
    StringList oligo_trainset_file_names = StringList();
    // Either is used depending on kernel and input file type
    vector<String> peptides;
    vector<AASequence> modified_peptides;
    Size number_of_peptides = 0;
    // The models
    SVMWrapper svm;
    OligoModel oligomodel = OligoModel();
    LibSVMEncoder encoder;
    // Usage depending on kernel
    svm_problem* prediction_data = nullptr;
    SVMData prediction_samples;
    // Restrictions
    String allowed_amino_acid_characters = "ACDEFGHIKLMNPQRSTVWY";
    Size max_number_of_peptides = getIntOption_("max_number_of_peptides");
    // Other options
    double total_gradient_time = 1.;
    bool separation_prediction = false;
    // Temporary results
    // TODO is transforming String to AASeq or the other way and doing the lookups
    // in only one map too inefficient? Would save a lot of if-statements.
    map<String, double> predicted_data;
    map<AASequence, double> predicted_modified_data;
    pair<double, double> temp_point;
    vector<double> predicted_retention_times;
    // Final results (Results for each PeptideHit from each PeptideID concatenated)
    vector<double> final_predicted_retention_times;
    vector<double> bootstrap_sd;
    vector<float> performance_retention_times;

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String outputfile_name_positive = getStringOption_("out_id:positive");
    String outputfile_name_negative = getStringOption_("out_id:negative");
    // for separation prediction, we require both files to be present!
    if (outputfile_name_positive != "" || outputfile_name_negative != "")
    {
      if (outputfile_name_positive != "" && outputfile_name_negative != "")
      {
        separation_prediction = true;
      }
      else
      {
        writeLog_("Both files for separation prediction required. Please specify the other one as well. Aborting!");
        return ILLEGAL_PARAMETERS;
      }
    }

    // either or
    String input_id = getStringOption_("in_id");
    String input_text = getStringOption_("in_text");
    if (input_text != "" && input_id != "")
    {
      writeLog_("Two input parameter files given, only one allowed! Use either -in_id:file or -in_text:file!");
      return ILLEGAL_PARAMETERS;
    }
    else if (input_text == "" && input_id == "")
    {
      writeLog_("No input file given. Aborting...");
      return ILLEGAL_PARAMETERS;
    }

    // OUTPUT
    // (can use both)
    String output_id = getStringOption_("out_id:file");
    String output_text = getStringOption_("out_text:file");
    if (output_text == "" && output_id == "" && !separation_prediction)
    {
      writeLog_("No output files given. Aborting...");
      return ILLEGAL_PARAMETERS;
    }

    svmfile_names = getStringList_("svm_model");
    oligo_model_file_names = getStringList_("in_oligo_params");
    oligo_trainset_file_names = getStringList_("in_oligo_trainset");
    total_gradient_time = getDoubleOption_("total_gradient_time");

    //-------------------------------------------------------------
    // Prepare reading input
    //-------------------------------------------------------------

    if (svmfile_names.size() <= 0)
    {
      writeLog_("No SVM model specified. Please train one using RTModel first."
                + String("\nAborting!"));
      return ILLEGAL_PARAMETERS;
    }
    if ((svmfile_names.size() > 1) && separation_prediction)
    {
        writeLog_("Calculating bootstrap statistics for peptide separation prediction is not supported."
                  + String("\nPlease use a single input model. Aborting!"));
        return ILLEGAL_PARAMETERS;
    }

    //Load first model file for actual predictions here to peak into the used SVM_TYPE
    String svmfile_name = svmfile_names[0];

    svm.loadModel(svmfile_name);
    int svm_type = svm.getIntParameter(SVMWrapper::SVM_TYPE);

    if ((svm_type == C_SVC || svm_type == NU_SVC) && !separation_prediction) {
      writeLog_("You cannot perform peptide separation prediction with a model trained for"
                + String("\npeptide retention time prediction. Aborting!"));
      return ILLEGAL_PARAMETERS;
    }
    if ((svm_type != C_SVC && svm_type != NU_SVC) && separation_prediction) {
      writeLog_("You cannot perform peptide retention time prediction with a model trained for\n"
                + String("peptide separation prediction. Aborting!"));
      return ILLEGAL_PARAMETERS;
    }

    //-------------------------------------------------------------
    // Reading peptide input
    //-------------------------------------------------------------
    if (input_text != "") {
      loadStrings_(input_text, peptides);
      if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
        for (Size i = 0; i < peptides.size(); ++i) {
          modified_peptides.push_back(AASequence::fromString(peptides[i]));
        }
        peptides.clear();
      }
    } else {
      String document_id;
      idXML_file.load(input_id, protein_identifications, identifications, document_id);
    }

    if (input_id != "") {
      for (Size i = 0; i < identifications.size(); i++) {
        temp_peptide_hits = identifications[i].getHits();
        for (Size j = 0; j < temp_peptide_hits.size(); j++) {
          if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
            modified_peptides.push_back(temp_peptide_hits[j].getSequence());
          } else {
            peptides.push_back(temp_peptide_hits[j].getSequence().toUnmodifiedString());
          }
        }
      }
    }

    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
      number_of_peptides = modified_peptides.size();
    } else {
      number_of_peptides = peptides.size();
    }

    //-------------------------------------------------------------
    // Reading additional Oligo model input
    //-------------------------------------------------------------
    // Since the POBK is not included in the libsvm we have to load
    // additional parameters from additional files.

    if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
      String in_params_name = oligo_model_file_names[0];
      if (in_params_name.empty()) {
        in_params_name = svmfile_name + "_additional_parameters";
        writeLog_("Warning: Using OLIGO kernel but in_oligo_params parameter is missing. Trying default filename: " +
                  in_params_name);
      }
      inputFileReadable_(in_params_name, "in_oligo_params");

      String in_trainset_name = oligo_trainset_file_names[0];
      if (in_trainset_name.empty()) {
        in_trainset_name = svmfile_name + "_samples";
        writeLog_(
                "Warning: Using OLIGO kernel but in_oligo_trainset parameter is missing. Trying default filename: " +
                in_trainset_name);
      }
      inputFileReadable_(in_trainset_name.c_str(), "in_oligo_trainset");

      ExitCodes parseOligoExit = readOligoModel_(in_params_name, in_trainset_name, separation_prediction, svm, oligomodel);
      if (parseOligoExit != 0)
      {
        return parseOligoExit;
      }
    }


    //-------------------------------------------------------------
    // Prediction on first (main) model
    //-------------------------------------------------------------

    vector<String>::iterator it_from = peptides.begin();
    vector<String>::iterator it_to = peptides.begin();
    vector<AASequence>::iterator it_from_mod = modified_peptides.begin();
    vector<AASequence>::iterator it_to_mod = modified_peptides.begin();
    Size counter = 0;

    while (counter < number_of_peptides) {
      vector<String> temp_peptides(0);
      vector<AASequence> temp_modified_peptides(0);
      vector<double> temp_rts(0);


      // Read in up to max_number_of_peptides
      Size temp_counter = 0;
      if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)
      {
        while ((temp_counter < max_number_of_peptides) && (it_to != peptides.end()))
        {
          ++it_to;
          ++temp_counter;
        }
        temp_peptides.insert(temp_peptides.end(), it_from, it_to);
        //temp_peptides.insert(temp_peptides.end(), peptides.begin(), peptides.end());
        temp_rts.resize(temp_peptides.size(), 0.0);

        prediction_data =
                encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(temp_peptides,
                                                                           temp_rts,
                                                                           allowed_amino_acid_characters,
                                                                           maximum_length);
        it_from = it_to;
      }
      else //OLIGO
      {
        while ((temp_counter < max_number_of_peptides) && (it_to_mod != modified_peptides.end()))
        {
          ++it_to_mod;
          ++temp_counter;
        }
        temp_modified_peptides.insert(temp_modified_peptides.end(), it_from_mod, it_to_mod);
        // temp_modified_peptides.insert(temp_modified_peptides.end(), modified_peptides.begin(), modified_peptides.end());
        temp_rts.resize(temp_modified_peptides.size(), 0.0);

        encoder.encodeProblemWithOligoBorderVectors(temp_modified_peptides,
                                                    oligomodel.k_mer_length,
                                                    allowed_amino_acid_characters,
                                                    oligomodel.border_length,
                                                    prediction_samples.sequences);
        prediction_samples.labels = temp_rts;
        it_from_mod = it_to_mod;
      }
      counter += temp_counter;

      if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
      {
        svm.predict(prediction_samples, predicted_retention_times);
        prediction_samples.labels.clear();
        prediction_samples.sequences.clear();
      }
      else
      {
        svm.predict(prediction_data, predicted_retention_times);
      }
      LibSVMEncoder::destroyProblem(prediction_data);


      for (Size i = 0; i < temp_counter; ++i)
      {
        double rt_store = predicted_retention_times[i] * total_gradient_time;

        // Inserts (sequence -> predicted_RT) into map
        // Duplicate entries will be predicted and inserted multiple times.
        // Last prediction is used but should be the same
        if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
        {
          if (output_id != "")
          {
            predicted_modified_data.insert(make_pair(temp_modified_peptides[i], rt_store));
          }
          else
          {
            predicted_data.insert(make_pair(temp_modified_peptides[i].toString(), rt_store));
          }
        }
        else /*(svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)*/
        {
          predicted_data.insert(make_pair(temp_peptides[i], rt_store));
        }
        final_predicted_retention_times.push_back(rt_store);
      }
      predicted_retention_times.clear();
    } // Main prediction finished


    //-------------------------------------------------------------
    // Reading and applying additional bootstrap models
    //-------------------------------------------------------------
    long t0 = 0;
    std::vector<double> t1(0);
    std::vector<double> t2(0);

    if (svmfile_names.size() > 2) {  // With Bootstrap statistics (at least 3 predictions needed for standard deviation)

      if (separation_prediction) {
        writeLog_("Calculating bootstrap statistics for peptide separation prediction is not supported."
                  + String("\nPlease use a single input model. Aborting!"));
        return ILLEGAL_PARAMETERS;
      }

      // Initialization of online data structures

      // t_k is the total of the values x_i to the k. t_x = sum_{i=1}^{n} x_i^n
      // for online updating of standard deviation.

      // basically the number of datapoints. we know that beforehand from the number of bootstraps.
      //TODO use RunningStatistics object (i.e. factor out from MassDevScorer).
      t0 = svmfile_names.size();
      t1 = final_predicted_retention_times;
      t2.resize(number_of_peptides);
      std::transform(t1.begin(), t1.end(), t2.begin(), computeSquare);

      // TODO factor out prediction loop body?
      //-------------------------------------------------------------
      // Prediction for every bootstrap model
      //-------------------------------------------------------------

      for (int k = 1; k < svmfile_names.size(); ++k)
      {

        vector<double> current_predicted_retention_times (0);
        SVMWrapper curr_svm;
        OligoModel curr_oligomodel;

        svmfile_name = svmfile_names[k];
        curr_svm.loadModel(svmfile_name);

        if (curr_svm.getIntParameter(SVMWrapper::SVM_TYPE) != svm_type) {
          writeLog_("Model from file " + svmfile_name + " was not trained with"
                    + String("\nthe same type of SVM. Aborting!"));
          return ILLEGAL_PARAMETERS;
        }

        // Since the POBK is not included in the libsvm we have to load
        // additional parameters from additional files.
        if (curr_svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
        {
          String in_params_name = "";
          if (k < oligo_model_file_names.size())
          {
            in_params_name = oligo_model_file_names[k];
          }
          else
          {
            writeLog_("Bootstrap models trained with Oligo kernel but not the same amount of oligo model"
                      + String("\nfiles as svm models present. Aborting!"));
            return ILLEGAL_PARAMETERS;
          }

          if (in_params_name.empty())
          {
            in_params_name = svmfile_name + "_additional_parameters";
            writeLog_(
                "Warning: Using OLIGO kernel but in_oligo_params parameter is missing. Trying default filename: " +
                in_params_name);
          }
          inputFileReadable_(in_params_name, "in_oligo_params");


          String in_trainset_name = "";
          if (k < oligo_trainset_file_names.size())
          {
            in_trainset_name = oligo_trainset_file_names[k];
          }
          else
          {
            writeLog_("Bootstrap models trained with Oligo kernel but not the same amount of oligo trainset"
                      + String("\nfiles as svm models present. Aborting!"));
            return ILLEGAL_PARAMETERS;
          }
          if (in_trainset_name.empty())
          {
            in_trainset_name = svmfile_name + "_samples";
            writeLog_(
                "Warning: Using OLIGO kernel but in_oligo_trainset parameter is missing. Trying default filename: " +
                in_trainset_name);
          }
          inputFileReadable_(in_trainset_name.c_str(), "in_oligo_trainset");

          ExitCodes parseOligoExit = readOligoModel_(in_params_name, in_trainset_name, separation_prediction, curr_svm,
                                                     curr_oligomodel);
          if (parseOligoExit != 0)
          {
            return parseOligoExit;
          }
        }

          vector<String>::iterator curr_it_from = peptides.begin();
          vector<String>::iterator curr_it_to = peptides.begin();
          vector<AASequence>::iterator curr_it_from_mod = modified_peptides.begin();
          vector<AASequence>::iterator curr_it_to_mod = modified_peptides.begin();
          Size curr_counter = 0;

          while (curr_counter < number_of_peptides) {
            vector<String> curr_temp_peptides(0);
            vector<AASequence> curr_temp_modified_peptides(0);
            vector<double> curr_temp_rts(0);

            // Read in up to max_number_of_peptides
            Size curr_temp_counter = 0;
            if (curr_svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)
            {
              while (curr_temp_counter < max_number_of_peptides && curr_it_to != peptides.end())
              {
                ++curr_it_to;
                ++curr_temp_counter;
              }
              curr_temp_peptides.insert(curr_temp_peptides.end(), curr_it_from, curr_it_to);
              curr_temp_rts.resize(curr_temp_peptides.size(), 0);

              prediction_data =
                      encoder.encodeLibSVMProblemWithCompositionAndLengthVectors(curr_temp_peptides,
                                                                                 curr_temp_rts,
                                                                                 allowed_amino_acid_characters,
                                                                                 maximum_length);
              curr_it_from = curr_it_to;
            }
            else //OLIGO
            {
              while ((curr_temp_counter < max_number_of_peptides) && (curr_it_to_mod != modified_peptides.end()))
              {
                ++curr_it_to_mod;
                ++curr_temp_counter;
              }
              curr_temp_modified_peptides.insert(curr_temp_modified_peptides.end(), curr_it_from_mod, curr_it_to_mod);
              curr_temp_rts.resize(curr_temp_modified_peptides.size(), 0);

              encoder.encodeProblemWithOligoBorderVectors(curr_temp_modified_peptides,
                                                          curr_oligomodel.k_mer_length,
                                                          allowed_amino_acid_characters,
                                                          curr_oligomodel.border_length,
                                                          prediction_samples.sequences);
              prediction_samples.labels = curr_temp_rts;
              curr_it_from_mod = curr_it_to_mod;
            }
            curr_counter += curr_temp_counter;

            if (curr_svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
            {
              curr_svm.predict(prediction_samples, predicted_retention_times);
              prediction_samples.labels.clear();
              prediction_samples.sequences.clear();
            }
            else
            {
              curr_svm.predict(prediction_data, predicted_retention_times);
            }
            LibSVMEncoder::destroyProblem(prediction_data);

            for (Size i = 0; i < curr_temp_counter; ++i) {
              if (curr_svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO)
              {
                if (output_id != "")
                {
                  predicted_modified_data.insert(make_pair(curr_temp_modified_peptides[i],
                                                           (predicted_retention_times[i] * total_gradient_time)));
                }
                else
                {
                  predicted_data.insert(make_pair(curr_temp_modified_peptides[i].toString(),
                                                  (predicted_retention_times[i] * total_gradient_time)));
                }
              }
              else /*(svm.getIntParameter(SVMWrapper::KERNEL_TYPE) != SVMWrapper::OLIGO)*/
              {
                predicted_data.insert(make_pair(curr_temp_peptides[i],
                                                (predicted_retention_times[i] * total_gradient_time)));
              }
              current_predicted_retention_times.push_back(predicted_retention_times[i] * total_gradient_time);
            }
            predicted_retention_times.clear();
          }

          //TODO maybe it is safer to have maps with the sequence instead of relying on the order to be the same.
          //Also would reduce space because of removing duplicates.

          // Update totals
          for (vector<double>::iterator
                   tot_it = t1.begin(),
                   sqtot_it = t2.begin(),
                   rt_it = current_predicted_retention_times.begin();
               tot_it != t1.end();
               ++tot_it, ++sqtot_it, ++rt_it)
          {
            *tot_it += *rt_it;
            *sqtot_it += *rt_it * *rt_it;
          }

          current_predicted_retention_times.clear();
        }
      } // End iterate through bootstrap models


    //-------------------------------------------------------------
    // Filling data structures
    //-------------------------------------------------------------

    //TODO debug output
    //std::ofstream tmp_outfile ("rt_devs.txt", std::ofstream::out | std::ofstream::trunc);
    // IdXML output
    if (input_id != "") {
      if (!separation_prediction) {
        for (Size i = 0; i < identifications.size(); i++) {
          temp_peptide_hits = identifications[i].getHits();

          for (Size j = 0; j < temp_peptide_hits.size(); j++) {
            double temp_rt = 0.;
            double temp_p_value = 0.;
            double temp_training_stddev = 0.;
            double temp_training_mean = 0.;
            double temp_bootstrap_mean = 0.;
            double temp_bootstrap_stddev = 0.;
            double temp_bootstrap_prob_correct = 0.;
            double temp_likelihood_correct = 0.;
            double temp_likelihood_false = 0.;
            double temp_bootstrap_prob_correct_alt = 0.;

            if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
              temp_rt = predicted_modified_data[temp_peptide_hits[j].getSequence()];
            } else {
              temp_rt = predicted_data[temp_peptide_hits[j].getSequence().toUnmodifiedString()];
            }

            temp_point.first = 0.;
            if (oligomodel.first_dim_rt) {
              temp_point.first = identifications[i].getMetaValue("first_dim_rt");
            } else if (identifications[i].hasRT()) {
              temp_point.first = identifications[i].getRT();
            }

            if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
              temp_point.second = temp_rt;
              temp_p_value = svm.getPValue(oligomodel.sigma_0, oligomodel.sigma_max, temp_point);
              if (!t1.empty()) //bootstrap data available
              {
                // experimental RT
                temp_training_mean = temp_point.first;
                // standard deviation based on RT-dependent fit of a 95% confidence interval of the trained model
                // loaded with SVMWrapper.
                temp_training_stddev = svm.getStdDevAtPoint(oligomodel.sigma_0, oligomodel.sigma_max, temp_point);
                // mean and SD of the bootstrap predictions
                // add indices to account for concatenation of all PeptideHits in all PeptideIDs
                temp_bootstrap_mean = t1[i+j] / t0;
                temp_bootstrap_stddev = (1.0/t0) * std::sqrt((t0 * t2[i+j]) - (t1[i+j] * t1[i+j]));
                // TODO remove?
                boost::math::normal_distribution<double> correct_dist = boost::math::normal_distribution<double>(
                  temp_bootstrap_mean,
                  450.0//temp_bootstrap_stddev. Hack: did EM in R
                );
                boost::math::normal_distribution<double> incorrect_dist = boost::math::normal_distribution<double>(
                    temp_bootstrap_mean,
                    2500.0//temp_bootstrap_stddev. Hack: did EM in R
                );
                temp_bootstrap_prob_correct_alt = boost::math::pdf(correct_dist,temp_point.first) /
                    (boost::math::pdf(correct_dist,temp_point.first) + boost::math::pdf(incorrect_dist,temp_point.first));

                // unnormalized expectation of likelihood coming from the correct distribution (from training data)
                temp_likelihood_correct = integrateGaussianPDFProduct_(temp_training_mean,
                                                                       temp_training_stddev,
                                                                       temp_bootstrap_mean,
                                                                       temp_bootstrap_stddev,
                                                                       0.,
                                                                       total_gradient_time);
                // unnormalized expectation of likelihood coming from a uniform distribution along the whole gradient
                temp_likelihood_false = integrateGaussianPDF_(temp_bootstrap_mean,
                                                              temp_bootstrap_stddev,
                                                              0.,
                                                              total_gradient_time) / total_gradient_time;
                // Bayes rule to get the probability of being a correct hit
                temp_bootstrap_prob_correct = temp_likelihood_correct / (temp_likelihood_correct + temp_likelihood_false);
              }
            }
            if (oligomodel.first_dim_rt) {
              if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
                temp_peptide_hits[j].setMetaValue("predicted_RT_p_value_first_dim", temp_p_value);
              }
              if (!t1.empty()) //bootstrap data available
              {
                temp_peptide_hits[j].setMetaValue("predicted_RT_prob_correct_first_dim", temp_bootstrap_prob_correct);
              }
              temp_peptide_hits[j].setMetaValue("predicted_RT_first_dim", temp_rt);
              performance_retention_times.push_back(identifications[i].getMetaValue("first_dim_rt"));
            } else {
              if (identifications[i].hasRT()) {
                performance_retention_times.push_back(identifications[i].getRT());
              } else {
                performance_retention_times.push_back(0.0f);
              }
              if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
                temp_peptide_hits[j].setMetaValue("predicted_RT_p_value", temp_p_value);
              }
              if (!t1.empty()) //bootstrap data available
              {
                temp_peptide_hits[j].setMetaValue("predicted_RT_prob_correct", temp_bootstrap_prob_correct);
                temp_peptide_hits[j].setMetaValue("predicted_RT_bootstrap_mean", temp_bootstrap_mean);
                temp_peptide_hits[j].setMetaValue("predicted_RT_bootstrap_std_dev", temp_bootstrap_stddev);
                temp_peptide_hits[j].setMetaValue("predicted_RT_training_std_dev", temp_training_stddev);
                temp_peptide_hits[j].setMetaValue("predicted_RT_prob_correct_alt", temp_bootstrap_prob_correct_alt);
              }
              temp_peptide_hits[j].setMetaValue("predicted_RT", temp_rt);
            }
          }
          identifications[i].setHits(temp_peptide_hits);
          if (getFlag_("out_id:rewrite_peptideidentification_rtmz")) {
            identifications[i].sort();
            Int charge = identifications[i].getHits().front().getCharge();
            double mz = identifications[i].getHits().front().getSequence().getMonoWeight(Residue::Full, charge) /
                        double(charge);
            double rt = identifications[i].getHits().front().getMetaValue("predicted_RT");

            identifications[i].setRT(rt);
            identifications[i].setMZ(mz);
          }

          identifications[i].setHits(temp_peptide_hits);
        }
      }
      else // separation prediction
      {
        vector<PeptideHit> hits_positive;
        vector<PeptideHit> hits_negative;

        PeptideIdentification temp_identification;

        for (Size i = 0; i < identifications.size(); i++) {
          hits_negative.clear();
          hits_positive.clear();

          temp_peptide_hits = identifications[i].getHits();
          for (vector<PeptideHit>::iterator it = temp_peptide_hits.begin();
               it != temp_peptide_hits.end();
               ++it) {
            if (svm.getIntParameter(SVMWrapper::KERNEL_TYPE) == SVMWrapper::OLIGO) {
              if (predicted_modified_data[it->getSequence()] > 0) {
                hits_positive.push_back(*it);
              } else {
                hits_negative.push_back(*it);
              }
            } else {
              if (predicted_data[it->getSequence().toUnmodifiedString()] > 0) {
                hits_positive.push_back(*it);
              } else {
                hits_negative.push_back(*it);
              }
            }
          }

          temp_identification.setMZ(identifications[i].getMZ());
          temp_identification.setRT(identifications[i].getRT());

          temp_identification = identifications[i];
          temp_identification.setHits(hits_positive);
          identifications_positive.push_back(temp_identification);
          temp_identification.setHits(hits_negative);
          identifications_negative.push_back(temp_identification);
        }
      }
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    if (separation_prediction)
    {
      idXML_file.store(outputfile_name_positive,
                       protein_identifications,
                       identifications_positive);
      idXML_file.store(outputfile_name_negative,
                       protein_identifications,
                       identifications_negative);
    }
    else
    {
      if (output_text != "") // text
      {
        writeStringLabelLines_(output_text, predicted_data);
        //TODO write out bootstrap statistics?
      }
      if (output_id != "") // idXML
      {
        idXML_file.store(output_id,
                         protein_identifications,
                         identifications);
        // TODO check if better solution
        if (!performance_retention_times.empty()) {
          writeDebug_("Linear correlation between predicted and measured rt is: "
                      + String(Math::pearsonCorrelationCoefficient(final_predicted_retention_times.begin(),
                                                                   final_predicted_retention_times.end(),
                                                                   performance_retention_times.begin(),
                                                                   performance_retention_times.end())), 1);
          writeDebug_("MSE between predicted and measured rt is: "
                      + String(Math::meanSquareError(final_predicted_retention_times.begin(),
                                                     final_predicted_retention_times.end(),
                                                     performance_retention_times.begin(),
                                                     performance_retention_times.end())), 1);
        }
      }
    }
    return EXECUTION_OK;
  } //end main_


private:
  inline double integrateGaussianPDFProduct_(double m1, double sd1, double m2, double sd2, double from, double to)
  {
    //TODO check for under/overflow and precision (e.g. with pi constant)
    double a = std::exp(-std::pow((m1 - m2),2.0) / (2.0 * (sd1*sd1 + sd2*sd2)));
    double erfdenom = (std::sqrt(2.0) * sd1 * sd2 * std::sqrt(sd1*sd1 + sd2*sd2));
    double b = -std::erf((from*sd1*sd1 - m2*sd1*sd1 + (from - m1)*sd2*sd2) / erfdenom);
    double c = std::erf((-m2*sd1*sd1 + sd1*sd1*to + sd2*sd2*(-m1 + to)) / erfdenom);
    double d = 2 * std::sqrt(2 * Constants::PI) * std::sqrt(sd1*sd1 + sd2*sd2);
    return a * (b + c) / d;
  }

  inline double integrateGaussianPDF_(double m, double sd, double from, double to)
  {
    double erfdenom = (std::sqrt(2.0) * sd);
    double a = -std::erf((m - to) / erfdenom);
    double b = std::erf((m - from) / erfdenom);
    return 0.5 * (a + b);
  }

};


int main(int argc, const char** argv)
{
  TOPPRTPredict tool;
  return tool.main(argc, argv);
}

/// @endcond
