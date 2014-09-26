// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: $
// $Authors: Hendrik Weisser, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>

#include <OpenMS/FORMAT/CsvFile.h>

///////////////////////////

START_TEST(TransformationModelInterpolated, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationModelInterpolated* ptr = 0;
TransformationModelInterpolated* nullPointer = 0;

TransformationModel::DataPoints dummy_data;
dummy_data.push_back(make_pair(0.0, 1.0));
dummy_data.push_back(make_pair(0.5, 4.0));
dummy_data.push_back(make_pair(1.0, 2.0));
dummy_data.push_back(make_pair(1.0, 4.0));

START_SECTION((TransformationModelInterpolated(const DataPoints &data, const Param &params)))
{
  ptr = new TransformationModelInterpolated(dummy_data, Param());
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~TransformationModelInterpolated()))
{
  delete ptr;
}
END_SECTION

START_SECTION((double evaluate(const double value) const))
{
  // NOTE: This test ensure first of all the compatibility of the interpolation
  //       to the original GSL based implementation of TransformationModelInterpolated

  /* load data generated by gsl interpolation
   */
  CsvFile data_points(OPENMS_GET_TEST_DATA_PATH("TransformationModelInterpolated_base_data.txt"), '\t');
  TransformationModel::DataPoints base_data;
  for (Size i = 0; i < data_points.size(); ++i)
  {
    StringList sl;
    if (data_points.getRow(i, sl))
    {
      base_data.push_back(make_pair(sl[0].toDouble(), sl[1].toDouble()));
    }
  }

  // create interpolations
  Param p_linear;
  TransformationModelInterpolated::getDefaultParameters(p_linear);
  p_linear.setValue("interpolation_type", "linear");
  TransformationModelInterpolated linear_interpolation(base_data, p_linear);

  Param p_cspline;
  TransformationModelInterpolated::getDefaultParameters(p_cspline);
  p_cspline.setValue("interpolation_type", "cspline");
  TransformationModelInterpolated cspline_interpolation(base_data, p_cspline);

  Param p_akima;
  TransformationModelInterpolated::getDefaultParameters(p_akima);
  p_akima.setValue("interpolation_type", "akima");
  TransformationModelInterpolated akima_interpolation(base_data, p_akima);


  // read gsl results
  CsvFile gsl_results(OPENMS_GET_TEST_DATA_PATH("TransformationModelInterpolated_gsl_data.txt"), '\t');
  std::vector<double> gsl_target_points;
  std::vector<double> gsl_cspline;
  std::vector<double> gsl_linear;
  std::vector<double> gsl_akima;

  for (Size i = 0; i < gsl_results.size(); ++i)
  {
    StringList sl;
    if (gsl_results.getRow(i, sl))
    {
      gsl_target_points.push_back(sl[0].toDouble());
      gsl_linear.push_back(sl[1].toDouble());
      gsl_akima.push_back(sl[2].toDouble());
      gsl_cspline.push_back(sl[3].toDouble());
    }
  }

  // test the interpolation
  for (size_t i = 0; i < gsl_target_points.size(); ++i)
  {
    const double x = gsl_target_points[i];
    TEST_REAL_SIMILAR(linear_interpolation.evaluate(x), gsl_linear[i])
    TEST_REAL_SIMILAR(cspline_interpolation.evaluate(x), gsl_cspline[i])
    TEST_REAL_SIMILAR(akima_interpolation.evaluate(x), gsl_akima[i])
  }
}
END_SECTION

START_SECTION(([EXTRA] TransformationModelInterpolated::evaluate() beyond the actual borders))
{
  // independent of the actual interpolation beyond the borders a linear extra-polation
  // based on TransformationModelLinear is performed
  ptr = new TransformationModelInterpolated(dummy_data, Param());

  // see TransformationModelLinear_test
  TEST_REAL_SIMILAR(ptr->evaluate(-0.5), 0.0);
  TEST_REAL_SIMILAR(ptr->evaluate(1.5), 4.0);
}

END_SECTION


START_SECTION((static void getDefaultParameters(Param & params)))
{
  Param p;
  TransformationModelInterpolated::getDefaultParameters(p);
  TEST_EQUAL(!p.empty(), true)
  TEST_EQUAL(p.getValue("interpolation_type"), "cspline")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST