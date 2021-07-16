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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/Utils/MZIndex.h>
///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <sstream>

using namespace OpenMS;
using namespace std;

// static_assert(OpenMS::Test::fulfills_rule_of_5<MSSpectrum>(), "Must fulfill rule of 5");
// static_assert(OpenMS::Test::fulfills_rule_of_6<MSSpectrum>(), "Must fulfill rule of 6");
// static_assert(OpenMS::Test::fulfills_fast_vector<MSSpectrum>(), "Must have fast vector semantics");
// static_assert(std::is_nothrow_move_constructible<MSSpectrum>::value, "Must have nothrow move constructible");

START_TEST(MSSpectrum, "$Id$")

    /////////////////////////////////////////////////////////////
    // Dummy peak data

    Peak1D p1;
    p1.setIntensity(1.0f);
    p1.setMZ(10.0);

    Peak1D p2;
    p2.setIntensity(2.0f);
    p2.setMZ(10.2);

    Peak1D p3;
    p3.setIntensity(3.0f);
    p3.setMZ(11.0);

    Peak1D p4;
    p4.setIntensity(3.0f);
    p4.setMZ(11.5);

    Peak1D p5;
    p5.setIntensity(3.0f);
    p5.setMZ(13.3);

    Peak1D p6;
    p6.setIntensity(3.0f);
    p6.setMZ(13.7);

    Peak1D p7;
    p7.setIntensity(3.0f);
    p7.setMZ(14.44444);

    Peak1D p8;
    p8.setIntensity(3.0f);
    p8.setMZ(15.00005);

    Peak1D p9;
    p9.setIntensity(3.0f);
    p9.setMZ(16.5);

    Peak1D p10;
    p10.setIntensity(3.0f);
    p10.setMZ(18.0);

    MSSpectrum s;
    s.push_back(p1);
    s.push_back(p2);
    s.push_back(p3);
    s.push_back(p4);
    s.push_back(p5);
    s.push_back(p6);
    s.push_back(p7);
    s.push_back(p8);
    s.push_back(p9);
    s.push_back(p10);

    MSExperiment e;
    e.addSpectrum(s);
    e.updateRanges();

    START_SECTION(MZIndex)
        MzIndex i(e,8);
        TEST_EQUAL(i.findNearest(0,0.0), 0);
        TEST_EQUAL(i.findNearest(0,10.0), 0);
        TEST_EQUAL(i.findNearest(0,11), 2);
        TEST_EQUAL(i.findNearest(0,12.5), 4);
        TEST_EQUAL(i.findNearest(0,18.0), 9);
        TEST_EQUAL(i.findNearest(0,19.0), 9);
    END_SECTION


END_TEST