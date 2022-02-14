/*
   Copyright 2022, Ansel Blumers

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#include "gtest/gtest.h"
#include "density/disChiSq.h"


namespace statanaly {


TEST( ChiSquare_Distribution_Tests, pdf ) {
    auto myChiSq = disChiSq(9);
    EXPECT_DOUBLE_EQ(0.10008447084954811, myChiSq.pdf(6));
};

TEST( ChiSquare_Distribution_Tests, cdf ) {
    auto myChiSq = disChiSq(9);
    EXPECT_NEAR(0.26008170790534629, myChiSq.cdf(6), 1e-7);
};

TEST( ChiSquare_Distribution_Tests, hash ) {
    // Distribution with the same parameters must be identical.

    disChiSq d1(3);
    disChiSq d2(3);
    disChiSq d3(5);

    EXPECT_TRUE(d1.hash() == d2.hash());
    EXPECT_TRUE(d1.hash() != d3.hash());
};


}