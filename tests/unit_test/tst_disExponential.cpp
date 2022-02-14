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
#include "density/disExponential.h"


namespace statanaly {

TEST( Exponential_Distribution, pdf ) {
    disExponential d1(1.);

    EXPECT_DOUBLE_EQ(0.36787944117144233, d1.pdf(1.));
}

TEST( Exponential_Distribution, cdf ) {
    disExponential d1(1.);

    EXPECT_DOUBLE_EQ(0.63212055882855767, d1.cdf(1.));
}


}