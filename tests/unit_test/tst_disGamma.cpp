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
#include "density/disGamma.h"


namespace statanaly {

TEST( Gamma_Distribution, pdf ) {
    constexpr double scale = 0.5, shape= 3.0;
    disGamma d1(scale, shape);
    
    EXPECT_DOUBLE_EQ( 0.2930502222197468846995, d1.pdf(2) );
}

TEST( Gamma_Distribution, cdf ) {
    constexpr double scale = 0.5, shape= 3.0;
    disGamma d1(scale, shape);

    EXPECT_NEAR( 0.7618966944464556561817, d1.cdf(2), 1e-5 );
}


}