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
#include "density/disErlang.h"


namespace statanaly {

TEST( Erlang_Distribution, pdf ) {
    disErlang d1(3,2);
    
    double x = 4;
    double expected = pow(2,3) / 2 * pow(x,2) / exp(2*x);
    EXPECT_DOUBLE_EQ( expected, d1.pdf(x) );
}

TEST( Erlang_Distribution, cdf ) {
    disErlang d1(3,2);

    double x = 4;
    double expected = lowerGamma(3,2*x)/2;
    EXPECT_DOUBLE_EQ( expected, d1.cdf(x) );
}

}