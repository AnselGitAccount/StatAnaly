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
#include "density/disUniform.h"

namespace statanaly {


TEST( Std_Uniform_Distribution, pdf ) {
    disStdUniform distr{};
    EXPECT_DOUBLE_EQ( 0 , distr.pdf(-1.) );
    EXPECT_DOUBLE_EQ( 1., distr.pdf(0.4) );
    EXPECT_DOUBLE_EQ( 0 , distr.pdf(3.1) );
}

TEST( Std_Uniform_Distribution, cdf ) {
    disStdUniform distr{};
    EXPECT_DOUBLE_EQ( 0  , distr.cdf(-1.) );
    EXPECT_DOUBLE_EQ( 0.4, distr.cdf(0.4) );
    EXPECT_DOUBLE_EQ( 1  , distr.cdf(3.1) );
}

TEST( Uniform_Distribution, pdf ) {
    disUniform distr(3, 5);
    EXPECT_DOUBLE_EQ( 0  , distr.pdf(-1. ) );
    EXPECT_DOUBLE_EQ( 0.5, distr.pdf(3.45) );
    EXPECT_DOUBLE_EQ( 0  , distr.pdf(6.12) );
}

TEST( Uniform_Distribution, cdf ) {
    disUniform distr(3, 5);
    EXPECT_DOUBLE_EQ( 0    , distr.cdf(1.  ) );
    EXPECT_DOUBLE_EQ( 0.225, distr.cdf(3.45) );
    EXPECT_DOUBLE_EQ( 1    , distr.cdf(6.12) );
}

TEST( Std_Uniform_Distribution, hash ) {
    disStdUniform d1{};
    disStdUniform d2{};

    EXPECT_TRUE(d1.hash() == d2.hash());
};

TEST( Uniform_Distribution, hash ) {
    disUniform d1(3, 8);
    disUniform d2(3, 8);
    disUniform d3(5, 8);

    EXPECT_TRUE(d1.hash() == d2.hash());
    EXPECT_TRUE(d1.hash() != d3.hash());
};


}