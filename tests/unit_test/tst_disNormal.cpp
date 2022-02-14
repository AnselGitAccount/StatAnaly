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
#include "density/disNormal.h"

namespace statanaly {


TEST( Normal_Distribution_Tests, pdf ) {
    auto myNorm = disNormal(-2., 0.25);
    EXPECT_DOUBLE_EQ(0.10798193302637613, myNorm.pdf(-1.));
};

TEST( Normal_Distribution_Tests, cdf ) {
    auto myNorm = disNormal(-2., 0.25);
    EXPECT_DOUBLE_EQ(0.97724986805182079 , myNorm.cdf(-1.));
};

TEST( Normal_Distribution_Tests, hash ) {
    disNormal d1(3, 8);
    disNormal d2(3, 8);
    disNormal d3(5, 8);

    EXPECT_TRUE(d1.hash() == d2.hash());
    EXPECT_TRUE(d1.hash() != d3.hash());
};

}