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
#include "density/disChi.h"


namespace statanaly {


TEST( Chi_Distribution_Tests, pdf ) {
    auto myChi = disChi(8);
    auto expe = (6*6*6*6*6*6*6)*exp(-18)/(8*6);
    EXPECT_DOUBLE_EQ(expe, myChi.pdf(6));
};

TEST( Chi_Distribution_Tests, cdf ) {
    auto myChi = disChi(9);
    auto expe = 11.63126723822420996976 / (11.63126723822420996976 + 4.611583432389593884654E-4);
    EXPECT_DOUBLE_EQ(expe, myChi.cdf(6));
};


}
