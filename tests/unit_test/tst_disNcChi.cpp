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
#include "density/disNcChi.h"


namespace statanaly {

TEST( Noncentral_Chi_Distribution, pdf ) {
    // https://www.wolframalpha.com/input/?i=besseli%5B0%2C6%5D*exp%28-6.5%29*2

    disNcChi d1(2,3);
    
    EXPECT_FLOAT_EQ(0.2021656851300834280130675704254515255, d1.pdf(2) );
}

TEST( Noncentral_Chi_Distribution, cdf ) {
    // https://www.wolframalpha.com/input/?i=1-marcumq%5B1%2C3%2C2%5D
    
    disNcChi d1(2,3);

    EXPECT_DOUBLE_EQ(0.1132792455976077429640956095105139916, d1.cdf(2) );
}

}