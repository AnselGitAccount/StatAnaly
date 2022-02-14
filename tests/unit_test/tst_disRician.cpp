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
#include "density/disRician.h"


namespace statanaly {

TEST( Rician_Distribution, pdf ) {
    // https://www.wolframalpha.com/input/?i=rice+distribution%282%2C3%29+pdf%284%29
    
    disRician d1(2,3);
    
    EXPECT_FLOAT_EQ(0.176667229995736222326, d1.pdf(4) );
}

TEST( Rician_Distribution, cdf ) {
    disRician d1(2,3);

    EXPECT_FLOAT_EQ(0.512532223392871104629, d1.cdf(4) );
}


}