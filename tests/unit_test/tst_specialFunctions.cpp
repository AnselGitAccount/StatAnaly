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
#include "density/specialFunc.h"


namespace statanaly {

TEST( Lower_Gamma_Function, computation ) {
    // Exact value is taken from https://keisan.casio.com/exec/system/1180573447

    EXPECT_DOUBLE_EQ( 5.745719328049896026173, lowerGamma(4,8) );
    EXPECT_DOUBLE_EQ( 257.713423595950588589 , lowerGamma(8,4) );
}

TEST( Upper_Gamma_Function, computation ) {
    EXPECT_FLOAT_EQ( 0.2542806719501039738266, upperGamma(4,8) );
    EXPECT_FLOAT_EQ( 4782.286576404049411411,  upperGamma(8,4) );
}


TEST( Regularized_Lower_Gamma_Function, computation ) {
    // Exact value is taken from https://keisan.casio.com/exec/system/1180573447

    constexpr double expected = 5.745719328049896026173/(5.745719328049896026173+0.2542806719501039738266);
    EXPECT_DOUBLE_EQ( expected, regLowerGamma(4,8) );

    constexpr double expected2 = 257.713423595950588589/(257.713423595950588589+4782.286576404049411411);
    EXPECT_DOUBLE_EQ( expected2, regLowerGamma(8,4) );
}

TEST( Regularized_Upper_Gamma_Function, computation ) {
    constexpr double expected = 0.2542806719501039738266/(5.745719328049896026173+0.2542806719501039738266);
    EXPECT_FLOAT_EQ( expected, regUpperGamma(4,8) );

    constexpr double expected2 = 4782.286576404049411411/(257.713423595950588589+4782.286576404049411411);
    EXPECT_DOUBLE_EQ( expected2, regUpperGamma(8,4) );
}

TEST( MarcumQ_Function, integer_M ) {
    // https://www.wolframalpha.com/input/?i=ScientificForm%28marcumq%5B3%2C1.3%2C1.5%5D%29

    EXPECT_FLOAT_EQ( 0.024753758803239678448239696602465, marcumQ(1, 13 , 15 ) );
    EXPECT_FLOAT_EQ( 0.580034435214153658566764755380761, marcumQ(1, 1.3, 1.5) );
    EXPECT_FLOAT_EQ( 0.944047836736561254878101527678606, marcumQ(3, 1.3, 1.5) );
}

TEST( MarcumQ_Function, non_integer_M ) {
    // https://www.wolframalpha.com/input/?i=ScientificForm%28marcumq%5B1.4%2C1.3%2C1.5%5D%29

    EXPECT_FLOAT_EQ( 0.692728870227774340018922937858990, marcumQ(1.4, 1.3, 1.5) );
}


}