#include "gtest/gtest.h"
#include "../../include/density/specialFunc.h"


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

}