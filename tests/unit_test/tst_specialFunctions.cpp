#include "gtest/gtest.h"
#include "../../include/density/specialFunc.h"


namespace statanaly {

TEST( Lower_Gamma_Function, computation ) {
    // Exact value is taken from https://keisan.casio.com/exec/system/1180573447

    EXPECT_FLOAT_EQ( 0.64664716763387308106, lowerGamma(3.0,2.0) );
}

TEST( Upper_Gamma_Function, computation ) {
    EXPECT_FLOAT_EQ( 1.35335283236612691894, upperGamma(3.0,2.0) );
}


TEST( Regularized_Lower_Gamma_Function, computation ) {
    // Exact value is taken from https://keisan.casio.com/exec/system/1180573447

    double expected = 0.64664716763387308106/(0.64664716763387308106+1.35335283236612691894);
    EXPECT_FLOAT_EQ( expected, regLowerGamma(3.0,2.0) );
}

TEST( Regularized_Upper_Gamma_Function, computation ) {
    double expected = 1.35335283236612691894/(0.64664716763387308106+1.35335283236612691894);
    EXPECT_FLOAT_EQ( expected, regUpperGamma(3.0,2.0) );
}

}