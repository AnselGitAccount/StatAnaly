#include "gtest/gtest.h"
#include "../../include/density/disGamma.h"


namespace statanaly {

TEST( Gamma_Distribution, pdf ) {
    disGamma d1(0.5, 9.0);
    
    EXPECT_DOUBLE_EQ( 0.059540362609726359, d1.pdf(2) );
}

TEST( Gamma_Distribution, cdf ) {
    disGamma d1(4,2);

    // EXPECT_DOUBLE_EQ( 0.89758361765043326, d1.cdf(10) );
}



}