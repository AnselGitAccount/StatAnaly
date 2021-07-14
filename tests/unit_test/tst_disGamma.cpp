#include "gtest/gtest.h"
#include "../../include/density/disGamma.h"


namespace statanaly {

TEST( Gamma_Distribution, pdf ) {
    constexpr double scale = 0.5, shape= 3.0;
    disGamma d1(scale, shape);
    
    EXPECT_DOUBLE_EQ( 0.2930502222197468846995, d1.pdf(2) );
}

TEST( Gamma_Distribution, cdf ) {
    constexpr double scale = 0.5, shape= 3.0;
    disGamma d1(scale, shape);

    EXPECT_NEAR( 0.7618966944464556561817, d1.cdf(2), 1e-5 );
}


}