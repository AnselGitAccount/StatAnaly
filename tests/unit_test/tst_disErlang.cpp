#include "gtest/gtest.h"
#include "density/disErlang.h"


namespace statanaly {

TEST( Erlang_Distribution, pdf ) {
    disErlang d1(3,2);
    
    double x = 4;
    double expected = pow(2,3) / 2 * pow(x,2) / exp(2*x);
    EXPECT_DOUBLE_EQ( expected, d1.pdf(x) );
}

TEST( Erlang_Distribution, cdf ) {
    disErlang d1(3,2);

    double x = 4;
    double expected = lowerGamma(3,2*x)/2;
    EXPECT_DOUBLE_EQ( expected, d1.cdf(x) );
}

}