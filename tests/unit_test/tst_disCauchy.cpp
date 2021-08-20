#include "gtest/gtest.h"
#include "density/disCauchy.h"


namespace statanaly {

TEST( Cauchy_Distribution, pdf ) {
    disCauchy d1(4,2);
    
    EXPECT_DOUBLE_EQ( 0.015915494309189534, d1.pdf(10) );
}

TEST( Cauchy_Distribution, cdf ) {
    disCauchy d1(4,2);

    EXPECT_DOUBLE_EQ( 0.89758361765043326, d1.cdf(10) );
}



}