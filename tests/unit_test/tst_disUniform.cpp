#include "gtest/gtest.h"
#include "../../include/density/disUniform.h"

namespace statanaly {


TEST( Std_Uniform_Distribution, pdf ) {
    disStdUniform distr{};
    EXPECT_DOUBLE_EQ( 0 , distr.pdf(-1.) );
    EXPECT_DOUBLE_EQ( 1., distr.pdf(0.4) );
    EXPECT_DOUBLE_EQ( 0 , distr.pdf(3.1) );
}

TEST( Std_Uniform_Distribution, cdf ) {
    disStdUniform distr{};
    EXPECT_DOUBLE_EQ( 0  , distr.cdf(-1.) );
    EXPECT_DOUBLE_EQ( 0.4, distr.cdf(0.4) );
    EXPECT_DOUBLE_EQ( 1  , distr.cdf(3.1) );
}

TEST( Uniform_Distribution, pdf ) {
    disUniform distr(3, 5);
    EXPECT_DOUBLE_EQ( 0  , distr.pdf(-1. ) );
    EXPECT_DOUBLE_EQ( 0.5, distr.pdf(3.45) );
    EXPECT_DOUBLE_EQ( 0  , distr.pdf(6.12) );
}

TEST( Uniform_Distribution, cdf ) {
    disUniform distr(3, 5);
    EXPECT_DOUBLE_EQ( 0    , distr.cdf(1.  ) );
    EXPECT_DOUBLE_EQ( 0.225, distr.cdf(3.45) );
    EXPECT_DOUBLE_EQ( 1    , distr.cdf(6.12) );
}


}