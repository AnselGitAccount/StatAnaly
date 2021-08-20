#include "gtest/gtest.h"
#include "density/disUniform.h"

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

TEST( Std_Uniform_Distribution, hash ) {
    disStdUniform d1{};
    disStdUniform d2{};

    EXPECT_TRUE(d1.hash() == d2.hash());
};

TEST( Uniform_Distribution, hash ) {
    disUniform d1(3, 8);
    disUniform d2(3, 8);
    disUniform d3(5, 8);

    EXPECT_TRUE(d1.hash() == d2.hash());
    EXPECT_TRUE(d1.hash() != d3.hash());
};


}