#include "gtest/gtest.h"
#include "../../include/density/disChiSq.h"


namespace statanaly {


TEST( ChiSquare_Distribution_Tests, pdf ) {
    auto myChiSq = disChiSq(9);
    EXPECT_DOUBLE_EQ(0.10008447084954811, myChiSq.pdf(6));
};

TEST( ChiSquare_Distribution_Tests, cdf ) {
    auto myChiSq = disChiSq(9);
    EXPECT_NEAR(0.26008170790534629, myChiSq.cdf(6), 1e-7);
};

TEST( ChiSquare_Distribution_Tests, hash ) {
    // Distribution with the same parameters must be identical.

    disChiSq d1(3);
    disChiSq d2(3);
    disChiSq d3(5);

    EXPECT_TRUE(d1.hash() == d2.hash());
    EXPECT_TRUE(d1.hash() != d3.hash());
};


}