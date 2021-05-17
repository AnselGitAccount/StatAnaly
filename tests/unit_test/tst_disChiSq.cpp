#include "gtest/gtest.h"
#include "../../include/density/disChiSq.h"


namespace statanaly {

auto tst_disChiSq_pdf(int x, int y) {
    auto myChiSq = disChiSq(x);
    return myChiSq.pdf(y);
};

TEST( ChiSquare_Distribution_Tests, pdf ) {
    EXPECT_DOUBLE_EQ(0.10008447084954811,
        tst_disChiSq_pdf(9, 6));
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