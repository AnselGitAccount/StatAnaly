#include "gtest/gtest.h"
#include "density/disNormal.h"

namespace statanaly {


TEST( Normal_Distribution_Tests, pdf ) {
    auto myNorm = disNormal(-2., 0.25);
    EXPECT_DOUBLE_EQ(0.10798193302637613, myNorm.pdf(-1.));
};

TEST( Normal_Distribution_Tests, cdf ) {
    auto myNorm = disNormal(-2., 0.25);
    EXPECT_DOUBLE_EQ(0.97724986805182079 , myNorm.cdf(-1.));
};

TEST( Normal_Distribution_Tests, hash ) {
    disNormal d1(3, 8);
    disNormal d2(3, 8);
    disNormal d3(5, 8);

    EXPECT_TRUE(d1.hash() == d2.hash());
    EXPECT_TRUE(d1.hash() != d3.hash());
};

}