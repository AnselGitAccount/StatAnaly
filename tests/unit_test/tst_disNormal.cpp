#include "gtest/gtest.h"
#include "../../include/density/disNormal.h"

namespace statanaly {

auto tst_disNormal_pdf(double x, double y, int z) {
    auto myNorm = disNormal(x, y);
    return myNorm.pdf(z);
};


TEST( Normal_Distribution_Tests, pdf ) {
    EXPECT_DOUBLE_EQ(0.10798193302637613,
        tst_disNormal_pdf(-2., 0.5, -1));
};


TEST( Normal_Distribution_Tests, hash ) {
    disNormal d1(3, 8);
    disNormal d2(3, 8);
    disNormal d3(5, 8);

    EXPECT_TRUE(d1.hash() == d2.hash());
    EXPECT_TRUE(d1.hash() != d3.hash());
};

}