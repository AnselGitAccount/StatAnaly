#include "gtest/gtest.h"
#include "../../include/density/disNormal.h"

namespace statanaly {

auto tst_disNormal_pdf(double x, double y, int z) {
    auto myNorm = disNormal(x, y);
    return myNorm.pdf(z);
};


TEST(Normal_Distribution_Tests, pdf) {
    EXPECT_DOUBLE_EQ(0.10798193302637613,
        tst_disNormal_pdf(-2., 0.5, -1));
};



}