#include "gtest/gtest.h"
#include "../../include/density/disChiSq.h"


namespace statanaly {

auto tst_disChiSq_pdf(int x, int y) {
    auto myChiSq = disChiSq(x);
    return myChiSq.pdf(y);
};



TEST(ChiSquare_Distribution_Tests, pdf) {
    EXPECT_DOUBLE_EQ(0.10008447084954811,
        tst_disChiSq_pdf(9, 6));
};




}