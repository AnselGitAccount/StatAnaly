#include "gtest/gtest.h"
#include "density/disExponential.h"


namespace statanaly {

TEST( Exponential_Distribution, pdf ) {
    disExponential d1(1.);

    EXPECT_DOUBLE_EQ(0.36787944117144233, d1.pdf(1.));
}

TEST( Exponential_Distribution, cdf ) {
    disExponential d1(1.);

    EXPECT_DOUBLE_EQ(0.63212055882855767, d1.cdf(1.));
}


}