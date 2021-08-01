#include "gtest/gtest.h"
#include "../../include/density/disNcChi.h"


namespace statanaly {

TEST( Noncentral_Chi_Distribution, pdf ) {
    // Expected value is calculated by hand.

    disNcChi d1(2,3,4);
    
    EXPECT_DOUBLE_EQ(0.0862211699813562513999409928676656829, d1.pdf(2) );
}

TEST( Noncentral_Chi_Distribution, cdf ) {
    disNcChi d1(2,3,4);

    EXPECT_DOUBLE_EQ(0.1132792455976077429640956095105139916, d1.cdf(2) );
}

}