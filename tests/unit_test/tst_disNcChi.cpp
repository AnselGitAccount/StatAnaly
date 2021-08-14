#include "gtest/gtest.h"
#include "../../include/density/disNcChi.h"


namespace statanaly {

TEST( Noncentral_Chi_Distribution, pdf ) {
    // https://www.wolframalpha.com/input/?i=besseli%5B0%2C6%5D*exp%28-6.5%29*2

    disNcChi d1(2,3);
    
    EXPECT_FLOAT_EQ(0.2021656851300834280130675704254515255, d1.pdf(2) );
}

TEST( Noncentral_Chi_Distribution, cdf ) {
    // https://www.wolframalpha.com/input/?i=1-marcumq%5B1%2C3%2C2%5D
    
    disNcChi d1(2,3);

    EXPECT_DOUBLE_EQ(0.1132792455976077429640956095105139916, d1.cdf(2) );
}

}