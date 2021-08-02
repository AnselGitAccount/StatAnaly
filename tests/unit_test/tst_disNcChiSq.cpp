#include "gtest/gtest.h"
#include "../../include/density/disNcChiSq.h"


namespace statanaly {

TEST( Noncentral_Chi_Sq_Distribution, pdf ) {
    // Expected value is calculated by hand.

    disNcChiSq d1(2,3,4);
    
    EXPECT_DOUBLE_EQ(0.0268863894563691405995707330958536593, d1.pdf(2) );
}

TEST( Noncentral_Chi_Sq_Distribution, cdf ) {
    disNcChiSq d1(2,3,4);

    EXPECT_DOUBLE_EQ(0.2522069424260390183594385903034531831, d1.cdf(2) );
}

}