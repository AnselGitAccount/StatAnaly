#include "gtest/gtest.h"
#include "../../include/density/disRician.h"


namespace statanaly {

TEST( Rician_Distribution, pdf ) {
    // https://www.wolframalpha.com/input/?i=rice+distribution%282%2C3%29+pdf%284%29
    
    disRician d1(2,3);
    
    EXPECT_FLOAT_EQ(0.176667229995736222326, d1.pdf(4) );
}

TEST( Rician_Distribution, cdf ) {
    disRician d1(2,3);

    EXPECT_FLOAT_EQ(0.512532223392871104629, d1.cdf(4) );
}

// TEST( Rician_Distribution, convert_to_Noncentral_Chi ) {

//     auto Ray = disRician(1.5);
//     auto con = static_cast<disNcChi>(Ray);
//     auto Chi = disNcChi(2,1.5);

//     EXPECT_DOUBLE_EQ(Ray.pdf(6), con.pdf(6));
//     EXPECT_DOUBLE_EQ(Chi.pdf(6), con.pdf(6));
// }

}