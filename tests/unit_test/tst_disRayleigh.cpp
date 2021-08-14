#include "gtest/gtest.h"
#include "../../include/density/disRayleigh.h"


namespace statanaly {

TEST( Rayleigh_Distribution, pdf ) {
    disRayleigh d1(4);
    
    EXPECT_DOUBLE_EQ(0.02746058351462963582922, d1.pdf(10) );
}

TEST( Rayleigh_Distribution, cdf ) {
    disRayleigh d1(4);

    EXPECT_DOUBLE_EQ(0.9560630663765925826733, d1.cdf(10) );
}


}