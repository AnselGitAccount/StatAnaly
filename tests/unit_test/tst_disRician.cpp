#include "gtest/gtest.h"
#include "density/disRician.h"


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


}