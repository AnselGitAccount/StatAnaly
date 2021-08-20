#include "gtest/gtest.h"
#include "density/disNcChiSq.h"


namespace statanaly {

TEST( Noncentral_Chi_Sq_Distribution, pdf ) {
    // https://www.wolframalpha.com/input/?i=PDF%5Bnoncentralchisquaredistribution%5B2%2C3%5D%2C+2%5D

    disNcChiSq d1(2,3);
    
    EXPECT_DOUBLE_EQ(0.1299236871288784895691920805701871785, d1.pdf(2) );
}

TEST( Noncentral_Chi_Sq_Distribution, cdf ) {
    // https://www.wolframalpha.com/input/?i=CDF%5Bnoncentralchisquaredistribution%5B2%2C3%5D%2C+2%5D

    disNcChiSq d1(2,3);

    EXPECT_DOUBLE_EQ(0.2522069424260390183594385903034531831, d1.cdf(2) );
}

}