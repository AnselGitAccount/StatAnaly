#include "gtest/gtest.h"
#include "density/disChi.h"


namespace statanaly {


TEST( Chi_Distribution_Tests, pdf ) {
    auto myChi = disChi(8);
    auto expe = (6*6*6*6*6*6*6)*exp(-18)/(8*6);
    EXPECT_DOUBLE_EQ(expe, myChi.pdf(6));
};

TEST( Chi_Distribution_Tests, cdf ) {
    auto myChi = disChi(9);
    auto expe = 11.63126723822420996976 / (11.63126723822420996976 + 4.611583432389593884654E-4);
    EXPECT_DOUBLE_EQ(expe, myChi.cdf(6));
};


}
