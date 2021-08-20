#include "gtest/gtest.h"
#include "dCompare.h"
#include "density/disChiSq.h"
#include "density/disNormal.h"
#include "density/disMixture.h"


namespace statanaly {

TEST(fl_equality, isEqual_fl_tol) {
    // Test for equality of floating-point number within a numerical value.

    double a = 1/203.0;
    double b = 1/203.1;

    // 1e-6 < a-b < 1e-5
    const auto tol1 = pow(static_cast<long double>(10), -static_cast<long double>(5));
    const auto tol2 = pow(static_cast<long double>(10), -static_cast<long double>(6));

    EXPECT_TRUE(isEqual_fl_tol(a,b,tol1));
    EXPECT_FALSE(isEqual_fl_tol(a,b,tol2));

    // Exact same numbers should have a tolerance of zero.
    EXPECT_TRUE(isEqual_fl_tol(b,b,0));
}


TEST(fl_equality, isEqual_fl_ulp) {
    // Test for equality of floating-point number upto N ulp.
    
    // a & b is 1 ULP apart.
    // diff(a,b) is 8.88178e-16

    double a = 5;
    double o = sqrt(a);
    double b = o*o;

    // a-b == 1ULP
    EXPECT_FALSE(isEqual_fl_ulp(a,b,0));
    EXPECT_TRUE(isEqual_fl_ulp(a,b,1));

    // The IEEE standard says that any comparison operation involving
    // a NAN must return false.
    EXPECT_FALSE(isEqual_fl_ulp(sqrt(-1.0),a,0));

    // Exact same numbers should have a tolerance of zero.
    EXPECT_TRUE(isEqual_fl_ulp(a,a,0));
}


TEST(Distribution_comparison, isEqual) {
    // Test for equality of distributions.
    // bit-wise equal.

    disNormal d1(1., 2.);
    disNormal d2(2., 2.);

    EXPECT_TRUE(isEqual(d1,d1));
    EXPECT_FALSE(isEqual(d1,d2));

    // Test for equality of distributio mixtures.

    disMixture mix1;
    mix1.insert(d1, 1);
    mix1.insert(d2, 1);

    disMixture mix2;
    mix2.insert(d1, 1);
    mix2.insert(d2, 2);

    EXPECT_TRUE(isEqual(mix1,mix1));
    EXPECT_FALSE(isEqual(mix1,mix2));
}

TEST(Distribution_comparison, isEqual_tol) {
    // Test for equality of distributions.
    // Every parameter must match within a numerical value.

    disNormal d1(1/203.0, 2.);
    disNormal d2(1/203.1, 2.);

    const auto tol1 = pow(static_cast<long double>(10), -static_cast<long double>(5));
    const auto tol2 = pow(static_cast<long double>(10), -static_cast<long double>(6));

    // 1e-6 < difference < 1e-5
    EXPECT_TRUE(isEqual_tol(d1,d2,tol1));
    EXPECT_FALSE(isEqual_tol(d1,d2,tol2));
}


TEST(Distribution_comparison, isEqual_ulp) {
    // Test for equality of distributions.
    // Every parameter must match within N ulp.

    disNormal d1(5., 2.);
    disNormal d2(sqrt(5)*sqrt(5), 2.);

    // difference == 1ULP
    EXPECT_FALSE(isEqual_ulp(d1,d2,0));
    EXPECT_TRUE(isEqual_ulp(d1,d2,1));
}

TEST(Distribution_mixture_comparison, isEqual_tol) {
    // Test for equality of distribution mixture.
    // Every distribution in the mixture must equal.
    // And the weights should also be equal.

    disNormal d1(1/203.0, 2.);
    disNormal d2(1/203.1, 2.);
    
    disMixture mix1;
    mix1.insert(d1, 1);
    mix1.insert(d1, 1);
    
    disMixture mix2;
    mix2.insert(d2, 1);
    mix2.insert(d2, 1);

    const auto tol1 = pow(static_cast<long double>(10), -static_cast<long double>(5));
    const auto tol2 = pow(static_cast<long double>(10), -static_cast<long double>(6));

    EXPECT_TRUE(isEqual_tol(mix1, mix2, tol1));
    EXPECT_FALSE(isEqual_tol(mix1, mix2, tol2));
}


TEST(Distribution_mixture_comparison, isEqual_ulp) {
    // Test for equality of distribution mixture.
    // Every distribution in the mixture must equal.
    // And the weights should also be equal.

    disNormal d1(5., 2.);
    disNormal d2(sqrt(5)*sqrt(5), 2.);
    
    disMixture mix1;
    mix1.insert(d1, 1);
    mix1.insert(d1, 1);
    
    disMixture mix2;
    mix2.insert(d2, 1);
    mix2.insert(d2, 1);

    // difference == 1ULP
    EXPECT_FALSE(isEqual_ulp(mix1,mix2,0));
    EXPECT_TRUE(isEqual_ulp(mix1,mix2,1));
}

}   // namespace statanaly