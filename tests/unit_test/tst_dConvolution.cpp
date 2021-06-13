#include "gtest/gtest.h"
#include "../../include/dConvolution.h"
#include "../../include/density/disNormal.h"
#include "../../include/density/disUniform.h"
#include "../../include/density/disIrwinHall.h"


namespace statanaly {

TEST( dConvolution, double_StdUniform ) {
    /* Convolute a standard uniform and a standard uniform --> IrwinHall */

    disStdUniform su1{}, su2{};
    probDensFunc* rsu = cnvl.go(su1, su2);
    disIrwinHall* rih = dynamic_cast<disIrwinHall*>(rsu);  // Cast to base is optional.

    EXPECT_TRUE( rih->hash() == disIrwinHall{2}.hash() );

    // Also accept pointer-to-base.
    probDensFunc* su3 = new disStdUniform{};
    probDensFunc* su4 = new disStdUniform{};
    probDensFunc* rsu2 = cnvl.go(*su3, *su4);

    EXPECT_TRUE( rih->hash() == rsu2->hash() );

    delete rih, su3, su4, rsu2;
};


TEST( dConvolution, StdUniform ) {
    /* Convolute a list of standard uniform distributions --> Irwinhall */

    disStdUniform su1{}, su2{}, su3{}, su4{};
    probDensFunc* rsu = convolve<disStdUniform>({su1, su2, su3, su4});

    EXPECT_TRUE( rsu->hash() == disIrwinHall{4}.hash() );

    delete rsu;
};


TEST( dConvolution, double_Normal ) {
    /* Convolute a normal and a normal --> normal */

    disNormal n1{1,2}, n2{2,1};
    probDensFunc* rn = cnvl.go(n1,n2);
    disNormal expe{3,3};

    EXPECT_FLOAT_EQ( rn->mean(), expe.mean() );
    EXPECT_FLOAT_EQ( rn->variance(), expe.variance() );

    delete rn;
};


TEST( dConvolution, Normal ) {
    /* Convolute a list of normal distribution --> normal */

    disNormal n1{1,2}, n2{2,1}, n3{4.3, 2.3}, n4{1.2, 0.4};
    probDensFunc* rn = convolve<disNormal>({n1, n2, n3, n4});
    disNormal expe{8.5, 5.7};

    EXPECT_FLOAT_EQ( rn->mean(), expe.mean() );
    EXPECT_FLOAT_EQ( rn->variance(), expe.variance() );

    delete rn;
};

}