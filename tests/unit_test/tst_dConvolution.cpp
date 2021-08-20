/* Test of convolution of probability distribution functions.
 *
 * (1) Sum of two Random Variables 
 *      R = X + Y
 *
 * (2) Sum of a list of Random Variabales
 *      R = A + B + C + D + ...
 */


#include "gtest/gtest.h"
#include "dConvolution.h"
#include "density/disNormal.h"
#include "density/disUniform.h"
#include "density/disIrwinHall.h"


namespace statanaly {

TEST( dConvolution, two_StdUniform ) {
    /* a standard uniform + a standard uniform --> IrwinHall */

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


TEST( dConvolution, list_StdUniform ) {
    /* Sum of a list of standard uniform RVs --> Irwinhall */

    disStdUniform su1{}, su2{}, su3{}, su4{};
    probDensFunc* rsu = convolve<disStdUniform>({su1, su2, su3, su4});

    EXPECT_TRUE( rsu->hash() == disIrwinHall{4}.hash() );

    delete rsu;
};


TEST( dConvolution, two_Normal ) {
    /* a normal + a normal --> normal */

    disNormal n1{1,2}, n2{2,1};
    probDensFunc* rn = cnvl.go(n1,n2);
    disNormal expe{3,3};

    EXPECT_FLOAT_EQ( rn->mean(), expe.mean() );
    EXPECT_FLOAT_EQ( rn->variance(), expe.variance() );

    delete rn;
};


TEST( dConvolution, list_Normal ) {
    /* Sum of a list of normal RVs --> normal */

    disNormal n1{1,2}, n2{2,1}, n3{4.3, 2.3}, n4{1.2, 0.4};
    probDensFunc* rn = convolve<disNormal>({n1, n2, n3, n4});
    disNormal expe{8.5, 5.7};

    EXPECT_FLOAT_EQ( rn->mean(), expe.mean() );
    EXPECT_FLOAT_EQ( rn->variance(), expe.variance() );

    delete rn;
};


TEST( dConvolution, two_Cauchy ) {
    /* a cauchy + a cauchy --> cauchy */

    disCauchy n1{1,2}, n2{2,1};
    probDensFunc* rn = cnvl.go(n1,n2);
    disCauchy* rc = dynamic_cast<disCauchy*>(rn);
    disCauchy expe{3,3};

    EXPECT_EQ( rc->ploc(), expe.ploc() );
    EXPECT_EQ( rc->pscale(), expe.pscale() );

    delete rn;
};


TEST( dConvolution, list_Cauchy ) {
    /* Sum of a list of Cauchy RVs --> Cauchy */

    disCauchy n1{1,2}, n2{2,1}, n3{4.3, 2.3}, n4{1.2, 0.4};
    probDensFunc* rn = convolve<disCauchy>({n1, n2, n3, n4});
    disCauchy* rc = dynamic_cast<disCauchy*>(rn);
    disCauchy expe{8.5, 5.7};

    EXPECT_EQ( rc->ploc(), expe.ploc() );
    EXPECT_EQ( rc->pscale(), expe.pscale() );

    delete rn;
};

TEST( dConvolution, two_Gamma ) {
    /* a Gamma + a Gamma --> Gamma */

    disGamma n1{1,2}, n2{1,1};
    probDensFunc* rn = cnvl.go(n1,n2);
    disGamma* rc = dynamic_cast<disGamma*>(rn);
    disGamma expe{1,3};

    EXPECT_EQ( rc->pscale(), expe.pscale() );
    EXPECT_EQ( rc->pshape(), expe.pshape() );

    delete rn;
};


TEST( dConvolution, list_Gamma ) {
    /* Sum of a list of Gamma RVs --> Gamma */

    disGamma n1{1,2}, n2{1,1}, n3{1., 2.3}, n4{1., 0.4};
    probDensFunc* rn = convolve<disGamma>({n1, n2, n3, n4});
    disGamma* rc = dynamic_cast<disGamma*>(rn);
    disGamma expe{1., 5.7};

    EXPECT_EQ( rc->pscale(), expe.pscale() );
    EXPECT_EQ( rc->pshape(), expe.pshape() );

    delete rn;
};


TEST( dConvolution, two_Exponential ) {
    /* a Exponential + a Exponential --> Exponential */

    disExponential n1{4}, n2{4};
    probDensFunc* rn = cnvl.go(n1,n2);
    disErlang* rc = dynamic_cast<disErlang*>(rn);
    disErlang expe{2,4};

    EXPECT_EQ( rc->pshape(), expe.pshape() );
    EXPECT_EQ( rc->prate(), expe.prate() );

    delete rn;
};


TEST( dConvolution, list_Exponential ) {
    /* Sum of a list of Exponential RVs --> Exponential */

    disExponential n1{4}, n2{4}, n3{4}, n4{4};
    probDensFunc* rn = convolve<disExponential>({n1, n2, n3, n4});
    disErlang* rc = dynamic_cast<disErlang*>(rn);
    disErlang expe{4,4};

    EXPECT_EQ( rc->pshape(), expe.pshape() );
    EXPECT_EQ( rc->prate(), expe.prate() );

    delete rn;
};

}