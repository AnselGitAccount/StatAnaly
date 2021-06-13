#include "gtest/gtest.h"
#include "../../include/density/disIrwinHall.h"


namespace statanaly {

TEST( Irwin_Hall_Distribution, pdf ) {
    /* n=1 follows a uniform distribution.
        Except that on the right edge x=1, Irwin_Hall's pdf gives zero.
        Is this the desired behavior?? */
    disIrwinHall d1(1);
    for (double x=-0.01; x<1.1; x+=0.01) {
        if (0<=x && x<1)
            EXPECT_DOUBLE_EQ( 1, d1.pdf(x) );
        else 
            EXPECT_DOUBLE_EQ( 0, d1.pdf(x) );
    }

    /* n=2 follows a triangular distribution. */
    disIrwinHall d2(2);
    for (double x=-0.01; x<2.1; x+=0.01) {
        if (0<=x && x<=1)
            EXPECT_DOUBLE_EQ( x, d2.pdf(x) );
        else if (1<x && x<=2)
            EXPECT_DOUBLE_EQ( 2-x, d2.pdf(x) );
        else 
            EXPECT_DOUBLE_EQ( 0, d2.pdf(x) );
    }
    
    /* n=3 */
    disIrwinHall d3(3);
    for (double x=-0.01; x<3.1; x+=0.01) {
        if (0<=x && x<=1)
            EXPECT_FLOAT_EQ( .5*x*x, d3.pdf(x) );
        else if (1<x && x<=2)
            EXPECT_FLOAT_EQ( .5*(-2*x*x+6*x-3), d3.pdf(x) );
        else if (2<x && x<=3)
            EXPECT_NEAR( .5*(x-3)*(x-3), d3.pdf(x), 1e-15 );
        else 
            EXPECT_NEAR( 0, d3.pdf(x), 1e-15 );
    }
    
    /* n=4 */
    disIrwinHall d4(4);
    for (double x=-0.01; x<4.1; x+=0.01) {
        if (0<=x && x<=1)
            EXPECT_FLOAT_EQ( 1./6*x*x*x, d4.pdf(x) );
        else if (1<x && x<=2)
            EXPECT_FLOAT_EQ( 1./6*(-3*x*x*x+12*x*x-12*x+4), d4.pdf(x) );
        else if (2<x && x<=3)
            EXPECT_FLOAT_EQ( 1./6*(3*x*x*x-24*x*x+60*x-44), d4.pdf(x) );
        else if (3<x && x<=4)
            EXPECT_NEAR( 1./6*(4-x)*(4-x)*(4-x), d4.pdf(x), 1e-14 );
        else 
            EXPECT_NEAR( 0, d4.pdf(x), 1e-14 );
    }
}

TEST( Irwin_Hall_Distribution, cdf ) {
    /* cdf is just pdf divided by n */
    disIrwinHall d1(1);
    for (double x=-0.01; x<1.1; x+=0.01) {
        if (0<=x && x<1)
            EXPECT_DOUBLE_EQ( 1, d1.cdf(x) );
        else 
            EXPECT_DOUBLE_EQ( 0, d1.cdf(x) );
    }

    disIrwinHall d2(2);
    for (double x=-0.01; x<2.1; x+=0.01) {
        if (0<=x && x<=1)
            EXPECT_DOUBLE_EQ( x/2, d2.cdf(x) );
        else if (1<x && x<=2)
            EXPECT_DOUBLE_EQ( (2-x)/2, d2.cdf(x) );
        else 
            EXPECT_DOUBLE_EQ( 0, d2.cdf(x) );
    }
    
    disIrwinHall d3(3);
    for (double x=-0.01; x<3.1; x+=0.01) {
        if (0<=x && x<=1)
            EXPECT_FLOAT_EQ( .5*x*x/3, d3.cdf(x) );
        else if (1<x && x<=2)
            EXPECT_FLOAT_EQ( .5*(-2*x*x+6*x-3)/3, d3.cdf(x) );
        else if (2<x && x<=3)
            EXPECT_NEAR( .5*(x-3)*(x-3)/3, d3.cdf(x), 1e-15 );
        else 
            EXPECT_NEAR( 0, d3.cdf(x), 1e-15 );
    }
    
    disIrwinHall d4(4);
    for (double x=-0.01; x<4.1; x+=0.01) {
        if (0<=x && x<=1)
            EXPECT_FLOAT_EQ( 1./6*x*x*x, d4.pdf(x) );
        else if (1<x && x<=2)
            EXPECT_FLOAT_EQ( 1./6*(-3*x*x*x+12*x*x-12*x+4), d4.pdf(x) );
        else if (2<x && x<=3)
            EXPECT_FLOAT_EQ( 1./6*(3*x*x*x-24*x*x+60*x-44), d4.pdf(x) );
        else if (3<x && x<=4)
            EXPECT_NEAR( 1./6*(4-x)*(4-x)*(4-x), d4.pdf(x), 1e-14 );
        else 
            EXPECT_NEAR( 0, d4.pdf(x), 1e-14 );
    }
}

}