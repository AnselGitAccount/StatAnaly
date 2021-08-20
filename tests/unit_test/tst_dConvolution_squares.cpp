/* Test of convolution of probability distribution functions.
 *
 * (1) Sum of Squares of two Random Variables 
 *      R = X^2 + Y^2
 *
 * (2) Sum of a list of Random Variabales
 *      R = A^2 + B^2 + C^2 + D^2 + ...
 * 
 * (3) Sum of Squares of two Random Variables, then take the square-root.
 *      R = sqrt(X^2 + Y^2)
 * 
 * (4) Sum of Squares of a list of Random Variables, then take the square-root.
 *      R = sqrt(A^2 + B^2 + C^2 + D^2 + ...)
 */


#include "gtest/gtest.h"
#include "../../include/dConvolution.h"
#include "../../include/density/disNormal.h"
#include "../../include/density/disRayleigh.h"
#include "../../include/density/disRician.h"

namespace statanaly {

TEST( dConvolution_Square, two_Normal_zero_mean_one_std ) {
    /* Let A be a RV of N(0,1),
       A**2 + A**2 --> Central Chi Square RV */

    disNormal a(0,1);
    disNormal b(0,1);
    disChiSq s = *static_cast<disChiSq*>(cnvlSq.go(a,b));

    disChiSq e(2);

    EXPECT_TRUE( e.isEqual_ulp(s,0) );
}

TEST( dConvolution_Square, list_Normal_zero_mean_one_std ) {
    /* Let {A,B,C...} be a RV of N(0,1),
       A**2 + B**2 + C**2 + ... --> Central Chi Square RV */

    disNormal a(0,1), b(0,1), c(0,1);
    disChiSq s = *static_cast<disChiSq*>(convolveSq<disNormal>({a,b,c}));

    disChiSq e(3);

    EXPECT_TRUE( e.isEqual_ulp(s,0) );
}

TEST( dConvolution_Square, two_Normal_one_std ) {
    /* Let A be a RV of N(a,1),
       Let B be a RV of N(b,1),
       A**2 + B**2 --> Noncentral Chi Square RV */

    disNormal a(3,1);
    disNormal b(4,1);
    disNcChiSq s = *static_cast<disNcChiSq*>(cnvlSq.go(a,b));

    disNcChiSq e(2,3*3+4*4);

    EXPECT_TRUE( e.isEqual_ulp(s,0) );
}

TEST( dConvolution_Square, list_Normal_one_std ) {
    /* Let A be a RV of N(a,1),
       Let B be a RV of N(b,1),
       Let C be ...
       A**2 + B**2 + C**2 + ... --> Noncentral Chi Square RV */

    disNormal a(3,1), b(4,1), c(5,1);
    disNcChiSq s = *static_cast<disNcChiSq*>(convolveSq<disNormal>({a,b,c}));

    disNcChiSq e(3, 3*3+4*4+5*5);

    EXPECT_TRUE( e.isEqual_ulp(s,0) );
}

TEST( dConvolution_SSqrt, list_Normal_zero_mean_one_variance ) {
    /* Let {A,B,C...} be RVs of N(0,1),
       sqrt(A**2 + B**2 + C**2 + ...) --> Central Chi RV */
    
    disNormal a(0,1), b(0,1), c(0,1);
    disChi s = *static_cast<disChi*>(convolveSSqrt<disNormal>({a,b,c}));
    
    disChi e(3);

    EXPECT_TRUE( e.isEqual_ulp(s,0) );
}

TEST( dConvolution_SSqrt, list_Normal_one_variance ) {
    /* Let A be a RV of N(a,1),
       Let B be a RV of N(b,1),
       Let C be ...
       sqrt(A**2 + B**2 + C**2 + ...) --> Noncentral Chi RV */
    
    disNormal a(4,1), b(3,1), c(2,1);
    disNcChi s = *static_cast<disNcChi*>(convolveSSqrt<disNormal>({a,b,c}));
    
    disNcChi e(3, sqrt(4*4+3*3+2*2));

    EXPECT_TRUE( e.isEqual_ulp(s,10) );
}

TEST( dConvolution_SSqrt, two_Normal_zero_mean ) {
    /* Let A be a RV of N(0,sig**2),
       sqrt(A**2 + A**2) --> Rayleigh RV */

    disNormal a(0,9);
    disNormal b(0,9);
    disRayleigh s = *static_cast<disRayleigh*>(cnvlSSqrt.go(a,b));

    disRayleigh e(sqrt(9));

    EXPECT_TRUE( e.isEqual_ulp(s,0) );
}

TEST( dConvolution_SSqrt, two_Normal ) {
    /* Let A be a RV of N(a,sig**2),
       Let B be a RV of N(b,sig**2),
       sqrt(A**2 + B**2) --> Rician RV */

    disNormal a(2,9);
    disNormal b(3,9);
    disRician s = *static_cast<disRician*>(cnvlSSqrt.go(a,b));

    disRician e(sqrt(2*2+3*3), sqrt(9));

    EXPECT_TRUE( e.isEqual_ulp(s,0) );
}

}