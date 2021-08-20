#include <iostream>
#include "gtest/gtest.h"

#include "density/disMixture.h"
#include "density/disNormal.h"
#include "density/disUniform.h"
#include "density/disIrwinHall.h"
#include "density/disChiSq.h"
#include "density/disCauchy.h"
#include "density/disGamma.h"
#include "density/disExponential.h"
#include "density/disChi.h"
#include "density/disRayleigh.h"
#include "density/disRician.h"
#include "density/disNcChi.h"
#include "density/disNcChiSq.h"
#include "dContainer.h"
#include "fl_comparison.h"
#include "dCompare.h"
#include "dConvolution.h"

using namespace statanaly;



int main(int argc, char** argv) {

    // auto d1 = disExponential(1.0);
    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << d1.pdf(1) << std::endl;
    // std::cout << 1/exp(1) << std::endl;
    // std::cout << d1.cdf(1) << std::endl;
    // std::cout << 1-1/exp(1) << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);

    // disNormal a(2,9);
    // disNormal b(3,9);
    // disRician s = *static_cast<disRician*>(cnvlSSqrt.go(a,b));

    // disRician e(sqrt(2*2+3*3), sqrt(9));
    
    // printf("%g %g %g \n", s.p_distance(), s.p_scale(), s.pdf(2));
    // printf("%g %g %g \n", e.p_distance(), s.p_scale(), e.pdf(2));
    // printf("%u\n", e.isEqual_ulp(s,0));
    // printf("%u\n", e.hash()==s.hash());
    
    
    disNormal a(3,1);
    disNormal b(4,1);
    disNcChiSq s = *static_cast<disNcChiSq*>(cnvlSq.go(a,b));

    disNcChiSq e(2,3*3+4*4);
    printf("%u %g %g \n", s.p_dof(), s.p_distance(), s.pdf(2));
    printf("%u %g %g \n", e.p_dof(), e.p_distance(), e.pdf(2));


    return 1;
}