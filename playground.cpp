#include <iostream>
#include "gtest/gtest.h"

#include "include/density/disMixture.h"
#include "include/density/disNormal.h"
#include "include/density/disUniform.h"
#include "include/density/disIrwinHall.h"
#include "include/density/disChiSq.h"
#include "include/density/disCauchy.h"
#include "include/density/disGamma.h"
#include "include/density/disExponential.h"
#include "include/dContainer.h"
#include "include/fl_comparison.h"
#include "include/dCompare.h"
#include "include/dConvolution.h"

using namespace statanaly;



int main(int argc, char** argv) {

    // auto d1 = disExponential(1.0);
    // std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << d1.pdf(1) << std::endl;
    // std::cout << 1/exp(1) << std::endl;
    // std::cout << d1.cdf(1) << std::endl;
    // std::cout << 1-1/exp(1) << std::endl;

    std::cout.precision(std::numeric_limits<double>::max_digits10);
    // std::cout << regLowerGamma(1.0,2.0) << std::endl;
    // std::cout << regUpperGamma(1.0,2.0) << std::endl;

    disGamma n1{1,2}, n2{1,1}, n3{1., 2.3}, n4{1., 0.4};
    probDensFunc* rn = convolve<disGamma>({n1, n2, n3, n4});
    disGamma* rc = dynamic_cast<disGamma*>(rn);
    disGamma expe{1., 5.7};

    if ( rc->pscale() == expe.pscale() ) {
        std::cout << "scale\n";
    }
    if ( rc->pshape() == expe.pshape() ) {
        std::cout << "shape\n";
    }

    return 1;
}