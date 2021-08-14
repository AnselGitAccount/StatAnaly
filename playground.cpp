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

    disNcChiSq d1(2,3);
    

    std::cout << d1.pdf(2) << std::endl;
    std::cout << d1.cdf(2) << std::endl;
    
    
    // printf("%g \n",regLowerGamma(-2.4, 1.125));

    return 1;
}