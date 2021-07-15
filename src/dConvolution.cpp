#include <iostream>
#include "../include/dConvolution.h"

namespace statanaly {

FnDispatcher<probDensFunc,probDensFunc,probDensFunc*> cnvl;


/* Register combinations for dispatches here. ----------- */
auto ConvolutionDoubleDispatcherInitialization = [](){
    cnvl.add<disStdUniform,disStdUniform,convolve>();
    cnvl.add<disNormal,disNormal,convolve>();
    cnvl.add<disCauchy,disCauchy,convolve>();
    cnvl.add<disGamma,disGamma,convolve>();
    return true;
}();



/* Callback functions for double dispatcher for Convolution. 
 * Because the argument types are concrete types, 
 *
 * Note: one can call those functions with concrete types directly 
 * (aka without callbacks). ---------------------------- */

// Convolving Standard Uniform Distribution
probDensFunc* convolve(disStdUniform& lhs, disStdUniform& rhs) {
    disIrwinHall* res = new disIrwinHall(2);
    return static_cast<probDensFunc*>(res);
};

// Convolving Normal Distribution
probDensFunc* convolve(disNormal& lhs, disNormal& rhs) {
    disNormal* res = new disNormal(lhs.mean()+rhs.mean(), lhs.variance()+rhs.variance());
    return static_cast<probDensFunc*>(res);
};

// Convolving Cauchy Distribution
probDensFunc* convolve(disCauchy& lhs, disCauchy& rhs) {
    disCauchy* res = new disCauchy(lhs.ploc()+rhs.ploc(), lhs.pscale()+rhs.pscale());
    return static_cast<probDensFunc*>(res);
};

// Convolving Gamma Distribution
probDensFunc* convolve(disGamma& lhs, disGamma& rhs) {
    // The scale parameters must be identical.
    if (lhs.pscale() != rhs.pscale()) 
        throw std::runtime_error("Convolving(Gamma,Gamma) requires Gamma distributions' scale parameters to be identical.");

    disGamma* res = new disGamma(lhs.pscale(), lhs.pshape()+rhs.pshape());
    return static_cast<probDensFunc*>(res);
};




/* Sum of more than 2 Independent Random Variables ------- */

// Convolving Standard Uniform Distribution
template<>
probDensFunc* convolve<disStdUniform> (std::initializer_list<disStdUniform> l) {
    disIrwinHall* res = new disIrwinHall(l.size());
    return static_cast<probDensFunc*>(res);
};

// Convolving Normal Distribution
template<>
probDensFunc* convolve<disNormal> (std::initializer_list<disNormal> l) {
    double m = 0, v = 0;
    for (auto& e : l) {
        m += e.mean();
        v += e.variance();
    }
    disNormal* res = new disNormal(m, v);

    return static_cast<probDensFunc*>(res);
};

// Convolving Cauchy Distribution
template<>
probDensFunc* convolve<disCauchy> (std::initializer_list<disCauchy> l) {
    double m = 0, v = 0;
    for (auto& e : l) {
        m += e.ploc();
        v += e.pscale();
    }
    disCauchy* res = new disCauchy(m, v);

    return static_cast<probDensFunc*>(res);
};

// Convolving Gamma Distribution
template<>
probDensFunc* convolve<disGamma> (std::initializer_list<disGamma> l) {
    double m = 0, n = l.begin()->pscale();
    for (auto& e : l) {
        m += e.pshape();
        if (n != e.pscale()) 
            throw std::runtime_error("Convolving(Gamma,Gamma) requires Gamma distributions' scale parameters to be identical.");
    }
    disGamma* res = new disGamma(n, m);

    return static_cast<probDensFunc*>(res);
};


}