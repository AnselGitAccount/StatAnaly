#include <iostream>
#include "../include/dConvolution.h"

namespace statanaly {

FnDispatcher<probDensFunc,probDensFunc,probDensFunc*> cnvl;
FnDispatcher<probDensFunc,probDensFunc,probDensFunc*> cnvlSSqrt;


/* Register combinations for dispatches here. ----------- */
auto ConvolutionDoubleDispatcherInitialization = [](){
    cnvl.add<disStdUniform,disStdUniform,convolve>();
    cnvl.add<disNormal,disNormal,convolve>();
    cnvl.add<disCauchy,disCauchy,convolve>();
    cnvl.add<disGamma,disGamma,convolve>();
    cnvl.add<disExponential,disExponential,convolve>();
    return true;
}();

auto ConvolutionSSqrtDoubleDispatcherInitialization = [](){
    cnvlSSqrt.add<disNormal,disNormal,convolveSSqrt>();
    return true;
}();


/* Callback functions for double dispatcher for Convolution. 
 * Because the argument types are concrete types, they can be called directly.
 *
 * Note: one can call those functions with concrete types directly 
 * (aka without callbacks). ---------------------------- */

probDensFunc* convolve(disStdUniform& lhs, disStdUniform& rhs) {
    probDensFunc* res = new disIrwinHall(2);
    return res;
};

probDensFunc* convolve(disNormal& lhs, disNormal& rhs) {
    probDensFunc* res = new disNormal(lhs.mean()+rhs.mean(), lhs.variance()+rhs.variance());
    return res;
};

probDensFunc* convolve(disCauchy& lhs, disCauchy& rhs) {
    probDensFunc* res = new disCauchy(lhs.ploc()+rhs.ploc(), lhs.pscale()+rhs.pscale());
    return res;
};

probDensFunc* convolve(disGamma& lhs, disGamma& rhs) {
    // The scale parameters must be identical.
    if (lhs.pscale() != rhs.pscale())
        throw std::invalid_argument("convolve(Gamma,Gamma) requires Gamma distributions' scale parameters to be identical.");

    probDensFunc* res = new disGamma(lhs.pscale(), lhs.pshape()+rhs.pshape());
    return res;
};

probDensFunc* convolve(disExponential& lhs, disExponential& rhs) {
    // The scale parameters must be identical.
    if (lhs.prate() != rhs.prate())
        throw std::invalid_argument("convolve(Exponential,Exponential) requires Exponential distributions' rate parameters to be identical.");

    probDensFunc* res = new disErlang(2, lhs.prate());
    return res;
};

probDensFunc* convolveSSqrt(disNormal& lhs, disNormal& rhs) {
    // lhs's and rhs's means must be identical.
    // lhs's and rhs's variances must be identical.
    if (lhs.stddev() != rhs.stddev())
        throw std::invalid_argument("convolveSSqrt(Normal,Normal) requires Normal distributions' scale parameters to be identical.");
    if (lhs.mean() != rhs.mean())
        throw std::invalid_argument("convolveSSqrt(Normal,Normal) requires Normal distributions' location parameters to be identical.");

    // If lhs's and rhs's means are zero, return Rayleigh distribution.
    // Else, return Rician distribution.
    probDensFunc* res = nullptr;
    if (lhs.mean()==0)  // && rhs.mean()==0 
        res = new disRayleigh(lhs.stddev());
    else {
        auto a = std::sqrt(lhs.mean()*lhs.mean() + rhs.mean()*rhs.mean());
        res = new disRician(a, lhs.stddev());
    }

    return res;
}


/* Sum of more than 2 Independent Random Variables ------- */

template<>
probDensFunc* convolve<disStdUniform> (std::initializer_list<disStdUniform> l) {
    probDensFunc* res = new disIrwinHall(l.size());
    return res;
};

template<>
probDensFunc* convolve<disNormal> (std::initializer_list<disNormal> l) {
    double m = 0, v = 0;
    for (auto& e : l) {
        m += e.mean();
        v += e.variance();
    }
    probDensFunc* res = new disNormal(m, v);

    return res;
};

template<>
probDensFunc* convolve<disCauchy> (std::initializer_list<disCauchy> l) {
    double m = 0, v = 0;
    for (auto& e : l) {
        m += e.ploc();
        v += e.pscale();
    }
    probDensFunc* res = new disCauchy(m, v);

    return res;
};

template<>
probDensFunc* convolve<disGamma> (std::initializer_list<disGamma> l) {
    double m = 0, n = l.begin()->pscale();
    for (auto& e : l) {
        m += e.pshape();
        if (n != e.pscale()) 
            throw std::invalid_argument("Convolve({Gamma_i}) requires Gamma distributions' scale parameters to be identical.");
    }
    probDensFunc* res = new disGamma(n, m);

    return res;
};

template<>
probDensFunc* convolve<disExponential> (std::initializer_list<disExponential> l) {
    double n = l.begin()->prate();
    for (auto& e : l) {
        if (n != e.prate()) 
            throw std::invalid_argument("Convolve({Exponential_i}) requires Exponential distributions' rate parameters to be identical.");
    }
    probDensFunc* res = new disErlang(l.size(), n);

    return res;
};

template<>
probDensFunc* convolveSSqrt<disNormal> (std::initializer_list<disNormal> l) {
    // If the Normal distributions' means are zero, then it is Central Chi.
    // Else, then it is Non-central Chi.
    double mu  = l.begin()->p_location();
    double sig = l.begin()->p_scale();
    long double a = 0;
    for (auto& e : l) {
        if (mu  != e.p_location()) 
            throw std::invalid_argument("ConvolveSSqrt({Normal_i}) requires Normal distributions' location parameters (ie, mean) to be zero.");
        if (1 != e.p_scale())
            throw std::invalid_argument("ConvolveSSqrt({Normal_i}) requires Normal distributions' scale parameters to be One");
        if (sig != e.p_scale())
            throw std::invalid_argument("ConvolveSSqrt({Normal_i}) requires Normal distributions' scale parameters (ie, std deviation) to be identical.");
        a += e.p_location()*e.p_location();
    }

    probDensFunc* res = nullptr;
    if (mu==0)
        res = new disChi(l.size());
    else 
        res = new disNcChi(l.size(), sqrt(a));

    return res;
};

template<>
probDensFunc* convolveSq<disNormal> (std::initializer_list<disNormal> l) {
    // If the Normal distributions' means are zero, then it is Central Chi Square.
    // Else, then it is Non-central Chi Square.
    double mu  = l.begin()->p_location();
    double sig = l.begin()->p_scale();
    long double a = 0;
    for (auto& e : l) {
        if (mu  != e.p_location()) 
            throw std::invalid_argument("ConvolveSq({Normal_i}) requires Normal distributions' location parameters (ie, mean) to be zero.");
        if (1 != e.p_scale())
            throw std::invalid_argument("ConvolveSq({Normal_i}) requires Normal distributions' scale parameters to be One");
        if (sig != e.p_scale())
            throw std::invalid_argument("ConvolveSq({Normal_i}) requires Normal distributions' scale parameters (ie, std deviation) to be identical.");
        a += e.p_location()*e.p_location();
    }

    probDensFunc* res = nullptr;
    if (mu==0)
        res = new disChiSq(l.size());
    else 
        res = new disNcChiSq(l.size(), a);

    return res;    
}

}   // namespace statanaly