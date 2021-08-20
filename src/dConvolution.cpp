#include <iostream>
#include "dConvolution.h"

namespace statanaly {

FnDispatcher<probDensFunc,probDensFunc,probDensFunc*> cnvl;
FnDispatcher<probDensFunc,probDensFunc,probDensFunc*> cnvlSq;
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

auto ConvolutionSqDoubleDispatcherInitialization = [](){
    cnvlSq.add<disNormal,disNormal,convolveSq>();
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

probDensFunc* convolve(disStdUniform& l, disStdUniform& r) {
    probDensFunc* res = new disIrwinHall(2);
    return res;
};

probDensFunc* convolve(disNormal& l, disNormal& r) {
    probDensFunc* res = new disNormal(l.mean()+r.mean(), l.variance()+r.variance());
    return res;
};

probDensFunc* convolve(disCauchy& l, disCauchy& r) {
    probDensFunc* res = new disCauchy(l.ploc()+r.ploc(), l.pscale()+r.pscale());
    return res;
};

probDensFunc* convolve(disGamma& l, disGamma& r) {
    // The scale parameters must be identical.
    if (l.pscale() != r.pscale())
        throw std::invalid_argument("convolve(Gamma,Gamma) requires Gamma distributions' scale parameters to be identical.");

    probDensFunc* res = new disGamma(l.pscale(), l.pshape()+r.pshape());
    return res;
};

probDensFunc* convolve(disExponential& l, disExponential& r) {
    // The scale parameters must be identical.
    if (l.prate() != r.prate())
        throw std::invalid_argument("convolve(Exponential,Exponential) requires Exponential distributions' rate parameters to be identical.");

    probDensFunc* res = new disErlang(2, l.prate());
    return res;
};

probDensFunc* convolveSSqrt(disNormal& l, disNormal& r) {
    // l's and r's variances must be identical.
    if (l.stddev() != r.stddev())
        throw std::invalid_argument("convolveSSqrt(Normal,Normal) requires Normal distributions' scale parameters to be identical.");

    // If l's and r's means are zero, return Rayleigh distribution.
    // Else, return Rician distribution.
    probDensFunc* res = nullptr;
    if (l.mean()==0 && r.mean()==0)
        res = new disRayleigh(l.stddev());
    else {
        auto a = std::sqrt(l.mean()*l.mean() + r.mean()*r.mean());
        res = new disRician(a, l.stddev());
    }

    return res;
};

probDensFunc* convolveSq(disNormal& l, disNormal& r) {
    // l's and r's variances must be equal one.
    if (l.stddev() != 1 || r.stddev() != 1)
        throw std::invalid_argument("convolveSq(Normal,Normal) requires Normal distributions' scale parameters to be one.");

    // If l's and r's means are zero, return central chi distribution.
    // Else, return noncentral chi distribution.
    probDensFunc* res = nullptr;
    if (l.mean()==0 && r.mean()==0)
        res = new disChiSq(2);
    else {
        auto a = l.mean()*l.mean() + r.mean()*r.mean();
        res = new disNcChiSq(2,a);
    }

    return res;
};


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
    double mu  = l.begin()->p_location();
    double sig = l.begin()->p_scale();
    long double a = 0;
    for (auto& e : l) {
        if (1 != e.p_scale())
            throw std::invalid_argument("ConvolveSSqrt({Normal_i}) requires Normal distributions' scale parameters to be One");
        if (sig != e.p_scale())
            throw std::invalid_argument("ConvolveSSqrt({Normal_i}) requires Normal distributions' scale parameters (ie, std deviation) to be identical.");
        a += e.p_location()*e.p_location();
    }

    // If the Normal distributions' means are zero, then it is Central Chi.
    // Else, then it is Non-central Chi.
    probDensFunc* res = nullptr;
    if (mu==0)
        res = new disChi(l.size());
    else 
        res = new disNcChi(l.size(), sqrt(a));

    return res;
};

template<>
probDensFunc* convolveSq<disNormal> (std::initializer_list<disNormal> l) {
    double mu  = l.begin()->p_location();
    double sig = l.begin()->p_scale();
    long double a = 0;
    for (auto& e : l) {
        if (1 != e.p_scale())
            throw std::invalid_argument("ConvolveSq({Normal_i}) requires Normal distributions' scale parameters to be One");
        if (sig != e.p_scale())
            throw std::invalid_argument("ConvolveSq({Normal_i}) requires Normal distributions' scale parameters (ie, std deviation) to be identical.");
        a += e.p_location()*e.p_location();
    }

    // If the Normal distributions' means are zero, then it is Central Chi Square.
    // Else, then it is Non-central Chi Square.
    probDensFunc* res = nullptr;
    if (mu==0)
        res = new disChiSq(l.size());
    else 
        res = new disNcChiSq(l.size(), a);

    return res;    
}

}   // namespace statanaly