#include <iostream>
#include "../include/dConvolution.h"

namespace statanaly {

FnDispatcher<probDensFunc,probDensFunc,probDensFunc*> cnvl;


/* Register combinations for dispatches here. ----------- */
auto ConvolutionDoubleDispatcherInitialization = [](){
    cnvl.add<disStdUniform,disStdUniform,convolve>();
    cnvl.add<disNormal,disNormal,convolve>();
    return true;
}();



/* Callback functions for double dispatcher for Convolution. 
   Because the argument types are concrete types, 
   one can call those functions with concrete types directly 
   (aka without callbacks). ---------------------------- */

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




/* Sum of >2 Independent Random Variables ---------------------------- */

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
}