#ifndef STATANALY_D_CONVOLUTION_H_
#define STATANALY_D_CONVOLUTION_H_

#include "double_dispatcher.h"
#include "density/disUniform.h"
#include "density/disNormal.h"
#include "density/disIrwinHall.h"
#include "density/disCauchy.h"
#include "density/disGamma.h"
#include "density/disExponential.h"

namespace statanaly {

// Create a Double Dispatcher for Convolution.
extern FnDispatcher<probDensFunc,probDensFunc,probDensFunc*> cnvl;


/* Sum of two or more Independent Random Variables is the 
    Convolution of their individual distributions --------------------- */


// Callback functions for double dispatcher for Convolution.
// Because the argument types are concrete types, one can call those functions with concrete types directly (aka without callbacks).

// Convolving Standard Uniform Distribution
probDensFunc* convolve(disStdUniform& lhs, disStdUniform& rhs);

// Convolving Normal Distribution
probDensFunc* convolve(disNormal& lhs, disNormal& rhs);

// Convolving Cauchy Distribution
probDensFunc* convolve(disCauchy& lhs, disCauchy& rhs);

// Convolving Gamma Distribution
probDensFunc* convolve(disGamma& lhs, disGamma& rhs);



/* Sum of >2 Independent Random Variables ---------------------------- */

// Note: initializer-list uses copy-semantics.
template<typename T>
probDensFunc* convolve(std::initializer_list<T> l) = delete;

// Convolving Standard Uniform Distribution
template<>
probDensFunc* convolve<disStdUniform> (std::initializer_list<disStdUniform> l);

// Convolving Normal Distribution
template<>
probDensFunc* convolve<disNormal> (std::initializer_list<disNormal> l);

// Convolving Cauchy Distribution
template<>
probDensFunc* convolve<disCauchy> (std::initializer_list<disCauchy> l);

// Convolving Gamma Distribution
template<>
probDensFunc* convolve<disGamma> (std::initializer_list<disGamma> l);


}   // namespace statanaly


#endif