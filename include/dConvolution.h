#ifndef STATANALY_D_CONVOLUTION_H_
#define STATANALY_D_CONVOLUTION_H_

#include "double_dispatcher.h"
#include "density/disUniform.h"
#include "density/disNormal.h"
#include "density/disIrwinHall.h"
#include "density/disCauchy.h"
#include "density/disGamma.h"
#include "density/disExponential.h"
#include "density/disErlang.h"
#include "density/disRayleigh.h"
#include "density/disRician.h"
#include "density/disChi.h"
#include "density/disChiSq.h"
#include "density/disNcChi.h"
#include "density/disNcChiSq.h"


/* Sum of two or more Independent Random Variables is the 
    Convolution of their individual distributions --------------------- */


namespace statanaly {

// Create a Double Dispatcher for Convolution.
extern FnDispatcher<probDistr,probDistr,probDistr*> cnvl;
extern FnDispatcher<probDistr,probDistr,probDistr*> cnvlSq;
extern FnDispatcher<probDistr,probDistr,probDistr*> cnvlSSqrt;

// Callback functions for double dispatcher for Convolution.
// Because the argument types are concrete types, one can call those functions with concrete types directly (aka without callbacks).

// Sum of two Standard Uniform RVs
// R = X + Y
probDistr* convolve(disStdUniform& lhs, disStdUniform& rhs);

// Sum of two Normal RVs
// R = X + Y
probDistr* convolve(disNormal& lhs, disNormal& rhs);

// Sum of two Cauchy RVs
// R = X + Y
probDistr* convolve(disCauchy& lhs, disCauchy& rhs);

// Sum of two Gamma RVs
// R = X + Y
probDistr* convolve(disGamma& lhs, disGamma& rhs);

// Sum of two Erlang RVs
// R = X + Y
probDistr* convolve(disExponential& lhs, disExponential& rhs);

// Sum of the square of two Normal RVs.
// The Normal RVs have zero mean, ie, N(0,sig^2).
// R = X^2 + Y^2
probDistr* convolveSq(disNormal& lhs, disNormal& rhs);

// Sum of the square of two Normal RVs; then take the sqrt of the sum.
// R = sqrt(X^2 + Y^2)
probDistr* convolveSSqrt(disNormal& lhs, disNormal& rhs);



// Disable un-implemented distributions.
// Note: initializer-list uses copy-semantics.
template<typename T>
probDistr* convolve(std::initializer_list<T> l) = delete;

// Sum of Standard Uniform RVs
template<>
probDistr* convolve<disStdUniform> (std::initializer_list<disStdUniform> l);

// Sum of Normal RVs
template<>
probDistr* convolve<disNormal> (std::initializer_list<disNormal> l);

// Sum of Cauchy RVs
template<>
probDistr* convolve<disCauchy> (std::initializer_list<disCauchy> l);

// Sum of Gamma RVs
template<>
probDistr* convolve<disGamma> (std::initializer_list<disGamma> l);

// Sum of Exponential RVs
template<>
probDistr* convolve<disExponential> (std::initializer_list<disExponential> l);



template<typename T>
probDistr* convolveSq(std::initializer_list<T> l) = delete;

// Sum of the square of Normal RVs.
// R = X**2 + Y**2 + Z**2 + ...
template<>
probDistr* convolveSq<disNormal> (std::initializer_list<disNormal> l);



template<typename T>
probDistr* convolveSSqrt(std::initializer_list<T> l) = delete;

// Sum of the square of Normal RVs; then take the sqrt of the sum.
// R = sqrt(X**2 + Y**2 + Z**2 + ...)
template<>
probDistr* convolveSSqrt<disNormal> (std::initializer_list<disNormal> l);


}   // namespace statanaly


#endif