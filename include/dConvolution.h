/*
   Copyright 2022, Ansel Blumers

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
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


/**
 * @file dConvolution.h
 * @brief Convolution of probability distributions.
 * 
 * This file contains:
 * 
 * 1. Callback functions for double dispatcher when convoluting two distributions.
 * 2. Functions for convoluting a collection of distribution of the same type. 
 * 
 * Because the argument types are concrete types, one can call 
 * those functions with concrete types directly.
 *
 * Remember: Sum of two or more Independent Random Variables is 
 * the Convolution of their individual distributions. 
 */

namespace statanaly {

// Create a Double Dispatcher for Convolution.
extern FnDispatcher<probDistr,probDistr,probDistr*> cnvl;
extern FnDispatcher<probDistr,probDistr,probDistr*> cnvlSq;
extern FnDispatcher<probDistr,probDistr,probDistr*> cnvlSSqrt;


/**
 * @brief Sum of two Standard Uniform RVs.
 * 
 * R = X + Y
 */
probDistr* convolve(disStdUniform& lhs, disStdUniform& rhs);

/**
 * @brief Sum of two Normal RVs.
 * 
 * R = X + Y
 */
probDistr* convolve(disNormal& lhs, disNormal& rhs);

/**
 * @brief Sum of two Cauchy RVs.
 * 
 * R = X + Y
 */
probDistr* convolve(disCauchy& lhs, disCauchy& rhs);

/**
 * @brief Sum of two Gamma RVs.
 * 
 * R = X + Y
 * Requires Gamma distributions' scale parameters to be identical.
 */
probDistr* convolve(disGamma& lhs, disGamma& rhs);

/**
 * @brief Sum of two Exponential RVs.
 * 
 * R = X + Y
 * Requires Exponential distributions' rate parameters to be identical.
 */
probDistr* convolve(disExponential& lhs, disExponential& rhs);

/**
 * @brief Sum of the square of two Normal RVs.
 * 
 * R = X^2 + Y^2
 * The Normal RVs have zero mean, ie, N(0,sig^2).
 */
probDistr* convolveSq(disNormal& lhs, disNormal& rhs);

/**
 * @brief Sum of the square of two Normal RVs, then take the sqrt of the sum.
 * 
 * R = sqrt(X^2 + Y^2)
 * Requires Normal distributions' scale parameters to be identical.
 */
probDistr* convolveSSqrt(disNormal& lhs, disNormal& rhs);



// Disable un-implemented distributions.
// Note: initializer-list uses copy-semantics.
template<typename T>
probDistr* convolve(std::initializer_list<T> l) = delete;

/**
 * @brief Sum of Standard Uniform RVs.
 * 
 * R = X + Y + Z + ...
 */
template<>
probDistr* convolve<disStdUniform> (std::initializer_list<disStdUniform> l);

/**
 * @brief Sum of Normal RVs.
 * 
 * R = X + Y + Z + ...
 */
template<>
probDistr* convolve<disNormal> (std::initializer_list<disNormal> l);

/**
 * @brief Sum of Cauchy RVs.
 * 
 * R = X + Y + Z + ...
 */
template<>
probDistr* convolve<disCauchy> (std::initializer_list<disCauchy> l);

/**
 * @brief Sum of Gamma RVs.
 * 
 * R = X + Y + Z + ...
 * Requires Gamma distributions' scale parameters to be identical.
 */
template<>
probDistr* convolve<disGamma> (std::initializer_list<disGamma> l);


/**
 * @brief Sum of Exponential RVs.
 * 
 * R = X + Y + Z + ...
 * Requires Exponential distributions' rate parameters to be identical.
 */
template<>
probDistr* convolve<disExponential> (std::initializer_list<disExponential> l);




// Disable un-implemented distributions.
template<typename T>
probDistr* convolveSq(std::initializer_list<T> l) = delete;


/**
 * @brief Sum of the square of Normal RVs.
 * 
 * R = X^2 + Y^2 + Z^2 + ...
 * 
 * Requires Normal distributions' scale parameters to be One.
 * Requires Normal distributions' scale parameters (ie, std deviation) to be identical.
 */
template<>
probDistr* convolveSq<disNormal> (std::initializer_list<disNormal> l);



// Disable un-implemented distributions.
template<typename T>
probDistr* convolveSSqrt(std::initializer_list<T> l) = delete;

/**
 * @brief Sum of the square of Normal RVs, then take the sqrt of the sum.
 * 
 * R = sqrt(X^2 + Y^2 + Z^2 + ...)
 * 
 * Requires Normal distributions' scale parameters to be One.
 * Requires Normal distributions' scale parameters (ie, std deviation) to be identical.
 */
template<>
probDistr* convolveSSqrt<disNormal> (std::initializer_list<disNormal> l);


}   // namespace statanaly


#endif