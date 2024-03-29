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
#ifndef STATANALY_D_COMPARE_H_
#define STATANALY_D_COMPARE_H_

#include "fl_comparison.h"
#include "density/probDistr.h"

namespace statanaly {    


/**
 * @brief Check equality by comparing two distributions' hash.
 * 
 * Compare two probabilities' hash values and check for equality.
 * The hash is computed by hashing distribution parameters and the distribution-type. 
 * This is the bit-wise equality.
 * 
 * @param a Probability distribution A.
 * @param b Probability distribution B.
 * @return true 
 * @return false 
 */
inline bool isEqual(const probDistr& a, const probDistr& b) {
    return a.hash() == b.hash();
}



/**
 * @brief Check equality by comparing two distributions' parameters upto a numerical tolerance.
 * 
 * Compare each of the corresponding pair of parameters in two distribution.
 * If any pair of parameters do not match within a numerical difference, then the distributions are not equal.
 * The tolerance is used to ensure the equality two floating-point numbers of equal apparent value, 
 * even though their bit-wise representation are different.
 * 
 * @param a Probability distribution A.
 * @param b Probability distribution B.
 * @param tol A numerical tolerance in double precision.
 * @return true 
 * @return false 
 */
inline bool isEqual_tol(const probDistr& a, const probDistr& b, const double tol) {
    // a & b must have the same concrete type.
    return a.isEqual_tol(b,tol);
}



/**
 * @brief Check equality by comparing two distributions' parameters upto a fixed ULP.
 * 
 * Compare each of the corresponding pair of parameters in two distribution.
 * If any pair of parameters do not match within a fixed ULP, then the distributions are not equal.
 * The tolerance is used to ensure the equality two floating-point numbers of equal apparent value, 
 * even though their bit-wise representation are different.
 *
 * @param a Probability distribution A.
 * @param b Probability distribution B.
 * @param ulp Number of ULPs.
 * @return true 
 * @return false 
 */
inline bool isEqual_ulp(const probDistr& a, const probDistr& b, const unsigned ulp) {
    return a.isEqual_ulp(b,ulp);
}


}   // namespace statanaly

#endif