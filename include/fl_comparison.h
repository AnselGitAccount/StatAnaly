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
#ifndef STATANALY_FL_COMPARISON_H_
#define STATANALY_FL_COMPARISON_H_

#include <cmath>
#include <stdint.h>

/**
 * @file fl_comparison.h
 * @brief Comparing floating point numbers.
 * 
 * A derivative of gtest. See gtest for detail.
 */

namespace statanaly {

/**
 * @brief The unsigned integer type that has the same size as the floating point number.
 * 
 * Use template specialization to prevent the user from using TypeWithSize<N> with incorrect size N.
 *
 * @tparam size Size of the representing floating point number.
 */
template <std::size_t size>
class TypeWithSize {
 public:
  typedef void UInt;
};

// The specialization for size 4.
template <>
class TypeWithSize<4> {
 public:
  // unsigned int has size 4 in both gcc and MSVC.
  //
  // As base/basictypes.h doesn't compile on Windows, we cannot use
  // uint32, uint64, and etc here.
  typedef int32_t Int;
  typedef uint32_t UInt;
};

// The specialization for size 8.
template <>
class TypeWithSize<8> {
 public:
  typedef int64_t Int;
  typedef uint64_t UInt;
};


// Taken from gtest:
// 
// Converts an integer from the sign-and-magnitude representation to
// the biased representation.  More precisely, let N be 2 to the
// power of (kBitCount - 1), an integer x is represented by the
// unsigned number x + N.
//
// For instance,
//
//   -N + 1 (the most negative number representable using
//          sign-and-magnitude) is represented by 1;
//   0      is represented by N; and
//   N - 1  (the biggest number representable using
//          sign-and-magnitude) is represented by 2N - 1.
//
// Read http://en.wikipedia.org/wiki/Signed_number_representations
// for more details on signed number representations.
template<typename RawType, typename Bits>
static Bits SignAndMagnitudeToBiased(const Bits &sam) {
    // # of bits in a number.
    static const size_t kBitCount = 8*sizeof(RawType);

    // The mask for the sign bit.
    static const Bits kSignBitMask = static_cast<Bits>(1) << (kBitCount - 1);

    if (kSignBitMask & sam) {
        // sam represents a negative number.
        return ~sam + 1;
    } else {
        // sam represents a positive number.
        return kSignBitMask | sam;
    }
};


// Given two numbers in the sign-and-magnitude representation,
// returns the distance between them as an unsigned number.
template<typename RawType, typename Bits>
static Bits DistanceBetweenSignAndMagnitudeNumbers(const Bits &sam1,
                                                   const Bits &sam2) {
    const Bits biased1 = SignAndMagnitudeToBiased<RawType>(sam1);
    const Bits biased2 = SignAndMagnitudeToBiased<RawType>(sam2);
    return (biased1 >= biased2) ? (biased1 - biased2) : (biased2 - biased1);
};



// The data type used to store the actual floating-point number.
template<typename RawType, typename Bits>
union FloatingPointUnion {
    RawType value_;  // The raw floating-point number.
    Bits bits_;      // The bits that represent the number.
};



// Returns true iff lhs is at most kMaxUlps ULP's away from
// rhs.  In particular, this function:
//
//   - returns false if either number is (or both are) NAN.
//   - treats really large numbers as almost equal to infinity.
//   - thinks +0.0 and -0.0 are 0 DLP's apart.
template<typename RawType, typename Bits>
bool AlmostEquals(const FloatingPointUnion<RawType,Bits>& lhs, 
                  const FloatingPointUnion<RawType,Bits>& rhs,
                  const uint32_t kMaxUlps) {
    // The IEEE standard says that any comparison operation involving
    // a NAN must return false.
    if (std::isnan(lhs.value_) || std::isnan(rhs.value_)) 
        return false;

    return DistanceBetweenSignAndMagnitudeNumbers<RawType,Bits>(lhs.bits_, rhs.bits_)
        <= kMaxUlps;
}


/**
 * @brief Floating-point comparator with ULP threshold.
 * 
 * Two floating-point numbers are equal if their difference is within "kMaxUlps" ULPs.
 * 
 * Default is kMaxUlps=4 ULPs. 
 * On Intel CPU all fp calculations are done with 80-bit precision.
 * Double has 64-bits. So 4 ULP should be enough for ordinary use.
 * 
 * @tparam RawType The raw floating-point type (either float or double).
 * @param a 
 * @param b 
 * @param kMaxUlps Number of ULPs as the threshold.
 */
template<typename RawType>
requires std::is_floating_point_v<RawType>
bool isEqual_fl_ulp(const RawType a, const RawType b, const uint32_t kMaxUlps) {

    // Template parameter:
    //
    //   RawType: the raw floating-point type (either float or double)
    //   Bits: the unsigned integer type that has the same size as the floating point number.

    typedef typename TypeWithSize<sizeof(RawType)>::UInt Bits;
    FloatingPointUnion<RawType,Bits> lhs{a}, rhs{b};

    return AlmostEquals<RawType,Bits>(lhs, rhs, kMaxUlps);
}


/**
 * @brief Floating-point comparator with a floating-point threshold.
 * 
 * Two floating-point numbers are equal if their difference is the floating-point thershold.
 * 
 * @tparam RawType The raw floating-point type (either float or double).
 * @tparam TolType The raw floating-point type (either float or double) for the threshold.
 * @param a 
 * @param b 
 * @param tol The floating-point theshold.
 */
template<class RawType, class TolType>
requires std::is_arithmetic_v<TolType> && std::is_floating_point_v<RawType>
bool isEqual_fl_tol(const RawType a, const RawType b, const TolType tol) {
    return ((a>=b) ? (a-b):(b-a)) <= tol;
}

}

#endif