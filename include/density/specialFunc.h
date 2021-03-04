#ifndef STATANALY_SPECIAL_FUNC_H_
#define STATANALY_SPECIAL_FUNC_H_

#include <type_traits>
#include <array>
#include <cstdint>
#include <cmath>
#include "../const_values.h"

namespace statanaly {

/* Factorial */
template<typename T> 
requires std::is_integral<T>::value
constexpr inline uint64_t genFactorial(T x) {
    uint64_t prod = 1;
    while (x>1) {prod *= x--;}
    return prod;
};

// Generate a Lookup Table at compile-time.
constexpr auto factorial = []{
    // uint64 can hold factorial(18).
    constexpr auto LUT_size = 19;
    std::array<uint64_t, LUT_size> arr{};
    for(int i=0; i<LUT_size; ++i) {
        arr[i] = genFactorial(i);
    }

    return arr;
}(); // () at end instantiates this lambda function.

static_assert(factorial[18] == 6402373705728000);



/* Gamma function, Integer */
/* When argument is an integer, 
   gamma function is equivalent to factorial:
   Gamma(n) = (n-1)!    */
constexpr auto gammaIntFunc = []{
    constexpr auto LUT_size = 20;
    std::array<uint64_t, LUT_size> arr{};

    arr[0] = 0;
    for(int i=1; i<LUT_size; ++i) {
        arr[i] = genFactorial(i-1);
    }

    return arr;
}(); // () at end instantiates this lambda function.
    
static_assert(gammaIntFunc[19] == 6402373705728000);


/* Gamma function, General (Integer + FloatingPoint) */
// http://www.numericana.com/answer/info/godfrey.htm

constexpr long double logGammaCoeffs(const long double x) {
    long double res = 
        0.99999999999999709182L               + 57.156235665862923517L    / (x+1)  \
        - 59.597960355475491248L    / (x+2)   + 14.136097974741747174L    / (x+3)  \
        - 0.49191381609762019978L   / (x+4)   + .33994649984811888699e-4L / (x+5)  \
        + .46523628927048575665e-4L / (x+6)   - .98374475304879564677e-4L / (x+7)  \
        + .15808870322491248884e-3L / (x+8)   - .21026444172410488319e-3L / (x+9)  \
        + .21743961811521264320e-3L / (x+10)  - .16431810653676389022e-3L / (x+11) \
        + .84418223983852743293e-4L / (x+12)  - .26190838401581408670e-4L / (x+13) \
        + .36899182659531622704e-5L / (x+14);
    return res;
}

template<typename T>
constexpr T logGamma(T x) {   
    x -= 1;
    T a = 607.L/128 + 0.5L;
    T t1 = (x + T(0.5))*std::log(x + T(a)) - (x + T(a));
    T t2 = T(SACV_LOG_SQRT_2PI) + T(std::log(logGammaCoeffs(x)));

    return t1 + t2;
}

static_assert(uint64_t(std::exp(logGamma(19.L))) < 6402373705728000);
static_assert(uint64_t(std::exp(logGamma(19.L))) > 6402373705727980);



} // stantanaly

#endif