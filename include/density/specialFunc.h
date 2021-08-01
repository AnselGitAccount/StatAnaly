#ifndef STATANALY_SPECIAL_FUNC_H_
#define STATANALY_SPECIAL_FUNC_H_

#include <type_traits>
#include <array>
#include <cstdint>
#include <cmath>
#include "../const_values.h"
#include <cstdio>
#include <stdexcept>

namespace statanaly {

/* Special Constants ----------------------------------
 */

// sqrt(pi/2)
constexpr double SQRT_PI_2 = M_PI_2 * M_2_SQRTPI * M_SQRT1_2;   


/* Factorial ------------------------------------------ 
 * Generate a Lookup Table at compile-time.
 */
template<typename T> 
requires std::is_integral<T>::value
constexpr inline uint64_t genFactorial(T x) {
    uint64_t prod = 1;
    while (x>1) {prod *= x--;}
    return prod;
};

constexpr auto factorial = []{
    // uint64 can hold upto factorial(18).
    std::array<uint64_t, 19> arr{};
    for(int i=0; i<19; ++i) {
        arr[i] = genFactorial(i);
    }

    return arr;
}(); // () at end instantiates this lambda function.

static_assert(factorial[18] == 6402373705728000);



/* Gamma function ------------------------------------ 
 * Integer parameters */
/* When argument is an integer, 
 * gamma function is equivalent to factorial:
 * Gamma(n) == (n-1)!    
 */
constexpr auto gammaIntFunc = []{
    std::array<uint64_t, 20> arr{};

    arr[0] = 0;
    for(int i=1; i<20; ++i) {
        arr[i] = genFactorial(i-1);
    }

    return arr;
}(); // () at end instantiates this lambda function.
    
static_assert(gammaIntFunc[19] == 6402373705728000);


/* Log Gamma function -----------------------------------
 * Takes in either Integer or FloatingPoint.
 * Reference: Godfrey  
 * http://www.numericana.com/answer/info/godfrey.htm 
 */

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
requires std::is_floating_point_v<T>
constexpr T logGamma(T x) {   
    x -= 1;
    T a = 607.L/128 + 0.5L;
    T t1 = (x + T(0.5))*std::log(x + T(a)) - (x + T(a));
    T t2 = T(SACV_LOG_SQRT_2PI) + T(std::log(logGammaCoeffs(x)));

    return t1 + t2;
}

static_assert(uint64_t(std::exp(logGamma(19.L))) < 6402373705728000);
static_assert(uint64_t(std::exp(logGamma(19.L))) > 6402373705727980);




/* erf function ---------------------------------------- 
 * Use erf defined in libm for type double.
 * Use erf defined in libstdc++ (cmath.h) for other types.
 * https://stackoverflow.com/questions/631629/erfx-and-math-h
 */

static_assert(std::erf(2.0) < 0.996);



/* Regularized Lower Gamma function --------------------
 * Power series expansion.
 * Source:
 *      - Wikipedia
 *      - AS245 (http://lib.stat.cmu.edu/apstat/245)
 *      - https://github.com/lh3/samtools/blob/master/bcftools/kfunc.c
 */

#define STATANALY_GAMMA_EPS 1e-14
#define STATANALY_GAMMA_TINY 1e-290

inline double _regLowerGamma(double s, double z) {
    double sum = 1, x = 1;
    for (int k = 1; k < 100; ++k) {
        x *= z / (s+k);
        sum += x;
		if (x/sum < STATANALY_GAMMA_EPS)
            break;
    }
    return exp(s * log(z) - z - logGamma(s+1.) + log(sum));
}


/* Regularized Upper Gamma function -------------------
 * Infinite fraction.
 * Source:
 *      - Numerical Recipes in C, 2nd Ed, section 5.2
 */

inline double _regUpperGamma(double s, double z) {
    double f = 1. + z - s;
    double C = f, D = 0;
    
    // Modified Lentz's algorithm to simplify the infinite fraction.
    for (int j = 1; j < 100; j++) {
        double a = j * (s - j);
        double b = (j<<1) + 1 + z - s;
        double d;
		D = b + a * D;
		if (D < STATANALY_GAMMA_TINY) 
            D = STATANALY_GAMMA_TINY;
		C = b + a / C;
		if (C < STATANALY_GAMMA_TINY) 
            C = STATANALY_GAMMA_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < STATANALY_GAMMA_EPS) 
            break;
    }

    return exp(s * log(z) - z - logGamma(s) - log(f));
}

// z is the bound of the integral.
// s is the power.
double regLowerGamma(double s, double z);
double regUpperGamma(double s, double z);


/* Lower and Upper Gamma function --------------------
 * Compute regularized Lower and Upper Gamma,
 * and then multiply by Complete Gamma function.
 */

// z is the bound of the integral.
// s is the power.
double upperGamma(double s, double z);
double lowerGamma(double s, double z);


/* Marcum Q-Function (For integer M) -----------------
 * Sum up ~100 terms for the series form.
 * M must be positive by definition
 */
template<class T, class U>
requires std::is_arithmetic_v<T> && std::is_arithmetic_v<U>
double marcumQ(const int m, const T a, const U b) {
    // M must be positive by definition
    if (m<0) throw std::invalid_argument("Marcum Q's order M must be positive.");

    const double p = double(a)*b;
    const double s = double(a)/b;
    double t = pow(s,1-m);
    double sum = 0;
    for(int k=1-m; k<100-m; k++) {
        sum += t * std::cyl_bessel_i(abs(k),p);
        t *= s;
    }
    return exp(-0.5*(a*a+b*b)) * sum;
}


/* Marcum Q-Function (For non-integer M) ------------- 
 */
template<class T, class U, class V>
requires std::is_arithmetic_v<T> && std::is_arithmetic_v<U> && std::is_floating_point_v<V>
double marcumQ(const V m, const T a, const U b) {
    // M must be positive by definition
    if (m<0) throw std::invalid_argument("Marcum Q's order M must be positive.");

    const double aa = a*a;
    const double bb = b*b;
    double t = 1.;
    double f = 1.;   // factorial
    double sum = 0;
    for(std::size_t k=0; k<100; k++) {
        sum += t / f * regLowerGamma(m+k,0.5*bb);
        t *= 0.5*aa;
        f *= k+1;
    }
    return 1. - exp(-0.5*aa) * sum;
}


} // namespace stantanaly

#endif