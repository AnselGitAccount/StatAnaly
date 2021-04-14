#ifndef STATANALY_UNIFORM_H_
#define STATANALY_UNIFORM_H_

#include "probDensFunc.h"


namespace statanaly {

/* Continuous Standard Uniform Distribution */
class disStdUniform : probDensFunc {
public:
    constexpr double pdf(const double x=0) override {
        if (0>x || 1<x) {return 0;}
        return 1;
    }
    
    double cdf(const double x) override {
        if (0>x) return 0;
        if (1<x) return 1;
        return x;
    }

    constexpr double mean() override {
        return 0.5;
    }

    constexpr double stddev() override {
        const double r = sqrt(1./12);
        return r;
    }

    constexpr double variance() override {
        const double r = 1./12;
        return r;
    }

    constexpr double skewness() override {
        return 0;
    }
};


/* Continuous Uniform distribution  */
template<class T>
requires std::is_arithmetic_v<T>
class disUniform : probDensFunc {
private:
    double a;
    double b;

public:
    disUniform(const T lower, const T upper) {
        a = lower;
        b = upper;
        assert( a!=b );
    }
    ~disUniform() = default;

    constexpr double pdf(const double x) override {
        if (a>x || b<x) {return 0;}
        return 1/(b-a);
    }

    constexpr double cdf(const double x) override {
        if (a>x) return 0;
        if (b<x) return 1;
        return (x-a)/(b-a);
    }

    constexpr double mean() override {
        return 0.5*(a+b);
    }

    constexpr double stddev() override {
        return (b-a)/sqrt(12.);
    }

    constexpr double variance() override {
        return (b-a)*(b-a)/12.;
    }

    constexpr double skewness() override {
        return 0;
    }
};


}


#endif 