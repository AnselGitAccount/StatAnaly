#ifndef STATANALY_UNIFORM_H_
#define STATANALY_UNIFORM_H_

#include "probDensFunc.h"


namespace statanaly {

/* Continuous Standard Uniform Distribution */
class disStdUniform : public probDensFunc {
public:
    double pdf(const double x=0) const override {
        if (0>x || 1<x) {return 0;}
        return 1;
    }
    
    double cdf(const double x=0) const override {
        if (0>x) return 0;
        if (1<x) return 1;
        return x;
    }

    double mean() const override {
        return 0.5;
    }

    double stddev() const override {
        const double r = sqrt(1./12);
        return r;
    }

    double variance() const override {
        constexpr double r = 1./12;
        return r;
    }

    double skewness() const override {
        return 0;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disStdUniform>(static_cast<disStdUniform const&>(*this));
    };

    disStdUniform* clone() const override {
        return new disStdUniform(*this);
    }

    const dFuncID id = dFuncID::STD_UNIFORM_DISTR;
};


/* Continuous Uniform distribution  */
class disUniform : public probDensFunc {
private:
    double a;
    double b;

public:
    template<class T>
    requires std::is_arithmetic_v<T>
    disUniform(const T lower, const T upper) {
        a = lower;
        b = upper;
        assert( a!=b );
    }
    disUniform() = delete;
    ~disUniform() = default;

    constexpr double pdf(const double x) const override  {
        if (a>x || b<x) {return 0;}
        return 1/(b-a);
    }

    constexpr double cdf(const double x) const override {
        if (a>x) return 0;
        if (b<x) return 1;
        return (x-a)/(b-a);
    }

    constexpr double mean() const override {
        return 0.5*(a+b);
    }

    constexpr double stddev() const override {
        return (b-a)/sqrt(12.);
    }

    constexpr double variance() const override {
        return (b-a)*(b-a)/12.;
    }

    constexpr double skewness() const override {
        return 0;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disUniform>(static_cast<disUniform const&>(*this));
    };

    disUniform* clone() const override {
        return new disUniform(*this);
    }

    const dFuncID id = dFuncID::UNIFORM_DISTR;
};


}


#endif 