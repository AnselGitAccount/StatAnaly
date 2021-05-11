#ifndef STATANALY_PROBDENSFUNC_H_
#define STATANALY_PROBDENSFUNC_H_

#include "specialFunc.h"


namespace statanaly {

enum class dFuncID {
    BASE_DISTR,
    NORMAL_DISTR,
    STD_UNIFORM_DISTR,
    UNIFORM_DISTR,
    CHISQ_DISTR,
    MIXTURE_DISTR,
    COUNT
};


class probDensFunc {
public:
    virtual ~probDensFunc() = default;

    virtual double pdf(const double=0) const    = 0;
    virtual double cdf(const double=0) const    = 0;
    virtual double mean() const                 = 0;
    virtual double stddev() const               = 0;
    virtual double variance() const             = 0;
    virtual double skewness() const             = 0;
    
    
    virtual std::unique_ptr<probDensFunc> cloneUnique() const = 0;
    virtual probDensFunc* clone() const = 0;    // Return type can be Covariant.

    const dFuncID id = dFuncID::BASE_DISTR;

protected:
    // This class cannot be instantiated. Only be inherited.
    probDensFunc() = default;
    probDensFunc(const probDensFunc&) = default;
    probDensFunc(probDensFunc&&) = default;
};

}

#endif