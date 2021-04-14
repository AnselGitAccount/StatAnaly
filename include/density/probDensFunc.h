#ifndef STATANALY_PROBDENSFUNC_H_
#define STATANALY_PROBDENSFUNC_H_

#include "specialFunc.h"

/* Abstract base class to probability distributions */
class probDensFunc {
public:
    virtual ~probDensFunc() = default;
    virtual double pdf(const double)    = 0;
    virtual double cdf(const double)    = 0;
    virtual double mean()               = 0;
    virtual double stddev()             = 0;
    virtual double variance()           = 0;
    virtual double skewness()           = 0;
};

#endif