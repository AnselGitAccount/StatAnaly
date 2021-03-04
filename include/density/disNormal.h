#ifndef STATANALY_DEN_NORMAL_H_
#define STATANALY_DEN_NORMAL_H_

#include "probDensFunc.h"


namespace statanaly {

template<class T> 
requires std::is_floating_point<T>::value
class disNormal : probDensFunc {
private:
    double mu;      // parameter -- mean
    double sig;     // parameter -- std deviation

public:
    disNormal(const T mean, const T stdDev){
        mu = mean;
        sig = stdDev;
    }
    ~disNormal() = default;

    double pdf(const double x) override {
        const double x_scaled = ((x-mu)/sig) * ((x-mu)/sig);
        const double r = exp(-log(sig) - SACV_LOG_SQRT_2PI - 0.5*x_scaled);
        return r;
    }

    double cdf(const double x) override {
        return 0;
    }

    double mean() override {
        return mu;
    }

    double variance() override {
        return mu*mu;
    }

    double skewness() override {
        return 0;
    }
    
};

} // namespace statanaly

/* TODO:
 * - Implement disNormal::cdf().
 */

#endif