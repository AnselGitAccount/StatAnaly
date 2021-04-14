#ifndef STATANALY_DIS_CHISQ_H_
#define STATANALY_DIS_CHISQ_H_

#include "probDensFunc.h"
// #include <bits/floatn-common.h>
// #include <bits/stdint-uintn.h>

namespace statanaly {

template<class T> 
requires std::is_integral<T>::value
class disChiSq : probDensFunc {
private:
    double k;    // parameter -- degree of freedom

public:
    disChiSq(const T k_){
        k = k_;
    }
    ~disChiSq() = default;

    double pdf(const double x) override {
        const double kd2 = k*0.5;
        const double lgf = logGamma(kd2);
        const double r   = exp((kd2-1.0)*log(x) - x*0.5 - kd2*M_LN2 - lgf);
        return r;
    }

    double cdf(const double x) override {
        return 0;
    }

    double mean() override {
        return k;
    }

    double stddev() override {
        return sqrt(k*2.);
    }

    double variance() override {
        return k*2.;
    }

    double skewness() override {
        return sqrt(8./k);
    }
    
};

}
/* TODO:
 * - Implement disChiSq::compute_cdf().
 */

#endif