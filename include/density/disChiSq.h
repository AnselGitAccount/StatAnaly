#ifndef STATANALY_DIS_CHISQ_H_
#define STATANALY_DIS_CHISQ_H_

#include "probDensFunc.h"
// #include <bits/floatn-common.h>
// #include <bits/stdint-uintn.h>

namespace statanaly {

class disChiSq : public probDensFunc {
private:
    double k;    // parameter -- degree of freedom

public:
    template<class T> 
    requires std::is_integral<T>::value
    disChiSq(const T k_){
        k = k_;
    }
    ~disChiSq() = default;

    double pdf(const double x) const override {
        const double kd2 = k*0.5;
        const double lgf = logGamma(kd2);
        const double r   = exp((kd2-1.0)*log(x) - x*0.5 - kd2*M_LN2 - lgf);
        return r;
    }

    double cdf(const double x) const override {
        return 0;
    }

    double mean() const override {
        return k;
    }

    double stddev() const override {
        return sqrt(k*2.);
    }

    double variance() const override {
        return k*2.;
    }

    double skewness() const override {
        return sqrt(8./k);
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disChiSq>(static_cast<disChiSq const&>(*this));
    };

    disChiSq* clone() const override {
        return new disChiSq(*this);
    }

    const dFuncID id = dFuncID::CHISQ_DISTR;
};

}
/* TODO:
 * - Implement disChiSq::compute_cdf().
 */

#endif