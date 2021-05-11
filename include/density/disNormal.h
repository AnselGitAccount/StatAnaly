#ifndef STATANALY_DEN_NORMAL_H_
#define STATANALY_DEN_NORMAL_H_

#include "probDensFunc.h"


namespace statanaly {

class disNormal : public probDensFunc {
private:
    double mu;      // parameter -- mean
    double sig;     // parameter -- std deviation

public:
    template<class T> 
    requires std::is_arithmetic_v<T>
    disNormal(const T mean, const T stdDev){
        mu = mean;
        sig = stdDev;
    }
    disNormal() = delete;
    ~disNormal() = default;

    double pdf(const double x) const override {
        const double x_scaled = ((x-mu)/sig) * ((x-mu)/sig);
        const double r = exp(-log(sig) - SACV_LOG_SQRT_2PI - 0.5*x_scaled);
        return r;
    }

    double cdf(const double x) const override {
        throw std::string("Need to implement");
        return 0;
    }

    double mean() const override {
        return mu;
    }

    double stddev() const override {
        return mu;
    }

    double variance() const override {
        return mu*mu;
    }

    double skewness() const override {
        return 0;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disNormal>(static_cast<disNormal const&>(*this));
    };

    disNormal* clone() const override {
        return new disNormal(*this);
    }

    const dFuncID id = dFuncID::NORMAL_DISTR;
};

} // namespace statanaly

/* TODO:
 * - Implement disNormal::cdf().
 */

#endif