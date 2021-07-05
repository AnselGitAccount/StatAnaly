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
    disNormal(const T mean, const T variance){
        mu = mean;
        sig = sqrt(variance);
    }
    disNormal() = delete;
    ~disNormal() = default;

    double pdf(const double x) const override {
        const double x_scaled = ((x-mu)/sig) * ((x-mu)/sig);
        const double r = exp(-log(sig) - SACV_LOG_SQRT_2PI - 0.5*x_scaled);
        return r;
    }

    double cdf(const double x) const override {
        const double scaledx = (x-mu)/sig;
        return 0.5 * (1. + std::erf(scaledx * M_SQRT1_2));
    }

    double mean() const override {
        return mu;
    }

    double stddev() const override {
        return sig;
    }

    double variance() const override {
        return sig*sig;
    }

    double skewness() const override {
        return 0;
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, mu);
        combine_hash(seed, sig);
        return seed;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disNormal>(static_cast<disNormal const&>(*this));
    };

    disNormal* clone() const override {
        return new disNormal(*this);
    }

    void print(std::ostream& output) const override {
        output << "Normal distribution -- mu = " << mu << "  sig = " << sig;
    }

    bool isEqual_tol(const probDensFunc& o, const double tol) const override {
        const disNormal& oo = dynamic_cast<const disNormal&>(o);
        bool r = true;
        r &= isEqual_fl_tol(mu, oo.mu, tol);
        r &= isEqual_fl_tol(sig, oo.sig, tol);
        return r;
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned ulp) const override {
        const disNormal& oo = dynamic_cast<const disNormal&>(o);
        bool r = true;
        r &= isEqual_fl_ulp(mu, oo.mu, ulp);
        r &= isEqual_fl_ulp(sig, oo.sig, ulp);
        return r;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::NORMAL_DISTR;
};

} // namespace statanaly


template<>
class std::hash<statanaly::disNormal> {
public:
    std::size_t operator() (const statanaly::disNormal& d) const {
        return d.hash();
    }
};


#endif