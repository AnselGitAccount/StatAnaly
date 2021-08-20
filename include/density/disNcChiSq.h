#ifndef STATANALY_DIS_NONCENTRAL_CHI_SQUARED_H_
#define STATANALY_DIS_NONCENTRAL_CHI_SQUARED_H_

#include "probDensFunc.h"


namespace statanaly {

/* Non-central Chi Squared Distribution */

class disNcChiSq : public probDensFunc {
private:
    unsigned k;     // degree of freedom
    double lambda;  // distance

public:
    template<class T, class U>
    requires std::is_integral_v<T> && std::is_arithmetic_v<U>
    disNcChiSq(const T dof, const U distance) {
        k = dof;
        lambda = distance;
    }
    disNcChiSq() = delete;
    ~disNcChiSq() = default;

    double pdf(const double x) const override {
        const double kp = 0.5*k - 1.;
        const double t = 0.5 * pow(x/lambda, 0.5*kp) * exp(-0.5*(x+lambda));
        return t * std::cyl_bessel_i(kp, std::sqrt(lambda*x));
    }

    double cdf(const double x) const override {
        return 1. - marcumQ(0.5*k, std::sqrt(lambda), std::sqrt(x));
    }

    double mean() const override {
        return k+lambda;
    }

    double stddev() const override {
        return std::sqrt(variance());
    }

    double variance() const override {
        return 2*(k+2*lambda);
    }

    double skewness() const override {
        return (k+3*lambda) / std::sqrt(0.125*(k+2*lambda)*(k+2*lambda)*(k+2*lambda));
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, k);
        combine_hash(seed, lambda);
        return seed;
    } 

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disNcChiSq>(static_cast<disNcChiSq const&>(*this));
    };

    disNcChiSq* clone() const override {
        return new disNcChiSq(*this);
    }

    void print(std::ostream& output) const override {
        output << "Non-central Chi-Sqaured distribution -- k = " << k << " lambda = " << lambda;
    }

    bool isEqual_tol(const probDensFunc& o, const double tol) const override {
        const disNcChiSq& oo = dynamic_cast<const disNcChiSq&>(o);
        bool r = true;
        r &= (k == oo.k);
        r &= isEqual_fl_tol(lambda, oo.lambda, tol);
        return r;
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned ulp) const override {
        const disNcChiSq& oo = dynamic_cast<const disNcChiSq&>(o);
        bool r = true;
        r &= (k == oo.k);
        r &= isEqual_fl_ulp(lambda, oo.lambda, ulp);
        return r;
    }

    auto p_dof() const noexcept {
        return k;
    }

    auto p_distance() const noexcept {
        return lambda;
    }
};

}   // namesapce statanaly


template<>
class std::hash<statanaly::disNcChiSq> {
public:
    std::size_t operator() (const statanaly::disNcChiSq& d) const {
        return d.hash();
    }
};


#endif