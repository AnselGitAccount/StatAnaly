#ifndef STATANALY_DIS_NONCENTRAL_CHI_H_
#define STATANALY_DIS_NONCENTRAL_CHI_H_

#include "probDensFunc.h"


namespace statanaly {

/* Non-Central Chi Distribution */

class disNcChi : public probDensFunc {
private:
    unsigned k;     // degree of freedom
    double lambda;  // related to mean
    double sigma;   // scale

public:
    template<class T, class U>
    requires std::is_integral_v<T> && std::is_arithmetic_v<U>
    disNcChi(const T dof, const U mean, const double scale = 1.) {
        k = dof;
        lambda = mean;
        sigma = scale;
    }
    ~disNcChi() = default;

    double pdf(const double x) const override {
        const double s2_inv = 1/(sigma*sigma);
        const double kh = 0.5*k;
        const double t = pow(lambda, 1.-kh) * s2_inv * pow(x,kh) * exp(-0.5*(x*x+lambda*lambda)*s2_inv);
        return t * std::cyl_bessel_i(kh-1., lambda*x*s2_inv);
    }

    double cdf(const double x) const override {
        return 1. - marcumQ(0.5*k,lambda,x);
    }

    double mean() const override {
        /* https://en.wikipedia.org/wiki/Noncentral_chi_distribution
         * https://math.stackexchange.com/questions/3187779/associated-laguerre-polynomials-of-half-integer-parameters
         * https://math.stackexchange.com/questions/2574995/expectation-noncentral-chi-distribution-reference-request
         * https://en.wikipedia.org/wiki/Confluent_hypergeometric_function
         */
        throw std::runtime_error("Non-central Chi distribution's mean is not implemented currently.");
        return 0;
    }

    double stddev() const override {
        throw std::runtime_error("Non-central Chi distribution's standard deviation is not implemented currently.");
        return 0;
    }

    double variance() const override {
        throw std::runtime_error("Non-central Chi distribution's variance is not implemented currently.");
        return 0;
    }

    double skewness() const override {
        throw std::runtime_error("Non-central Chi distribution's skewness is not implemented currently.");
        return 0;
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, k);
        combine_hash(seed, lambda);
        combine_hash(seed, sigma);
        return seed;
    } 

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disNcChi>(static_cast<disNcChi const&>(*this));
    };

    disNcChi* clone() const override {
        return new disNcChi(*this);
    }

    void print(std::ostream& output) const override {
        output << "Non-central Chi distribution -- k = " << k << " lambda = " << lambda << " sigma = " << sigma;
    }

    bool isEqual_tol(const probDensFunc& o, const double tol) const override {
        const disNcChi& oo = dynamic_cast<const disNcChi&>(o);
        bool r = true;
        r &= (k == oo.k);
        r &= isEqual_fl_tol(lambda, oo.lambda, tol);
        r &= isEqual_fl_tol(sigma, oo.sigma, tol);
        return r;
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned ulp) const override {
        const disNcChi& oo = dynamic_cast<const disNcChi&>(o);
        bool r = true;
        r &= (k == oo.k);
        r &= isEqual_fl_ulp(lambda, oo.lambda, ulp);
        r &= isEqual_fl_ulp(sigma, oo.sigma, ulp);
        return r;
    }
};

}   // namesapce statanaly


template<>
class std::hash<statanaly::disNcChi> {
public:
    std::size_t operator() (const statanaly::disNcChi& d) const {
        return d.hash();
    }
};


#endif