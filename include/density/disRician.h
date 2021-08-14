#ifndef STATANALY_DIS_RICIAN_H_
#define STATANALY_DIS_RICIAN_H_

#include "probDensFunc.h"
#include "disNcChi.h"


namespace statanaly {

class disRician : public probDensFunc {
private:
    double nu;      // distance
    double sigma;   // spread

public:
    template<class T, class U>
    requires std::is_arithmetic_v<T> && std::is_arithmetic_v<U> 
    disRician(const T distance, const U spread) {
        nu = std::abs(distance);
        sigma = spread;
    }
    disRician() = delete;
    ~disRician() = default;

    double pdf(const double x) const override {
        const double s2inv = 1/(sigma*sigma);
        const double x_e = exp(-(x*x+nu*nu)*0.5*s2inv) * std::cyl_bessel_i(0,x*nu*s2inv);
        return x * s2inv * x_e;
    }

    double cdf(const double x) const override {
        const double s_inv = 1/sigma;
        return 1 - marcumQ(1, nu*s_inv, x*s_inv);
    }

    double mean() const override {
        const double x = -0.5*nu*nu/(sigma*sigma);
        const double lague = exp(x/2) * 
            ((1-x)*std::cyl_bessel_i(0,-0.5*x) - x*std::cyl_bessel_i(1,-0.5*x));
        return sigma * SQRT_PI_2 * lague;
    }

    double stddev() const override {
        return std::sqrt(variance());
    }

    double variance() const override {
        const double x = -0.5*nu*nu/(sigma*sigma);
        const double lague = exp(x/2) * 
            ((1-x)*std::cyl_bessel_i(0,-0.5*x) - x*std::cyl_bessel_i(1,-0.5*x));
        return 2*sigma*sigma + nu*nu - M_PI_2*sigma*sigma*lague*lague;
    }

    double skewness() const override {
        throw std::runtime_error("Rician distribution's skewness is too complicated.");
        return 0;
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, nu);
        combine_hash(seed, sigma);
        return seed;
    } 

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disRician>(static_cast<disRician const&>(*this));
    };

    disRician* clone() const override {
        return new disRician(*this);
    }

    void print(std::ostream& output) const override {
        output << "Rician distribution -- nu = " << nu << " sigma = " << sigma;
    }

    bool isEqual_tol(const probDensFunc& o, const double tol) const override {
        const disRician& oo = dynamic_cast<const disRician&>(o);
        bool r = true;
        r &= isEqual_fl_tol(nu, oo.nu, tol);
        r &= isEqual_fl_tol(sigma, oo.sigma, tol);
        return r;
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned ulp) const override {
        const disRician& oo = dynamic_cast<const disRician&>(o);
        bool r = true;
        r &= isEqual_fl_ulp(nu, oo.nu, ulp);
        r &= isEqual_fl_ulp(sigma, oo.sigma, ulp);
        return r;
    }
};

}   // namespace statanaly


template<>
class std::hash<statanaly::disRician> {
public:
    std::size_t operator() (const statanaly::disRician& d) const {
        return d.hash();
    }
};


#endif