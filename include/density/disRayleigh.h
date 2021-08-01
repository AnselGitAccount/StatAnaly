#ifndef STATANALY_DIS_RAYLEIGH_H_
#define STATANALY_DIS_RAYLEIGH_H_

#include "probDensFunc.h"
#include "disChi.h"


namespace statanaly {

class disRayleigh : public probDensFunc {
private:
    double sigma;   // scale

public:
    template<class T>
    requires std::is_arithmetic_v<T>
    disRayleigh(const T scale) {
        sigma = scale;
    }
    disRayleigh() = delete;
    ~disRayleigh() = default;

    double pdf(const double x) const override {
        const double inv_ss = 1/(sigma*sigma);
        const double r = x*inv_ss * std::exp(-0.5*x*x*inv_ss);
        return r;
    }

    double cdf(const double x) const override {
        const double s = std::exp(-0.5*x*x/(sigma*sigma));
        return 1-s;
    }

    double mean() const override {
        constexpr double s = std::sqrt(M_PI/2);
        return sigma*s; 
    }

    double stddev() const override {
        constexpr double s = std::sqrt(2-M_PI_2);
        return s*sigma;
    }

    double variance() const override {
        constexpr double s = 2-M_PI_2;
        return s*sigma*sigma;
    }

    constexpr double skewness() const override {
        constexpr double s = 2.*std::sqrt(M_PI)*(M_PI-3) / std::sqrt((4-M_PI)*(4-M_PI)*(4-M_PI));
        return s;
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, sigma);
        return seed;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disRayleigh>(static_cast<disRayleigh const&>(*this));
    };

    disRayleigh* clone() const override {
        return new disRayleigh(*this);
    }

    void print(std::ostream& output) const override {
        output << "Rayleigh distribution -- sigma = " << sigma;
    }

    bool isEqual_tol(const probDensFunc& o, const double tol) const override {
        const disRayleigh& oo = dynamic_cast<const disRayleigh&>(o);
        return isEqual_fl_tol(sigma, oo.sigma, tol);
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned ulp) const override {
        const disRayleigh& oo = dynamic_cast<const disRayleigh&>(o);
        return isEqual_fl_ulp(sigma, oo.sigma, ulp);
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::RAYLEIGH_DISTR;

    // Rayleigh distribution is equivalent to Central Chi distribution 
    // with the same sigma and degree-of-freedom=2.
    explicit operator disChi() const {
        return disChi(2,sigma);
    }
};

}   // namespace statanaly


template<>
class std::hash<statanaly::disRayleigh> {
public:
    std::size_t operator() (const statanaly::disRayleigh& d) const {
        return d.hash();
    }
};


#endif