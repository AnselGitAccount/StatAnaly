#ifndef STATANALY_DIS_CHI_H_
#define STATANALY_DIS_CHI_H_

#include "probDensFunc.h"


namespace statanaly {

/* Central Chi Distribution */

class disChi : public probDensFunc {
private:
    unsigned k;     // degree of freedom

public:
    template<class T>
    requires std::is_integral_v<T>
    disChi(const T dof) {
        k = dof;
    }
    ~disChi() = default;

    double pdf(const double x) const override {
        const double ko2 = pow(2,k/2.-1) * std::tgamma(k/2.);
        return pow(x,k-1) * exp(-x*x/2) / ko2;
    }

    double cdf(const double x) const override {
        return regLowerGamma(k/2.0, x*x/2);
    }

    double mean() const override {
        return M_SQRT2 * std::tgamma((k+1)/2.) / std::tgamma(k/2.);
    }

    double stddev() const override {
        return std::sqrt(variance());
    }

    double variance() const override {
        const double mu = mean();
        return k - mu*mu;
    }

    double skewness() const override {
        const double mu = mean();
        const double st = variance();
        return mu * (1 - 2*(k-mu*mu)) / (st*st*st);
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, k);
        return seed;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disChi>(static_cast<disChi const&>(*this));
    };

    disChi* clone() const override {
        return new disChi(*this);
    }

    void print(std::ostream& output) const override {
        output << "Central Chi distribution -- k = " << k;
    }

    bool isEqual_tol(const probDensFunc& o, const double tol=0) const override {
        const disChi& oo = dynamic_cast<const disChi&>(o);
        bool r = true;
        r &= k==oo.k;
        return r;
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned ulp=0) const override {
        const disChi& oo = dynamic_cast<const disChi&>(o);
        bool r = true;
        r &= k==oo.k;
        return r;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::CHI_DISTR;
};

}   // namespace statanaly


template<>
class std::hash<statanaly::disChi> {
public:
    std::size_t operator() (const statanaly::disChi& d) const {
        return d.hash();
    }
};

#endif