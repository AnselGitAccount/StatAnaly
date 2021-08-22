#ifndef STATANALY_DIS_ERLANG_H_
#define STATANALY_DIS_ERLANG_H_

#include "probDistr.h"

namespace statanaly {

/**
 * @brief Erlang distribution
 * 
 * @param k Shape
 * @param lambda Rate
 * 
 */

class disErlang : public probDistr {
private:
    unsigned k;     // shape
    double lambda;  // rate

public:
    template<typename T, typename P>
    requires std::is_integral_v<T> && std::is_arithmetic_v<P>
    disErlang(const T shape, const P rate) {
        k = shape;
        lambda = rate;
    }
    disErlang() = delete;


    double pdf(const double x) const override {
        double r = pow(lambda,k) / factorial[k-1] * pow(x,k-1) / exp(lambda*x);
        return r;
    }

    double cdf(const double x) const override {
        return lowerGamma(k,lambda*x) / factorial[k-1];
    }

    double mean() const override {
        return k/lambda;
    }

    double stddev() const override {
        return sqrt(variance());
    }

    double variance() const override {
        return k/lambda/lambda;
    }

    double skewness() const override {
        return 2./sqrt(k);
    }

    auto pshape() const noexcept {return k;}
    auto prate() const noexcept {return lambda;}

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, k);
        combine_hash(seed, lambda);
        return seed;
    }

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disErlang>(static_cast<disErlang const&>(*this));
    };

    disErlang* clone() const override {
        return new disErlang(*this);
    }

    void print(std::ostream& output) const override {
        output << "Erlang distribution -- k = " << k << " lambda = " << lambda;
    }

    bool isEqual_tol(const probDistr& o, const double tol) const override {
        const disErlang& oo = dynamic_cast<const disErlang&>(o);
        bool r = true;
        r &= k == oo.k;
        r &= isEqual_fl_tol(lambda, oo.lambda, tol);
        return r;
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp) const override {
        const disErlang& oo = dynamic_cast<const disErlang&>(o);
        bool r = true;
        r &= k == oo.k;
        r &= isEqual_fl_ulp(lambda, oo.lambda, ulp);
        return r;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::ERLANG_DISTR;
};

}   // namespace statanaly


/**
 * @brief STL hasher overload
 * 
 * @tparam Erlang distribution
 */

template<>
class std::hash<statanaly::disErlang> {
public:
    std::size_t operator() (const statanaly::disErlang& d) const {
        return d.hash();
    }
};

#endif