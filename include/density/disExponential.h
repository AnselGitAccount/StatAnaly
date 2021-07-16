#ifndef STATANALY_DIS_EXPONENTIAL_H_
#define STATANALY_DIS_EXPONENTIAL_H_

#include "probDensFunc.h"

namespace statanaly {

class disExponential : public probDensFunc {
private:
    double lambda;      // rate

public:
    template<class T> 
    requires std::is_arithmetic_v<T>
    disExponential(const T rate){
        lambda = rate;
    }
    disExponential() = delete;


    double pdf(const double x) const override {
        return lambda * exp(-lambda*x);
    }

    double cdf(const double x) const override {
        return 1. - exp(-lambda*x);
    }

    double mean() const override {
        return 1./lambda;
    }

    double stddev() const override {
        return 1./lambda;
    }

    double variance() const override {
        return 1./(lambda*lambda);
    }

    double skewness() const override {
        return 2;
    }

    auto prate() const noexcept {return lambda;}

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, lambda);
        return seed;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disExponential>(static_cast<disExponential const&>(*this));
    };

    disExponential* clone() const override {
        return new disExponential(*this);
    }

    void print(std::ostream& output) const override {
        output << "Exponential distribution -- lambda = " << lambda;
    }

    bool isEqual_tol(const probDensFunc& o, const double tol) const override {
        const disExponential& oo = dynamic_cast<const disExponential&>(o);
        return isEqual_fl_tol(lambda, oo.lambda, tol);
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned ulp) const override {
        const disExponential& oo = dynamic_cast<const disExponential&>(o);
        return isEqual_fl_ulp(lambda, oo.lambda, ulp);
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::EXPONENTIAL_DISTR;
};

} // namespace statanaly


template<>
class std::hash<statanaly::disExponential> {
public:
    std::size_t operator() (const statanaly::disExponential& d) const {
        return d.hash();
    }
};


#endif