#ifndef STATANALY_DIS_CHISQ_H_
#define STATANALY_DIS_CHISQ_H_

#include "probDensFunc.h"
// #include <bits/floatn-common.h>
// #include <bits/stdint-uintn.h>

namespace statanaly {

class disChiSq : public probDensFunc {
private:
    unsigned k;    // parameter -- degree of freedom

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

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, k);
        return seed;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disChiSq>(static_cast<disChiSq const&>(*this));
    };

    disChiSq* clone() const override {
        return new disChiSq(*this);
    }

    void print(std::ostream& output) const override {
        output << "Chi Square distribution -- k = " << k;
    }

    bool isEqual_tol(const probDensFunc& o, const double tol=0) const override {
        const disChiSq& oo = dynamic_cast<const disChiSq&>(o);
        return k==oo.k;
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned nlp=0) const override {
        const disChiSq& oo = dynamic_cast<const disChiSq&>(o);
        return k==oo.k;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::CHISQ_DISTR;
};

}   // namespace statanaly


template<>
class std::hash<statanaly::disChiSq> {
public:
    std::size_t operator() (const statanaly::disChiSq& d) const {
        return d.hash();
    }
};


/* TODO:
 * - Implement disChiSq::compute_cdf().
 */

#endif