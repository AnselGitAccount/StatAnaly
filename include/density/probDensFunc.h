#ifndef STATANALY_PROBDENSFUNC_H_
#define STATANALY_PROBDENSFUNC_H_

#include "specialFunc.h"
#include "fl_comparison.h"
#include "hasher.h"
#include <memory>

namespace statanaly {

enum class dFuncID {
    BASE_DISTR,
    NORMAL_DISTR,
    STD_UNIFORM_DISTR,
    UNIFORM_DISTR,
    CHI_DISTR,
    CHISQ_DISTR,
    MIXTURE_DISTR,
    IRWIN_HALL,
    CAUCHY_DISTR,
    GAMMA_DISTR,
    EXPONENTIAL_DISTR,
    ERLANG_DISTR,
    RAYLEIGH_DISTR,
    COUNT
};


class probDensFunc {
public:
    virtual ~probDensFunc() = default;

    virtual double pdf(const double=0) const    = 0;
    virtual double cdf(const double=0) const    = 0;
    virtual double mean() const                 = 0;
    virtual double stddev() const               = 0;
    virtual double variance() const             = 0;
    virtual double skewness() const             = 0;

    virtual std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, id);
        return seed;
    }
    
    virtual std::unique_ptr<probDensFunc> cloneUnique() const = 0;
    virtual probDensFunc* clone() const = 0;    // Return type can be Covariant.

    // Virtual friend idiom -- The friendship will get pass down to derived classes.
    friend std::ostream& operator << (std::ostream&, const probDensFunc&);
    virtual void print(std::ostream&) const = 0;

    virtual bool isEqual_tol(const probDensFunc&, const double) const = 0;
    virtual bool isEqual_ulp(const probDensFunc&, const unsigned) const = 0;

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::BASE_DISTR;

protected:
    // This class cannot be instantiated. Only be inherited.
    probDensFunc() = default;
    probDensFunc(const probDensFunc&) = default;
    probDensFunc(probDensFunc&&) = default;
};


// Human friendly text.
std::ostream& operator << (std::ostream& output, const probDensFunc& distr);

}

#endif