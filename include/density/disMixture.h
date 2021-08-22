#ifndef STATANALY_DIS_MIXTURE_H_
#define STATANALY_DIS_MIXTURE_H_

#include "probDistr.h"
#include "dContainer.h"
#include <vector>
#include <map>

namespace statanaly {

/**
 * @brief Mixture of Distributions
 * 
 * A Mixture distribution is a distribution that contains a collection of distributions.
 * Each distribution has an associated weight.
 * PDF of a mixture is the weighted sum of pdf of each component.
 * CDF of a mixture is the weighted sum of cdf of components. 
 * Mean of a mixture is the weighted sum of mean of components. 
 * Variance of a mixture is computed via the Law of Total Variance.
 * 
 * @param ctr Container for a collection of distributions.
 */

class disMixture : public probDistr {
    using weightType = double;
    
    dCtr ctr;

public:
    disMixture() = default;

    /** Copy constructor: deep-copy, do the same as clone(). */
    disMixture(const disMixture& o) {
        // clone the container.
        ctr = o.ctr;
    }

    /** Copy assignment: deep-copy */
    disMixture& operator = (const disMixture& o) {
        // clone the container
        ctr = o.ctr;
        return *this;
    };
    
    /** Move constructor.
     * Move the named distribution.
     */
    disMixture(disMixture&& o) {
        ctr = std::move(o.ctr);
    }

    /** Move assignment */
    disMixture& operator = (disMixture&& o) {
        // Move the named distribution.
        ctr = std::move(o.ctr);
        return *this;
    }

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disMixture>(static_cast<disMixture const&>(*this));
    };

    disMixture* clone() const override {
        return new disMixture(*this);
    };

    /** Insert a distribution and its weight */
    template<typename F, typename W>
    requires std::is_arithmetic_v<W>
    void insert(F&& distr, W weight) {
        // Forward to container's insert.
        ctr.insert( std::forward<F>(distr), weight);
    }

    /** Find a distribution in the mixture.
     * Check each component's hash (ie, type and parameters). 
     */
    template<typename F>
    inline auto find(F&& distr) const {
        return ctr.find( std::forward<F>(distr) );
    }

    inline const auto& get() const {return ctr.get();}
    inline const auto end() const {return ctr.end();}
    inline void clear() {ctr.clear();}

    /** pdf of a mixture is the weighted sum of pdf of each component. */
    double pdf(const double x) const override {
        double res = 0;
        for (const auto& [d, ws] : ctr.get()) {
            const double w = ws.second;
            res += d->pdf(x) * w;
        }
        return res;
    }

    /** cdf of a mixture is the weighted sum of cdf of each component. */
    double cdf(const double x) const override {
        double res = 0;
        for (const auto& [d, ws] : ctr.get()) {
            const double w = ws.second;
            res += d->cdf(x) * w;
        }
        return res;
    }

    /** mean of a mixture is the weighted sum of mean of each component. */
    double mean() const override {
        double res = 0;
        for (const auto & [d, ws] : ctr.get()) {
            const double w = ws.second;
            res += d->mean() * w;
        }
        return res;
    }

    double stddev() const override {
        double res = std::sqrt(variance());
        return res;
    }

    /** Variance of a mixture is computed via the Law of Total Variance. */
    double variance() const override {
        double tmp = 0;
        double wmu = 0;
        for (const auto & [d, ws] : ctr.get()) {
            const double w = ws.second;
            const double m = d->mean();
            wmu += m * w;
            tmp += (d->variance() + m * m) * w;
        }
        return tmp - wmu * wmu;
    }

    double skewness() const override {
        double tmp = 0;
        for (const auto & [d, ws] : ctr.get()) {
            const double w = ws.second;
            const double s = d->stddev();
            const double m = d->mean();
            tmp += w*s*s*s*d->skewness() + 3*w*m*s*s + w*m*m*m;
        }
        const double s = stddev();
        const double m = mean();
        tmp += - 3*m*s*s - m*m*m;
        tmp /= s*s*s;
        return tmp;
    }

    /** See dContainer hash() */
    inline std::size_t hash() const noexcept {
        // disMixture only contains a dContainer.
        // So the hash should be that of the dContainer.

        return ctr.hash();
    }

    void print(std::ostream& output) const override {
        output << "Mixture distribution : \n";
        for (const auto& [d,ws] : ctr.get()) {
            output << "    - " << *d << " @ weight = " << ws.second << "\n";
        }
    }

    bool isEqual_tol(const probDistr& o, const double tol) const override {
        const disMixture& oo = dynamic_cast<const disMixture&>(o);
        return ctr.isEqual_tol(oo.ctr, tol);
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp) const override {
        const disMixture& oo = dynamic_cast<const disMixture&>(o);
        return ctr.isEqual_ulp(oo.ctr, ulp);
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::MIXTURE_DISTR;
};

} // namespace 


/**
 * @brief STL hasher overload
 * 
 * @tparam Mixture distribution
 */

template<>
class std::hash<statanaly::disMixture> {
public:
    std::size_t operator() (const statanaly::disMixture& d) const {
        return d.hash();
    }
};


#endif