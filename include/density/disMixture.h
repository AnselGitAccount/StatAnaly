#ifndef STATANALY_DIS_MIXTURE_H_
#define STATANALY_DIS_MIXTURE_H_

#include "probDensFunc.h"
#include "dContainer.h"
#include <vector>
#include <map>

namespace statanaly {

class disMixture : public probDensFunc {
    using weightType = double;
    
    dCtr ctr;

public:
    disMixture() = default;

    // Copy constructor: deep-copy, do the same as clone().
    disMixture(const disMixture& o) {
        // clone the container.
        ctr = o.ctr;
    }

    // Copy assignment: deep-copy
    disMixture& operator = (const disMixture& o) {
        // clone the container
        ctr = o.ctr;
        return *this;
    };
    
    // Move constructor
    disMixture(disMixture&& o) {
        // Move the named distribution.
        ctr = std::move(o.ctr);
    }

    // Move assignment
    disMixture& operator = (disMixture&& o) {
        // Move the named distribution.
        ctr = std::move(o.ctr);
        return *this;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disMixture>(static_cast<disMixture const&>(*this));
    };

    disMixture* clone() const override {
        return new disMixture(*this);
    };

    // Forward to container's insert.
    template<typename F, typename W>
    requires std::is_arithmetic_v<W>
    void insert(F&& distr, W weight) {
        ctr.insert( std::forward<F>(distr), weight);
    }

    // Find a distribution that match distr's hash (ie, type and parameters).
    template<typename F>
    inline auto find(F&& distr) const {
        return ctr.find( std::forward<F>(distr) );
    }

    inline const auto& get() const {return ctr.get();}
    inline const auto end() const {return ctr.end();}
    inline void clear() {ctr.clear();}

    double pdf(const double x) const override {
        // pdf of a mixture is the weighted sum of pdf of components.
        double res = 0;
        for (const auto& [d, ws] : ctr.get()) {
            const double w = ws.second;
            res += d->pdf(x) * w;
        }
        return res;
    }

    double cdf(const double x) const override {
        // cdf of a mixture is the weighted sum of cdf of components.
        double res = 0;
        for (const auto& [d, ws] : ctr.get()) {
            const double w = ws.second;
            res += d->cdf(x) * w;
        }
        return res;
    }

    double mean() const override {
        // mean of a mixture is the weighted sum of mean of components.
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

    double variance() const override {
        // Law of Total Variance
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

    inline std::size_t hash() const noexcept {
        // disMixture only contains a dContainer.
        // So the hash should be that of a dContainer.

        return ctr.hash();
    }

    void print(std::ostream& output) const override {
        output << "Mixture distribution : \n";
        for (const auto& [d,ws] : ctr.get()) {
            output << "    - " << *d << " @ weight = " << ws.second << "\n";
        }
    }

    bool isEqual_tol(const probDensFunc& o, const double tol) const override {
        const disMixture& oo = dynamic_cast<const disMixture&>(o);
        return ctr.isEqual_tol(oo.ctr, tol);
    }

    bool isEqual_ulp(const probDensFunc& o, const unsigned ulp) const override {
        const disMixture& oo = dynamic_cast<const disMixture&>(o);
        return ctr.isEqual_ulp(oo.ctr, ulp);
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::MIXTURE_DISTR;
};

} // namespace 


template<>
class std::hash<statanaly::disMixture> {
public:
    std::size_t operator() (const statanaly::disMixture& d) const {
        return d.hash();
    }
};


#endif