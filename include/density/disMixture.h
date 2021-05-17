#ifndef STATANALY_DIS_MIXTURE_H_
#define STATANALY_DIS_MIXTURE_H_

#include "probDensFunc.h"
#include "../dContainer.h"
#include <vector>
#include <map>

namespace statanaly {

class disMixture : public probDensFunc {
    using weightType = double;
    
    // Use a container object.
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
        double r = 0;
        // for (const auto& [key, fs] : ingred) {
        //     for (const auto& f : fs) {
        //         r += f.pdf(x);
        //     }
        // }
        return r;
    }

    double cdf(const double x) const override {
        double r = 0;
        // for (const auto& [key, fs] : ingred) {
        //     for (const auto& f : fs) {
        //         r += f.cdf(x);
        //     }
        // }
        return r;
    }

    double mean() const override {
        throw std::string("Need to implement");
        return 0;
    }

    double stddev() const override {
        throw std::string("Need to implement");
        return 0;
    }

    double variance() const override {
        throw std::string("Need to implement");
        return 0;
    }

    double skewness() const override {
        throw std::string("Need to implement");
        return 0;
    }

    inline std::size_t hash() const noexcept {
        // Mixture distribution is different when the contents are different.
        // Insertion order of the ingredients does not matter.

        std::vector<std::size_t> hashes;
        for (const auto& [d,w] : ctr.get()) {
            hashes.push_back(d->hash());
        }
        std::sort(hashes.begin(), hashes.end());

        std::size_t seed = 0;
        combine_hash(seed, char(id));
        for (const auto h : hashes) {combine_hash(seed, h);}

        return seed;
    }

    void print(std::ostream& output) const override {
        output << "Mixture distribution : \n";
        for (const auto& [d,ws] : ctr.get()) {
            output << "    - " << *d << " @ weight = " << ws.second << "\n";
        }
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