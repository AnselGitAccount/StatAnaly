#ifndef STATANALY_DIS_MIXTURE_H_
#define STATANALY_DIS_MIXTURE_H_

#include "probDensFunc.h"

#include <vector>
#include <map>

namespace statanaly {

class disMixture : public probDensFunc {
    using weightType = double;
    std::vector<weightType> weights;
    
    // Stores the list of named distribution functions 
    // that made up the mixture distribution.
    std::vector<std::pair<std::unique_ptr<probDensFunc>, weightType>> ingreds;

public:
    disMixture() = default;
    
    // Copy constructor: deep-copy, do the same as clone().
    disMixture(const disMixture& o) {
        // clone the named distribution.
        for (const auto& [d,w] : o.ingreds) {
            ingreds.push_back(std::make_pair(d->cloneUnique(),w));
        }
        weights = o.weights;
    }

    // Copy assignment: deep-copy
    disMixture& operator = (const disMixture& o) {
        // clone the named distribution.
        for (const auto& [d,w] : o.ingreds) {
            ingreds.push_back(std::make_pair(d->cloneUnique(),w));
        }
        weights = o.weights;
        return *this;
    };
    
    // Move constructor
    disMixture(disMixture&& o) {
        // Move the named distribution.
        ingreds = std::move(o.ingreds);
        weights = std::move(o.weights);
    }

    // Move assignment
    disMixture& operator = (disMixture&& o) {
        // Move the named distribution.
        ingreds = std::move(o.ingreds);
        weights = std::move(o.weights);
        return *this;
    }

    std::unique_ptr<probDensFunc> cloneUnique() const override {
        return std::make_unique<disMixture>(static_cast<disMixture const&>(*this));
    };

    disMixture* clone() const override {
        return new disMixture(*this);
    };

    // disMixture makes a deep copy of the distribution functions that is passing in.
    template<typename F, typename W>
    requires std::is_arithmetic_v<W>
    void insert(const F* distr, W weight) {
        // Make a deep-copy
        ingreds.push_back(
            std::make_pair(distr->cloneUnique(), static_cast<weightType>(weight)) );

        // Rescale the weights so that they sum up to one.
        rescale();
    }

    template<typename F, typename W>
    requires std::is_arithmetic_v<W>
    void insert(const std::unique_ptr<F>& distr, W weight) {
        ingreds.push_back( 
            std::make_pair(distr->cloneUnique(), static_cast<weightType>(weight)) );    

        rescale();
    }

    void rescale() {
        weightType sum = 0;
        for (const auto& [d,w] : ingreds) {sum += w;}

        weights.resize(ingreds.size());
        for (std::size_t i=0; i<ingreds.size(); i++) {    
            weights[i] = (ingreds[i].second) / sum;
        }
    }

    void clear() {
        ingreds.clear();
        weights.clear();
    }

    // // Map the content to a 
    // std::vector<int> remap() {
    //     std::vector<int> cnts(static_cast<std::size_t>(dFuncID::COUNT), 0);
    //     for (const auto& [d,w] : ingreds) {
    //         switch (d->id) {
    //             case dFuncID::BASE_DISTR :
    //                 cnts[ static_cast<std::size_t>(dFuncID::BASE_DISTR) ]++;
    //                 break;
    //             case dFuncID::NORMAL_DISTR :
    //                 cnts[ static_cast<std::size_t>(dFuncID::NORMAL_DISTR) ]++;
    //                 break;
    //             case dFuncID::STD_UNIFORM_DISTR :
    //                 cnts[ static_cast<std::size_t>(dFuncID::STD_UNIFORM_DISTR) ]++;
    //                 break;
    //             case dFuncID::UNIFORM_DISTR :
    //                 cnts[ static_cast<std::size_t>(dFuncID::UNIFORM_DISTR) ]++;
    //                 break;
    //             case dFuncID::CHISQ_DISTR :
    //                 cnts[ static_cast<std::size_t>(dFuncID::CHISQ_DISTR) ]++;
    //                 break;
    //             case dFuncID::MIXTURE_DISTR :
    //                 cnts[ static_cast<std::size_t>(dFuncID::MIXTURE_DISTR) ]++;
    //                 break;
    //         }
    //     }
    //     return cnts;
    // }

    // // Encode the distribution into a string.
    // std::string encode() {
    //     return std::string{};
    // }

    // void print() {
    //     auto cnts = encode();

    //     // Print to screen
    //     for (std::size_t i=0; i<cnts.size(); i++) {printf("%d ",cnts[i]);}
    //     printf("\n");

    //     // TODO: overload << operators for every distribution.
    //     // TODO: Encode each distribution.
    // }

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

    const dFuncID id = dFuncID::MIXTURE_DISTR;
};




} // namespace 

#endif