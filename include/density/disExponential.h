/*
   Copyright 2022, Ansel Blumers

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#ifndef STATANALY_DIS_EXPONENTIAL_H_
#define STATANALY_DIS_EXPONENTIAL_H_

#include "probDistr.h"

namespace statanaly {


/**
 * @brief Exponential Distribution
 * 
 * @param lambda Rate
 */

class disExponential : public probDistr {
private:
    double lambda;

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

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disExponential>(static_cast<disExponential const&>(*this));
    };

    disExponential* clone() const override {
        return new disExponential(*this);
    }

    void print(std::ostream& output) const override {
        output << "Exponential distribution -- lambda = " << lambda;
    }

    bool isEqual_tol(const probDistr& o, const double tol) const override {
        const disExponential& oo = dynamic_cast<const disExponential&>(o);
        return isEqual_fl_tol(lambda, oo.lambda, tol);
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp) const override {
        const disExponential& oo = dynamic_cast<const disExponential&>(o);
        return isEqual_fl_ulp(lambda, oo.lambda, ulp);
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::EXPONENTIAL_DISTR;
};

} // namespace statanaly


/**
 * @brief STL hasher overload
 * 
 * @tparam Exponential distribution
 */

template<>
class std::hash<statanaly::disExponential> {
public:
    std::size_t operator() (const statanaly::disExponential& d) const {
        return d.hash();
    }
};


#endif