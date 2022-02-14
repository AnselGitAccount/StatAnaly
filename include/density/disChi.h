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
#ifndef STATANALY_DIS_CHI_H_
#define STATANALY_DIS_CHI_H_

#include "probDistr.h"

namespace statanaly {

/**
 * @brief Central Chi Distribution.
 * 
 * @param k Degree of freedom
 */

class disChi : public probDistr {
private:

    unsigned k;

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

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disChi>(static_cast<disChi const&>(*this));
    };

    disChi* clone() const override {
        return new disChi(*this);
    }

    void print(std::ostream& output) const override {
        output << "Central Chi distribution -- k = " << k;
    }

    bool isEqual_tol(const probDistr& o, const double tol=0) const override {
        const disChi& oo = dynamic_cast<const disChi&>(o);
        bool r = true;
        r &= k==oo.k;
        return r;
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp=0) const override {
        const disChi& oo = dynamic_cast<const disChi&>(o);
        bool r = true;
        r &= k==oo.k;
        return r;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::CHI_DISTR;

    unsigned p_dof() const noexcept{
        return k;
    }
};

}   // namespace statanaly


/**
 * @brief STL hasher overload.
 * 
 * @tparam Central Chi distribution
 */

template<>
class std::hash<statanaly::disChi> {
public:
    std::size_t operator() (const statanaly::disChi& d) const {
        return d.hash();
    }
};

#endif