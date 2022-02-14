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
#ifndef STATANALY_DIS_NONCENTRAL_CHI_SQUARED_H_
#define STATANALY_DIS_NONCENTRAL_CHI_SQUARED_H_

#include "probDistr.h"


namespace statanaly {


/**
 * @brief Non-central Chi Squared Distribution
 * 
 * @param k Degree of freedom
 * @param lambda Distance
 */

class disNcChiSq : public probDistr {
private:

    unsigned k;
    double lambda;

public:
    template<class T, class U>
    requires std::is_integral_v<T> && std::is_arithmetic_v<U>
    disNcChiSq(const T dof, const U distance) {
        k = dof;
        lambda = distance;
    }
    disNcChiSq() = delete;
    ~disNcChiSq() = default;

    double pdf(const double x) const override {
        const double kp = 0.5*k - 1.;
        const double t = 0.5 * pow(x/lambda, 0.5*kp) * exp(-0.5*(x+lambda));
        return t * std::cyl_bessel_i(kp, std::sqrt(lambda*x));
    }

    double cdf(const double x) const override {
        return 1. - marcumQ(0.5*k, std::sqrt(lambda), std::sqrt(x));
    }

    double mean() const override {
        return k+lambda;
    }

    double stddev() const override {
        return std::sqrt(variance());
    }

    double variance() const override {
        return 2*(k+2*lambda);
    }

    double skewness() const override {
        return (k+3*lambda) / std::sqrt(0.125*(k+2*lambda)*(k+2*lambda)*(k+2*lambda));
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, k);
        combine_hash(seed, lambda);
        return seed;
    } 

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disNcChiSq>(static_cast<disNcChiSq const&>(*this));
    };

    disNcChiSq* clone() const override {
        return new disNcChiSq(*this);
    }

    void print(std::ostream& output) const override {
        output << "Non-central Chi-Sqaured distribution -- k = " << k << " lambda = " << lambda;
    }

    bool isEqual_tol(const probDistr& o, const double tol) const override {
        const disNcChiSq& oo = dynamic_cast<const disNcChiSq&>(o);
        bool r = true;
        r &= (k == oo.k);
        r &= isEqual_fl_tol(lambda, oo.lambda, tol);
        return r;
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp) const override {
        const disNcChiSq& oo = dynamic_cast<const disNcChiSq&>(o);
        bool r = true;
        r &= (k == oo.k);
        r &= isEqual_fl_ulp(lambda, oo.lambda, ulp);
        return r;
    }

    auto p_dof() const noexcept {
        return k;
    }

    auto p_distance() const noexcept {
        return lambda;
    }
};

}   // namesapce statanaly


/**
 * @brief STL hasher overload
 * 
 * @tparam Non-Central Chi Squared distribution
 */

template<>
class std::hash<statanaly::disNcChiSq> {
public:
    std::size_t operator() (const statanaly::disNcChiSq& d) const {
        return d.hash();
    }
};


#endif