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
#ifndef STATANALY_DIS_GAMMA_H_
#define STATANALY_DIS_GAMMA_H_

#include "probDistr.h"


namespace statanaly {

/**
 * @brief Gamma Distribution
 * 
 * @param theta Scale
 * @param alpha Shape
 */

class disGamma : public probDistr {
private:
    double theta;
    double alpha;

public:
    template<class T>
    requires std::is_arithmetic_v<T>
    disGamma(const T scale, const T shape) {
        theta = scale;
        alpha = shape;
    }
    disGamma() = delete;
    
    double pdf (const double x) const override {
        return pow(x,alpha-1)/pow(theta,alpha) * exp(-x/theta) / std::tgamma(alpha);
    }

    double cdf (const double x) const override {
        return regLowerGamma(alpha, x/theta);
    }

    double mean() const override {
        return alpha*theta;
    }

    double stddev() const override {
        return std::sqrt(variance());
    }

    double variance() const override {
        return alpha*theta*theta;
    }

    double skewness() const override {
        return 2./sqrt(alpha);
    }

    auto pscale() const noexcept {return theta;}
    auto pshape() const noexcept {return alpha;}

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, theta);
        combine_hash(seed, alpha);
        return seed;
    }

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disGamma>(static_cast<disGamma const&>(*this));
    };

    disGamma* clone() const override {
        return new disGamma(*this);
    }

    void print(std::ostream& output) const override {
        output << "Gamma distribution -- theta = " << theta << " alpha = " << alpha;
    }

    bool isEqual_tol(const probDistr& o, const double tol=0) const override {
        const disGamma& oo = dynamic_cast<const disGamma&>(o);
        bool r = true;
        r &= isEqual_fl_tol(theta, oo.theta, tol);
        r &= isEqual_fl_tol(alpha, oo.alpha, tol);
        return r;
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp=0) const override {
        const disGamma& oo = dynamic_cast<const disGamma&>(o);
        bool r = true;
        r &= isEqual_fl_ulp(theta, oo.theta, ulp);
        r &= isEqual_fl_ulp(alpha, oo.alpha, ulp);
        return r;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::GAMMA_DISTR;
};

}   // namespace statanaly


/**
 * @brief STL hasher overload
 * 
 * @tparam Gamma distribution
 */

template<>
class std::hash<statanaly::disGamma> {
public:
    std::size_t operator() (const statanaly::disGamma& d) const {
        return d.hash();
    }
};

#endif