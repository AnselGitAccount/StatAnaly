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
#ifndef STATANALY_DIS_CHISQ_H_
#define STATANALY_DIS_CHISQ_H_

#include "probDistr.h"

namespace statanaly {

/**
 * @brief Central Chi Square distribution
 * 
 * @param k Degree of freedom
 * 
 */

class disChiSq : public probDistr {
private:
    unsigned k;

public:
    template<class T> 
    requires std::is_integral_v<T>
    disChiSq(const T dof){
        k = dof;
    }
    ~disChiSq() = default;

    double pdf(const double x) const override {
        const double kd2 = k*0.5;
        const double lgf = logGamma(kd2);
        const double r   = exp((kd2-1.0)*log(x) - x*0.5 - kd2*M_LN2 - lgf);
        return r;
    }

    double cdf(const double x) const override {
        const double lgf = regLowerGamma(k/2.0, x/2.0);
        return lgf;
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

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disChiSq>(static_cast<disChiSq const&>(*this));
    };

    disChiSq* clone() const override {
        return new disChiSq(*this);
    }

    void print(std::ostream& output) const override {
        output << "Central Chi Square distribution -- k = " << k;
    }

    bool isEqual_tol(const probDistr& o, const double tol=0) const override {
        const disChiSq& oo = dynamic_cast<const disChiSq&>(o);
        return k==oo.k;
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp=0) const override {
        const disChiSq& oo = dynamic_cast<const disChiSq&>(o);
        return k==oo.k;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::CHISQ_DISTR;

    unsigned p_dof() const noexcept{
        return k;
    }
};

}   // namespace statanaly


/**
 * @brief STL hasher overload.
 * 
 * @tparam Central Chi Square distribution
 */

template<>
class std::hash<statanaly::disChiSq> {
public:
    std::size_t operator() (const statanaly::disChiSq& d) const {
        return d.hash();
    }
};

#endif