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
#ifndef STATANALY_DIS_CAUCHY_H_
#define STATANALY_DIS_CAUCHY_H_

#include "probDistr.h"


namespace statanaly {

/**
 * @brief Cauchy distribution.
 * 
 * @param t location
 * @param s scale
 */

class disCauchy : public probDistr {
private:
    double t;
    double s;

public:
    template<class T>
    requires std::is_arithmetic_v<T>
    disCauchy(const T loc, const T scale) {
        t = loc;
        s = scale;
    }
    disCauchy() = delete;

    double pdf(const double x) const override {
        const double scaledx = (x-t)/s;
        const double tmp = s*M_PIf64*(1.0+scaledx*scaledx);
        return 1.0/tmp;
    }

    double cdf(const double x) const override {
        const double scaledx = (x-t)/s;
        return 0.5 + atan(scaledx)*M_1_PIf64;
    }

    double mean() const override {
        throw std::runtime_error("Mean of Cauchy distribution is undefined.");
        return 0;
    }

    double stddev() const override {
        throw std::runtime_error("Standard Deviation of Cauchy distribution is undefined.");
        return 0;
    }

    double variance() const override {
        throw std::runtime_error("Variance of Cauchy distribution is undefined.");
        return 0;
    }

    double skewness() const override {
        throw std::runtime_error("Skewness of Cauchy distribution is undefined.");
        return 0;        
    }
    
    auto ploc() const noexcept {return t;}
    auto pscale() const noexcept {return s;}

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, t);
        combine_hash(seed, s);
        return seed;
    }

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disCauchy>(static_cast<disCauchy const&>(*this));
    }

    disCauchy* clone() const override {
        return new disCauchy(*this);
    }

    void print(std::ostream& output) const override {
        output << "Cauchy distribution -- s = " << s << "  t = " << t;
    }

    bool isEqual_tol(const probDistr& o, const double tol) const override {
        const disCauchy& oo = dynamic_cast<const disCauchy&>(o);
        bool r = true;
        r &= isEqual_fl_tol(t, oo.t, tol);
        r &= isEqual_fl_tol(s, oo.s, tol);
        return r;
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp) const override {
        const disCauchy& oo = dynamic_cast<const disCauchy&>(o);
        bool r = true;
        r &= isEqual_fl_ulp(t, oo.t, ulp);
        r &= isEqual_fl_ulp(s, oo.s, ulp);
        return r;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::CAUCHY_DISTR;
};

}   // namespace statanaly


/**
 * @brief STL hasher overload.
 * 
 * @tparam Cauchy distribution
 */

template<>
class std::hash<statanaly::disCauchy> {
public:
    std::size_t operator() (const statanaly::disCauchy& d) const {
        return d.hash();
    }
};

#endif