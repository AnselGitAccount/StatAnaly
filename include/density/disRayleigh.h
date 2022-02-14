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
#ifndef STATANALY_DIS_RAYLEIGH_H_
#define STATANALY_DIS_RAYLEIGH_H_

#include "probDistr.h"
#include "disChi.h"


namespace statanaly {

/**
 * @brief Rayleigh Distribution
 * 
 * @param sigma Scale
 */

class disRayleigh : public probDistr {
private:

    double sigma;

public:

    template<class T>
    requires std::is_arithmetic_v<T>
    disRayleigh(const T scale) {
        sigma = scale;
    }
    disRayleigh() = delete;
    ~disRayleigh() = default;

    double pdf(const double x) const override {
        const double inv_ss = 1/(sigma*sigma);
        const double r = x*inv_ss * std::exp(-0.5*x*x*inv_ss);
        return r;
    }

    double cdf(const double x) const override {
        const double s = std::exp(-0.5*x*x/(sigma*sigma));
        return 1-s;
    }

    double mean() const override {
        constexpr double s = std::sqrt(M_PI/2);
        return sigma*s; 
    }

    double stddev() const override {
        constexpr double s = std::sqrt(2-M_PI_2);
        return s*sigma;
    }

    double variance() const override {
        constexpr double s = 2-M_PI_2;
        return s*sigma*sigma;
    }

    constexpr double skewness() const override {
        constexpr double s = 2.*std::sqrt(M_PI)*(M_PI-3) / std::sqrt((4-M_PI)*(4-M_PI)*(4-M_PI));
        return s;
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, sigma);
        return seed;
    }

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disRayleigh>(static_cast<disRayleigh const&>(*this));
    };

    disRayleigh* clone() const override {
        return new disRayleigh(*this);
    }

    void print(std::ostream& output) const override {
        output << "Rayleigh distribution -- sigma = " << sigma;
    }

    bool isEqual_tol(const probDistr& o, const double tol) const override {
        const disRayleigh& oo = dynamic_cast<const disRayleigh&>(o);
        return isEqual_fl_tol(sigma, oo.sigma, tol);
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp) const override {
        const disRayleigh& oo = dynamic_cast<const disRayleigh&>(o);
        return isEqual_fl_ulp(sigma, oo.sigma, ulp);
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::RAYLEIGH_DISTR;

    double p_scale() const {
        return sigma;
    }
};

}   // namespace statanaly


/**
 * @brief STL hasher overload
 * 
 * @tparam Rayleigh distribution
 */

template<>
class std::hash<statanaly::disRayleigh> {
public:
    std::size_t operator() (const statanaly::disRayleigh& d) const {
        return d.hash();
    }
};


#endif