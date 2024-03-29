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
#ifndef STATANALY_UNIFORM_H_
#define STATANALY_UNIFORM_H_

#include "probDistr.h"

namespace statanaly {

/**
 * @brief Standard Uniform Distribution.
 * 
 * No class parameters.
 */

class disStdUniform : public probDistr {
public:
    double pdf(const double x=0) const override {
        if (0>x || 1<x) {return 0;}
        return 1;
    }
    
    double cdf(const double x=0) const override {
        if (0>x) return 0;
        if (1<x) return 1;
        return x;
    }

    double mean() const override {
        return 0.5;
    }

    double stddev() const override {
        const double r = sqrt(1./12);
        return r;
    }

    double variance() const override {
        constexpr double r = 1./12;
        return r;
    }

    double skewness() const override {
        return 0;
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        return seed;
    }

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disStdUniform>(static_cast<disStdUniform const&>(*this));
    };

    disStdUniform* clone() const override {
        return new disStdUniform(*this);
    }

    void print(std::ostream& output) const override {
        output << "Std Uniform distribution -- a = 0  b = 1";
    }

    constexpr bool isEqual_tol(const probDistr& o, const double tol=0) const override {
        return true;
    }
    constexpr bool isEqual_ulp(const probDistr& o, const unsigned ulp=0) const override {
        return true;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::STD_UNIFORM_DISTR;
};


/**
 * @brief Uniform Distribution.
 * 
 * @param a Lower limit.
 * @param b Upper limit.
 */

class disUniform : public probDistr {
private:

    double a;
    double b;

public:
    template<class T>
    requires std::is_arithmetic_v<T>
    disUniform(const T lower, const T upper) {
        a = lower;
        b = upper;
        if(a == b) {
            throw std::runtime_error("Uniform distribution must have valid boundary.");
        }
    }
    disUniform() = delete;
    ~disUniform() = default;

    constexpr double pdf(const double x) const override  {
        if (a>x || b<x) {return 0;}
        return 1/(b-a);
    }

    constexpr double cdf(const double x) const override {
        if (a>x) return 0;
        if (b<x) return 1;
        return (x-a)/(b-a);
    }

    constexpr double mean() const override {
        return 0.5*(a+b);
    }

    constexpr double stddev() const override {
        return (b-a)/sqrt(12.);
    }

    constexpr double variance() const override {
        return (b-a)*(b-a)/12.;
    }

    constexpr double skewness() const override {
        return 0;
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, a);
        combine_hash(seed, b);
        return seed;
    }

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disUniform>(static_cast<disUniform const&>(*this));
    };

    disUniform* clone() const override {
        return new disUniform(*this);
    }

    void print(std::ostream& output) const override {
        output << "Uniform distribution -- a = " << a << "  b = " << b;
    }

    bool isEqual_tol(const probDistr& o, const double tol) const override {
        const disUniform& oo = dynamic_cast<const disUniform&>(o);
        bool r = true;
        r &= isEqual_fl_tol(a, oo.a, tol);
        r &= isEqual_fl_tol(b, oo.b, tol);
        return r;
    }

    bool isEqual_ulp(const probDistr& o, const unsigned ulp) const override {
        const disUniform& oo = dynamic_cast<const disUniform&>(o);
        bool r = true;
        r &= isEqual_fl_ulp(a, oo.a, ulp);
        r &= isEqual_fl_ulp(b, oo.b, ulp);
        return r;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::UNIFORM_DISTR;
};
}   // namespace statanaly


/**
 * @brief STL hasher overload.
 * 
 * @tparam Standard Normal distribution.
 */

template<>
class std::hash<statanaly::disStdUniform> {
public:
    std::size_t operator() (const statanaly::disStdUniform& d) const {
        return d.hash();
    }
};

/**
 * @brief STL hasher overload.
 * 
 * @tparam Normal distribution.
 */

template<>
class std::hash<statanaly::disUniform> {
public:
    std::size_t operator() (const statanaly::disUniform& d) const {
        return d.hash();
    }
};


#endif 