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
#ifndef STATANALY_IRWIN_HALL_H_
#define STATANALY_IRWIN_HALL_H_

#include "probDistr.h"


namespace statanaly {


/**
 * @brief Central Chi Distribution
 * 
 * Irwin–Hall random variable is defined as the sum of a number of independnet random variables.
 * 
 * @param n Num of IDD of Uniform distributions.
 */

class disIrwinHall : public probDistr {
private:

    unsigned n;

public:
    template<class T>
    requires std::is_integral_v<T>
    disIrwinHall(const T a) {
        if (a<=0) 
            throw std::runtime_error("Irwin Hall distribution takes a positive non-zero parameter.\n");
        n = a;
    }
    disIrwinHall() = delete;
    ~disIrwinHall() = default;

    double pdf(const double x=0) const override {
        if (x<0) return 0;

        double sum = 0;
        int sign = -1;
        for (std::size_t k=0; k<=std::size_t(x); k++) {
            sign = -sign;
            sum += sign * pow(x-k,n-1) / factorial[k] / factorial[n-k];
        }
        double res = sum * n;
        return res;
    }
    
    double cdf(const double x=0) const override {
        if (x<0) return 0;

        double sum = 0;
        int sign = -1;
        for (std::size_t k=0; k<=std::size_t(x); k++) {
            sign = -sign;
            sum += sign * pow(x-k,n-1) / factorial[k] / factorial[n-k];
        }
        return sum;
    }

    double mean() const override {
        return n/2.;
    }

    double stddev() const override {
        return std::sqrt(variance());
    }

    double variance() const override {
        return n/12.;
    }

    double skewness() const override {
        return 0;
    }

    inline std::size_t hash() const noexcept {
        std::size_t seed = 0;
        combine_hash(seed, char(id));
        combine_hash(seed, n);
        return seed;
    }

    std::unique_ptr<probDistr> cloneUnique() const override {
        return std::make_unique<disIrwinHall>(static_cast<disIrwinHall const&>(*this));
    };

    disIrwinHall* clone() const override {
        return new disIrwinHall(*this);
    }

    void print(std::ostream& output) const override {
        output << "Irwin Hall distribution -- n = " << n;
    }

    bool isEqual_tol(const probDistr& o, const double tol=0) const override {
        const disIrwinHall& oo = dynamic_cast<const disIrwinHall&>(o);
        return n==oo.n;
    }
    bool isEqual_ulp(const probDistr& o, const unsigned ulp=0) const override {
        const disIrwinHall& oo = dynamic_cast<const disIrwinHall&>(o);
        return n==oo.n;
    }

    virtual dFuncID getID() const {return id;};
    const dFuncID id = dFuncID::IRWIN_HALL;
};

} // namespace statanaly


/**
 * @brief STL hasher overload
 * 
 * @tparam Irwin Hall distribution
 */

template<>
class std::hash<statanaly::disIrwinHall> {
public:
    std::size_t operator() (const statanaly::disIrwinHall& d) const {
        return d.hash();
    }
};

#endif