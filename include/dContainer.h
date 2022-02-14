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
#ifndef STATANALY_D_CONTAINER_H_
#define STATANALY_D_CONTAINER_H_

#include <unordered_map>
#include <vector>
#include "density/probDistr.h"


namespace statanaly {

/**
 * @brief Container class storing a collection of distributions.
 * 
 * This is the container for Mixture distribution. It stores a list of named distributions and their weight.
 * When a distribution and its weight are inserted, the weights of existing distributions in the container are rebalanced.
 * For example:
 * 
 * Inserting (Normal, weight=1) to 
 * 
 *      {(Uniform, weight=0.5) (Normal, weight=0.5)}.
 * 
 * After insertion and weight rebalancing:
 * 
 *      {(Uniform, weight=1/3) (Normal, weight=1/3) (Normal, weight=1/3)}
 * 
 * @param ingreds A collection of distribution and their weights.
 * @see disMixture
 */
class dCtr {
    using weightType = double;
    
    /**
     * @brief Internal storage for the collection.
     * 
     * The first weight value is the unnormalized weight.
     * The second weight value is the normalized weight.
     */
    std::unordered_map<probDistr*, std::pair<weightType,weightType>> ingreds;

public:

    dCtr() = default;
    ~dCtr() {
        for (auto& [d,ws] : ingreds) {delete d;}
    };

    /** Copy constructor: deep-copy, do the same as clone(). */
    dCtr(const dCtr& o) {
        // clone the named distribution.
        for (const auto& [d,w] : o.ingreds) {
            ingreds.emplace( d->clone(), w );
        }
    }

    /** Copy assignment: deep-copy */
    dCtr& operator = (const dCtr& o) {
        // clone the named distribution.
        for (const auto& [d,w] : o.ingreds) {
            ingreds.emplace( d->clone(), w );
        }
        return *this;
    };
    
    /** Move constructor */
    dCtr(dCtr&& o) {
        // Move the named distribution.
        ingreds = std::move(o.ingreds);
    }

    /** Move assignment */
    dCtr& operator = (dCtr&& o) {
        // Move the named distribution.
        ingreds = std::move(o.ingreds);
        return *this;
    }

    std::unique_ptr<dCtr> cloneUnique() const {
        return std::make_unique<dCtr>(static_cast<dCtr const&>(*this));
    };

    dCtr* clone() const {
        return new dCtr(*this);
    };

    /** Rescale weight distribution after named distribution insertion and deletion. */
    void rescale() {
        weightType sum = 0;
        for (const auto& [d,ws] : ingreds) {sum += ws.first;}
        for (auto& [d,ws] : ingreds) {ws.second = (ws.first) / sum;}
    }

    /** Insert a distribution and its weight */
    template<typename F, typename W>
    requires std::is_arithmetic_v<W>
    void insert(F&& distr, W weight) {
        // Make a deep-copy
        auto tmp = std::make_pair<weightType,weightType>(static_cast<weightType>(weight), 0);
        ingreds.emplace( distr.clone(), tmp );

        // Rescale the weights so that they sum up to one.
        rescale();
    }
    
    /** Insert a distribution and set its weight to one. */
    template<typename F>
    void insert(F&& distr) {
        // Default weight to 1.
        insert( std::forward<F>(distr), 1);
    }

    /** Find a distribution that match another distribution's hash (ie, type and parameters). */
    template<typename F>
    inline auto find(F&& distr) const {
        for (auto it = ingreds.begin(); it != ingreds.end(); it++) {
            if (it->first->hash() == distr.hash()) {return it;}
        }
        return ingreds.end();
    }

    inline const auto& get() const {return ingreds;}
    inline const auto end() const {return ingreds.end();}
    inline void clear() {ingreds.clear();}

    /**
     * @brief Computing the hash of this class instance.
     * 
     * Each containing distribution has a hash based on parameters of that distribution and distribution type.
     * Two containers are equal when the hashes of containing distribution are equal.
     * Note that insertion order of the ingredients does NOT matter.
     * 
     * @return std::size_t 
     */
    inline std::size_t hash() const noexcept {
        std::vector<std::size_t> hashes;
        for (const auto& [d,w] : ingreds) {
            // Contents -- distribution parameters and distribution weight.
            std::size_t h = 0;
            combine_hash(h, d->hash());
            combine_hash(h, w.first);
            hashes.push_back(h);
        }
        std::sort(hashes.begin(), hashes.end());

        std::size_t seed = 0;
        combine_hash(seed, char(dFuncID::COUNT));   // id is the last element.
        for (const auto h : hashes) {combine_hash(seed, h);}

        return seed;
    }

    /**
     * @brief Catelog the containing distributions by their type.
     * 
     * Any insertion order of the same set of probability distributions should return the same container.
     * So their should return the same catelog.
     * 
     * map[0] is all of the base distributions;
     * map[1] is all of the normal distributions;
     * map[2] is all of the std uniform distributions;
     * ...
     * 
     * @return std::vector<std::vector<const probDistr*>> 
     */
    std::vector<std::vector<const probDistr*>> catelog() const {
        std::vector<std::vector<const probDistr*>> map(
            static_cast<std::size_t>(dFuncID::COUNT), 
            std::vector<const probDistr*>{});

        for (const auto& [d,w] : ingreds) {
            map[ static_cast<std::size_t>(d->getID()) ].push_back(d);
        }

        // Sort distributions based on their hash values, for each type of distribution.
        for (auto& ds : map) {
            std::sort(ds.begin(), ds.end(), [](const auto* a, const auto* b)
                { return a->hash() < b->hash(); });
        }

        return map;
    }

    /**
     * @brief Print the detail of the container to STDOUT in a human-friendly format.
     * 
     * Print distributions in the order as they appear in the catelog. 
     * So it might be different from insertion order.
     * 
     * @return std::ostream& 
     */
    friend std::ostream& operator << (std::ostream&, const dCtr&);
    void print(std::ostream& output) const {
        std::vector<std::vector<const probDistr*>> map = catelog();

        output << "Container content : \n";
        for (const std::vector<const probDistr*>& ds : map) {
            if (ds.size()==0) continue;
            output << "  " << ds.size() << " following type of distribution ...\n";
            for (const probDistr* d : ds) {
                const auto it = ingreds.find(const_cast<probDistr*>(d));
                output << "    - " << *d << " @ weight = " << it->second.second << "\n";
            }
        }
        
    }

    // Every distribution has to be less than tol.
    bool isEqual_tol(const dCtr& o, const double tol) const {
        if (ingreds.size() != o.ingreds.size())
            return false;

        for (auto it = ingreds.begin(), oit = o.ingreds.begin(); 
            it != ingreds.end(), oit != o.ingreds.end(); 
            it++, oit++) {
            const bool is = it->first->isEqual_tol(*oit->first, tol);
            if (!is) return false;
        }

        return true;
    }

    // Every distribution has to be less than N ulp.
    bool isEqual_ulp(const dCtr& o, const unsigned ulp) const {
        if (ingreds.size() != o.ingreds.size())
            return false;

        for (auto it = ingreds.begin(), oit = o.ingreds.begin(); 
            it != ingreds.end(), oit != o.ingreds.end(); 
            it++, oit++) {
            const bool is = it->first->isEqual_ulp(*oit->first, ulp);
            if (!is) return false;
        }
        
        return true;
    }

};

std::ostream& operator << (std::ostream&, const dCtr&);

}   // namespace

/**
 * @brief STL hasher overload
 * 
 * @tparam Distribution container.
 */
template<>
class std::hash<statanaly::dCtr> {
public:
    std::size_t operator() (const statanaly::dCtr& d) const {
        return d.hash();
    }
};



#endif