#ifndef STATANALY_D_CONTAINER_H_
#define STATANALY_D_CONTAINER_H_

#include <unordered_map>
#include <vector>
#include "density/probDensFunc.h"


namespace statanaly {

/* A container to store a collection of various distributions. */

class dCtr {
    using weightType = double;
    
    // Stores the list of named distribution functions that are the contents.
    // The second weight value is the weight distribution.
    std::unordered_map<probDensFunc*, std::pair<weightType,weightType>> ingreds;

public:

    dCtr() = default;
    ~dCtr() {
        for (auto& [d,ws] : ingreds) {delete d;}
    };

    // Copy constructor: deep-copy, do the same as clone().
    dCtr(const dCtr& o) {
        // clone the named distribution.
        for (const auto& [d,w] : o.ingreds) {
            ingreds.emplace( d->clone(), w );
        }
    }

    // Copy assignment: deep-copy
    dCtr& operator = (const dCtr& o) {
        // clone the named distribution.
        for (const auto& [d,w] : o.ingreds) {
            ingreds.emplace( d->clone(), w );
        }
        return *this;
    };
    
    // Move constructor
    dCtr(dCtr&& o) {
        // Move the named distribution.
        ingreds = std::move(o.ingreds);
    }

    // Move assignment
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

    // Rescale weight distribution after named distribution insertion and deletion.
    void rescale() {
        weightType sum = 0;
        for (const auto& [d,ws] : ingreds) {sum += ws.first;}
        for (auto& [d,ws] : ingreds) {ws.second = (ws.first) / sum;}
    }

    template<typename F, typename W>
    requires std::is_arithmetic_v<W>
    void insert(F&& distr, W weight) {
        // Make a deep-copy
        auto tmp = std::make_pair<weightType,weightType>(static_cast<weightType>(weight), 0);
        ingreds.emplace( distr.clone(), tmp );

        // Rescale the weights so that they sum up to one.
        rescale();
    }
    
    template<typename F>
    void insert(F&& distr) {
        // Default weight to 1.
        insert( std::forward<F>(distr), 1);   // redirect
    }

    // Find a distribution that match distr's hash (ie, type and parameters).
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

    inline std::size_t hash() const noexcept {
        // Container is different when the contents are different.
        // Insertion order of the ingredients does NOT matter.

        std::vector<std::size_t> hashes;
        for (const auto& [d,w] : ingreds) {
            hashes.push_back(d->hash());
        }
        std::sort(hashes.begin(), hashes.end());

        std::size_t seed = 0;
        combine_hash(seed, char(dFuncID::COUNT));   // id is the last element.
        for (const auto h : hashes) {combine_hash(seed, h);}

        return seed;
    }

    std::vector<std::vector<const probDensFunc*>> mapping() const {
        /* map[0] is all of the base distributions;
         * map[1] is all of the normal distributions;
         * map[2] is all of the std uniform distributions;
         * ... */
        std::vector<std::vector<const probDensFunc*>> map(
            static_cast<std::size_t>(dFuncID::COUNT), 
            std::vector<const probDensFunc*>{});

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

    // Print detail about this container.
    // Print order might be different from insertion order.
    friend std::ostream& operator << (std::ostream&, const dCtr&);
    void print(std::ostream& output) const {
        // first, get a mapping.
        std::vector<std::vector<const probDensFunc*>> map = mapping();

        output << "Container content : \n";
        for (const std::vector<const probDensFunc*>& ds : map) {
            if (ds.size()==0) continue;
            output << "  " << ds.size() << " following type of distribution ...\n";
            for (const probDensFunc* d : ds) {
                const auto it = ingreds.find(const_cast<probDensFunc*>(d));
                output << "    - " << *d << " @ weight = " << it->second.second << "\n";
            }
        }
        
    }
};

std::ostream& operator << (std::ostream&, const dCtr&);

}   // namespace


template<>
class std::hash<statanaly::dCtr> {
public:
    std::size_t operator() (const statanaly::dCtr& d) const {
        return d.hash();
    }
};



#endif