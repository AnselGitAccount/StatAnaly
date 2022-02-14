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
#ifndef STATANALY_HASHER_H_
#define STATANALY_HASHER_H_

#include <tuple>

namespace statanaly {

template<class T>
inline void combine_hash(std::size_t& seed, const T& val) noexcept {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

} // namespace statanaly


/**
 * @brief STL hasher specialized for std::tuple.
 * 
 * std::tuple can have any number of arguments.
 * https://www.variadic.xyz/2018/01/15/hashing-stdpair-and-stdtuple/
 * 
 * @tparam TupleArgs 
 */
template<class... TupleArgs>
class std::hash<std::tuple<TupleArgs...>> {

    template<size_t Idx, typename... TupleTypes>
    inline typename std::enable_if_t<Idx == sizeof...(TupleTypes), void>
    hash_tup(size_t& seed, const std::tuple<TupleTypes...>& tup) const {}

    template<size_t Idx, typename... TupleTypes>
    inline typename std::enable_if_t<Idx <  sizeof...(TupleTypes), void>
    hash_tup(size_t& seed, const std::tuple<TupleTypes...>& tup) const {
        statanaly::combine_hash(seed, std::get<Idx>(tup));
        hash_tup<Idx + 1>(seed, tup);
    }

public:
    auto operator()(const std::tuple<TupleArgs...>& tupleValue) const {
        // need Idx to extract element from std::tuple.
        size_t seed = 0;
        hash_tup<0>(seed, tupleValue);
        return seed;
    }
};



#endif