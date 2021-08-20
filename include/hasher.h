#ifndef STATANALY_HASHER_H_
#define STATANALY_HASHER_H_

#include <tuple>

namespace statanaly {

template<class T>
inline void combine_hash(std::size_t& seed, const T& val) noexcept {
    seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

} // namespace statanaly


// Specialized std::hash for std::tuple.
// The std::tuple can have any number of arguments.
// https://www.variadic.xyz/2018/01/15/hashing-stdpair-and-stdtuple/
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