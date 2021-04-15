#include "gtest/gtest.h"
#include "../../include/rand_num_gen.h"
#include <unordered_set>

namespace statanaly {

TEST(rng_unix, get_random_numbers) {
    // Test whether the generated random numbers are unique.
    std::size_t N = 100;
    rng_unix mygenerator;
    auto rns = mygenerator.getRN<std::uint64_t>(N);
    std::unordered_set<std::uint64_t> bag;
    for (std::size_t i=0; i<N; i++) {
        EXPECT_FALSE( bag.contains(rns[i]) );
        bag.insert(rns[i]);
    }
}

}