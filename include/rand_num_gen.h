/* Random Number Generator
 * 
 * For Linux systems, get random numbers from /dev/urandom.
 */

#ifndef STATANALY_RAND_NUM_GEN_H
#define STATANALY_RAND_NUM_GEN_H

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <fcntl.h>  // O_RDONLY
#include <cassert>
#include <unistd.h>   // read
#include <cstring>
#include <memory>
#include <random>
#include <algorithm>
#include <iterator>

namespace statanaly {


/**
 * @brief Pseudo RNG
 * 
 * Draw from /dev/urandom
 * More "direct" than libstd.
 * Non-blocking.
 * For Unix-like OS only.
 */
class rng_unix {
    int fd = 0;

    void openUrandom() {
        fd = open("/dev/urandom", O_RDONLY);
        assert(fd>=0);
    }

public:
    rng_unix() {
        openUrandom();
    }
    ~rng_unix() {
        if (fd>=0) close(fd);
    }

    // Disable copy operations.
    rng_unix(const rng_unix&) = delete;
    rng_unix& operator = (const rng_unix&) = delete;
    
    // Move operations only
    rng_unix(rng_unix&& other) {
        fd = std::exchange(other.fd, -1);
    }
    rng_unix& operator = (rng_unix&& other) {
        fd = std::exchange(other.fd, -1);
        return *this;
    }

    inline int getFD() noexcept{
        return fd;
    }

    template<class T>
    auto getRN() {
        T buf;
        read(fd, &buf, sizeof(T));
        return buf;
    }

    template<class T>
    auto getRN(uint num) {
        std::unique_ptr<T[]> buf = std::make_unique<T[]>(num);
        read(fd, buf.get(), sizeof(T)*num);
        return buf;
    }
};



/**
 * @brief Seed STL mt19937 RN generators.
 * 
 * Seeding Mersenne Twister with enough entropy for its state-size, i.e., 624x4 bytes
 * 
 * @param eng RN generator engine.
 */
void setProperSeed(std::mt19937& eng) {
    // use "default" token because it is portable.
    std::random_device rd;
    std::random_device::result_type rnd_nums[std::mt19937::state_size];
    std::generate(std::begin(rnd_nums), std::end(rnd_nums), std::ref(rd));
    std::seed_seq seeds(std::begin(rnd_nums), std::end(rnd_nums));
    
    eng.seed(seeds);
};


} // namespace 

#endif