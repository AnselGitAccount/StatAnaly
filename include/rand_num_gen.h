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

/* Pseudo RNG -- Mersenne Twister
 * Seed-able and properly seeded,
 * Portable,
 * Platform agnostic */
class rng {
    std::random_device rd;  // use "default" token because it is portable.
    std::mt19937 eng;   // typedef to std::mersenne_twister_engine

public:
    rng() {
        setProperSeed();
    }
    rng(std::size_t x) {
        setSeed(x);
    }

    // Disable copy operations -- std::random_device prohibits copying.
    rng(const rng&) = delete;
    rng& operator = (const rng&) = delete;

    // Disable move operations. 
    rng(rng&&) = delete;
    rng& operator = (rng&&) = delete;


    void setProperSeed() {
        // Seeding Mersenne Twister with enough entropy for its state-size, i.e., 624x4 bytes
        std::random_device::result_type rnd_nums[std::mt19937::state_size];
        std::generate(std::begin(rnd_nums), std::end(rnd_nums), std::ref(rd));
        std::seed_seq seeds(std::begin(rnd_nums), std::end(rnd_nums));
        
        eng.seed(seeds);
    };

    void setSeed(std::size_t x) {
        // Simple seeding, for RN reproducibility.
        std::seed_seq seeds{x};
        eng.seed(seeds);
    }
    
};


/* Pseudo RNG that is non-blocking. 
 * For Linux kernels only. */
class rng_linux {
    int fd = 0;

    void openUrandom() {
        fd = open("/dev/urandom", O_RDONLY);
        assert(fd>=0);
    }

public:
    rng_linux() {
        openUrandom();
    }
    ~rng_linux() {
        if (fd>=0) close(fd);
    }

    // Disable copy operations.
    rng_linux(const rng_linux&) = delete;
    rng_linux& operator = (const rng_linux&) = delete;
    
    // Move operations only
    rng_linux(rng_linux&& other) {
        fd = std::exchange(other.fd, -1);
    }
    rng_linux& operator = (rng_linux&& other) {
        fd = std::exchange(other.fd, -1);
        return *this;
    }

    inline int getFD() noexcept{
        return fd;
    }

    template<class T>
    T getRN() {
        T buf;
        read(fd, &buf, sizeof(T));
        return buf;
    }

    template<class T>
    std::unique_ptr<T[]> getRN(uint num) {
        std::unique_ptr<T[]> buf = std::make_unique<T[]>(num);
        read(fd, buf.get(), sizeof(T)*num);
        return buf;
    }
};

} // namespace 


/* TODO:
 *  Test the good-ness of rng and rng_linux.
 */

#endif