#ifndef AARAND_AARAND_HPP
#define AARAND_AARAND_HPP

#include <cmath>
#include <limits>
#include <stdexcept>

/**
 * @file aarand.hpp
 *
 * @brief Collection of random distribution functions.
 */

namespace aarand {

/**
 * @tparam T Floating point type.
 * @tparam Engine A random number generator class with `operator()`, `min()` (static) and `max()` (static) methods,
 * where the `result_type` is an unsigned integer value.
 *
 * @param eng Instance of an RNG class like `std::mt19937_64`.
 *
 * @return Draw from a standard uniform distribution.
 */
template<typename T = double, class Engine>
T standard_uniform(Engine& eng) {
    typedef typename Engine::result_type R;
    static_assert(std::numeric_limits<R>::is_integer);

    // Can't be bothered to figure out whether the range fits into 'R' for signed values.
    // So instead, we just require unsigned integers, where the range will always fit.
    static_assert(!std::numeric_limits<R>::is_signed); 

    // Stolen from Boost, see https://www.boost.org/doc/libs/1_67_0/boost/random/uniform_01.hpp
    // The +1 probably doesn't matter for 64-bit generators, but is helpful for engines with 
    // fewer output bits, to reduce the (small) probability of sampling 1's.
    constexpr double factor = 1.0 / (static_cast<T>(Engine::max() - Engine::min()) + 1.0);

    // Note that it still might be possible to get a result = 1, depending on
    // the numerical precision used to compute the product; hence the loop.
    double result;
    do {
        result = static_cast<T>(eng() - Engine::min()) * factor;
    } while (result == 1.0);
    return result;
}

/**
 * @tparam T Floating point type.
 * @tparam Engine A random number generator class with `operator()`, `min()` (static) and `max()` (static) methods,
 * where the `result_type` is an unsigned integer value.
 *
 * @param eng Instance of an RNG class like `std::mt19937_64`.
 *
 * @return A pair of independent draws from a standard normal distribution with mean 0 and variance 1.
 */
template<typename T = double, class Engine>
std::pair<T, T> standard_normal(Engine& eng) {
    constexpr double pi = 3.14159265358979323846;

    // Box-Muller gives us two random values at a time.
    double constant = std::sqrt(-2 * std::log(standard_uniform<T>(eng)));
    double angle = 2 * pi * standard_uniform<T>(eng);
    return std::make_pair(constant * std::sin(angle), constant * std::cos(angle));
}

/**
 * @tparam T Floating point type.
 * @tparam Engine A random number generator class with `operator()`, `min()` (static) and `max()` (static) methods,
 * where the `result_type` is an unsigned integer value.
 *
 * @param eng Instance of an RNG class like `std::mt19937_64`.
 *
 * @return Draw from a standard exponential distribution.
 */
template<typename T = double, class Engine>
T standard_exponential(Engine& eng) {
    T val;
    do {
        val = standard_uniform<T>(eng);
    } while (val == 0);
    return -std::log(val);
}

/**
 * @tparam T Integer type.
 * @tparam Engine A random number generator class with `operator()`, `min()` (static) and `max()` (static) methods,
 * where the `result_type` is an unsigned integer value.
 *
 * @param eng Instance of an RNG class like `std::mt19937_64`.
 * @param bound Positive integer specifying the upper bound of the discrete distribution.
 *
 * @return Draw from a discrete uniform distribution in `[0, bound)`.
 */
template<typename T = int, class Engine>
T discrete_uniform(Engine& eng, T bound) {
    typedef typename Engine::result_type R;
    static_assert(std::numeric_limits<R>::is_integer);
    static_assert(!std::numeric_limits<R>::is_signed); // don't want to figure out how to store the range.

    constexpr R range = Engine::max() - Engine::min();
    if (bound > range) {
        throw std::runtime_error("'bound' should be less than the RNG range");
    }

    static_assert(std::numeric_limits<T>::is_integer);
    if (bound <= 0) {
        throw std::runtime_error("'bound' should be a positive integer");
    }

    // The limit is necessary to provide uniformity in the presence of the
    // modulus. The idea is to re-sample if we get a draw above the limit.
    // Technically this can have problems as bound approaches range, in which
    // case we might end up discarding a lot of the sample space... but this
    // is unlikely to happen in practice, so whatever. Note that the +1 is
    // necessary because range is inclusive but bound is not.
    const R limit = range - (range % bound + 1);
    
    // In addition, we don't have to deal with the crap about combining draws
    // to get enough entropy, which is 90% of the Boost implementation.
    R draw;
    do {
        draw = eng() - Engine::min();
    } while (draw > limit);

    return draw % bound;
}

/**
 * @tparam In Random-access iterator or pointer.
 * @tparam Engine A random number generator class with `operator()`, `min()` (static) and `max()` (static) methods,
 * where the `result_type` is an unsigned integer value.
 *
 * @param values Iterator or pointer to an array of values to shuffle.
 * @param n Number of values in the array pointed to by `values`.
 * @param eng Instance of an RNG class like `std::mt19937_64`.
 *
 * @return Contents of `values` are randomly permuted in place using the Fisher-Yates algorithm.
 */
template<class In, class Engine>
void shuffle(In values, size_t n, Engine& eng) {
    if (n) {
        using std::swap;  
        for (size_t i = 0; i < n - 1; ++i) {
            auto chosen = discrete_uniform(eng, n - i);
            swap(*(values + i), *(values + i + chosen));
        }
    }
    return;
}

/**
 * @tparam In Random-access iterator or pointer for the inputs.
 * @tparam Out Random-access iterator or pointer for the outputs.
 * @tparam Engine A random number generator class with `operator()`, `min()` (static) and `max()` (static) methods,
 * where the `result_type` is an unsigned integer value.
 *
 * @param values Iterator or pointer to an array of values to sample from.
 * @param n Number of values in the array pointed to by `values`.
 * @param s Number of values to sample.
 * @param output Iterator or pointer to an array of length `s`, to store the sampled values. 
 * @param eng Instance of an RNG class like `std::mt19937_64`.
 *
 * @return `output` is filled with `s` sampled values from `values`.
 *
 * If `s > n`, `values` is copied into the first `n` elements of `output` and the remaining values of `output` are undefined.
 */
template<class In, class Out, class Engine>
void sample(In values, size_t n, size_t s, Out output, Engine& eng) {
    for (size_t i = 0; i < n && s; ++i, ++values) {
        const double threshold = static_cast<double>(s)/(n - i);
        if (threshold >= 1 || standard_uniform(eng) <= threshold) {
            *output = *values;
            ++output;
            --s;
        }
    }
}

/**
 * @tparam Out Random-access iterator or pointer for the outputs.
 * @tparam Engine A random number generator class with `operator()`, `min()` (static) and `max()` (static) methods,
 * where the `result_type` is an unsigned integer value.
 *
 * @param bound Upper bound of the indices to sample from.
 * @param s Number of values to sample.
 * @param output Iterator or pointer to an array of length `s`, to store the sampled values. 
 * @param eng Instance of an RNG class like `std::mt19937_64`.
 *
 * @return `output` is filled with `s` sampled values from the sequence of integers in `{0, 1, ..., bound - 1}`.
 *
 * If `s > bound`, the first `n` elements of `output` will contain the sequence of integers from `0` to `bound - 1`.
 * The remaining values of `output` are undefined.
 */
template<class Out, class Engine>
void sample(size_t bound, size_t s, Out output, Engine& eng) {
    for (size_t i = 0; i < bound && s; ++i) {
        const double threshold = static_cast<double>(s)/(bound - i);
        if (threshold >= 1 || standard_uniform(eng) <= threshold) {
            *output = i;
            ++output;
            --s;
        }
    }
}

}

#endif
