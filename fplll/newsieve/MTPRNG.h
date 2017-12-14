/**
This file provides an interface for a PRNG in a possibly multi-threaded setting.
Notably, we use a master seed to initialize #threads many child RNGs that operate independently.
This is done to make sure that (at least if we don't change basis), the set of vectors that are
sampled is at least somewhat consistent among different runs.

This file also contains some basic rejection sampling routines for sampling Gaussians on Z.
*/

#ifndef MTPRNG_H
#define MTPRNG_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "Typedefs.h"

namespace GaussSieve
{

// wrapper around (a vector of) random number engines of type Engine
// This is used to unify the single and multi-threaded case
template <class Engine, bool MT, class Sseq> class MTPRNG;

/**

 These functions sample from a discrete Gaussian distribution with parameter s and center c
 on (the 1-dimensional lattice given by) the integers Z.

 We cutoff the Gaussian at s*cutoff. This means that the distribution is discrete on a subset of
 Z with output probability for x being proportional to exp(-pi(x-c)^2/s^2). Note the scaling by pi
 in the exponent.

 For reasons of numerical stability, center should not be very large in absolute value (it is
 possible to reduce to |center|<1 anyway), such that  center +/- cutoff * s does not overflow.
 Z must be an integral POD type (e.g. short, int, long, long long).
 center needs to be representable as an (exact) double.
 We do NOT support mpz_t here! The output will take the role of coefficients wrt a given basis.

 We only support double for the floating point numbers. For sieving algorithms, there is no really
 good reason for now to support different precisions, as sampling does not dominate anyway.
 Furthermore, the algorithm is not very sensitive to the quality of the samples.

 Note: if one ever wants to have higher precision, one also needs to adjust the PRNGs to actually
 output high precision. (std::... only supports double)

 engine is supposed to be a random number engine (as defined by the STL).

 The variant sample_z_gaussian_VMD takes
 s2pi = s^2 / pi and maxdeviation = cutoff * s as parameters.

 The implementation is just a simple rejection sampling.
 */

template <class Z, class Engine>
Z sample_z_gaussian(double s, double const center, Engine &engine, double const cutoff);

template <class Z, class Engine>
Z sample_z_gaussian_VMD(double const s2pi, double const center, Engine &engine,
                        double const maxdeviation);
template <class Engine> int sample_uniform(int max_val, Engine &engine);

/**
The class MTPRNG is just a wrapper around a PRNG engine to facilitate switching to multi-threaded.
Due to the fact that we support multi-threading, MTPRNG<Engine,true,.> is a wrapper around
a vector of Engines, whereas, MTPRNG<Engine,false,.> is a wrapper around Engine.
reseed seeds *all* engines. Use rnd(thread-id) to obtain the underlying Engine.
Thread-safety: Init and reseed are not thread-safe. Concurrent calls to rnd are fine.
You may concurrently call rnd and use init to increase (but not decrease) the number of threads.

Randomness path:
The global (master) seed that is input to the constructor resp. to reseed is used to create
20x32 = 640 bit per-thread-seeds for each actual Engine.

The output from theses engine(s) is then accessed using rnd(thread-number) with
0<=thread-number < num_threads.
Note that rnd(thread-number) returns a reference to the engine. To get actual random data,
use rnd(thread-number)() or feed the engine to a distribution, e.g.
std::uniform_int_distribution<int> six_sided_die(1,6)
int result = six_sided_die(rnd(thread-id));

The usage syntax is the same for the single-threaded case to ensure consistency.
(in particular, rnd takes a thread-id, which should be zero). The randomness path is the same,
so for a given master seed, rnd(0) does not depend on whether the single- or multi-threaded variant
is used.
Note that for obtaining the per-thread seeds from the master seeds, we use a fixed Mersenne twister
engine and not the engine given as template parameter.
*/

// MTPRNG for MT = true, i.e. multi-threaded case.
template <class Engine, class Sseq> class MTPRNG<Engine, true, Sseq>
{
public:
  // clang-format off
  // constructs an uninitialized MTPRNG
  // note that _seq is actually changed in seeded(_seq) !
  explicit MTPRNG(Sseq &_seq) : seeder(_seq), engines(0), num_threads(0)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing (yet-uninitialized) MT RNG engines.")
  }
  // clang-format on
  inline void reseed(Sseq &_seq);

  /**
  will make sure at least _num_threads engines are actually running, starting new ones as desired.
  Will never reseed/restart already running engines. Reducing the number of threads and increasing
  it back saves the random state (unless we reseed).
  */
  inline void init(unsigned int const _num_threads);
  inline Engine &rnd(unsigned int const thread)
  {
// clang-format off
    // Note that if init is used to decrease num_threads, then
    // rnd(thread-number) for thread-number that was was made invalid by the decrease
    // would actually work until reseed is called.
    // This is considered wrong usage of the class and not allowed.
// clang-format on
#ifdef DEBUG_SIEVE_MTPRNG_THREAD_RANGE
    assert(thread < num_threads);
#endif
    return engines[thread];
  }

private:
  // seeded with initial seq and consecutively used to seed the children PRNGs.
  std::mt19937_64 seeder;
  std::vector<Engine> engines;
  // number of initialized engines. May differ from size of the vector.
  // In particular, num_threads = 0 means uninitialized.
  unsigned int num_threads;
  Sseq seq;
  /**
     number of 32bit values to use as seeds for the underlying engine(s). Technically, we could use
     state_size if the engine provides it, but not even all default engines do.
  */
  static unsigned int constexpr seed_length = 20;
};

// singlethreaded case of MTPRNG: just wrapper around Engine
template <class Engine, class Sseq> class MTPRNG<Engine, false, Sseq>
{
public:
  explicit MTPRNG(Sseq &_seq) : engine()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing ST RNG Engine.")
    reseed(_seq);
  }
  inline void reseed(Sseq &_seq);

  // does nothing. The (ignored) argument corresponds is to ensure a consistent interface to
  // to the multithreaded case, where it is the number of threads.
  FORCE_INLINE inline void init(unsigned int const = 1)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing Single-Threaded RNG Engines.")
  }

  // Argument is number of thread. It is ignored.
  // Note: non-const constexpr (in C++14) is exactly what we need.
  CPP14CONSTEXPR inline Engine &rnd(unsigned int const = 0) { return engine; }

private:
  Engine engine;
  /**
  number of 32bit values to use as seeds for the underlying engine(s). Technically, we could
  use state_size if the engine provides it, but not even all default engines do.
    */
  static unsigned int constexpr seed_length = 20;
};  // End of MTPRNG

/**
  out-of class definitions:
*/

// clang-format off
template <class Engine, class Sseq>
inline void MTPRNG<Engine, true, Sseq>::reseed(Sseq &_seq)
// clang-format on
{
  seeder.seed(_seq);
  unsigned int old_threads = num_threads;
  num_threads              = 0;
  init(old_threads);  // will restart all engines, because num_threads = 0;
}

template <class Engine, class Sseq>
inline void MTPRNG<Engine, true, Sseq>::init(unsigned int const new_num_threads)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("(Re-)initializing Multithreaded RNG Engines")
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("setting number of threads as" << new_num_threads)

  if (new_num_threads <= num_threads)
  {
    // no need to initalize. We do not actually change num_threads until reseed is called.
    return;
  }
  engines.resize(new_num_threads);
  engines.shrink_to_fit();
  uint32_t per_engine_seed[seed_length];
  // else initialize remaining threads
  for (unsigned int i = num_threads; i < new_num_threads; ++i)
  {

    for (unsigned int j = 0; j < seed_length; ++j)
    {
      per_engine_seed[j] = seeder();
    }
    std::seed_seq per_engine_see_seq(per_engine_seed, per_engine_seed + seed_length);
    engines[i].seed(per_engine_see_seq);
  }
  num_threads = new_num_threads;
}

// clang-format off
template <class Engine, class Sseq>
inline void MTPRNG<Engine, false, Sseq>::reseed(Sseq &_seq)
// clang-format on
{
  std::mt19937_64 seeder(_seq);
  uint32_t per_engine_seed[seed_length];
  for (unsigned int j = 0; j < seed_length; ++j)
  {
    per_engine_seed[j] = seeder();
  }
  std::seed_seq derived_seed_seq(per_engine_seed, per_engine_seed + seed_length);
  engine.seed(derived_seed_seq);
}

// implementation of the samplers. Just plain old rejection sampling.
template <class Z, class Engine>
inline Z sample_z_gaussian(double s, double const center, Engine &engine, double const cutoff)
{
  static_assert(std::is_integral<Z>::value,
                "Return type for sample_z_gaussian must be POD integral type.");

  // maximum deviation of the Gaussian from the center. Note that maxdev may be 1.
  Z const maxdev = static_cast<Z>(std::ceil(s * cutoff));

  // uniform distribution on the set of possible outputs.
  std::uniform_int_distribution<Z> uniform_in_range(std::floor(center - maxdev),
                                                    std::ceil(center + maxdev));

  // defaults to value from [0,1), used in rejection sampling to make the decision.
  std::uniform_real_distribution<double> rejection_test(0.0, 1.0);

  // closest int to center, i.e. most likely value.
  Z const closest_int = std::round(center);

  // negative squared distance to most likely value.
  // This is used to scale up the Gaussian weight function s.t. it is 1 at the most likely value.
  double const adj = -(center - closest_int) * (center - closest_int);

  s = s * s / pi;  // overwriting (local copy of) s. We only care about s^2/pi anyway.

  // does not really work in a portable way, so we do not bother
  // std::fenv_t env;
  // feholdexcept( &env); //This disables all floating-point exceptions.

  // use rejection sampling
  while (true)
  {
    Z result    = uniform_in_range(engine);  // sample uniform result.
    double dist = result - center;
    // compute Gaussian weight. std::fma(dist,dist,adj) computes dist^2 + adj = (result-center)^2  -
    // MIN{(result-center)^2 | result integral}.
    //(fma = fused-multiply-add)
    // Recall that s was overwritten to be s^2/pi.

    // Note that the argument of the exp-function might be a tiny positive value due to numeric
    // error
    //(Even if result==closest_int, adj = ROUND((closest_int-center)^2), the computation of
    // std::fma(dist,dist,adj) does not round the intermediate dist^2, leading to a non-zero
    // argument)
    // In particular, it is conceivable that floating point underruns occur in the std::fma - call.
    // Furthermore, if cutoff is large or if s<<1 (in this case, the issue is the rounding when we
    // determined the range), the argument to exp can be extremely small, leading to further
    // potential underruns.
    // We do not care about this for now...

    if (rejection_test(engine) < std::exp(-std::fma(dist, dist, adj) / s))
    {
      // std::feclearexcept(FE_UNDERFLOW);
      // std::feupdateenv(&env);
      return result;
    }
  }
}

// Version taking in s^2/pi and (absolute) maximum deviation. Works just as above.
template <class Z, class Engine>
inline Z sample_z_gaussian_VMD(double const s2pi, double const center, Engine &engine,
                               double const maxdeviation)
{
  // Note : The following allows to access / modify floating point exceptions and modes.
  //#pragma STDC FENV_ACCESS on
  // This is too compiler/implementation-specific and does not work most of the time...

  static_assert(std::is_integral<Z>::value,
                "Return type for sample_z_gaussian must be POD integral type.");
  std::uniform_int_distribution<Z> uniform_in_range(std::floor(center - maxdeviation),
                                                    std::ceil(center + maxdeviation));
  std::uniform_real_distribution<double> rejection_test(0.0, 1.0);
  Z const closest_int = std::round(center);
  double const adj    = -(center - closest_int) * (center - closest_int);

  // std::fenv_t env;
  // feholdexcept( &env); //This disables all floating-point exceptions.

  while (true)
  {
    Z result    = uniform_in_range(engine);  // sample uniform result.
    double dist = result - center;

    if (rejection_test(engine) < std::exp(-std::fma(dist, dist, adj) / s2pi))
    {
      // std::feclearexcept(FE_UNDERFLOW);
      // std::feupdateenv(&env);
      return result;
    }
  }
}

// Samples uniformly at random from the interval [0, max_val]
// clang-format off
template <class Engine>
inline unsigned int sample_uniform(unsigned int max_val, Engine &engine)
// clang-format on
{
  std::uniform_int_distribution<int> uniform_in_range(0, max_val);
  return uniform_in_range(engine);
}

}  // end of namespace GaussSieve

#endif
