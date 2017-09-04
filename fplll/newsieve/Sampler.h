/**
 This file provided the interface for the lattice point sampler.
 To this end, this file provides a virtual interface template class Sampler, from which the actual
 samplers are derived.

 Dependencies:
 The sampler instance contains a pointer back to the sieve main object, creating circular
 dependencies.
 For this reason, parts of the implementation are in Sampler_impl.h or Sampler.cpp (TODO: Decide
 on whether we want header-only)
*/

#ifndef SAMPLER_H
#define SAMPLER_H

#include "MTPRNG.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include <iostream>
#include <random>
#include <type_traits>
#include "LatticeBases.h"
//#include <cfenv>
#include "DebugAll.h"

namespace GaussSieve
{

// forward declarations
template <class SieveTraits, bool MT> class Sieve;

// declared in this file, as specialization for each possible MT.
template <class SieveTraits, bool MT, class Engine, class Sseq> class Sampler;

// printing
template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::ostream &operator<<(std::ostream &os,
                                Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr);

// reading (may also be used by constructor from istream)
template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::istream &operator>>(std::istream &is,
                                Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr);

enum class SamplerType
{
  user_defined     = 0,
  elliptic_sampler = 1,
  shi_sampler      = 2,
  gauss_sampler    = 3
};

/**
This class is a generic Sampler. It serves as interface and all other sampler are derived from it.

Sseq is supposed to satisfy the C++ concept "SeedSequence". The standard library has std::seed_seq
as a canonical example.
Engine is supposed to satisfy the C++ concept of a "Random Number Engine". <random> provides several
of those, e.g. std::mt19937_64.
*/

// because operator<<< looks bad and it screws up the declarations by inserting line breaks.
// clang-format off

template <class SieveTraits, bool MT, class Engine, class Sseq> class Sampler
{
public:
  using GaussSampler_ReturnType = typename SieveTraits::GaussSampler_ReturnType;
  friend std::ostream & operator<< <SieveTraits, MT, Engine, Sseq>
            (std::ostream &os, Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr);
  friend std::istream & operator>> <SieveTraits, MT, Engine, Sseq>
            (std::istream &is, Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr);

  explicit Sampler<SieveTraits, MT, Engine, Sseq>(Sseq &initial_seed)
          : engine(initial_seed), sieveptr(nullptr)
    {
      DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing Sampler (general).")
    }

  // We call init first, then custom_init (via init).
  inline void init(Sieve<SieveTraits, MT> *const sieve,
    SieveLatticeBasis<SieveTraits,MT> const & input_basis );
  virtual ~Sampler() = 0;  // needs to be virtual

  /**
  run-time type information.

  This may be used to determine how to interpret a dump file. Defaults to user-defined.
  Other values mean that the GaussSieve dumping routine is aware of the type, simplifying the
  syntax for dumping / reading.

  TODO: Dumping / reading is not yet implemented.
  */
  virtual SamplerType sampler_type() const { return SamplerType::user_defined; };

  virtual GaussSampler_ReturnType sample(int thread = 0) = 0;  // thread is the index of the calling
                                                               // thread (we need to keep separate
                                                               // PRNGs for each thread)

  // TODO : Allow sampling in subspaces, updating basis.

private:
  // called before any points are sampled. This function is called from init after sieveptr is set.
  virtual void custom_init(SieveLatticeBasis<SieveTraits,MT> const & input_basis)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("No custom initialization for sampler requested.")
  };
  // dummy implementation of << operator.
  virtual std::ostream &dump_to_stream(std::ostream &os) { return os; };
  // dummy implementation of >> operator.
  virtual std::istream &read_from_stream(std::istream &is) { return is; };

protected:
  MTPRNG<Engine, MT, Sseq> engine;  // or engines
  Sieve<SieveTraits, MT> *sieveptr;  // pointer to parent sieve. Set in init();
};

// clang-format on

}  // end namespace

#endif