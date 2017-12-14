// clang-format status: NOT OK (reason: templates)

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

#include "DefaultIncludes.h"
#include "LatticeBases.h"
#include "MTPRNG.h"
#include "SieveUtility.h"
#include "Typedefs.h"

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
  user_defined        = 0,
  elliptic_sampler    = 1,  // not implemented
  GPV_sampler         = 2,  // GPVSampler.h
  uniform_sampler     = 3,  // UniformSampler.h
  GPVExtended_sampler = 4   // GPVSamplerExtended.h
};

/**
This class is a generic Sampler. It serves as interface and all other sampler are derived from it.

Sseq is supposed to satisfy the C++ concept "SeedSequence". The standard library has std::seed_seq
as a canonical example.
Engine is supposed to satisfy the C++ concept of a "Random Number Engine". <random> provides several
of those, e.g. std::mt19937_64.

IMPORTANT:  The main sieve can take a user-provided sampler, which may be of any possibly
            user-defined class derived from Sampler.
            This user-defined sampler object will be constructed *before* the sieve and may outlive
            the sieve. For that reason, some initializations are defered from the constructor to
            init() / custom_init, which get called after the sampler is associated with the sieve.
            If the macro DEBUG_SIEVE_STANDALONE_SAMPLER is set, we may even use the sampler without
            any associated sieve. (This latter case need not work with user-defined samplers)
*/

// because operator<<< looks bad and it screws up the declarations by inserting line breaks.
// in fact, clang-format parses it wrongly, leading to operator<< and operator>> getting formatted
// differently...
// clang-format off

template <class SieveTraits, bool MT, class Engine, class Sseq> class Sampler
{
public:
  using GaussSampler_ReturnType = typename SieveTraits::GaussSampler_ReturnType;
  friend std::ostream &operator<< <SieveTraits, MT, Engine, Sseq>(
      std::ostream &os, Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr);
  friend std::istream &operator>> <SieveTraits, MT, Engine, Sseq>(
      std::istream &is, Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr);

  explicit Sampler<SieveTraits, MT, Engine, Sseq>(Sseq &initial_seed)
      : engine(initial_seed), sieveptr(nullptr)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing Sampler (general).")
  }

  // We call init first, then custom_init (via init).
  inline void init(Sieve<SieveTraits, MT> *const sieve,
                   SieveLatticeBasis<SieveTraits, MT> const &input_basis);
  virtual ~Sampler() = 0;  // needs to be virtual

#ifdef PROGRESSIVE
  // set_progressive_rank is virtual, so a child might overwrite it, e.g. to update internal data
  // structures or to output diagnostics whenever the progressive rank changes.
  virtual void set_progressive_rank(uint_fast16_t const new_progressive_rank)
  {
    progressive_rank = new_progressive_rank;
  }
  uint_fast16_t get_progressive_rank() const { return progressive_rank; }
#endif

  /**
  run-time type information.

  This may be used to determine how to interpret a dump file. Defaults to user-defined.
  Other values mean that the GaussSieve dumping routine is aware of the type, simplifying the
  syntax for dumping / reading.

  TODO: Dumping / reading is not yet implemented.
  */
  virtual SamplerType sampler_type() const { return SamplerType::user_defined; }

  // thread is the index of the calling thread (we need to keep separate PRNGs for each thread)
  virtual GaussSampler_ReturnType sample(int const thread = 0) = 0;

  // TODO : Allow sampling in subspaces, updating basis.

private:
  // called before any points are sampled. This function is called from init after sieveptr is set.
  virtual void custom_init(SieveLatticeBasis<SieveTraits,MT> const &input_basis)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("No custom initialization for sampler requested.")
  }
  // dummy implementation of << operator.
  virtual std::ostream &dump_to_stream(std::ostream &os) { return os; }
  // dummy implementation of >> operator.
  virtual std::istream &read_from_stream(std::istream &is) { return is; }

protected:
  MTPRNG<Engine, MT, Sseq> engine;   // or engines
  Sieve<SieveTraits, MT> *sieveptr;  // pointer to parent sieve. Set in init();

#ifdef PROGRESSIVE
  uint_fast16_t progressive_rank;  // In case we sample from a sublattice.
#endif
};

}  // end namespace

#endif
