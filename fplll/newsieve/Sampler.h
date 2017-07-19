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

// forward declarations
#include "MTPRNG.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include <iostream>
#include <random>
#include <type_traits>
//#include <cfenv>

// forward declarations

template <class ET, bool MT, int nfixed> class Sieve;

template <class ET, bool MT, class Engine, class Sseq, int nfixed> class Sampler;

// printing
template <class ET, bool MT, class Engine, class Sseq, int nfixed>
inline std::ostream &operator<<(std::ostream &os,
                                Sampler<ET, MT, Engine, Sseq, nfixed> *const samplerptr);

// reading (may also be used by constructor from istream)
template <class ET, bool MT, class Engine, class Sseq, int nfixed>
inline std::istream &operator>>(std::istream &is,
                                Sampler<ET, MT, Engine, Sseq, nfixed> *const samplerptr);

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

// clang-format off

template <class ET, bool MT, class Engine, class Sseq, int nfixed> class Sampler
{
public:
  using GaussSampler_ReturnType = typename GaussSieve::GaussSampler_ReturnType<ET, MT, nfixed>;
  friend std::ostream & operator<< <ET, MT, Engine, Sseq, nfixed>
            (std::ostream &os, Sampler<ET, MT, Engine, Sseq, nfixed> *const samplerptr);
  friend std::istream & operator>> <ET, MT, Engine, Sseq, nfixed>
            (std::istream &is, Sampler<ET, MT, Engine, Sseq, nfixed> *const samplerptr);

  explicit Sampler<ET, MT, Engine, Sseq, nfixed>(Sseq &initial_seed)
      : engine(initial_seed), sieveptr(nullptr) {}

  // We call init first, then custom_init (via init).
  void init(Sieve<ET, MT, nfixed> *const sieve);
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
  virtual void custom_init(){};
  // dummy implementation of << operator.
  virtual std::ostream &dump_to_stream(std::ostream &os) { return os; };
  // dummy implementation of >> operator.
  virtual std::istream &read_from_stream(std::istream &is) { return is; };

protected:
  MTPRNG<Engine, MT, Sseq> engine;  // or engines
  Sieve<ET, MT, nfixed> *sieveptr;  // pointer to parent sieve. Set in init();
};

// clang-format on

// clang-format off


#include "SieveGauss.h"


//actually needed, even though destructor is pure virtual as the base class destructor is eventually called implicitly.
template <class ET,bool MT, class Engine, class Sseq, int nfixed>
Sampler<ET,MT, Engine,Sseq,nfixed>::~Sampler() {}

template <class ET,bool MT, class Engine, class Sseq, int nfixed> void Sampler<ET,MT,Engine,Sseq,nfixed>::init(Sieve<ET,MT,nfixed> * const sieve)
{
    sieveptr = sieve;
//    std::cout << "Initializing RNGS engines" << std::endl << flush;
    engine.init(sieve->get_num_threads());
//    cout << "Done. Starting custom initialization of specific sampler" << endl << flush;
    custom_init();
//    cout << "Finished custom initialization" << endl << flush;
}

template<class ET,bool MT, class Engine, class Sseq, int nfixed>
inline std::ostream & operator<<(std::ostream &os,Sampler<ET,MT,Engine,Sseq,nfixed>* const samplerptr)
{
  return samplerptr->dump_to_stream(os);
}

template<class ET,bool MT, class Engine, class Sseq, int nfixed>
inline std::istream & operator>>(std::istream &is,Sampler<ET,MT,Engine,Sseq,nfixed>* const samplerptr)
{
  return samplerptr->read_from_stream(is);
}

#endif

//clang-format on
