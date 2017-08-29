/**
 Implementation for some methods of the sampler. These are in a separate file because they need to
 access the header of the main Sieve class
*/

#ifndef SAMPLER_IMPL_H
#define SAMPLER_IMPL_H

#include "Typedefs.h"
#include "DebugAll.h"
#include "Sampler.h"
#include "SieveGauss.h"
#include <iostream>
#include <random>
#include <type_traits>

namespace GaussSieve
{

// actually needed, even though destructor is pure virtual as the base class destructor is
// eventually called implicitly.
template <class SieveTraits, bool MT, class Engine, class Sseq>
Sampler<SieveTraits, MT, Engine, Sseq>::~Sampler()
{
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
void Sampler<SieveTraits, MT, Engine, Sseq>::init(Sieve<SieveTraits, MT> *const sieve)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing Sampler:")
  sieveptr = sieve;
  //    std::cout << "Initializing RNGS engines" << std::endl << std::flush;
  engine.init(sieve->get_num_threads());
  //    cout << "Done. Starting custom initialization of specific sampler" << endl << flush;
  custom_init();
  //    cout << "Finished custom initialization" << endl << flush;
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::ostream &operator<<(std::ostream &os,
                                Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr)
{
  return samplerptr->dump_to_stream(os); // virtual dispatch
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::istream &operator>>(std::istream &is,
                                Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr)
{
  assert(false); // The line below looks wrong.
  return samplerptr->read_from_stream(is); //virtual dispatch
}

template class MTPRNG<std::mt19937_64, false, std::seed_seq>;
// template class MTPRNG<std::mt19937,true,  std::seed_seq>;
// template class Sampler<fplll::Z_NR<long>, false, std::mt19937_64, std::seed_seq, -1>;
template class Sampler<DefaultSieveTraits, false, std::mt19937_64, std::seed_seq>;

// template class Sampler<Z_NR<long>, true,  std::mt19937,std::seed_seq>;
}

#endif
