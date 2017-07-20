/**
 Implementation for some methods of the sampler. These are in a separate file because they need to
 access the header of the main Sieve class
*/

#ifndef SAMPLER_IMPL_H
#define SAMPLER_IMPL_H

#include "Sampler.h"
#include "SieveGauss.h"
#include <type_traits>
#include <iostream>
#include <random>

// actually needed, even though destructor is pure virtual as the base class destructor is
// eventually called implicitly.
template <class ET, bool MT, class Engine, class Sseq, int nfixed>
Sampler<ET, MT, Engine, Sseq, nfixed>::~Sampler()
{
}

template <class ET, bool MT, class Engine, class Sseq, int nfixed>
void Sampler<ET, MT, Engine, Sseq, nfixed>::init(Sieve<ET, MT, nfixed> *const sieve)
{
  sieveptr = sieve;
  //    std::cout << "Initializing RNGS engines" << std::endl << std::flush;
  engine.init(sieve->get_num_threads());
  //    cout << "Done. Starting custom initialization of specific sampler" << endl << flush;
  custom_init();
  //    cout << "Finished custom initialization" << endl << flush;
}

template <class ET, bool MT, class Engine, class Sseq, int nfixed>
inline std::ostream &operator<<(std::ostream &os,
                                Sampler<ET, MT, Engine, Sseq, nfixed> *const samplerptr)
{
  return samplerptr->dump_to_stream(os);
}

template <class ET, bool MT, class Engine, class Sseq, int nfixed>
inline std::istream &operator>>(std::istream &is,
                                Sampler<ET, MT, Engine, Sseq, nfixed> *const samplerptr)
{
  return samplerptr->read_from_stream(is);
}

template class MTPRNG<std::mt19937_64, false, std::seed_seq>;
// template class MTPRNG<std::mt19937,true,  std::seed_seq>;
template class Sampler<fplll::Z_NR<long>, false, std::mt19937_64, std::seed_seq, -1>;
// template class Sampler<Z_NR<long>, true,  std::mt19937,std::seed_seq>;

#endif
