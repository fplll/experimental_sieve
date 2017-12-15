// clang-format status: NOT OK (nested ifdefs)

/**
 Implementation for some methods of the sampler. These are in a separate file because they need to
 access the header of the main Sieve class
*/

#ifndef SAMPLER_IMPL_H
#define SAMPLER_IMPL_H

#include "DefaultIncludes.h"
#include "LatticeBases.h"
#include "Sampler.h"
#include "Typedefs.h"

namespace GaussSieve
{

// actually needed, even though destructor is pure virtual as the base class destructor is
// eventually called implicitly.
template <class SieveTraits, bool MT, class Engine, class Sseq>
Sampler<SieveTraits, MT, Engine, Sseq>::~Sampler()
{
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
void Sampler<SieveTraits, MT, Engine, Sseq>::init(
    Sieve<SieveTraits, MT> *const sieve, SieveLatticeBasis<SieveTraits, MT> const &input_basis)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing Sampler:")
  sieveptr = sieve;
#ifndef DEBUG_SIEVE_STANDALONE_SAMPLER
  assert(sieveptr != nullptr);
#else
  assert(sieveptr == nullptr);
#endif

#ifdef DEBUG_SIEVE_STANDALONE_SAMPLER
  engine.init(1);
#else
  engine.init(sieve->get_num_threads());
#endif

#ifdef PROGRESSIVE
#ifndef DEBUG_SIEVE_STANDALONE_SAMPLER
  progressive_rank = sieveptr->get_progressive_rank();
  std::cout << "initialized progressive rank to " << progressive_rank << std::endl;
#else
  // PROGRESSIVE && !STANDALONE :
  progressive_rank = input_basis.lattice_rank;
#endif
#endif

  custom_init(input_basis);
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Finished Initializing Sampler.")
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::ostream &operator<<(std::ostream &os,
                                Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr)
{
  std::cerr << "Not implemented yet. (Generic Sampler,streamout)" << std::endl << std::flush;
  return os;
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::istream &operator>>(std::istream &is,
                                Sampler<SieveTraits, MT, Engine, Sseq> *const samplerptr)
{
  std::cerr << "Not implemented yet. (Generic Sampler,streamin)" << std::endl << std::flush;
  return is;
}

// explicit instantiation of template
// template class MTPRNG<std::mt19937_64, false, std::seed_seq>;

}  // end namespace GaussSieve

#endif  // include guards
