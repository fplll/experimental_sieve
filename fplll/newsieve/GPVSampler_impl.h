#ifndef GPV_SAMPLER_IMPL_H
#define GPV_SAMPLER_IMPL_H

#include "DefaultIncludes.h"
#include "Sampler.h"
#include "ShiSampler.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include <random>
#include <vector>
#include <math.h>
#include "LatticeBases.h"

namespace GaussSieve
{
  
  template <class SieveTraits, bool MT, class Engine, class Sseq>
  void GPVSampler<SieveTraits, MT, Engine, Sseq>::custom_init(SieveLatticeBasis<SieveTraits,MT> const & input_basis)
  {
  }
  
  template <class SieveTraits, bool MT, class Engine, class Sseq>
  typename SieveTraits::GaussSampler_ReturnType
  GPVSampler<SieveTraits, MT, Engine, Sseq>::sample(int const thread)
  {
    assert(initialized);
#ifdef DEBUG_SIEVE_STANDALONE_SAMPLER
    assert(sieveptr==nullptr);
#else
    assert(sieveptr!=nullptr);
#endif
    
    typename SieveTraits::PlainPoint vec;
    vec.fill_with_zero();
    
  }
  
}
