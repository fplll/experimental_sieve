#ifndef UNIFORM_SAMPLER_H
#define UNIFORM_SAMPLER_H

#include "DefaultIncludes.h"
#include "Sampler.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"
#include <random>
#include <vector>
#include "LatticeBases.h"


namespace GaussSieve
{
  template <class SieveTraits, bool MT> class Sieve;
  
  
  template <class SieveTraits, bool MT, class Engine, class Sseq> class UniformSampler;
  
  template <class SieveTraits, bool MT, class Engine, class Sseq>
  class UniformSampler final : public Sampler<SieveTraits, MT, Engine, Sseq>
  
public:
  using DimensionType = typename SieveTraits::DimensionType;
  using EntryType     = typename SieveTraits::EntryType;
  using RetType       = typename SieveTraits::GaussSampler_ReturnType;
}


#endif