#ifndef UNIFORM_SAMPLER_IMPL_H
#define UNIFORM_SAMPLER_IMPL_H

#include "DefaultIncludes.h"
#include "LatticeBases.h"
#include "Sampler.h"
#include "UniformSampler.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"

namespace GaussSieve
{

template <class SieveTraits, bool MT, class Engine, class Sseq>
void UniformSampler<SieveTraits, MT, Engine, Sseq>::custom_init(
    SieveLatticeBasis<SieveTraits, MT> const &input_basis)
{
  assert(!initialized);
#ifndef DEBUG_SIEVE_STANDALONE_SAMPLER
  assert(sieveptr != nullptr);
#else
  assert(sieveptr == nullptr);
#endif

  if (static_init_plainpoint != nullptr)
  {
    assert(false);
  }
  if (static_init_rettype != nullptr)
  {
    assert(false);
  }

  dim          = input_basis.ambient_dimension;
  lattice_rank = input_basis.lattice_rank;
  basis.resize(lattice_rank);

  static_init_rettype    = new StaticInitializer<RetType>(MaybeFixed<SieveTraits::get_nfixed>{dim});
  static_init_plainpoint = new StaticInitializer<typename SieveTraits::PlainPoint>(
      MaybeFixed<SieveTraits::get_nfixed>{dim});

  // copy the basis
  for (uint_fast16_t i = 0; i < lattice_rank; ++i)
  {
    basis[i] = input_basis.get_basis_vector(i).make_copy();
  }
  initialized = true;
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
typename SieveTraits::GaussSampler_ReturnType
UniformSampler<SieveTraits, MT, Engine, Sseq>::sample(int const thread)
{
  assert(initialized);
#ifdef DEBUG_SIEVE_STANDALONE_SAMPLER
  assert(sieveptr == nullptr);
#else
  assert(sieveptr != nullptr);
#endif

  typename SieveTraits::PlainPoint vec;
  vec.fill_with_zero();

  for (unsigned int i = 0; i < sparcity; ++i)
  {
    vec += basis[sample_uniform<Engine>(dim - 1, engine.rnd(thread))];
    vec -= basis[sample_uniform<Engine>(dim - 1, engine.rnd(thread))];
  }

  typename SieveTraits::GaussSampler_ReturnType ret;

  // TODO : Fix conversion here.
  ret = make_from_any_vector<typename SieveTraits::GaussSampler_ReturnType>(vec, dim);
  return ret;
}

}  // end namespace GaussSieve

#endif
