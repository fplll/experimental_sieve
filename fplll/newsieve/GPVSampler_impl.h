/**
Implementation file for the GPV sampler.

*/

#ifndef GPV_SAMPLER_IMPL_H
#define GPV_SAMPLER_IMPL_H

#include "DefaultIncludes.h"
#include "Sampler.h"
#include "GPVSampler.h"
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
  assert(!initialized);
#ifndef DEBUG_SIEVE_STANDALONE_SAMPLER
  assert(sieveptr!=nullptr);
#else
  assert(sieveptr==nullptr);
#endif

  dim           = input_basis.ambient_dimension;
  lattice_rank  = input_basis.lattice_rank;
  mu_matrix     = input_basis.get_mu_matrix();

  // vectors of length lattice_rank
  s2pi.resize(lattice_rank);
  maxdeviations.resize(lattice_rank);
  basis.resize(lattice_rank);

  auto const maxbistar2 = input_basis.get_maxbistar2();

  for (uint_fast16_t i = 0; i < lattice_rank; ++i)
  {
    //
    // Note: input_basis.get_g(i,i) might be of type mpz_class.
    // We can not convert after division, because double / mpz_class always gives 0.
    // (and also, double / mpz_class has an expression template as return type, which screws up
    // convert_to_double)

    double res = maxbistar2 / convert_to_double(input_basis.get_g(i,i));

    // the scaling of 1.0 is somewhat arbitrary but works ok
    //
    s2pi[i] = 1.0*res / GaussSieve::pi; // We rescale by pi to avoid doing this during sampling.
    maxdeviations[i] = sqrt(res) * cutoff;

    basis[i] = input_basis.get_basis_vector(i).make_copy();
  }

  using RetType = typename SieveTraits::GaussSampler_ReturnType;

  if(static_init_plainpoint!=nullptr)
  {
    assert(false);
  }
  if(static_init_rettype!=nullptr)
  {
    assert(false);
  }

  static_init_rettype   = new StaticInitializer<RetType>(MaybeFixed<SieveTraits::get_nfixed>{dim});
  static_init_plainpoint= new StaticInitializer<typename SieveTraits::PlainPoint>(MaybeFixed<SieveTraits::get_nfixed>{dim});

  initialized = true;
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


  // shifts expressed in coordinates wrt the Gram-Schmidt basis.
  std::vector<double> shifts(lattice_rank, 0.0);
  
  
  // Note: This is a while - loop, because --j will cause trouble on unsigned j.
  // (With signed j, the correct for loop would be for(int j = lattice_rank-1 ; j>=0;--j) )
  //std::cout << vec << std::endl;
  while ( vec.is_zero())
  {
    #ifdef PROGRESSIVE
    uint_fast16_t j = this->get_progressive_rank();
    #else
    uint_fast16_t j = lattice_rank;
    #endif

    //std::cout << "j = " << j <<std::endl;
    while(j>0)
    {
      --j;
      long const newcoeff = sample_z_gaussian_VMD<long, Engine>(
        s2pi[j], shifts[j], engine.rnd(), maxdeviations[j]);  // coefficient of b_j in vec.

      vec += basis[j] * newcoeff;

      for (uint_fast16_t i = 0; i < j; ++i)  // adjust shifts
      {
        shifts[i] -= newcoeff * (mu_matrix[j][i]);
      }
    }
  }

  typename SieveTraits::GaussSampler_ReturnType ret;
  ret = make_from_any_vector<typename SieveTraits::GaussSampler_ReturnType>(vec, dim);
  return ret;
}
} // end namespace

#endif
