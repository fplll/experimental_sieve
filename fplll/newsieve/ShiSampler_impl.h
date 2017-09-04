/**
Implementation file for ShiSampler.

TODO: Change internal representation of basis.
*/

#ifndef SHI_SAMPLER_IMPL_H
#define SHI_SAMPLER_IMPL_H

#include "Sampler.h"
#include "ShiSampler.h"
//#include "SieveGauss.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include <random>
#include <vector>
#include <math.h>
#include "LatticeBases.h"
#include "DebugAll.h"

namespace GaussSieve
{

template <class SieveTraits, bool MT, class Engine, class Sseq>
void ShiSampler<SieveTraits, MT, Engine, Sseq>::custom_init(SieveLatticeBasis<SieveTraits,MT> const & input_basis)
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

    s2pi[i] = res / GaussSieve::pi; // We rescale to avoid doing this during sampling.
    maxdeviations[i] = sqrt(res) * cutoff;

    basis[i] = input_basis.get_basis_vector(i).make_copy();

//    tmp.set_z(g(i, i));
//    tmp2.div(maxbistar2, tmp);  // s'_i^2 = max GS length^2 / lenght^2 of basis vector
//    s2pi[i] = tmp2.get_d() / GaussSieve::pi;
//    tmp2.sqrt(tmp2);
//    maxdeviations[i] = tmp2.get_d() * cutoff;
  }
  using RetType = typename SieveTraits::GaussSampler_ReturnType;
  bool s = RetType::class_init(MaybeFixed<SieveTraits::get_nfixed>{dim});
  assert(s); // TODO: Clean up and throw exception instead.
  s = SieveTraits::PlainPoint::class_init(MaybeFixed<SieveTraits::get_nfixed>{dim});
  assert(s);
  initialized = true;
}


template <class SieveTraits, bool MT, class Engine, class Sseq>
typename SieveTraits::GaussSampler_ReturnType
ShiSampler<SieveTraits, MT, Engine, Sseq>::sample(int thread)
{
  assert(initialized);
#ifdef DEBUG_SIEVE_STANDALONE_SAMPLER
  assert(sieveptr==nullptr);
#else
  assert(sieveptr!=nullptr);
#endif

  typename SieveTraits::PlainPoint vec;
  vec.fill_with_zero();

  // shift, expressed in coordinates wrt the Gram-Schmidt basis.
  std::vector<double> shifts(lattice_rank, 0.0);

  // Note: This is a while - loop, because --j will cause trouble on unsigned j.
  // (With signed j, the correct for loop would be for(int j = lattice_rank-1 ; j>=0;--j) )
  {
    uint_fast16_t j = lattice_rank;
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