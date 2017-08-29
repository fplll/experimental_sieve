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

namespace GaussSieve
{

template <class SieveTraits, bool MT, class Engine, class Sseq>
void ShiSampler<SieveTraits, MT, Engine, Sseq>::custom_init(SieveLatticeBasis<SieveTraits,MT> const & input_basis)
{
  assert(!initialized);
//  SieveLatticeBasis<SieveTraits,MT> const & current_basis = sieveptr->get_original_basis();

  dim           = input_basis.ambient_dimension;
  lattice_rank  = input_basis.lattice_rank;
  mu_matrix     = input_basis.get_mu_matrix();

  // vectors of length lattice_rank
  s2pi.resize(lattice_rank);
  maxdeviations.resize(lattice_rank);
  basis.resize(lattice_rank);

  auto const maxbistar2 = input_basis.get_maxbistar2();

//  g = GSO.get_g_matrix();

//  fplll::FP_NR<double> maxbistar2 = GSO.get_max_bstar();

//  fplll::FP_NR<double> tmp;
//  fplll::FP_NR<double> tmp2;
  for (uint_fast16_t i = 0; i < lattice_rank; ++i)
  {
    //
    // Note: input_basis.get_g(i,i) might be of type mpz_class.
    // In this case, maxbistar / ... is NOT of type mpz_class, but rather of type
    // __gmp_expr<...>, which is an expression template. This is how mpz_class implements lazy
    // evaluation. The static_cast essentially turns of this lazy evaluation, because otherwise the
    // argument type of convert_to_double (which is specialized for mpz_class) does not work.
    // Alternatively, we could fiddle around with convert_to_double.

//    double res = convert_to_double( static_cast<decltype(input_basis.get_g(i,i))> (
//      maxbistar2 / input_basis.get_g(i,i)
//      ));

    double res = maxbistar2 / convert_to_double(input_basis.get_g(i,i));

    std::cout << "Setting i" << i << "as" << maxbistar2 / input_basis.get_g(i,i);
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
//  auto it = helper_current_basis.cend();
//  for (unsigned int i = 0; i < lattice_rank; ++i)
//  {
//    it = helper_current_basis.cend();
//    helper_current_basis.emplace(it, current_basis[i]);
//  }
}


template <class SieveTraits, bool MT, class Engine, class Sseq>
typename SieveTraits::GaussSampler_ReturnType
ShiSampler<SieveTraits, MT, Engine, Sseq>::sample(int thread)
{
  assert(initialized);

//  assert(sieveptr!=nullptr);

  //    MyLatticePoint<ET,nfixed> *vec = new MyLatticePoint<ET,nfixed>(dim);
  typename SieveTraits::PlainPoint vec;
  vec.fill_with_zero();
//  MyLatticePoint<ET, nfixed> vec;
//  vec.fill_with_zero();
  // vec->NumVect<ET>::fill(0); //current vector built up so far. //Note: We treat vec as a NumVect
  // until the end, because we don't want to normalize intermediate results.

  // shift, expressed in coordinates wrt the Gram-Schmidt basis.
  std::vector<double> shifts(lattice_rank, 0.0);
  for (uint_fast16_t j = lattice_rank - 1; j >= 0; --j) // Note: This would loop until underflow.
  {
    long const newcoeff = sample_z_gaussian_VMD<long, Engine>(
        s2pi[j], shifts[j], engine.rnd(), maxdeviations[j]);  // coefficient of b_j in vec.
//    ET newcoeffET;
//    newcoeffET = newcoeff;
    // vec+= current_basis[j].get_underlying_row(); //build up vector

    vec += basis[j] * newcoeff;

//    vec = add(vec, scalar_mult(helper_current_basis[j], newcoeffET));

    //        vec->NumVect<ET>::addmul_si(current_basis[j].get_underlying_row(), newcoeff);
    for (uint_fast16_t i = 0; i < j; ++i)  // adjust shifts
    {
      shifts[i] -= newcoeff * (mu_matrix[j][i]);
    }

    if(j==0) break;

  }
  typename SieveTraits::GaussSampler_ReturnType ret;
  ret = make_from_any_vector<typename SieveTraits::GaussSampler_ReturnType>(vec, dim);
  //    vec->normalize();
  return ret;
  //    return static_cast<CompressedPoint<ET,MT,-1> > (vec); //Note: This makes copies, which are
  //    unneccessary...
}
}

#endif
