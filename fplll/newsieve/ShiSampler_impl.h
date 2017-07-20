/**
Implementation file for ShiSampler.

TODO: Change internal representation of basis.
*/

#ifndef SHI_SAMPLER_IMPL_H
#define SHI_SAMPLER_IMPL_H

#include "Sampler.h"
#include "ShiSampler.h"
#include "SieveGauss.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include <vector>
#include <random>

template <class ET, bool MT, class Engine, class Sseq, int nfixed>
void ShiSampler<ET, MT, Engine, Sseq, nfixed>::custom_init()
{
  current_basis = sieveptr->get_original_basis();
  dim           = static_cast<Dimension<nfixed>>(sieveptr->get_ambient_dimension());
  lattice_rank  = sieveptr->get_lattice_rank();
  fplll::Matrix<ET> u, u_inv, g;  // intentionally uninitialized.

  fplll::MatGSO<ET, fplll::FP_NR<double>> GSO(current_basis, u, u_inv,
                                              fplll::MatGSOInterfaceFlags::GSO_INT_GRAM);
  GSO.update_gso();  // todo: raise exception in case of error.

  mu = GSO.get_mu_matrix();

  s2pi.resize(lattice_rank);
  maxdeviations.resize(lattice_rank);

  g = GSO.get_g_matrix();

  fplll::FP_NR<double> maxbistar2 = GSO.get_max_bstar();

  fplll::FP_NR<double> tmp;
  fplll::FP_NR<double> tmp2;
  for (unsigned int i = 0; i < lattice_rank; ++i)
  {
    tmp.set_z(g(i, i));
    tmp2.div(maxbistar2, tmp);  // s'_i^2 = max GS length^2 / lenght^2 of basis vector
    s2pi[i] = tmp2.get_d() / GaussSieve::pi;
    tmp2.sqrt(tmp2);
    maxdeviations[i] = tmp2.get_d() * cutoff;
  }

  auto it = helper_current_basis.cend();
  for (unsigned int i = 0; i < lattice_rank; ++i)
  {
    it = helper_current_basis.cend();
    helper_current_basis.emplace(it, current_basis[i]);
  }
}

template <class ET, bool MT, class Engine, class Sseq, int nfixed>
ShiSampler<ET, MT, Engine, Sseq, nfixed>::~ShiSampler()
{
}

template <class ET, bool MT, class Engine, class Sseq, int nfixed>
typename GaussSieve::GaussSampler_ReturnType<ET, MT, nfixed>
ShiSampler<ET, MT, Engine, Sseq, nfixed>::sample(int thread)
{
  assert(sieveptr != nullptr);
  //    MyLatticePoint<ET,nfixed> *vec = new MyLatticePoint<ET,nfixed>(dim);
  MyLatticePoint<ET, nfixed> vec;
  vec.fill_with_zero();
  // vec->NumVect<ET>::fill(0); //current vector built up so far. //Note: We treat vec as a NumVect
  // until the end, because we don't want to normalize intermediate results.

  // shift, expressed in coordinates wrt the Gram-Schmidt basis.
  std::vector<double> shifts(lattice_rank, 0.0);
  for (int j = lattice_rank - 1; j >= 0; --j)
  {
    long const newcoeff = GaussSieve::sample_z_gaussian_VMD<long, Engine>(
        s2pi[j], shifts[j], engine.rnd(), maxdeviations[j]);  // coefficient of b_j in vec.
    ET newcoeffET;
    newcoeffET = newcoeff;
    // vec+= current_basis[j].get_underlying_row(); //build up vector

    vec = add(vec, scalar_mult(helper_current_basis[j], newcoeffET));

    //        vec->NumVect<ET>::addmul_si(current_basis[j].get_underlying_row(), newcoeff);
    for (int i = 0; i < j; ++i)  // adjust shifts
    {
      shifts[i] -= newcoeff * (mu[j][i].get_d());
    }
  }
  //    vec->normalize();
  return vec;
  //    return static_cast<CompressedPoint<ET,MT,-1> > (vec); //Note: This makes copies, which are
  //    unneccessary...
}

#endif
