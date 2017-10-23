#ifndef TEST_BASIS_UTILS_H
#define TEST_BASIS_UTILS_H

#include "../LatticeBases.h"
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"
#include <iostream>
#include "../GlobalStaticData.h"

bool test_basis_utils()
{
  int constexpr dim = 20;
  fplll::ZZ_mat<mpz_t> B;
  B.resize(dim, dim);
  srand (1);

  //generates GM lattice
  B.gen_qary_prime(1, 10*dim);

  fplll::lll_reduction(B, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::LM_WRAPPER);

  std::cout << "Generating random basis: Result is" << std::endl;
  std::cout << B << std::endl << std::flush;

  using Traits = GaussSieve::DefaultSieveTraits<mpz_class, false, -1>;
  typename Traits::GlobalStaticDataInitializer init_args(dim);


  GaussSieve::SieveLatticeBasis<Traits,false> sieve_basis(B,init_args);

  double maxbistar2 = sieve_basis.get_maxbistar2();

  std::cout << "maxbistar = " << maxbistar2 << std::endl;

  std::vector<std::vector<double>> mu = sieve_basis.get_mu_matrix();
  std::cout << "Mu-matrix:" << std::endl;
  for(uint_fast16_t i=0; i < dim ; ++i)
  {
    for(uint_fast16_t j=0;j<dim;++j)
    {
      std::cout << mu[i][j] << " ";

      if (j>i)
      {
        assert(mu[i][j] == sieve_basis.get_mu(i,j));
      }
    }
    std::cout << std::endl;
  }

  auto g = sieve_basis.get_g_matrix();
  std::cout << "g-matrix:" << std::endl;
  for(uint_fast16_t i=0; i < dim ; ++i)
  {
    for(uint_fast16_t j=0;j<dim;++j)
    {
      std::cout << g[i][j] << " ";
      if (j>=i)
      {
        assert(g[i][j] == sieve_basis.get_g(i,j));
      }
    }
    std::cout << std::endl;
  }
  for(uint_fast16_t i=0; i<dim;++i)
  {
    std::cout << sieve_basis.get_basis_vector(i) << std::endl;
  }

  std::cout << "Minkowski bound" << sieve_basis.get_minkowski_bound() << std::endl;

  return true;
}


#endif
