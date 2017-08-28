#ifndef TEST_SHI_SAMPLER_H
#define TEST_SHI_SAMPLER_H

//#include "../ShiSampler.h"
#include <iostream>
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"

bool test_shi_sampler()
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

  return true;
}


#endif
