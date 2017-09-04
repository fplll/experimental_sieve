#ifndef TEST_SHI_SAMPLER_H
#define TEST_SHI_SAMPLER_H

#include <iostream>
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"
#include "../ShiSampler.h"
#include <random>
#include "../LatticeBases.h"

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
  using Sampler = GaussSieve::ShiSampler<Traits,false,std::mt19937,std::seed_seq>;
  GaussSieve::SieveLatticeBasis<Traits,false> sieve_basis(B); // convert to SieveLatticeBasis

  std::seed_seq sseq {1,2,3};

  Sampler sampler(sseq);
  sampler.init(nullptr, sieve_basis);
  for(int i=0; i<100 ; ++i)
  {
    std::cout << sampler.sample() << std::endl << std::flush;
  }


  return true;
}


#endif
