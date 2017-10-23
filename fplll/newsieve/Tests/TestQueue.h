#ifndef TEST_QUEUE_H
#define TEST_QUEUE_H

#include <iostream>
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"
#include "../ShiSampler.h"
#include <random>
#include "../LatticeBases.h"
#include "../GaussQueue.h"

bool test_queue()
{
  int constexpr dim = 20;
  fplll::ZZ_mat<mpz_t> B;
  B.resize(dim, dim);
  srand (1);

  //generates GM lattice
  B.gen_qary_prime(1, 10*dim);

  fplll::lll_reduction(B, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::LM_WRAPPER);
  using Traits = GaussSieve::DefaultSieveTraits<mpz_class, false, -1>;
  typename Traits::GlobalStaticDataInitializer init_arg (dim);
  using Sampler = GaussSieve::ShiSampler<Traits,false,std::mt19937,std::seed_seq>;
  GaussSieve::SieveLatticeBasis<Traits,false> sieve_basis(B,init_arg); // convert to SieveLatticeBasis


  // sampler2 is a free-standing sampler for B.
  std::seed_seq sseq {1,2,3,5};
  Sampler sampler2(sseq);
  sampler2.init(nullptr, sieve_basis);

  // gauss_queue contains another sampler.
  GaussSieve::GaussQueue<Traits,false> gauss_queue(nullptr,init_arg);
  gauss_queue.sampler->init(nullptr, sieve_basis);

  std::cout << "taking 20 elements from the queue (sampling)" << std::endl;
  for(int i=0;i<20;++i)
  {
  std::cout << gauss_queue.true_pop() << std::endl;
  }

  std::cout << "sampling single element" << std::endl;

  auto LP = sampler2.sample(); // LP is non-zero for the above sseq.
  std::cout << LP << std::endl;
  auto LP2 = LP.make_copy();
  gauss_queue.push(std::move(LP));
  assert(gauss_queue.size()==1);
  LP = gauss_queue.true_pop();
  assert(gauss_queue.size()==0);
  assert(LP==LP2);

  return true;
}

#endif
