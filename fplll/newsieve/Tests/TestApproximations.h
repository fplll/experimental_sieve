#ifndef TEST_APPROXIMATIONS_H
#define TEST_APPROXIMATIONS_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "../EMVApproximation.h"
#include "gmpxx.h"
#include <array>

#include "../ApproximatedPoint.h"
#include "../EMVApproximation.h"

bool test_approximations()
{
  using namespace GaussSieve;
  long x =0L;
  long y =1L;
  using LeveledLong = MakeLeveledScalar<long>;
  LeveledLong z = y;
  assert(x<z);
  assert(z.access<0>() == 1);
  LeveledLong z2 = 56L;
  assert(2*z2 == 112);
  assert(z2*2 == 112);
  using Approx = AddApproximation<LeveledLong,EMVScalar>;
  Approx a = 32L;
  Approx b = 34L;
  assert(a==32L);
  assert(a<b);


  int constexpr dimfixed = 25;
  using LatP = ExactLatticePoint<mpz_class,dimfixed>;
  using LeveledLatP = MakeLeveledLatticePoint<LatP>;
  StaticInitializerArg<MaybeFixed<dimfixed>> init_arg{ MaybeFixed<dimfixed>{dimfixed}};
  StaticInitializer<LeveledLatP> init1 (init_arg);

  std::array<mpz_class,dimfixed> arr;
  for(int i=0;i<dimfixed;++i)
  {
    arr[i] = 400 *i * i;
  }
  LeveledLatP latp = make_from_any_vector<LatP>(arr,MaybeFixed<dimfixed>{dimfixed});

  /*
  using ExactLP = GaussSieve::ExactLatticePoint<mpz_class, dimfixed>;
  using Approx  = GaussSieve::EMVApproximation<dimfixed>;
  using CombinedLP = GaussSieve::VectorWithApproximation<ExactLP,Approx>;
  using GaussSieve::MaybeFixed;
  GaussSieve::StaticInitializerArg<MaybeFixed<dimfixed>> init_arg{ MaybeFixed<dimfixed>{dimfixed} };
  GaussSieve::StaticInitializer<CombinedLP> init1 (init_arg);

  std::array<mpz_class,dimfixed> arr;
  for(int i=0;i<dimfixed;++i)
  {
    arr[i] = 400 *i * i;
  }

  ExactLP latp_exact = GaussSieve::make_from_any_vector<ExactLP>(arr,MaybeFixed<dimfixed>{dimfixed});
  CombinedLP latp_comb (std::move(latp_exact));
  std::cout << latp_comb << std::endl;
//
  std::cout << compute_sc_product(latp_comb.access_exact(),latp_comb.access_exact() ) << std::endl;
  std::cout << compute_sc_product_approx(latp_comb.access_approx(), latp_comb.access_approx() ) << std::endl;
  std::cout << latp_comb.access_exact().get_norm2() << std::endl;
  mpz_class X = latp_comb.get_norm2(); std::cout << "Delayed:" << X << std::endl;
//
//  std::cout <<std::endl<<std::endl << "--------------" << std::endl;

*/

  return true;
}








#endif
