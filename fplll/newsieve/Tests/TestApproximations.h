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

bool test_approximations()
{
  /*
  int constexpr dimfixed = 25;

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
