#ifndef TEST_APPROXIMATIONS_H
#define TEST_APPROXIMATIONS_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "../EMVApproximation.h"
#include "gmpxx.h"

#include "../ApproximatedPoint.h"

bool test_approximations()
{
  int constexpr dimfixed = 25;

  using ExactLP = GaussSieve::ExactLatticePoint<mpz_class, dimfixed>;
  using Approx  = GaussSieve::EMVApproximation<dimfixed>;
  using CombinedLP = GaussSieve::VectorWithApproximation<ExactLP,Approx>;
  using GaussSieve::MaybeFixed;
  GaussSieve::StaticInitializerArg<MaybeFixed<dimfixed>> init_arg{ MaybeFixed<dimfixed>{dimfixed} };
  GaussSieve::StaticInitializer<CombinedLP> init1 (init_arg);

  return true;
}








#endif
