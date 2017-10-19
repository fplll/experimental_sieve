#ifndef TEST_BITAPPROX_H
#define TEST_BITAPPROX_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "../BitApproximation.h"
#include <limits>
#include "../ExactLatticePoint.h"
#include "gmpxx.h"

bool test_bit_approx()
{
  int constexpr dim = 25;
  int constexpr dimfixed=25;
  
  using LP = GaussSieve::ExactLatticePoint<mpz_class, dimfixed>;
  using GaussSieve::MaybeFixed;
  
  
  return true;
}


#endif