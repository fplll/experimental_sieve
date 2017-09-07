#ifndef TEST_EMV_H
#define TEST_EMV_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "../EMVApproximation.h"
#include <limits>
#include "../ExactLatticePoint.h"
#include "gmpxx.h"

bool test_emv()
{
  double xd = 0.6;
  int xi = 2352626;
  long xl = 2352627;
  mpz_class xm = 236;
  int xi2 = -65536;

  using GaussSieve::EMVScalar;
  std::cout << "Size of approximated Norm2s is" << std::numeric_limits<typename GaussSieve::EMVApproximationTraits::ApproxNorm2Type>::digits << " bits." << std::endl;
  std::cout << "Size of approximated entries is" << std::numeric_limits<typename GaussSieve::EMVApproximationTraits::ApproxEntryType>::digits << " bits." << std::endl;

  EMVScalar ad(xd);
  EMVScalar ai = static_cast<EMVScalar> (xi);
  EMVScalar ai22(-1,16);
  EMVScalar ai2 = static_cast<EMVScalar> (xi2);
  EMVScalar am = static_cast<EMVScalar>(xm);
  EMVScalar al(xl);

  std::cout << "Exact:" << xd << " approximated by " << ad << std::endl;
  std::cout << "Exact:" << xi << " approximated by " << ai << std::endl;
  std::cout << "Exact:" << xl << " approximated by " << al << std::endl;
  std::cout << "Exact:" << xm << " approximated by " << am << std::endl;
  std::cout << "Exact:" << xi2 << " approximated by " << ai2 << std::endl;

  mpz_class mpz1 = 1024;
  mpz_class mpz2 = 1025;
  mpz_class mpz3 =-1024;
  mpz_class mpz4 =-1025;

  EMVScalar app1(mpz1);
  EMVScalar app2(mpz2);
  EMVScalar app3(mpz3);
  EMVScalar app4(mpz4);

  assert(app1 < app2);
  assert(app2 > app1);

  assert(app3 < app1);
  assert(app4 < app1);
  assert(app3 < app2);
  assert(app4 < app2);
  assert(app4 < app3);
  assert(app3 > app4);
  assert(app1 > 0);
  assert(app1 > 1023.5);
  assert(app1 < 1025);
  assert(app1 < 1024.01);

  double d1 = 32;
  double d2 = -33;
  double d3 = 0.24;
  double d4 = 31.999;

  long l1 = 32;
  long l2 = -33;
  long l3 = 0;

  mpz_class m1 = 32;
  mpz_class m2 = -33;
  mpz_class m3 = 0;

  #define PRINT_VAL_WITH_EXP(arg) std::cout << "Value = " << arg << " exponent =" << EMVScalar::get_exponent(arg) << std::endl;

  PRINT_VAL_WITH_EXP(d1);
  PRINT_VAL_WITH_EXP(d2);
  PRINT_VAL_WITH_EXP(d3);
  PRINT_VAL_WITH_EXP(d4);
  PRINT_VAL_WITH_EXP(l1);
  PRINT_VAL_WITH_EXP(l2);
  PRINT_VAL_WITH_EXP(l3);
  PRINT_VAL_WITH_EXP(m1);
  PRINT_VAL_WITH_EXP(m2);
  PRINT_VAL_WITH_EXP(m3);

  int constexpr dim = 25;
  int constexpr dimfixed=25;

  using LP = GaussSieve::ExactLatticePoint<mpz_class, dimfixed>;
  using GaussSieve::MaybeFixed;
  using GaussSieve::EMVApproximation;
  GaussSieve::StaticInitializer<LP> init1 (MaybeFixed<dimfixed>{dim});
  GaussSieve::StaticInitializer<EMVApproximation<dimfixed>> init2 (MaybeFixed<dimfixed>{dim});

  std::array<mpz_class,dim> arr;
  for(int i=0;i<dim;++i)
  {
    arr[i] = 400 *i * i;
  }

  LP latp = GaussSieve::make_from_any_vector<LP>(arr,MaybeFixed<dimfixed>{dim});

  EMVApproximation<dimfixed> emv(latp);

  std::cout << emv;


  return true;
}


#endif
