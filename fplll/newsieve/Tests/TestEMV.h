#ifndef TEST_EMV_H
#define TEST_EMV_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "../EMVApproximation.h"
#include <limits>

bool test_emv()
{
  double xd = 0.6;
  int xi = 2352626;
  long xl = 2352627;
  mpz_class xm = 236;
  int xi2 = -65536;

  using GaussSieve::EMVScalarProduct;
  std::cout << "Size of approximated Norm2s is" << std::numeric_limits<typename GaussSieve::EMVApproximationTraits::ApproxNorm2Type>::digits << " bits." << std::endl;
  std::cout << "Size of approximated entries is" << std::numeric_limits<typename GaussSieve::EMVApproximationTraits::ApproxEntryType>::digits << " bits." << std::endl;

  EMVScalarProduct ad(xd);
  EMVScalarProduct ai = static_cast<EMVScalarProduct> (xi);
  EMVScalarProduct ai22(-1,16);
  EMVScalarProduct ai2 = static_cast<EMVScalarProduct> (xi2);
  EMVScalarProduct am = static_cast<EMVScalarProduct>(xm);
  EMVScalarProduct al(xl);

  std::cout << "Exact:" << xd << " approximated by " << ad << std::endl;
  std::cout << "Exact:" << xi << " approximated by " << ai << std::endl;
  std::cout << "Exact:" << xl << " approximated by " << al << std::endl;
  std::cout << "Exact:" << xm << " approximated by " << am << std::endl;
  std::cout << "Exact:" << xi2 << " approximated by " << ai2 << std::endl;

  mpz_class mpz1 = 1024;
  mpz_class mpz2 = 1025;
  mpz_class mpz3 =-1024;
  mpz_class mpz4 =-1025;

  EMVScalarProduct app1(mpz1);
  EMVScalarProduct app2(mpz2);
  EMVScalarProduct app3(mpz3);
  EMVScalarProduct app4(mpz4);

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

  return true;
}


#endif
