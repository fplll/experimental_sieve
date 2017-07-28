// clang-format off

#define USE_REGULAR_QUEUE //priority queue not working atm.

#define DEBUG_SIEVE_SILENT_ALL

// very verbose...
//#define DEBUG_SIEVE_LP_TRACEGENERIC

//#include "SieveGauss_main.h"
#include <iostream>

#include "Tests/TestTraitClasses.h"
#include "Tests/TestPlainLatticePoint.h"
#include "Tests/TestMaybeFixedInt.h"
#include "Tests/TestExactLatticePoint.h"
#include "Tests/TestMTPRNG.h"

int main(int argc, char **argv)
{
  if (test_trait_classes())
  {
    std::cout << "Trait utility macros work as expected." << std::endl;
  }
  if (test_maybe_fixed_int())
  {
    std::cout << "Tests for possibly-compile-time-fixed ints worked." << std::endl;
  }
  if (test_plain_LP())
  {
    std::cout << "Plain Lattice Points work as expected." << std::endl;
  }
  if (test_exact_LP() )
  {
    std::cout << "Exact Lattice Points work as expected." << std::endl;
  }
  if (test_mtprng() )
  {
    std::cout << "MTPRNG work as expected." << std::endl;
  }

  return 0; // indicating success.

}

