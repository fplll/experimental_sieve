#define USE_REGULAR_QUEUE //priority queue not working atm.

#define DEBUG_SIEVE_SILENT_ALL
#define DEBUG_SIEVE_INITIALIZATIONS
#define DEBUG_SIEVE_STANDALONE_MODULES_ALL


// very verbose...
//#define DEBUG_SIEVE_LP_TRACEGENERIC

//#include "SieveGauss_main.h"
#include <iostream>

#include "Tests/TestTraitClasses.h"
#include "Tests/TestPlainLatticePoint.h"
#include "Tests/TestMaybeFixedInt.h"
#include "Tests/TestExactLatticePoint.h"
#include "Tests/TestMTPRNG.h"
#include "Tests/TestBasisUtils.h"
#include "Tests/TestShiSampler.h"
#include "Tests/TestQueue.h"
#include "Tests/TestList.h"

// clang-reorder-guard
#define ASSADGDGSDGKHSDLHEFOIEGFOIGFOSGFVOSGVFSAOPIGFASPOFGAPOFG

#include "ShiSampler_impl.h"
#include "Sampler_impl.h"
#include "GaussQueue_impl.h"

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

  if (test_basis_utils() )
  {
    std::cout << "Basis utilities work as expected." << std::endl;
  }

  if (test_shi_sampler() )
  {
    std::cout << "Sampler works as expected." << std::endl;
  }

  if (test_queue() )
  {
    std::cout << "Queue works as expected." << std::endl;
  }

  if (test_list() )
    std::cout << "List works as expected." << std::endl;
  {

  }

  return 0; // indicating success.

}

