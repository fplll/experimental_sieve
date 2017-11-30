#define USE_REGULAR_QUEUE //priority queue not working atm.

#define DEBUG_SIEVE_SILENT_ALL
#define DEBUG_SIEVE_INITIALIZATIONS
#define DEBUG_SIEVE_STANDALONE_MODULES_ALL

//#define TEST_ALL
#define TEST_EMV
//#define TEST_LAZY
//#define TEST_APPROXIMATIONS
//#define TEST_PLAIN_LATTICE_POINT
//#define TEST_EXACT_LATTICE_POINT
#define TEST_BITAPPROX
//#define TEST_RELEVANT_COORDS

#ifdef TEST_ALL
  #define TEST_TRAIT_CLASSES
  #define TEST_PLAIN_LATTICE_POINT
  #define TEST_MAYBE_FIXED_INT
  #define TEST_EXACT_LATTICE_POINT
  #define TEST_MTPRNG
  #define TEST_BASIS_UTILS
  #define TEST_SHI_SAMPLER
  #define TEST_QUEUE
  #define TEST_LIST
  #define TEST_EMV
  #define TEST_BITAPPROX
  #define TEST_APPROXIMATIONS
  #define TEST_LAZY
#endif

// very verbose...
//#define DEBUG_SIEVE_LP_TRACEGENERIC

//#include "SieveGauss_main.h"
#include <iostream>

#ifdef TEST_TRAIT_CLASSES
  #include "Tests/TestTraitClasses.h"
#endif

#ifdef TEST_PLAIN_LATTICE_POINT
  #include "Tests/TestPlainLatticePoint.h"
#endif

#ifdef TEST_MAYBE_FIXED_INT
  #include "Tests/TestMaybeFixedInt.h"
#endif

#ifdef TEST_EXACT_LATTICE_POINT
  #include "Tests/TestExactLatticePoint.h"
#endif

#ifdef TEST_MTPRNG
  #include "Tests/TestMTPRNG.h"
#endif

#ifdef TEST_BASIS_UTILS
  #include "Tests/TestBasisUtils.h"
#endif

#ifdef TEST_SHI_SAMPLER
  #include "Tests/TestShiSampler.h"
#endif

#ifdef TEST_QUEUE
  #include "Tests/TestQueue.h"
#endif

#ifdef TEST_LIST
  #include "Tests/TestList.h"
#endif

#ifdef TEST_EMV
  #include "Tests/TestEMV.h"
#endif

#ifdef TEST_APPROXIMATIONS
  #include "Tests/TestApproximations.h"
#endif

#ifdef TEST_LAZY
  #include "Tests/TestLazy.h"
#endif

#ifdef TEST_BITAPPROX
  #include "Tests/TestBitApprox.h"
#endif

#ifdef TEST_RELEVANT_COORDS
  #include "Tests/TestRelevantCoords.h"
#endif

// clang-reorder-guard
#define ASSADGDGSDGKHSDLHEFOIEGFOIGFOSGFVOSGVFSAOPIGFASPOFGAPOFG

#ifdef TEST_SHI_SAMPLER
  #include "ShiSampler_impl.h"
  #include "Sampler_impl.h"
#endif

#ifdef TEST_QUEUE
  #include "GaussQueue_impl.h"
#endif

int main(int argc, char **argv)
{
#ifdef TEST_TRAIT_CLASSES
  if (test_trait_classes())
  {
    std::cout << "Trait utility macros work as expected." << std::endl;
  }
#endif
#ifdef TEST_MAYBE_FIXED_INT
  if (test_maybe_fixed_int())
  {
    std::cout << "Tests for possibly-compile-time-fixed ints worked." << std::endl;
  }
#endif
#ifdef TEST_PLAIN_LATTICE_POINT
  if (test_plain_LP())
  {
    std::cout << "Plain Lattice Points work as expected." << std::endl;
  }
#endif
#ifdef TEST_EXACT_LATTICE_POINT
  if (test_exact_LP() )
  {
    std::cout << "Exact Lattice Points work as expected." << std::endl;
  }
#endif
#ifdef TEST_MTPRNG
  if (test_mtprng() )
  {
    std::cout << "MTPRNG work as expected." << std::endl;
  }
#endif
#ifdef TEST_BASIS_UTILS
  if (test_basis_utils() )
  {
    std::cout << "Basis utilities work as expected." << std::endl;
  }
#endif

#ifdef TEST_SHI_SAMPLER
  if (test_shi_sampler() )
  {
    std::cout << "Sampler works as expected." << std::endl;
  }
#endif

#ifdef TEST_QUEUE
  if (test_queue() )
  {
    std::cout << "Queue works as expected." << std::endl;
  }
#endif

#ifdef TEST_LIST
  if (test_list() )
  {
    std::cout << "List works as expected." << std::endl;
  }
#endif

#ifdef TEST_EMV
  if (test_emv() )
  {
    std::cout << "EMV Approximation works as expected." << std::endl;
  }
#endif

#ifdef TEST_APPROXIMATIONS
  if( test_approximations() )
  {
    std::cout << "Approximations combiner works as expected" << std::endl;
  }
#endif

#ifdef TEST_LAZY
  if (test_lazy() )
  {
    std::cout << "Lazy Evaluations work as expected" << std::endl;
  }
#endif

#ifdef TEST_BITAPPROX
  if (test_bit_approx())
  {
    std::cout << "Bit Approximation works as expected" <<std::endl;
  }
#endif

#ifdef TEST_RELEVANT_COORDS
  if(test_relevant_coords())
  {

  }
#endif
  return 0; // indicating success.

}
