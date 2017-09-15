#ifndef TEST_LAZY_H
#define TEST_LAZY_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "gmpxx.h"
#include "../ExactLatticePoint.h"
#include "../EMVApproximation.h"

#include "../Lazy.h"

bool test_lazy()
{
  unsigned int constexpr dim = 10;
  using namespace GaussSieve;
  using ELP = ExactLatticePoint<long,dim>;
  using Approx = EMVApproximation<dim>;

  using Combiner = LazyEval::Lazy_FromExactAndApprox<ELP,Approx>;
  using Wrapper  = LazyEval::SieveLazyEval<ELP,Approx,Combiner,ELP,Approx>;

  StaticInitializer<ELP> init (dim);

  long A[dim] = {12,113,262,151,141,723,-1,25,35,0};
  ELP exact_point = make_from_any_vector<ELP>(A,dim);
  std::cout << exact_point;
  Approx approx_point = static_cast<Approx>(exact_point);

//  Wrapper w(exact_point,approx_point);


  return true;
}








#endif
