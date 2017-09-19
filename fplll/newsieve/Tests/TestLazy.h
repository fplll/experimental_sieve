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

//  using Combiner = LazyEval::Lazy_VectorFromExactAndApprox<ELP,Approx>;
//  using VectorWrapper  = LazyEval::SieveLazyEval<ELP,Approx,Combiner,ELP,Approx>;

  StaticInitializer<ELP> init_ELP (dim);
  StaticInitializer<Approx> init_Approx (dim);

  long A[dim] = {12,113,262,151,141,723,-1,25,35,0};
  ELP exact_point = make_from_any_vector<ELP>(A,dim);
  std::cout << exact_point << std::endl;
  Approx approx_point = static_cast<Approx>(exact_point);

  using LazyWrapperEV = LazyEval::LazyWrapExactVector<ELP,Approx>;
  using LazyWrapperBV = LazyEval::LazyWrapExactAndApproxVector<ELP,Approx>;
  LazyWrapperEV vector_wrapper(exact_point);
  LazyWrapperBV vector_wrapper2(exact_point,approx_point);

  std::cout << vector_wrapper.eval_exact() << std::endl;
  std::cout << vector_wrapper.eval_approx() << std::endl;

  std::cout << vector_wrapper2.eval_exact() << std::endl;
  std::cout << vector_wrapper2.eval_approx() << std::endl;

  using VecIdentity = LazyEval::Lazy_Identity<ELP,Approx,LazyWrapperEV>;

  std::cout << VecIdentity::eval_exact(std::tuple<typename LazyWrapperEV::TreeType>(vector_wrapper.args)) << std::endl;

  return true;
}








#endif
