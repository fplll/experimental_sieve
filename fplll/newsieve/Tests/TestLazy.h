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

  using Combiner = LazyEval::Lazy_VectorFromExactAndApprox<ELP,Approx>;
  using VectorWrapper  = LazyEval::SieveLazyEval<ELP,Approx,Combiner,ELP,Approx>;

  StaticInitializer<ELP> init_ELP (dim);
  StaticInitializer<Approx> init_Approx (dim);

  long A[dim] = {12,113,262,151,141,723,-1,25,35,0};
  ELP exact_point = make_from_any_vector<ELP>(A,dim);
  std::cout << exact_point << std::endl;
  Approx approx_point = static_cast<Approx>(exact_point);
  auto both = std::tie(exact_point,approx_point);
  std::cout << Combiner::eval_exact(both) << std::endl;
  VectorWrapper w(exact_point,approx_point);

  using ScProdFun = LazyEval::Lazy_ScalarProduct<ELP,Approx,VectorWrapper,VectorWrapper>;
  using ScProdWrapper = LazyEval::SieveLazyEval<ELP,Approx, ScProdFun, VectorWrapper, VectorWrapper>;

  std::cout << ScProdWrapper(w,w).eval_exact() << std::endl;
//  std::cout << ScProdWrapper(w,w).eval_approx() << std::endl;






  return true;
}








#endif
