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
#include <typeinfo>

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
  long B[dim] = {34, 12,666, 50, 14,-20,50,61,89,-1000};
  ELP exact_point = make_from_any_vector<ELP>(A,dim);
  ELP second_exact_point = make_from_any_vector<ELP>(B,dim);
  Approx approx_point = static_cast<Approx>(exact_point);
  Approx second_approx_point = static_cast<Approx>(second_exact_point);

  std::cout << "-- Testing Lazy evalutation module with this point --" << std::endl;

  std::cout << exact_point << std::endl;
  std::cout << approx_point << std::endl;


  std::cout << std::endl<< std::endl << "-- Wrappers: --" << std::endl;

  using LazyWrapperEV = LazyEval::LazyWrapExactVector<ELP,Approx>;
  using LazyWrapperBV = LazyEval::LazyWrapExactAndApproxVector<ELP,Approx>;
  LazyWrapperEV vector_wrapper(exact_point);
  LazyWrapperEV second_vector_wrapper(second_exact_point);
  LazyWrapperBV vector_wrapper2(exact_point,approx_point);
  LazyWrapperBV second_vector_wrapper2(second_exact_point, second_approx_point);


  std::cout << vector_wrapper.eval_exact() << std::endl;
  std::cout << vector_wrapper.eval_approx() << std::endl;

  std::cout << second_vector_wrapper2.eval_exact() << std::endl;
  std::cout << second_vector_wrapper2.eval_approx() << std::endl << std::flush;

  std::cout <<std::endl << std::endl << "-- Direct calling of identity function: --" << std::endl;

  using VecIdentity = LazyEval::Lazy_Identity<ELP,Approx,LazyWrapperEV>;
  using VecIdentity2 = LazyEval::Lazy_Identity<ELP,Approx,LazyWrapperBV>;

  std::cout << VecIdentity::eval_exact(std::tuple<typename LazyWrapperEV::TreeType>(vector_wrapper.args)) << std::endl;
  std::cout << VecIdentity2::eval_exact(std::tuple<typename LazyWrapperBV::TreeType>(vector_wrapper2.args)) << std::endl<<std::flush;

  std::cout << std::endl << std::endl << "-- Calling though delayed evaluation object: --" << std::endl;

  using DelayedIdentity = LazyEval::SieveLazyEval<VecIdentity>;
  using DelayedIdentity2 = LazyEval::SieveLazyEval<VecIdentity2>;

  DelayedIdentity delayed_identity{  std::tuple<typename LazyWrapperEV::TreeType>(vector_wrapper.args )   };
  DelayedIdentity2 delayed_identity2{ std::tuple<typename LazyWrapperBV::TreeType>(vector_wrapper2.args)   };

  std::cout << delayed_identity.eval_exact() << std::endl;
  std::cout << delayed_identity.eval_approx() << std::endl;
  std::cout << delayed_identity2.eval_exact() << std::endl;
  std::cout << delayed_identity2.eval_approx() << std::endl << std::flush;

  std::cout << std::endl << std::endl << "-- Creating delayed scalar product object: --" << std::endl;

  using ScalarProdFun = LazyEval::Lazy_ScalarProduct<ELP,Approx,LazyWrapperEV, LazyWrapperBV>;
  using DelayedScProd = LazyEval::SieveLazyEval<ScalarProdFun>;

  DelayedScProd delayed_sc_prod{ std::tuple<typename LazyWrapperEV::TreeType,typename LazyWrapperBV::TreeType>(vector_wrapper.args,second_vector_wrapper2.args)};
//  std::cout << "BLAH" << std::endl;
  std::cout << "Exact Scalar Product is:" << delayed_sc_prod.eval_exact() << std::endl;
  std::cout << "Approximate Scalar Product is:"<< delayed_sc_prod.eval_approx() << std::endl;

  bool is_smaller = (delayed_sc_prod<170000);
  std::cout << "ScProduct is smaller than 170000:" << is_smaller << std::endl;
  bool is_larger  = (delayed_sc_prod>175000);
  std::cout << "ScProduct is larger than 175000:"  << is_larger << std::endl;

  return true;
}








#endif
