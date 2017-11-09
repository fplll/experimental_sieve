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
  using namespace GaussSieve::LazyEval;
  using ELP = ExactLatticePoint<long,dim>;
  using Approx = EMVApproximation<dim>;

  using ExactScalar = long;
  using ApproxScalar = EMVScalar;
  using CombinedScalar = ObjectWithApproximation<ExactScalar,ApproxScalar>;
  using CombinedVector = ObjectWithApproximation<ELP,Approx>;

//  using Combiner = LazyEval::Lazy_VectorFromExactAndApprox<ELP,Approx>;
//  using VectorWrapper  = LazyEval::SieveLazyEval<ELP,Approx,Combiner,ELP,Approx>;

  StaticInitializer<ELP> init_ELP (dim);
  StaticInitializer<Approx> init_Approx (dim);
  long A[dim] = {12,113,262,151,141,723,-1,25,35,0};
  long B[dim] = {34, 12,666, 50, 14,-20,50,61,89,-1000};
  long C = 8000000L;

  std::cout << std::endl << "-- Creating vectors, approximations and combination -- " << std::endl << std::flush;

  ELP exact_point1 = make_from_any_vector<ELP>(A,dim);
  ELP exact_point2 = make_from_any_vector<ELP>(B,dim);
  Approx approx_point1 = static_cast<Approx>(exact_point1);
  Approx approx_point2(exact_point2);
  ExactScalar exact_scalar1 = C;
  ApproxScalar approx_scalar1 (exact_scalar1);
//  ExactScalar exact_scalar2 = -24623682362L;
  CombinedScalar combined_scalar(C);
  CombinedVector combined_vector1(exact_point1.make_copy());

  // swapping twice to check for std::move-correctness.
  ELP exact_vector_swp (std::move(combined_vector1));
  combined_vector1 = static_cast<CombinedVector>(std::move(exact_vector_swp));

  Approx approx_vector_tmp(exact_point2);
  CombinedVector combined_vector2(exact_point2.make_copy(),std::move(approx_vector_tmp));

  std::cout << exact_point1 << approx_point1 << std::endl;
  std::cout << exact_point2 << approx_point2 << std::endl;
  std::cout << combined_vector1.access<0>() << combined_vector1.access<1>() << std::endl;

  std::cout << exact_scalar1 << " approx. by " << approx_scalar1 << std::endl;
  std::cout << combined_scalar.access_exact() << " approx. by " << combined_scalar.access_approx() << std::endl;

  std::cout << std::endl << "-- Wrappers --" << std::endl << std::flush;


  using LazyWrapS  = LazyWrapCR<CombinedScalar>;
  using LazyWrapV  = LazyWrapCR<CombinedVector,1>;


  LazyWrapS wrap_scalar(combined_scalar);
  LazyWrapV wrap_vector1(combined_vector1);

  std::cout << "Vector: " << wrap_vector1.eval<0>() << wrap_vector1.eval<1>() << std::endl;
  std::cout << "Scalar: " << wrap_scalar.eval<0>() << " approx. by " << wrap_scalar.eval<1>() << std::endl;


  // rvalue versions:

  using LazyWrapS2  = LazyWrapRV<CombinedScalar>;
  using LazyWrapV2  = LazyWrapRV<CombinedVector,1>;

  //copy everything:

  CombinedScalar combined_scalar_m(combined_scalar);
  CombinedVector combined_vector_m1( combined_vector1.access<0>().make_copy() );

  // rvalue wrappers.

  LazyWrapS2 wrap_scalar_m(combined_scalar_m);
  LazyWrapV2 wrap_vector_m1(combined_vector_m1);

  // Note that from a specification POV, the arguments to the constructors are now in a possibly
  // invalid state, so we won't use them anymore.

  std::cout << std::endl << "-- Direct calling of identity function: --" << std::endl << std::flush;

  using IdentityFunScalar = Lazy_Identity<CombinedScalar>;
  using IdentityFunVector = Lazy_Identity<CombinedVector>;

  std::cout << IdentityFunScalar::call<0>(exact_scalar1) << " approx. by " << IdentityFunScalar::call<1>(approx_scalar1) << std::endl;
  std::cout << IdentityFunVector::call<0>(exact_point1 ) << IdentityFunVector::call<1>(approx_point1)  << std::endl;
  std::cout << IdentityFunScalar::call<0>(combined_scalar.access<0>())  << "approx. by " << IdentityFunScalar::call<1>(combined_scalar.access<1>()) << std::endl;
  std::cout << IdentityFunVector::call<0>(combined_vector1.access<0>()) << IdentityFunVector::call<1>(combined_vector1.access<1>()) << std::endl;

  std::cout << std::endl << "-- Lazyly calling identity function: --" << std::endl << std::flush;

  using IDNode1 = SieveLazyEval<IdentityFunVector,LazyWrapV2>;
  using IDNode2 = SieveLazyEval<IdentityFunVector,IDNode1>;
  IDNode1 lazy_id1{ std::move(wrap_vector_m1) };
  IDNode2 lazy_id2{ std::move(lazy_id1) };
  std::cout << "Lazy Eval:" << lazy_id2.eval<0>() << lazy_id2.eval<1>() << std::endl;



//  using IdentityFunES = Lazy_Identity<ELP,Approx,ScalarOrVector::scalar_type>;
//  using IdentityFunBS = Lazy_Identity<ELP,Approx,ScalarOrVector::scalar_type>;
//  using IdentityFunCS = Lazy_Identity<ELP,Approx,ScalarOrVector::scalar_type>;
//  using IdentityFunEV = Lazy_Identity<ELP,Approx,ScalarOrVector::vector_type>;
//  using IdentityFunBV = Lazy_Identity<ELP,Approx,ScalarOrVector::vector_type>;
//  using IdentityFunCV = Lazy_Identity<ELP,Approx,ScalarOrVector::vector_type>;

//  std::cout << IdentityFunBS::eval_exact (std::tuple<typename LazyWrapBS::TreeType>(wrap_scalar1.args)) << " approx. "
//            << IdentityFunBS::eval_approx(std::tuple<typename LazyWrapBS::TreeType>(wrap_scalar1.args)) << std::endl;
//  std::cout << IdentityFunES::eval_exact (std::tuple<typename LazyWrapES::TreeType>(wrap_scalar2.args)) << " approx. "
//            << IdentityFunES::eval_approx(std::tuple<typename LazyWrapES::TreeType>(wrap_scalar2.args)) << std::endl;
//  std::cout << IdentityFunCS::eval_exact (std::tuple<typename LazyWrapCS::TreeType>(wrap_scalar3.args)) << " approx. "
//            << IdentityFunCS::eval_approx(std::tuple<typename LazyWrapCS::TreeType>(wrap_scalar3.args)) << std::endl;
//  std::cout << IdentityFunEV::eval_exact (std::tuple<typename LazyWrapEV::TreeType>(wrap_vector1.args))
//            << IdentityFunEV::eval_approx(std::tuple<typename LazyWrapEV::TreeType>(wrap_vector1.args)) << std::endl;
//  std::cout << IdentityFunBV::eval_exact (std::tuple<typename LazyWrapBV::TreeType>(wrap_vector2.args))
//            << IdentityFunBV::eval_approx(std::tuple<typename LazyWrapBV::TreeType>(wrap_vector2.args)) << std::endl;
//  std::cout << IdentityFunCV::eval_exact (std::tuple<typename LazyWrapCV::TreeType>(wrap_vector3.args))
//            << IdentityFunCV::eval_approx(std::tuple<typename LazyWrapCV::TreeType>(wrap_vector3.args)) << std::endl;



//
//  using VecIdentity = LazyEval::Lazy_Identity<ELP,Approx,LazyWrapperEV>;
//  using VecIdentity2 = LazyEval::Lazy_Identity<ELP,Approx,LazyWrapperBV>;
//
//  std::cout << VecIdentity::eval_exact(std::tuple<typename LazyWrapperEV::TreeType>(vector_wrapper.args)) << std::endl;
//  std::cout << VecIdentity2::eval_exact(std::tuple<typename LazyWrapperBV::TreeType>(vector_wrapper2.args)) << std::endl<<std::flush;
//
//  std::cout << std::endl << std::endl << "-- Calling though delayed evaluation object: --" << std::endl;
//
//  using DelayedIdentity = LazyEval::SieveLazyEval<VecIdentity>;
//  using DelayedIdentity2 = LazyEval::SieveLazyEval<VecIdentity2>;
//
//  DelayedIdentity delayed_identity{  std::tuple<typename LazyWrapperEV::TreeType>(vector_wrapper.args )   };
//  DelayedIdentity2 delayed_identity2{ std::tuple<typename LazyWrapperBV::TreeType>(vector_wrapper2.args)   };
//
//  std::cout << delayed_identity.eval_exact() << std::endl;
//  std::cout << delayed_identity.eval_approx() << std::endl;
//  std::cout << delayed_identity2.eval_exact() << std::endl;
//  std::cout << delayed_identity2.eval_approx() << std::endl << std::flush;
//
//  std::cout << std::endl << std::endl << "-- Creating delayed scalar product object: --" << std::endl;
//
//  using ScalarProdFun = LazyEval::Lazy_ScalarProduct<ELP,Approx,LazyWrapperEV, LazyWrapperBV>;
//  using DelayedScProd = LazyEval::SieveLazyEval<ScalarProdFun>;
//
//  DelayedScProd delayed_sc_prod{ std::tuple<typename LazyWrapperEV::TreeType,typename LazyWrapperBV::TreeType>(vector_wrapper.args,second_vector_wrapper2.args)};
////  std::cout << "BLAH" << std::endl;
//  std::cout << "Exact Scalar Product is:" << delayed_sc_prod.eval_exact() << std::endl;
//  std::cout << "Approximate Scalar Product is:"<< delayed_sc_prod.eval_approx() << std::endl;
//
//  bool is_smaller = (delayed_sc_prod<170000);
//  std::cout << "ScProduct is smaller than 170000:" << is_smaller << std::endl;
//  bool is_larger  = (delayed_sc_prod>175000);
//  std::cout << "ScProduct is larger than 175000:"  << is_larger << std::endl;


  return true;
}








#endif
