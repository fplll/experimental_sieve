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


  using LazyWrapS_CR  = LazyWrapCR<CombinedScalar>;
  using LazyWrapV_CR  = LazyWrapCR<CombinedVector,1>;


  LazyWrapS_CR wrap_scalar(combined_scalar);
  LazyWrapV_CR wrap_vector1(combined_vector1);

  std::cout << "Vector: " << wrap_vector1.eval<0>() << wrap_vector1.eval<1>() << std::endl;
  std::cout << "Scalar: " << wrap_scalar.eval<0>() << " approx. by " << wrap_scalar.eval<1>() << std::endl;


  // rvalue versions:

  using LazyWrapS_RV  = LazyWrapRV<CombinedScalar>;
  using LazyWrapV_RV  = LazyWrapRV<CombinedVector,1>;

  //copy everything:

  CombinedScalar combined_scalar_m(combined_scalar);
  CombinedVector combined_vector_m1( combined_vector1.access<0>().make_copy() );

  // rvalue wrappers.

  LazyWrapS_RV wrap_scalar_m(combined_scalar_m);
  LazyWrapV_RV wrap_vector_m1(combined_vector_m1);

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


  using IDNode1 = SieveLazyEval<IdentityFunVector,LazyWrapV_RV>;
  using IDNode2 = SieveLazyEval<IdentityFunVector,IDNode1>;
  IDNode1 lazy_id1{ std::move(wrap_vector_m1) };
  IDNode2 lazy_id2{ std::move(lazy_id1) };
  std::cout << "Lazy Eval:" << lazy_id2.eval<0>() << lazy_id2.eval<1>() << std::endl;

  using IDNode1S = SieveLazyEval<IdentityFunScalar, LazyWrapS_CR>;
  using IDNode2S = SieveLazyEval<IdentityFunScalar, IDNode1S>;
  IDNode1S lazy_id1_s{ wrap_scalar };
  IDNode1S lazy_id1_s_copy {lazy_id1_s};
  IDNode2S lazy_id2_s{ lazy_id1_s  };
  IDNode1S lazy_id3_s{ wrap_scalar };
  std::cout << "Compare:" << std::endl;
  std::cout << ( lazy_id2_s < lazy_id3_s ) << std::endl;
  std::cout << (abs(lazy_id2_s)< lazy_id2_s) << std::endl;
  std::cout << (abs(lazy_id2_s) <= lazy_id2_s) << std::endl;
  std::cout << (lazy_id2_s + lazy_id3_s).eval<0>() << std::endl;
  return true;
}








#endif
