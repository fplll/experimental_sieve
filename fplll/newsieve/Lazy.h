#ifndef GAUSS_SIEVE_LAZY_H
#define GAUSS_SIEVE_LAZY_H

#include "DebugAll.h"
#include <type_traits>
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include <iostream>
#include <utility>
#include <tuple>



namespace GaussSieve{
namespace LazyEval{     // sub-namespace to inject free functions like abs


template<class ELP, class Approximation, class LazyFunction, class... Args> class SieveLazyEval;

template<class ELP, class Approximation> class Lazy_ScalarProduct;
template<class ELP, class Approximation> class Lazy_FromExact;
template<class ELP, class Approximation> class Lazy_FromExactAndApprox;



template<class ELP, class Approximation> class Lazy_FromExactAndApprox
{
  public:
  using ExactScalarProductType = typename GetScalarProductStorageType<ELP>::type;
  using ApproxScalarProductType = typename Approximation::ScalarProductType;
  using DataType = std::tuple<ExactScalarProductType const &, ApproxScalarProductType const &>;

  static ExactScalarProductType eval_exact( DataType const & data) { return std::get<0>(data); }

  static ApproxScalarProductType eval_approx( DataType const &data) { return std::get<1>(data); }
};



template<class ELP, class Approximation, class LazyFunction, class... Args> class SieveLazyEval
{

  public:
  using ExactScalarProductType = typename GetScalarProductStorageType<ELP>::type;
  using ApproxScalarProductType = typename Approximation::ScalarProductType;

  constexpr explicit SieveLazyEval(Args const &... args) : function_operands(args...) {};

  auto eval_exact() const // If we ever use C++14, just remove the following line!
  -> decltype( LazyFunction::eval_exact(std::declval< std::tuple<Args const&...> const >() ) )
  {
    return LazyFunction::eval_exact(function_operands);
  }
  auto eval_approx() const
  //-> decltype( LazyFunction::eval_approx(std::declval<SieveLazyEval>().function_operands) )
  {
    return LazyFunction::eval_approx(function_operands);
  }

  operator ExactScalarProductType() const { return static_cast<ExactScalarProductType>( eval_exact() ); }
  bool operator< ( ExactScalarProductType const & rhs) const { return eval_exact() < rhs; }

  template<class LF2, class... Args2>
  bool operator< ( SieveLazyEval<ELP,Approximation,LF2,Args2...> const & rhs) const
  {
    return eval_exact() < rhs.eval_exact();
  }

  public:
  std::tuple<Args const &...> const function_operands;
//  std::tuple<Args const &...> get_ops() const { return function_operands; }
};





}} //end namespaces





#endif
