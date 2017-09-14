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
  using ExactScalarProductType = typename GetScPType<ELP>::type;
  using ApproxScalarProductType = typename Approximation::ScalarProductType;

  std::tuple<ExactScalarProductType const &, ApproxScalarProductType const &> const data;

};



template<class ELP, class Approximation, class LazyFunction, class... Args> class SieveLazyEval
{

  public:
  using ExactScalarProductType = typename GetScPType<ELP>::type;
  using ApproxScalarProductType = typename Approximation::ScalarProductType;

  constexpr explicit SieveLazyEval(Args const &... args) : function_operands(args...) {};

  auto eval_exact() const // If we ever use C++14, just remove the following line!
  -> decltype( LazyFunction::eval_exact(std::declval<SieveLazyEval>().function_operands) )
  {
    return LazyFunction::eval_exact(function_operands);
  }
  auto eval_approx() const
  -> decltype( LazyFunction::eval_approx(std::declval<SieveLazyEval>().function_operands) )
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

  private:
  std::tuple<Args const &...> const function_operands;
};





}} //end namespaces





#endif
