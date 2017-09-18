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

template<class ELP, class Approximation> // Lattice Point, Approximation
struct GetType
{
  using ExactVectorType = ELP;
  using ExactScalarType = typename GetScalarProductStorageType<ELP>::type;
  using ApproxVectorType = Approximation;
  using ApproxScalarType = typename Approximation::ScalarProductType;
};


template<class ELP, class Approximation, class LazyFunction, class... Args> class SieveLazyEval;

template<class ELP, class Approximation> class Lazy_ScalarProduct;
template<class ELP, class Approximation> class Lazy_FromExact;
template<class ELP, class Approximation> class Lazy_ScalarFromExactAndApprox;



template<class ELP, class Approximation> class Lazy_ScalarFromExactAndApprox
{
  public:
  using ExactScalarType =  typename GetType<ELP,Approximation>::ExactScalarType;
  using ApproxScalarType = typename GetType<ELP,Approximation>::ApproxScalarType;
  using DataType = std::tuple<ExactScalarType const &, ApproxScalarType const &>;

  static ExactScalarType const & eval_exact( DataType const & data) { return std::get<0>(data); }
  static ApproxScalarType const & eval_approx( DataType const &data) { return std::get<1>(data); }
};

template<class ELP, class Approximation> class Lazy_VectorFromExactAndApprox
{
  public:
  using ExactVectorType =  typename GetType<ELP,Approximation>::ExactVectorType;
  using ApproxVectorType = typename GetType<ELP,Approximation>::ApproxVectorType;
  using DataType = std::tuple<ExactVectorType const &, ApproxVectorType const &>;
  static ExactVectorType const & eval_exact( DataType const & data ) { return std::get<0>(data); }
  static ApproxVectorType const & eval_approx( DataType const & data ) { return std::get<1>(data); }
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
  -> decltype( LazyFunction::eval_approx(std::declval< std::tuple<Args const&...> const >() ) )
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
