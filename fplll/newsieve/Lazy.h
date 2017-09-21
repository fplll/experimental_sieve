#ifndef GAUSS_SIEVE_LAZY_H
#define GAUSS_SIEVE_LAZY_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include <tuple>


namespace GaussSieve{
namespace LazyEval{     // sub-namespace to inject free functions like abs


#define BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation) \
  using ExactVectorType = ELP; \
  using ExactScalarType = typename GetScalarProductStorageType<ELP>::type; \
  using ApproxVectorType = Approximation; \
  using ApproxScalarType = typename Approximation::ScalarProductType;

enum struct ScalarOrVector{ scalar_type, vector_type };


template<class ELP, class Approximation, class LazyFunction, class... Args> class SieveLazyEval;

template<class ELP, class Approximation, class LHS, class RHS> class Lazy_ScalarProduct;
template<class ELP, class Approximation, class Arg> class Lazy_Identity;


/*
  Encapsulates evaluation of a function F(x1,...), where x_1,... are of types Args...
  and the function F is encapsulated by a LazyFunction class.
  Available LazyFunction classes have names Lazy_F and adhere to a certain format.

  Args... is a tuple of types that are usually of the form SieveLazyEval again or are leaves of the
  expression tree.

  supports:
  eval_*
  assignment / conversion to exact / approximate
  operator<
*/

template<class ELP, class Approximation> class LazyWrapExactScalar
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation)
  using TreeType = ExactScalarType const &;
  using ExactEvalType = ExactScalarType const &;
  using ApproxEvalType = ApproxScalarType;
  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::scalar_type;

  constexpr LazyWrapExactScalar(ExactScalarType const &exact_scalar): args(exact_scalar) {};

  constexpr ExactEvalType eval_exact() { return args; }
  //[[deprecated]]
  constexpr ApproxEvalType eval_approx() { return static_cast< ApproxScalarType > (args); }

  TreeType const args;
};

template<class ELP, class Approximation> class LazyWrapExactAndApproxScalar
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation)
  using TreeType = std::tuple<ExactScalarType const &,ApproxScalarType const &>;
  using ExactEvalType = ExactScalarType const &;
  using ApproxEvalType = ApproxScalarType const &;

  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::scalar_type;

  constexpr LazyWrapExactAndApproxScalar(ExactScalarType const &exact_scalar, ApproxScalarType const &approx_scalar)
  : args(std::tie(exact_scalar,approx_scalar)  ) {}
  constexpr ExactEvalType eval_exact() const { return std::get<0>(args); }
  constexpr ApproxEvalType eval_approx() const { return std::get<1>(args); }

  TreeType const args;
};

template<class ELP, class Approximation> class LazyWrapExactVector
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation)
  using TreeType = ExactVectorType const &;
  using ExactEvalType = ExactVectorType const &;
  using ApproxEvalType = ApproxVectorType;

  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::vector_type;

  constexpr LazyWrapExactVector(ExactVectorType const &exact_vector): args(exact_vector) {};
  constexpr ExactEvalType eval_exact() const { return args; }
  //[[deprecated]]
  constexpr ApproxEvalType eval_approx() const { return static_cast< ApproxVectorType>(args); }


  TreeType const args;
};

template<class ELP, class Approximation> class LazyWrapExactAndApproxVector
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation)
  using TreeType = std::tuple<ExactVectorType const &, ApproxVectorType const &>;
  using ExactEvalType  = ExactVectorType const &;
  using ApproxEvalType = ApproxScalarType const &;

  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::vector_type;
  constexpr LazyWrapExactAndApproxVector(ExactVectorType const &exact_vector, ApproxVectorType const &approx_vector)
  :args(std::tie(exact_vector,approx_vector)) {}
  constexpr ExactVectorType const & eval_exact() { return std::get<0>(args); }
  constexpr ApproxVectorType const & eval_approx() { return std::get<1>(args); }

  TreeType const args;
};

template<class ELP, class Approximation, class Arg> class Lazy_Identity
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation)
  using ArgTree = typename Arg::TreeType;
  static constexpr int nargs = 1;
  static constexpr ScalarOrVector scalar_or_vector = Arg::scalar_or_vector;
  using ExactEvalType = typename Arg::ExactEvalType;
  using ApproxEvalType= typename Arg::ApproxEvalType;

  // C++14 decltype(auto) would really be helpful here...
  inline static ExactEvalType eval_exact( std::tuple<ArgTree const> const & arg)
  {
    return Arg(std::get<0>(arg)).eval_exact();
  }
  inline static ApproxEvalType eval_approx(std::tuple<ArgTree const> const & arg)
  {
    return Arg(std::get<0>(arg)).eval_approx();
  }
};

template<class ELP, class Approximation, class LHS, class RHS> class Lazy_ScalarProduct
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation)
  using ArgTreeLeft = typename LHS::TreeType;
  using ArgTreeRight= typename RHS::TreeType;
  static constexpr int nargs = 2;
  static_assert(LHS::scalar_or_vector == ScalarOrVector::vector_type,"Left hand side is no vector");
  static_assert(RHS::scalar_or_vector == ScalarOrVector::vector_type,"Right hand side is no vector");
  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::scalar_type;
  using ExactEvalType  = ExactScalarType;
  using ApproxEvalType = ApproxScalarType;

  inline static ExactEvalType eval_exact(std::tuple<ArgTreeLeft const, ArgTreeRight const> const & arg)
  {
    return compute_sc_product_exact(
      LHS(std::get<0>(arg)).eval_exact(),
      RHS(std::get<1>(arg)).eval_exact()
      );
  }
  inline static ApproxEvalType eval_approx(std::tuple<ArgTreeLeft const, ArgTreeRight const> const & arg)
  {
    return compute_sc_product_exact(
    LHS(std::get<0>(arg)).eval_approx(),
    RHS(std::get<1>(arg)).eval_approx()
    );
  }
};

template<class ELP, class Approximation, class LazyFunction, class... Args> class SieveLazyEval
{
  static_assert(sizeof...(Args) == LazyFunction::nargs, "Wrong number of arguments");
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation)
  using TreeType = std::tuple< typename Args::TreeType const...  >;
  using ExactEvalType = typename LazyFunction::ExactEvalType;
  using ApproxEvalType =typename LazyFunction::ApproxEvalType;
  static constexpr ScalarOrVector scalar_or_vector = LazyFunction::scalar_or_vector;

  TreeType const args;

  constexpr explicit SieveLazyEval(TreeType const & fn_args) : args(fn_args) {};

  ExactEvalType eval_exact() const { return LazyFunction::eval_exact(args); }
  ApproxEvalType eval_approx() const { return LazyFunction::eval_approx(args); }

  operator ExactEvalType() const { return eval_exact(); }
  operator ApproxEvalType() const { return eval_approx(); }

  bool operator< ( ExactScalarType const & rhs) const
  {
    static_assert(scalar_or_vector == ScalarOrVector::scalar_type,"comparing vector with scalar.");
    return eval_exact() < rhs;
  }

  template<class LazyFunctionRHS , class... ArgsRHS>
  bool operator< ( SieveLazyEval<ELP,Approximation,LazyFunctionRHS,ArgsRHS...> const & rhs) const
  {
    return eval_exact() < rhs.eval_exact();
  }
};





}} //end namespaces





#endif
