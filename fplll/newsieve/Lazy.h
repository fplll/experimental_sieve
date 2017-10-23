#ifndef GAUSS_SIEVE_LAZY_H
#define GAUSS_SIEVE_LAZY_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include <tuple>

/**
  This file provides functionality for lazy evaluation via expression templates.
  Notably, it provides a functionality, by which one can use expressions like
  A + B, ScalarProduct(A,B)
  whose return type is a "delayed addition of A,B" resp. "delayed scalar product of A,B".
  No actual computation is performed (at runtime). The encapsulation into a
  "delayed computation object" should be performed at compile time as to give no overhead
  (in theory, modulo bugs)
  Only when such a delayed object is _assigned_ to something is the computation performed.
  More generally, computation is performed whenever an appropriate conversion is (explicitly or
  implicitly requested). A particular case is comparison via <=, which triggers computation.

  The allowed data types for A,B are "vectors" or "scalars" or delayed expressions returning such.
  (this allows building up complex expression trees). A and B are supposed to support both an
  exact and an approximate mode.
  Consequently, all involved classes are template-parameterized by ELP and Approximation, where
  ELP is an (exact) lattice point class and
  Approximation is a class for approximations to a lattice point.
  The "scalar" types are obtained as the result types of scalar products of these classes.

  Note that currently, the implementation is not "eidetic".

  Note further that (even recursive) expression trees only ever store references to the *original*
  vectors / scalars involved in their creation. In particular, it is possible to store an
  delayed expression as long as all original objects are in scope.
  i.e. auto Res = A + B + C may actually be valid code, as long as references to A,B,C are valid.
  (This is in stark contrast to the way e.g. libgmp implements expression templates.)
*/


namespace GaussSieve{
namespace LazyEval{     // sub-namespace to inject free functions like abs

#define LAZY_FUNCTION \
  using ELP_t = ELP;    \
  using Approximation_t = Approximation; \
  using IsLazyFunction = std::true_type

// ELP is an exact lattic point class.
// Approximation is an approximation class.
// This #define just serves to bring the appropriate typedefs into scope.
// Last semicolon intentionally missing.
#define BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation) \
  using ExactVectorType = ELP; \
  using ExactScalarType = typename Get_ScalarProductStorageType<ELP>::type; \
  using ApproxVectorType = Approximation; \
  using ApproxScalarType = typename Approximation::ScalarProductType

// to differentiate between vectors and scalars. Mostly used for static_asserts to catch bugs.
enum struct ScalarOrVector{ scalar_type, vector_type };

// forward declarations:
template<class LazyFunction> class SieveLazyEval;
template<class ELP, class Approximation, class LHS, class RHS> class Lazy_ScalarProduct;
template<class ELP, class Approximation, class Arg> class Lazy_Identity;

/**
  SieveLazyEval is the main class of this module:

  SieveLazyEval<ELP,Approximation, LazyFunction, Args>
  encapsulates evaluation of a function F(x1,...), where x_1,... are of types Args...
  and the function F is encapsulated by a LazyFunction class.

  Available LazyFunction classes have names Lazy_F and adhere to a certain format (see below)

  Args... is a tuple of types that are usually of the form SieveLazyEval again or are leaves of the
  expression tree.

  args of type TreeType models the expression tree. Note that the template argument is a variadic
  template parameter (i.e. SieveLazyEval can take an arbitrary number of class arguments).
  By contrast, TreeType is a single class, defined as a tuple<Arg1, Arg2,...>. The individual
  Argi's are either again such TreeTypes or (for leaves) const-references to a scalar/vector.
  Note that this differs from a standard alternative implementation, where each Argi would be a
  reference to a Leaf / a reference to a TreeType (i.e. TreeType contains the whole expression tree
  *without* having to resolve internal references).
  This allows users to store expression templates (sometimes, at least).
  The downside may be some overhead due to copying of (sub-)trees, depending on compiler
  optimizations

  eval_exact and eval_approx cause exact resp. approximate evaluation.
  The returned types are announced as ExactEvalType resp. ApproxEvalType.
  scalar_or_vector indicates whether the result is a scalar or a vector.

  We support:
  eval_* to obtain exact/approximate evaluation
  assignment to / conversion to the appropriate type (same as eval_*)
  comparison operators

  Note that leaves are of different types, but expose nearly the same syntax.
*/

CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(IsLazyNode, std::true_type, Has_IsLazyNode);

template<class LazyFunction> class SieveLazyEval
{
  static_assert(LazyFunction::IsLazyFunction::value,"No Lazy Function");
//  static_assert(sizeof...(Args) == LazyFunction::nargs, "Wrong number of arguments");
  public:
  using ELP = typename LazyFunction::ELP_t;
  using Approximation = typename LazyFunction::Approximation_t;
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
//  using TreeType = std::tuple< typename Args::TreeType const...  >;
  using TreeType      = typename LazyFunction::TreeType;
  using ExactEvalType = typename LazyFunction::ExactEvalType;
  using ApproxEvalType =typename LazyFunction::ApproxEvalType;
  using IsLazyNode = std::true_type;
  static constexpr ScalarOrVector scalar_or_vector = LazyFunction::scalar_or_vector;

  TreeType const args;

  constexpr explicit SieveLazyEval(TreeType const & fn_args) : args(fn_args) {};
  SieveLazyEval() = delete;
  constexpr explicit SieveLazyEval(SieveLazyEval const &other) = default;
  constexpr explicit SieveLazyEval(SieveLazyEval && other) = default;
  constexpr SieveLazyEval& operator=(SieveLazyEval const &other) = default;
  constexpr SieveLazyEval& operator=(SieveLazyEval &&other) = default;
  // TODO:
//  template<class... Args>  constexpr explicit SieveLazyEval(Args const &... args)
//    : SieveLazyEval(std::tuple) {}

  inline ExactEvalType eval_exact() const { return LazyFunction::eval_exact(args); }
  inline ApproxEvalType eval_approx() const { return LazyFunction::eval_approx(args); }

  inline explicit operator ExactEvalType() const { return eval_exact(); }
  inline explicit operator ApproxEvalType() const { return eval_approx(); }

  inline bool operator< ( ExactScalarType const & rhs) const
  {
    static_assert(scalar_or_vector == ScalarOrVector::scalar_type,"comparing vector with scalar.");
    return eval_exact() < rhs;
  }

  template<class LazyRHS, TEMPL_RESTRICT_DECL2( Has_IsLazyNode<typename std::decay<LazyRHS>::type> )>
  inline bool operator< ( LazyRHS const & rhs) const
  {
    using RHSType = typename std::decay<LazyRHS>::type;
    static_assert(std::is_same<typename RHSType::ELP,ELP>::value, ""  );
    static_assert(std::is_same<typename RHSType::Approximation, Approximation>::value,"");
    return (eval_approx() < rhs.eval_approx() ) && (eval_exact() < rhs.eval_exact() );
  }

  inline bool operator> ( ExactScalarType const & rhs) const
  {
    static_assert(scalar_or_vector == ScalarOrVector::scalar_type,"comparing vector with scalar.");
    return eval_exact() > rhs;
  }

  template<class LazyRHS, TEMPL_RESTRICT_DECL2( Has_IsLazyNode<typename std::decay<LazyRHS>::type> )>
  inline bool operator> ( LazyRHS const & rhs) const
  {
    using RHSType = typename std::decay<LazyRHS>::type;
    static_assert(std::is_same<typename RHSType::ELP,ELP>::value, ""  );
    static_assert(std::is_same<typename RHSType::Approximation, Approximation>::value,"");
    return (eval_approx() > rhs.eval_approx() ) && (eval_exact() > rhs.eval_exact() );
  }



};

/**
  LazyWrap* model the leaves of the expression trees.
  These classes are just wrappers around a reference to a *
  and expose the same interface as SieveLazyEval.

  Note that the functionality is a bit more limited than that of SieveLazyEval (due to lazyness of
  the coder). As a workaround, you may use SieveLazyEval with an IdentityFunction.
*/

/**
  Wraps around (a reference to) an exact scalar.
*/

template<class ELP, class Approximation> class LazyWrapExactScalar
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  using TreeType = ExactScalarType const &;
  using IsLazyNode = std::true_type;
  using ExactEvalType = ExactScalarType const &;
  using ApproxEvalType = ApproxScalarType; // no reference here! We create a new object
  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::scalar_type;

  constexpr LazyWrapExactScalar(ExactScalarType const &exact_scalar): args(exact_scalar) {};

  inline constexpr ExactEvalType eval_exact() { return args; }
  //[[deprecated]]
  inline constexpr ApproxEvalType eval_approx() { return static_cast< ApproxScalarType > (args); }

  TreeType const args;
};

/**
  Wraps around a pair of Exact/Approx Scalar.
*/

template<class ELP, class Approximation> class LazyWrapExactAndApproxScalar
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  using TreeType = std::tuple<ExactScalarType const &,ApproxScalarType const &>;
  using IsLazyNode = std::true_type;
  using ExactEvalType = ExactScalarType const &;
  using ApproxEvalType = ApproxScalarType const &;

  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::scalar_type;

  constexpr explicit LazyWrapExactAndApproxScalar(ExactScalarType const &exact_scalar, ApproxScalarType const &approx_scalar)
  : args(std::tie(exact_scalar,approx_scalar)  ) {}
  constexpr LazyWrapExactAndApproxScalar(TreeType const &arg_wrap):args(arg_wrap) {};
  // These constexprs might require C++14
  inline CPP14CONSTEXPR ExactEvalType eval_exact() const { return std::get<0>(args); }
  inline CPP14CONSTEXPR ApproxEvalType eval_approx() const { return std::get<1>(args); }

  TreeType const args;
};

/**
  Wraps around an exact vector
*/

template<class ELP, class Approximation> class LazyWrapExactVector
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  using TreeType = ExactVectorType const &;
  using IsLazyNode = std::true_type;
  using ExactEvalType = ExactVectorType const &;
  using ApproxEvalType = ApproxVectorType; // No reference here!

  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::vector_type;

  constexpr LazyWrapExactVector(ExactVectorType const &exact_vector): args(exact_vector) {};
  constexpr ExactEvalType eval_exact() const { return args; }
  //[[deprecated]]
  constexpr ApproxEvalType eval_approx() const { return static_cast< ApproxVectorType>(args); }


  TreeType const args;
};

/**
  Wraps around a exact/approximate pair for a vector.
*/

template<class ELP, class Approximation> class LazyWrapExactAndApproxVector
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  using TreeType = std::tuple<ExactVectorType const &, ApproxVectorType const &>;
  using IsLazyNode = std::true_type;
  using ExactEvalType  = ExactVectorType const &;
  using ApproxEvalType = ApproxVectorType const &;

  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::vector_type;
  constexpr explicit LazyWrapExactAndApproxVector(ExactVectorType const &exact_vector, ApproxVectorType const &approx_vector)
  :args(std::tie(exact_vector,approx_vector)) {}
  constexpr LazyWrapExactAndApproxVector(TreeType const &arg_wrap): args(arg_wrap){};
  CPP14CONSTEXPR ExactVectorType const & eval_exact() { return std::get<0>(args); }
  CPP14CONSTEXPR ApproxVectorType const & eval_approx() { return std::get<1>(args); }

  TreeType const args;
};

/**
  Layz_F functions.
  These classes provide a unified interface to call various functions F on
  expression trees.
  The interface is as follows:
  Lazy_F<ELP,Approximation,Arg1,Arg2,...> has 2+nargs template arguments of class type, where Argi
  is the type of the ith argument (either SieveLazyEval or LazyWrap*)

  The class is static (i.e. cannot be instantiated)
  It provides the following public typedefs / static constants:
  static constexpr int nargs : number of arguments.
  static constexpr ScalarOrVector scalar_or_vector : Result of evaluation
  typename ExactEvalType : Result Type of Exact Evaluation
  typename ApproxEvalType: Result Type of Approximate Evaluation
  static functions eval_exact and eval_approx
    These take a single argument of Type std::tuple<arg_tree1,arg_tree2,...>
    where arg_treei is Argi::TreeType
    The functions perform the actual evaluation of F by unpacking the tuple,
    evaluating the subtrees and then actually performing F.
*/

/**
  Identity function. F(x) = x. Can be used to promote a
  LazyWrap* to a SieveLazyEval<ELP,Approxiomation,Lazy_Identity<...>,...>
  */

template<class ELP, class Approximation, class Arg> class Lazy_Identity
{
  static_assert(Arg::IsLazyNode::value,"Invalid arg");
  public:
  LAZY_FUNCTION;
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  using ArgTree = typename Arg::TreeType;
  using TreeType= std::tuple<ArgTree const>;
  static constexpr int nargs = 1;
  static constexpr ScalarOrVector scalar_or_vector = Arg::scalar_or_vector;
  using ExactEvalType = typename Arg::ExactEvalType;
  using ApproxEvalType= typename Arg::ApproxEvalType;
  Lazy_Identity(...) = delete;
  // C++14 decltype(auto) would really be helpful here...
  CPP14CONSTEXPR inline static ExactEvalType eval_exact( TreeType const & arg)
  {
    return Arg(std::get<0>(arg)).eval_exact();
  }
  CPP14CONSTEXPR inline static ApproxEvalType eval_approx(TreeType const & arg)
  {
    return Arg(std::get<0>(arg)).eval_approx();
  }
};

/**
  Scalar Product function. Delegates to
    compute_sc_product_exact resp. compute_sc_product_approx
  of its arguments.
*/

template<class ELP, class Approximation, class LHS, class RHS> class Lazy_ScalarProduct
{
  static_assert(LHS::IsLazyNode::value, "Left hand argument invalid");
  static_assert(RHS::IsLazyNode::value, "Right hand argument invalid");
  public:
  LAZY_FUNCTION;
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  using ArgTreeLeft = typename LHS::TreeType;
  using ArgTreeRight= typename RHS::TreeType;
  using TreeType = std::tuple<ArgTreeLeft const, ArgTreeRight const>;
  static constexpr int nargs = 2;
  static_assert(LHS::scalar_or_vector == ScalarOrVector::vector_type,"Left hand side is no vector");
  static_assert(RHS::scalar_or_vector == ScalarOrVector::vector_type,"Right hand side is no vector");
  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::scalar_type;
  using ExactEvalType  = ExactScalarType;
  using ApproxEvalType = ApproxScalarType;

  inline static ExactEvalType eval_exact(TreeType const & arg)
  {
    return compute_sc_product_exact(
      LHS(std::get<0>(arg)).eval_exact(),
      RHS(std::get<1>(arg)).eval_exact()
      );
  }
  inline static ApproxEvalType eval_approx(TreeType const & arg)
  {
    return compute_sc_product_exact(
    LHS(std::get<0>(arg)).eval_approx(),
    RHS(std::get<1>(arg)).eval_approx()
    );
  }
};

template<class ELP, class Approximation, class Arg> class Lazy_Norm2
{
  static_assert(Arg::IsLazyNode::value, "Arg is invalid");
  public:
  LAZY_FUNCTION;
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  using ArgTree = typename Arg::TreeType;
  using TreeType = std::tuple<ArgTree const>;
  static constexpr int nargs = 1;
  static_assert(Arg::scalar_or_vector == ScalarOrVector::vector_type,"Can only take norm2 of vectors");
  static constexpr ScalarOrVector scalar_or_vector = ScalarOrVector::scalar_type;
  using ExactEvalType = ExactScalarType;
  using ApproxEvalType = ApproxScalarType;

  inline static ExactEvalType eval_exact(TreeType const & arg)
  {
    return static_cast<ExactEvalType>(Arg(std::get<0>(arg)).eval_exact().get_norm2());
  }
  inline static ApproxEvalType eval_approx(TreeType const & arg)
  {
    return static_cast<ApproxEvalType>(Arg(std::get<0>(arg)).eval_approx().get_approx_norm() );
  }
};

}} //end namespaces

#endif
