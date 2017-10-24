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

CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(IsLazyNode, std::true_type, Has_IsLazyNode);

// forward declarations:
template<class LazyFunction, class... Args> class SieveLazyEval;
template<class ELP, class Approximation, class LHS, class RHS> class Lazy_ScalarProduct;
//template<class ELP, class Approximation, class Arg> class Lazy_Identity;
template<class ELP, class Approximation, ScalarOrVector s_or_v> class Lazy_Identity;


// combines and *stores* a scalar/vector with an approximation thereof.
template<class ELP, class Approximation> struct ScalarWithApproximation;
template<class ELP, class Approximation> struct VectorWithApproximation;

/**
  This class stores a scalar with an approximation to that scalar.
  Note that this is for storage only. We have no arithmetic etc.
  Such functionality is defered to ApproximatedPoint.h
*/
template<class ELP, class Approximation>
struct ScalarWithApproximation
{
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  using ExactType = ExactScalarType;
  using ApproxType = ApproxScalarType;

  constexpr explicit ScalarWithApproximation(ExactScalarType const &exact, ApproxScalarType const &approx)
    :exact_scalar(exact), approx_scalar(approx) {};
  constexpr explicit ScalarWithApproximation(ExactScalarType const &exact)
    :exact_scalar(exact), approx_scalar(static_cast<ApproxScalarType>(exact)){};

  constexpr explicit operator ExactScalarType() const { return exact_scalar;}
  constexpr explicit operator ApproxScalarType() const { return approx_scalar; }

  constexpr ExactScalarType const & access_exact() const { return exact_scalar; }
  constexpr ExactScalarType & access_exact() { return exact_scalar; }
  constexpr ApproxScalarType const & access_approx() const { return approx_scalar;}
  constexpr ApproxScalarType & access_approx() { return approx_scalar; }

  // const-ness restriction is for debug purposes, mostly.
  ExactScalarType  exact_scalar;
  ApproxScalarType approx_scalar;


//  template<class LazyObject, TEMPL_RESTRICT_DECL2(Has_IsLazyNode<typename std::decay<LazyObject>::type> )>
//  constexpr explicit ScalarWithApproximation(LazyObject &&lazy_object)
//    :exact_scalar(lazy_object.eval_exact()),approx_scalar(lazy_object.eval_approx()) //todo: move semantics
//  {
//    using Arg = typename std::decay<LazyObject>::type;
//    static_assert(Arg::scalar_or_vector == ScalarOrVector::scalar_type,"Trying to assign a vector to a scalar");
//  }
};

template<class ELP, class Approximation>
struct VectorWithApproximation
{
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  using ExactType = ELP;
  using ApproxType= Approximation;
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  constexpr explicit VectorWithApproximation(ExactVectorType &&exact, ApproxVectorType &&approx)
    :exact_vector(std::move(exact)), approx_vector(std::move(approx)) {}
  // Note: This also captures the mixed cases of && / const &.
  constexpr explicit VectorWithApproximation(ExactVectorType const &exact, ApproxVectorType const &approx) = delete;

  constexpr explicit VectorWithApproximation(ExactVectorType const &exact) =delete;
  constexpr explicit VectorWithApproximation(ExactVectorType &&exact)
    :exact_vector(std::move(exact)), approx_vector(static_cast<ApproxVectorType>(exact)) {}

//  constexpr explicit operator ExactVectorType const & () const & { return exact_vector;}
  constexpr explicit operator ExactVectorType() const & { return exact_vector;}
  constexpr explicit operator ExactVectorType() && { return std::move(exact_vector); }
  constexpr explicit operator ApproxVectorType() const & { return approx_vector;}
  constexpr explicit operator ApproxVectorType() && { return std::move(approx_vector); }
  constexpr ExactVectorType const & access_exact() const { return exact_vector; }
  constexpr ExactVectorType & access_exact() { return exact_vector; }
  constexpr ApproxVectorType const & access_approx() const { return approx_vector;}
  constexpr ApproxVectorType & access_approx() { return approx_vector; }

  ExactVectorType exact_vector;
  ApproxVectorType approx_vector;


//  template<class LazyObject, TEMPL_RESTRICT_DECL2(Has_IsLazyNode<typename std::decay<LazyObject>::type>)>
//  constexpr explicit VectorWithApproximation(LazyObject &&lazy_object)
//    :exact_vector(lazy_object.eval_exact()), approx_vector(lazy_object.eval_approx()) // move semantics are important here!
//  {
//    using Arg = typename std::decay<LazyObject>::type;
//    static_assert(Arg::scalar_or_vector == ScalarOrVector::vector_type,"Trying to assign a scalar to a vector");
//  }
};

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

template<class LazyFunction, class... Args> class SieveLazyEval
{
  static_assert(LazyFunction::IsLazyFunction::value,"No Lazy Function");
  static_assert(sizeof...(Args) == LazyFunction::nargs, "Wrong number of arguments");
  static_assert(MyConjunction< Has_IsLazyNode<Args>... >::value,"Some argument is wrong.");

  public:
  using ELP = typename LazyFunction::ELP_t;
  using Approximation = typename LazyFunction::Approximation_t;
  // assert all ELP's / Approx's of Arg... are the same?
  // Actually, this may be too restrictive.
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
//  using TreeType = std::tuple< typename Args::TreeType const...  >;
//  using TreeType      = typename LazyFunction::TreeType;
  using ExactEvalType = typename LazyFunction::ExactEvalType;
  using ApproxEvalType =typename LazyFunction::ApproxEvalType;
  using IsLazyNode = std::true_type;
  static constexpr ScalarOrVector scalar_or_vector = LazyFunction::scalar_or_vector;
  //TODO: ensure consistency.

  using TreeType = std::tuple<Args...>;
  TreeType args; //const?
  SieveLazyEval() = delete;
  constexpr explicit SieveLazyEval(SieveLazyEval const &other) = default;
  constexpr explicit SieveLazyEval(SieveLazyEval && other) = default;
  constexpr SieveLazyEval& operator=(SieveLazyEval const &other) = default;
  constexpr SieveLazyEval& operator=(SieveLazyEval &&other) = default;

  //  TreeType const args;

  constexpr explicit SieveLazyEval(TreeType const & fn_args) : args(fn_args)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy function wrapper object (via tuple) for function " << LazyFunction::fun_name() << std::endl;
#endif
  }
  constexpr explicit SieveLazyEval(TreeType && fn_args) : args(std::move(fn_args))
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy function wrapper object (via tuple-move) for function " << LazyFunction::fun_name() << std::endl;
#endif
  }
  // TODO: Add direct init via Args...

  template<std::size_t... iarg>
  inline ExactEvalType do_eval_exact(MyIndexSeq<iarg...>) const
  {
    static_assert(sizeof...(iarg) == sizeof...(Args),"Something is very wrong");
#ifdef DEBUG_SIEVE_LAZY_TRACE_EVALS
    std::cout << "Calling function exactly. Function is " << LazyFunction::fun_name() << std::endl;
#endif
    return LazyFunction::call_exact( std::get<iarg>(args).eval_exact()... );
  }

  inline ExactEvalType eval_exact() const
  {
    return do_eval_exact(MyMakeIndexSeq<sizeof...(Args)>{} );
  }

  inline explicit operator ExactEvalType() const { return eval_exact(); }
//  inline explicit operator ApproxEvalType() const { return eval_approx(); }

//  inline bool operator< ( ExactScalarType const & rhs) const
//  {
//    static_assert(scalar_or_vector == ScalarOrVector::scalar_type,"comparing vector with scalar.");
//    return eval_exact() < rhs;
//  }
//
//  template<class LazyRHS, TEMPL_RESTRICT_DECL2( Has_IsLazyNode<typename std::decay<LazyRHS>::type> )>
//  inline bool operator< ( LazyRHS const & rhs) const
//  {
//    using RHSType = typename std::decay<LazyRHS>::type;
//    static_assert(std::is_same<typename RHSType::ELP,ELP>::value, ""  );
//    static_assert(std::is_same<typename RHSType::Approximation, Approximation>::value,"");
//    return (eval_approx() < rhs.eval_approx() ) && (eval_exact() < rhs.eval_exact() );
//  }
//
//  inline bool operator> ( ExactScalarType const & rhs) const
//  {
//    static_assert(scalar_or_vector == ScalarOrVector::scalar_type,"comparing vector with scalar.");
//    return eval_exact() > rhs;
//  }
//
//  template<class LazyRHS, TEMPL_RESTRICT_DECL2( Has_IsLazyNode<typename std::decay<LazyRHS>::type> )>
//  inline bool operator> ( LazyRHS const & rhs) const
//  {
//    using RHSType = typename std::decay<LazyRHS>::type;
//    static_assert(std::is_same<typename RHSType::ELP,ELP>::value, ""  );
//    static_assert(std::is_same<typename RHSType::Approximation, Approximation>::value,"");
//    return (eval_approx() > rhs.eval_approx() ) && (eval_exact() > rhs.eval_exact() );
//  }
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

template<class ExactType, class ApproxType>
class LazyWrapExactCR
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType const &;
  using ApproxEvalType = ApproxType;
  using MayInvalidateExact = std::false_type;
  using MayInvalidateApprox= std::false_type;
  constexpr LazyWrapExactCR(ExactType const &init_exact) : exact_value(init_exact)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for exact values." << std::endl;
#endif
  };
  LazyWrapExactCR(ExactType const &&) = delete;

  inline constexpr ExactEvalType  eval_exact()  const { return exact_value; }
  inline constexpr ApproxEvalType eval_approx() const
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_APPROX
    std::cout << "Computing approximation inside wrapper, as requested" << std::endl;
#endif
    return static_cast<ApproxType>(exact_value);
  }

  ExactType const & exact_value;
};

template<class ExactType, class ApproxType>
class LazyWrapExactRV
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType &&; // !!!
  using ApproxEvalType= ApproxType;
  using MayInvalidateExact = std::true_type;
  using MayInvalidateApprox= std::false_type;
  constexpr LazyWrapExactRV(ExactType & init_exact) : exact_value(init_exact)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for exact values, move version" << std::endl;
#endif
  }
  LazyWrapExactRV(ExactType const &&) = delete;

  inline ExactType && eval_exact() { return std::move(exact_value); }
  inline ApproxType eval_approx() const
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_APPROX
    std::cout << "Computing approximation inside wrapper, as requested (move version)" << std::endl;
    return static_cast<ApproxType>(exact_value);
#endif
  }
  ExactType & exact_value;
};

template<class ExactType, class ApproxType>
class LazyWrapBothCR
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType const &;
  using ApproxEvalType = ApproxType const &;
  using MayInvalidateExact = std::false_type;
  using MayInvalidateApprox = std::false_type;
  constexpr LazyWrapBothCR(ExactType const &init_exact, ApproxType const &init_approx)
    :exact_value(init_exact),approx_value(init_approx)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for exact value, precomputed approximation" << std::endl;
#endif
  }
  LazyWrapBothCR(ExactType const &&,ApproxType const &) = delete;
  LazyWrapBothCR(ExactType const &, ApproxType const &&)= delete;
  LazyWrapBothCR(ExactType const &&,ApproxType const &&)= delete;

  inline ExactEvalType eval_exact() const { return exact_value; }
  inline ApproxEvalType eval_approx() const { return approx_value; }

  ExactType const & exact_value;
  ApproxType const & approx_value;
};

template<class ExactType, class ApproxType>
class LazyWrapBothRV
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType &&;
  using ApproxEvalType = ApproxType &&;
  using MayInvalidateExact = std::true_type;
  using MayInvalidateApprox= std::true_type;
  constexpr LazyWrapBothRV(ExactType &init_exact,ApproxType &init_approx)
    :exact_value(init_exact), approx_value(init_approx)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for exact value, precomputed approximation, MOVE" << std::endl;
#endif
  }
  LazyWrapBothRV(ExactType const &&, ApproxType const &&) = delete;
  LazyWrapBothRV(ExactType const &,  ApproxType const &&) = delete;
  LazyWrapBothRV(ExactType const &&, ApproxType const &)  = delete;

  inline ExactType&& eval_exact() { return std::move(exact_value); }
  inline ApproxType&& eval_approx() { return std::move(approx_value); }
  ExactType & exact_value;
  ApproxType & approx_value;
};

template<class CombinedType>
class LazyWrapCombinedCR
{
  public:
  using ExactType = typename CombinedType::ExactType;
  using ApproxType= typename CombinedType::ApproxType;
  using IsLazyNode= std::true_type;
  using IsLazyLeaf= std::true_type;
  using ExactEvalType = ExactType const &;
  using ApproxEvalType= ApproxType const &;
  using MayInvalidateExact = std::false_type;
  using MayInvalidateApprox= std::false_type;
  constexpr LazyWrapCombinedCR(CombinedType const &init_combined):combined_value(init_combined)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for combined value." << std::endl;
#endif
  }
  LazyWrapCombinedCR(CombinedType const &&init_combined) = delete;

  inline ExactType const & eval_exact() const { return combined_value.access_exact(); }
  inline ApproxType const & eval_approx() const { return combined_value.access_approx(); }

  CombinedType const & combined_value;
};

template<class CombinedType>
class LazyWrapCombinedRV
{
  public:
  using ExactType = typename CombinedType::ExactType;
  using ApproxType= typename CombinedType::ApproxType;
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType  = ExactType  &&;
  using ApproxEvalType = ApproxType &&;
  using MayInvalidateExact = std::true_type;
  using MayInvalidateApprox= std::true_type;
  constexpr LazyWrapCombinedRV(CombinedType &init_combined) : combined_value(init_combined)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creaye Lazy Wrapper for combined value, MOVE version." << std::endl;
#endif
  }
  LazyWrapCombinedRV(CombinedType const &&) = delete;
  inline ExactEvalType eval_exact() { return combined_value.access_exact(); }
  inline ApproxEvalType eval_approx() { return combined_value.access_approx(); }
  CombinedType & combined_value;
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

template<class ELP, class Approximation, ScalarOrVector s_or_v> class Lazy_Identity
{
//  static_assert(Arg::IsLazyNode::value,"Invalid arg");
  public:
  LAZY_FUNCTION;
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
//  using ArgTree = typename Arg::TreeType;
//  using TreeType= std::tuple<ArgTree const>;
  static constexpr int nargs = 1;
  static constexpr ScalarOrVector scalar_or_vector = s_or_v;
  using ExactEvalType = typename std::conditional<s_or_v==ScalarOrVector::scalar_type,ExactScalarType,ExactVectorType>::type;
  using ApproxEvalType= typename std::conditional<s_or_v==ScalarOrVector::scalar_type,ApproxScalarType,ApproxVectorType>::type;
  Lazy_Identity(...) = delete;

  static std::string fun_name() {return "Identity function";}

  constexpr inline static ExactEvalType  call_exact ( ExactEvalType const &arg) {return arg;}
  constexpr inline static ApproxEvalType call_approx( ApproxEvalType const &arg) { return arg; }

//  CPP14CONSTEXPR inline static ExactEvalType eval_exact( TreeType const & arg)
//  {
//    return Arg(std::get<0>(arg)).eval_exact();
//  }
//  CPP14CONSTEXPR inline static ApproxEvalType eval_approx(TreeType const & arg)
//  {
//    return Arg(std::get<0>(arg)).eval_approx();
//  }
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

  static std::string fun_name() {return "Scalar Product";};
  inline static ExactEvalType eval_exact(TreeType const & arg)
  {
    return compute_sc_product_exact(
      LHS(std::get<0>(arg)).eval_exact(),
      RHS(std::get<1>(arg)).eval_exact()
      );
  }
  inline static ApproxEvalType eval_approx(TreeType const & arg)
  {
    return compute_sc_product_approx(
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
  static std::string fun_name() {return "Norm2";};

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
