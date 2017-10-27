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

// We might remove this:
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
// DEPRECATED:
//[[deprecated]]
enum struct [[deprecated]] ScalarOrVector{ scalar_type, vector_type };

CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(IsLazyNode, std::true_type, Has_IsLazyNode);

#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
  #define CONSTEXPR_IN_NON_DEBUG_TC
#else
  #define CONSTEXPR_IN_NON_DEBUG_TC constexpr
#endif


// forward declarations:
//template<class LazyFunction, class... Args> class SieveLazyEval;

/**
  This class stores a object with an approximation to that object.
  Note that this is for storage only. We have no arithmetic etc.
  Such functionality is defered to ApproximatedPoint.h
*/

template<class ExactClass, class ApproximationClass>
struct ObjectWithApproximation
{
  using ExactType = ExactClass;
  using ApproxType= ApproximationClass;
  ExactType  exact_object;
  ApproxType approx_object;
  static constexpr unsigned int ApproxLevel = ApproxLevelOf<ExactClass>::value + 1;

  // This would cause ambiguous overloads.
  static_assert(!(std::is_same<ExactType,ApproxType>::value),"Can not approximate by itself currently");

  constexpr      explicit ObjectWithApproximation(ExactType const& exact,ApproxType const &approx)
    :exact_object(exact),approx_object(approx){}
  CPP14CONSTEXPR explicit ObjectWithApproximation(ExactType     && exact,ApproxType const &approx)
    :exact_object(std::move(exact)),approx_object(approx){}
  CPP14CONSTEXPR explicit ObjectWithApproximation(ExactType const& exact,ApproxType      &&approx)
    :exact_object(exact),approx_object(std::move(approx)){}
  CPP14CONSTEXPR explicit ObjectWithApproximation(ExactType     && exact,ApproxType      &&approx)
    :exact_object(std::move(exact)),approx_object(std::move(approx)){}
  constexpr      explicit ObjectWithApproximation(ExactType const &exact)
    :exact_object(exact), approx_object(exact) {}
  CPP14CONSTEXPR explicit ObjectWithApproximation(ExactType && exact)
    :exact_object(std::move(exact)), approx_object(exact_object) {}
  constexpr      explicit operator ExactType()  const & { return exact_object;}
  CPP14CONSTEXPR explicit operator ExactType()  &&      { return std::move(exact_object);}
  constexpr      explicit operator ApproxType() const & { return approx_object;}
  CPP14CONSTEXPR explicit operator ApproxType() &&      { return std::move(approx_object);}

  constexpr      ExactType  const & access_exact()  const { return exact_object; }
  CPP14CONSTEXPR ExactType        & access_exact()        { return exact_object; }
  constexpr      ApproxType const & access_approx() const { return approx_object; }
  CPP14CONSTEXPR ApproxType       & access_approx()       { return approx_object; }

  //TODO
  //
////  template<class LazyObject, TEMPL_RESTRICT_DECL2(Has_IsLazyNode<typename std::decay<LazyObject>::type> )>
////  constexpr explicit ScalarWithApproximation(LazyObject &&lazy_object)
////    :exact_scalar(lazy_object.eval_exact()),approx_scalar(lazy_object.eval_approx()) //todo: move semantics
////  {
////    using Arg = typename std::decay<LazyObject>::type;
////    static_assert(Arg::scalar_or_vector == ScalarOrVector::scalar_type,"Trying to assign a vector to a scalar");
////  }
//};
};

// TODO: Update documentation to reflect recent changes.

/**
  SieveLazyEval is the main class of this module:

  SieveLazyEval<LazyFunction, Args>
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

template<class LazyFunction, unsigned int approx_level, class... Args> class SieveLazyEval
{
  static_assert(LazyFunction::IsLazyFunction::value,"No Lazy Function");
  static_assert(sizeof...(Args) == LazyFunction::nargs, "Wrong number of arguments");
  static_assert(MyConjunction< Has_IsLazyNode<Args>... >::value,"Some argument is wrong.");
  static_assert(approx_level >0, "Approximation level is 0");

  public:
//  using ELP = typename LazyFunction::ELP_t;
//  using Approximation = typename LazyFunction::Approximation_t;
  // assert all ELP's / Approx's of Arg... are the same?
  // Actually, this may be too restrictive.
//  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
//  using TreeType = std::tuple< typename Args::TreeType const...  >;
//  using TreeType      = typename LazyFunction::TreeType;
  using ExactEvalType =  typename LazyFunction::ExactEvalType;
  using ApproxEvalType = typename LazyFunction::ApproxEvalType;
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::false_type;
//  static constexpr ScalarOrVector scalar_or_vector = LazyFunction::scalar_or_vector;
  //TODO: ensure consistency.
  using MayInvalidateExact = MyDisjunction< typename Args::MayInvalidateExact...>; //OR of Args
  using MayInvalidateApprox= MyDisjunction< typename Args::MayInvalidateApprox...>;//OR of Args
  static constexpr unsigned int ApproxLevel = approx_level;

  using TreeType = std::tuple<Args...>; //TODO: const-correctness
  TreeType args; //const?
  SieveLazyEval() = delete; static_assert( sizeof...(Args)>0,"" ); //We might remove this.
  constexpr explicit SieveLazyEval(SieveLazyEval const &other) = default;
  constexpr explicit SieveLazyEval(SieveLazyEval && other) = default;
  SieveLazyEval& operator=(SieveLazyEval const &other) = default;
  SieveLazyEval& operator=(SieveLazyEval &&other) = default;

  //  TreeType const args;

  CONSTEXPR_IN_NON_DEBUG_TC explicit SieveLazyEval(TreeType const & fn_args) : args(fn_args)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy function wrapper object (via tuple) for function " << LazyFunction::fun_name() << std::endl;
#endif
  }
  CONSTEXPR_IN_NON_DEBUG_TC explicit SieveLazyEval(TreeType && fn_args) : args(std::move(fn_args))
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy function wrapper object (via tuple-move) for function " << LazyFunction::fun_name() << std::endl;
#endif
  }
  // TODO: Add direct init via Args...

  template<std::size_t... iarg>
  inline ExactEvalType do_eval_exact(MyIndexSeq<iarg...>)
  {
    static_assert(sizeof...(iarg) == sizeof...(Args),"Something is very wrong");
#ifdef DEBUG_SIEVE_LAZY_TRACE_EVALS
    std::cout << "Calling function exactly. Function is " << LazyFunction::fun_name() << std::endl;
#endif
    return LazyFunction::call_exact( std::get<iarg>(args).eval_exact()... );
  }

  template<std::size_t... iarg>
  inline ExactEvalType do_eval_approx(MyIndexSeq<iarg...>)
  {
    static_assert(sizeof...(iarg) == sizeof...(Args),"Something is very wrong");
#ifdef DEBUG_SIEVE_LAZY_TRACE_EVALS
    std::cout << "Calling function approximately. Function is " << LazyFunction::fun_name() << std::endl;
#endif
    return LazyFunction::call_approx( std::get<iarg>(args).eval_approx()... );
  }


  inline ExactEvalType eval_exact()
  {
    return do_eval_exact(MyMakeIndexSeq<sizeof...(Args)>{} );
  }

  inline ApproxEvalType eval_approx()
  {
    return do_eval_approx(MyMakeIndexSeq<sizeof...(Args)>{} );
  }

  inline explicit operator ExactEvalType() { return eval_exact(); }
  inline explicit operator ApproxEvalType() { return eval_approx(); }

  inline bool operator< ( ExactEvalType const & rhs)
  {
    return eval_exact() < rhs;
  }

  inline bool operator> ( ExactEvalType const & rhs)
  {
    return eval_exact() > rhs;
  }

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
  using ExactEvalType = ExactType;
  using ApproxEvalType = ApproxType;
  using MayInvalidateExact = std::false_type;
  using MayInvalidateApprox= std::false_type;
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ExactType>::value + 1;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapExactCR(ExactType const &init_exact) : exact_value(init_exact)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for exact values." << std::endl;
#endif
  };
  LazyWrapExactCR(ExactType const &&) = delete;

  inline constexpr ExactType const & eval_exact()  const { return exact_value; }
  inline CONSTEXPR_IN_NON_DEBUG_TC ApproxType eval_approx() const
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
  using ExactEvalType = ExactType; // Note that we actually return a rvalue-ref!
  using ApproxEvalType= ApproxType;
  using MayInvalidateExact = std::true_type;
  using MayInvalidateApprox= std::false_type;
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ExactType>::value + 1;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapExactRV(ExactType & init_exact) : exact_value(init_exact)
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
  using ExactEvalType = ExactType;
  using ApproxEvalType = ApproxType;
  using MayInvalidateExact = std::false_type;
  using MayInvalidateApprox = std::false_type;
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ExactType>::value + 1;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapBothCR(ExactType const &init_exact, ApproxType const &init_approx)
    :exact_value(init_exact),approx_value(init_approx)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for exact value, precomputed approximation" << std::endl;
#endif
  }
  LazyWrapBothCR(ExactType const &&,ApproxType const &) = delete;
  LazyWrapBothCR(ExactType const &, ApproxType const &&)= delete;
  LazyWrapBothCR(ExactType const &&,ApproxType const &&)= delete;

  inline ExactType  const & eval_exact() const { return exact_value; }
  inline ApproxType const & eval_approx() const { return approx_value; }

  ExactType const & exact_value;
  ApproxType const & approx_value;
};

template<class ExactType, class ApproxType>
class LazyWrapBothRV
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType;
  using ApproxEvalType = ApproxType;
  using MayInvalidateExact = std::true_type;
  using MayInvalidateApprox= std::true_type;
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ExactType>::value + 1;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapBothRV(ExactType &init_exact,ApproxType &init_approx)
    :exact_value(init_exact), approx_value(init_approx)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for exact value, precomputed approximation, MOVE" << std::endl;
#endif
  }
  LazyWrapBothRV(ExactType const &&, ApproxType const &&) = delete;
  LazyWrapBothRV(ExactType const &,  ApproxType const &&) = delete;
  LazyWrapBothRV(ExactType const &&, ApproxType const &)  = delete;

  inline ExactType&&  eval_exact()  { return std::move(exact_value); }
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
  using ExactEvalType = ExactType;
  using ApproxEvalType= ApproxType;
  using MayInvalidateExact = std::false_type;
  using MayInvalidateApprox= std::false_type;
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ExactType>::value + 1;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapCombinedCR(CombinedType const &init_combined):combined_value(init_combined)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy Wrapper for combined value." << std::endl;
#endif
  }
  LazyWrapCombinedCR(CombinedType const &&init_combined) = delete;

  inline ExactType  const & eval_exact()  const { return combined_value.access_exact(); }
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
  using ExactEvalType  = ExactType;
  using ApproxEvalType = ApproxType;
  using MayInvalidateExact = std::true_type;
  using MayInvalidateApprox= std::true_type;
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ExactType>::value + 1;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapCombinedRV(CombinedType &init_combined) : combined_value(init_combined)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creaye Lazy Wrapper for combined value, MOVE version." << std::endl;
#endif
  }
  LazyWrapCombinedRV(CombinedType const &&) = delete;
  inline ExactType  && eval_exact()  { return std::move(combined_value.access_exact()); }
  inline ApproxType && eval_approx() { return std::move(combined_value.access_approx()); }
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

template<class ExactType, class ApproxType>
class Lazy_Identity
{
  public:
  static int constexpr nargs = 1;
  static std::string fun_name() {return "Identity function";}
  using IsLazyFunction = std::true_type;
  using ExactEvalType = typename std::decay<ExactType>::type;
  using ApproxEvalType= typename std::decay<ApproxType>::type;
  static_assert(!(std::is_same<ExactEvalType,ApproxEvalType>::value),"");
  template<class Arg, TEMPL_RESTRICT_DECL2(std::is_same<ExactEvalType,typename std::decay<Arg>::type>)>
  inline static Arg && call_exact(Arg &&exact) { return std::forward<Arg>(exact);}
  template<class Arg, TEMPL_RESTRICT_DECL2(std::is_same<ApproxEvalType,typename std::decay<Arg>::type>)>
  inline static Arg && call_approx(Arg &&approx) { return std::forward<Arg>(approx);}
};


/**
  Scalar Product function.
  Delegates to compute_sc_product_exact resp. compute_sc_product_approx
*/

template<class ELP, class Approximation> class Lazy_ScalarProduct
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
//  using TreeType = std::tuple<ArgTreeLeft const, ArgTreeRight const>;
  static constexpr int nargs = 2;
  using IsLazyFunction = std::true_type;
  using ExactEvalType  = ExactScalarType;
  using ApproxEvalType = ApproxScalarType;
  static std::string fun_name() {return "Scalar Product";}
  template<class LHS, class RHS>
  inline static ExactScalarType call_exact(LHS &&lhs, RHS &&rhs)
  {
    static_assert(std::is_same<typename std::decay<LHS>::type,ExactVectorType>::value,"LHS wrong type.");
    static_assert(std::is_same<typename std::decay<RHS>::type,ExactVectorType>::value,"RHS wrong type.");
    return compute_sc_product(std::forward<LHS>(lhs),std::forward<RHS>(rhs));
  }
  template<class LHS, class RHS>
  inline static ApproxScalarType call_approx(LHS &&lhs, RHS &&rhs)
  {
    static_assert(std::is_same<typename std::decay<LHS>::type, ApproxVectorType>::value,"LHS wrong type.");
    static_assert(std::is_same<typename std::decay<RHS>::type, ApproxVectorType>::value,"RHS wrong type.");
    return compute_sc_product_approx(std::forward<LHS>(lhs),std::forward<RHS>(rhs));
  }
};

template<class ELP, class Approximation> class Lazy_Norm2
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  static constexpr int nargs = 1;
  using IsLazyFunction = std::true_type;
  using ExactEvalType  = ExactScalarType;
  using ApproxEvalType = ApproxScalarType;
  static std::string fun_name() {return "Norm2";};

  inline static ExactScalarType call_exact(ExactVectorType const &arg)
  {
    return arg.get_norm2();
  }
  inline static ApproxScalarType call_approx(ApproxVectorType const &arg)
  {
    return arg.get_approx_norm2();
  }
};

}} //end namespaces

#undef CONSTEXPR_IN_NON_DEBUG_TC

#endif
