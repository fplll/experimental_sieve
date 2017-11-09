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
  "delayed computation object" should be performed mostly at compile time as to give no overhead
  (in theory, modulo bugs, lack of compile optimization)
  Only when such a delayed object is _assigned_ to something is the computation performed.
  More generally, computation is performed whenever an appropriate conversion is (explicitly or
  implicitly requested). A particular case is comparison via <,<=, etc. which trigger computation.

  The allowed data types for A,B may be delayed expressions themselves.
  (this allows building up complex expression trees). A and B are supposed to support both an
  exact and an approximate mode.

  Note that the implementation is not "eidetic", i.e. results are not stored.
  The use-case is that we usually only trigger evaluation (at most) once.

  Note further that (even recursive) expression trees only ever store references to the *original*
  objects involved in their creation and not to internal nodes of the expression trees.

  This makes it possible to store a delayed expression as long as all original objects are in scope.
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
  using ExactScalarType = Get_ScalarProductStorageType<ELP>; \
  using ApproxVectorType = Approximation; \
  using ApproxScalarType = typename Approximation::ScalarProductType

// to differentiate between vectors and scalars. Mostly used for static_asserts to catch bugs.
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
  This class defines an interface for objects that store an exact value and an approximation.
  This is for storage only, we have no arithmetic etc.
  This class itself is only used for testing.

  Such functionality is defered to ApproximatedPoint.h.
  Classes encapsulating such a combination of exact/approx. object should adhere to this interface.
  (as far as meaningful)
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

  /*
  constexpr      ExactType    eval_exact() const &  { return exact_object;}
  CPP14CONSTEXPR ExactType&&  eval_exact() &&       { return std::move(exact_object);}
  constexpr      ApproxType   eval_approx() const & { return approx_object;}
  CPP14CONSTEXPR ApproxType&& eval_approx() &&      { return std::move(approx_object);}
  */

  constexpr      ExactType  const & access_exact()  const { return exact_object; }
  CPP14CONSTEXPR ExactType        & access_exact()        { return exact_object; }
  constexpr      ApproxType const & access_approx() const { return approx_object; }
  CPP14CONSTEXPR ApproxType       & access_approx()       { return approx_object; }

  template<class Arg, TEMPL_RESTRICT_DECL(ApproxLevelOf<Arg>::value > ApproxLevel)>
  CPP14CONSTEXPR explicit ObjectWithApproximation(Arg &&arg)
    :ObjectWithApproximation( std::forward<Arg>(arg).eval_exact() ) {}
  template<class Arg, TEMPL_RESTRICT_DECL( ApproxLevelOf<Arg>::value < ApproxLevelOf<ExactClass>::value)>
  constexpr explicit operator Arg() const & { return static_cast<Arg>(exact_object); }
  template<class Arg, TEMPL_RESTRICT_DECL( ApproxLevelOf<Arg>::value < ApproxLevelOf<ExactClass>::value)>
  CPP14CONSTEXPR explicit operator Arg() && { return static_cast<Arg>(std::move(exact_object)); }

};

/**
  SieveLazyEval is the main class of this module:

  SieveLazyEval<LazyFunction, Args>
  encapsulates evaluation of a function F(x1,...), where x_1,... are of types Args...
  and the function F is encapsulated by a LazyFunction class.

  Available LazyFunction classes have names Lazy_F and adhere to a certain format (see below)
  Note that LazyFunction models a function that takes its arguments either by value, rvalue-ref or const lvalue-ref.
  (i.e. no non-const lvalue-ref)

  Args... is a tuple of types that are usually of the form SieveLazyEval again or are leaves of the
  expression tree.

  Note that leaves are of different (wrapper) types, but expose nearly the same syntax.
  (See below for how the leaves work)

  Every node of the tree, including SieveLazyEval<LazyFunction, Args> supports
  two "modes" of evaluation: eval_exact and eval_approx, which trigger either
  exact or approximate evaluation.

  In order for templates to recognize them. each node has the tag IsLazyNode set to std::true_type.
  ApproxLevel denotes the number of (non-bit) approximations attached to the (result) object
  represented. Typically ApproxLevel == 1. We emphasize that this is completely unrelated to the
  depth of the evaluation tree.


  We support:
  Assignment / conversion operators to the corresponding types (These trigger evaluations).
  Comparison operators (These trigger approximate and *possibly* exact evaluations)
  NOTE: Comparisons short-circuit if the approximate comparison is false.
        In particular, A < B and B > A are not equivalent:
          We might sometimes err and return false (even thogh the exact values compare true),
          so the direction of comparison determines on which side the error is.
          Furthermore, efficiency is optimized if a false result occurs more often
*/


template<class LazyFunction, class... Args> class SieveLazyEval
{
  static_assert(LazyFunction::IsLazyFunction::value,"No Lazy Function");
  static_assert(sizeof...(Args) == LazyFunction::nargs, "Wrong number of arguments");
  static_assert(mystd::conjunction< Has_IsLazyNode<Args>... >::value,"Some argument is wrong.");

  public:
  // cv-unqualified return types of eval_*
  using ExactEvalType =  typename LazyFunction::ExactEvalType;
  using ApproxEvalType = typename LazyFunction::ApproxEvalType;
  // tags used for various static_assert's and template overload selection:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::false_type; // not a leaf of the expression tree.
  static constexpr unsigned int ApproxLevel = LazyFunction::ApproxLevel;
  static_assert(ApproxLevel >0, "Approximation level is 0.");

  // for now:
  static_assert(mystd::conjunction<mystd::bool_constant<ApproxLevel == ApproxLevelOf<Args>::value>...>::value,"All arguments must have the same approximation level");

  // EvalOnce means that calling eval_* might actually invalidate the data stored to / refered to.
  // If this is set, we
  // - store a non-const tuple of arguments
  // - optionally (implicitly) disable copying of SieveLazyEval
  // (via recursively disabling the automatic generation of the copy constructors / copy assignment)
  // (Note that disabling automagically generated constructors based on traits is quite non-trivial
  // in C++ and subject to many gotchas... The issue being that if you template these functions, the
  // compiler may no longer recognize them as "special".)
  // - eval_* may only be called on rvalue.
  // - may only call eval_* once. More precisely, calling eval_exact() forbids any future calls to
  //   either eval_* and calling eval_approx forbids calls to eval_approx.
  // These restrictions might not be enforced.
  // (Note that these restriction could be slightly relaxed, separating the restrictions, but it
  //  seems not worth the hassle.)
  // This happens if the leaves of the trees encode pass by rvalue-semantics.

  // EvalOnce is either std::true_type or std::false_type, EvalOnce_v is either true or false.
  using EvalOnce = MyDisjunction< typename Args::EvalOnce...>; //OR of Args
  static constexpr bool EvalOnce_v = EvalOnce::value;

  using TreeType = std::tuple<MaybeConst<!(Args::EvalOnce_v),Args>...>; // const-ness depends on EvalOnce.
  MaybeConst<!EvalOnce_v,TreeType> args;
  // This restriction can easily be removed. The issue is just that we want no default-constructor.
  // (at least for nargs > 0)
  SieveLazyEval() = delete; static_assert( sizeof...(Args)>0,"0-ary functions not supported yet." );

  template<class... FunArgs>
  CONSTEXPR_IN_NON_DEBUG_TC explicit SieveLazyEval(FunArgs&&... fun_args)
    : args(std::forward<FunArgs>(fun_args)...)
  {
    // Args and FunArgs must be the same (up to const and reference-ness)
    static_assert(sizeof...(Args) == sizeof...(FunArgs),"wrong number of arguments to constructor.");
    static_assert(mystd::conjunction< std::is_same<Args,mystd::decay_t<FunArgs>>...>::value,"Wrong type of arguments to constructor.");

    // If Args_i::EvalOnce is set, then fun_args_i must be passed via rvalue-reference (enforcing move!).
    // If fun_args_i is passed via rvalue-ref, then (via the way forwarding reference capture works)
    // FunArgs_i itself is a non-refence type.
    // This means that we must have (FunArgs_i is non reference || Args_i::EvalOnce_v == false) for each i
    static_assert(mystd::conjunction<
      mystd::bool_constant<
        (std::is_reference<FunArgs>::value == false) || (Args::EvalOnce_v == false)
      >...
    >::value,"Use (possibly explicit) move semantics for Lazy objects that may encapsulate RValues." );
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy function wrapper object for function " << LazyFunction::fun_name() << std::endl;
#endif
  }

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

// MyMakeIndexSeq generates a(n empty) struct MyIndexSeq<0,1,2,...>. This is used as a dummy paramter
// to select the correct version of do_eval_approx and allows to actually un-std::tuple the argument.
// (This is the least roundabout way of doing it and corresponds to the implementation of std::apply)
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
};


// General conditions on LHS, RHS to perform the comparisons via eval_*:
// Clearly, these definitions only make sense if we have approximations, so we require at least one
// argument to have an positive ApproxLevel.
// Also, we may NOT use the functions below to compare lattice points:
// (we may use it on their norm2's)
// Note that lattice_point_1 < lattice_point_2 is defined to compute
// lattice_point_1.get_norm2_exact() < lattice_point_2.get_norm2_exact() in LatticePointConcept.h.
//  In particular, comparing lattice points directly bypasses all approximations
// (get_norm2_exact() goes to approximation-level 0). This is to ensure that the comparison function
// on lattice points actually is an ordering, which approximation errors might violate.
//
// Comparing a lattice point with a non-lattice point should also not call the functions below.
// (Indeed, it is completely unclear what this should be)

// Be aware that these functions are inside namespace GaussSieve::LazyEval,
// so they are only considered if at least one argument is in that namespace as well.
// In particular, the "Disable this template for lattice points" is redundant.

#define SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS                                  \
      (MyNOR<IsALatticePoint<LHS>,IsALatticePoint<RHS>>::value                   \
  &&  ((ApproxLevelOf<LHS>::value > 0) || (ApproxLevelOf<RHS>::value >0)))

// general comparison for < : We bring down the arguments to the same approximation level,
// then we compare approximately. If that results in true, we actually do an exact check as well.
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value > ApproxLevelOf<RHS>::value))>
inline bool operator< (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs).eval_exact() < std::forward<RHS>(rhs);
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value < ApproxLevelOf<RHS>::value))>
inline bool operator< (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs) < std::forward<RHS>(rhs).eval_exact();
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value == ApproxLevelOf<RHS>::value))>
inline bool operator< (LHS && lhs, RHS && rhs)
{
// Note that C++ short-circuits && (unless someone overloads && for the return type of <)
  return (std::forward<LHS>(lhs).eval_approx() < std::forward<RHS>(rhs).eval_approx() )
      && (std::forward<LHS>(lhs).eval_exact()  < std::forward<RHS>(rhs).eval_exact()  );
}

// same for <=
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value > ApproxLevelOf<RHS>::value))>
inline bool operator<= (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs).eval_exact() <= std::forward<RHS>(rhs);
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value < ApproxLevelOf<RHS>::value))>
inline bool operator<= (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs) <= std::forward<RHS>(rhs).eval_exact();
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value == ApproxLevelOf<RHS>::value))>
inline bool operator<= (LHS && lhs, RHS && rhs)
{
  return (std::forward<LHS>(lhs).eval_approx() <= std::forward<RHS>(rhs).eval_approx() )
      && (std::forward<LHS>(lhs).eval_exact()  <= std::forward<RHS>(rhs).eval_exact()  );
}

// same for >
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value > ApproxLevelOf<RHS>::value))>
inline bool operator> (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs).eval_exact() > std::forward<RHS>(rhs);
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value < ApproxLevelOf<RHS>::value))>
inline bool operator> (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs) > std::forward<RHS>(rhs).eval_exact();
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value == ApproxLevelOf<RHS>::value))>
inline bool operator> (LHS && lhs, RHS && rhs)
{
  return (std::forward<LHS>(lhs).eval_approx() > std::forward<RHS>(rhs).eval_approx() )
      && (std::forward<LHS>(lhs).eval_exact()  > std::forward<RHS>(rhs).eval_exact()  );
}

// and for >=
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value > ApproxLevelOf<RHS>::value))>
inline bool operator>= (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs).eval_exact() >= std::forward<RHS>(rhs);
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value < ApproxLevelOf<RHS>::value))>
inline bool operator>= (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs) >= std::forward<RHS>(rhs).eval_exact();
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value == ApproxLevelOf<RHS>::value))>
inline bool operator>= (LHS && lhs, RHS && rhs)
{
  return (std::forward<LHS>(lhs).eval_approx() >= std::forward<RHS>(rhs).eval_approx() )
      && (std::forward<LHS>(lhs).eval_exact()  >= std::forward<RHS>(rhs).eval_exact()  );
}

// we also do this for ==

template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value > ApproxLevelOf<RHS>::value))>
inline bool operator== (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs).eval_exact() == std::forward<RHS>(rhs);
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value < ApproxLevelOf<RHS>::value))>
inline bool operator== (LHS && lhs, RHS && rhs)
{
  return std::forward<LHS>(lhs) == std::forward<RHS>(rhs).eval_exact();
}
template<class LHS, class RHS, TEMPL_RESTRICT_DECL ( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS
  && (ApproxLevelOf<LHS>::value == ApproxLevelOf<RHS>::value))>
inline bool operator== (LHS && lhs, RHS && rhs)
{
//  return (std::forward<LHS>(lhs).eval_approx() == std::forward<RHS>(rhs).eval_approx() )
//      && (std::forward<LHS>(lhs).eval_exact()  == std::forward<RHS>(rhs).eval_exact()  );
  return std::forward<LHS>(lhs).eval_exact() == std::forward<RHS>(rhs).eval_exact();
}

// != is defined in terms of ==

template<class LHS, class RHS, TEMPL_RESTRICT_DECL( SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS)>
inline bool operator != (LHS && lhs, RHS && rhs)
{
  return !(std::forward<LHS>(lhs)==std::forward<RHS>(rhs));
}

#undef SIEVE_GAUSS_LAZY_COMPARISON_CONDITIONS


/**
  LazyWrap* model the leaves of the expression trees.
  These classes are just wrappers around a reference to a *
  and expose the same interface as SieveLazyEval.

  Note that the functionality may be a bit more limited than that of SieveLazyEval (due to lazyness
  of the coder). As a workaround, you may use SieveLazyEval with an IdentityFunction.

  There are several different wrapper that differ in the following aspects:
  A wrapper may either wrap
  - only an Exact object. In this case, an approximation is computed on demand. This should probably
    never be used.
  - an exact and an approximate object.
  - a single combined object (such as ObjectWithApproximation defined above) that encapsulates both.
    This object is required to have access_exact() and access_approx() members that return references

  Furthermore, there is an additional differentiation regarding calling conventions.
  We support (up to) four different calling conventions:
  - Const-Ref (CR): The wrapper is initialized with an lvalue and stores a const-reference.
                    upon eval_*, the lvalue used to initialize must still be valid.
                    The function used in eval_ must pass this by value of const-ref

  - Delayed-move(RV): This wrapper models delayed (explicit) rvalue-passing.
                      The wrapper is initialized with an lvalue.
                      Upon eval_*, the lvalue used to initialize must still be valid.
                      The function used in eval_ will be called on an rvalue.
                      (This means that the wrapper (lazily) takes ownership of the wrapped object,
                      and usage of eval_* might modify the object.
                      If the evaluation tree contains such an object, EvalOnce is set.)

  Not implemented yet (not needed so far)
  - Move :  This models immediate rvalue-passing. The wrapper is initialized with an rvalue and
            directly move-stores the wrapped object. The function used in eval_* will be called on
            an rvalue.
            This also sets EvalOnce

  Not implemented yes (not needed so far)
  - Capture:  This mode stores the values directly (not as reference). The wrapper may be initialized
              with either rvalues or lvalues. The argument will be passed to eval_* by value or const-ref.
*/



/**
  Exact object, Const-Ref
*/

template<class ExactType, class ApproxType>
class LazyWrapExactCR
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType;
  using ApproxEvalType = ApproxType;
  using EvalOnce = std::false_type;
  constexpr static bool EvalOnce_v = false;
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

/**
  Exact object, delayed-move
*/

template<class ExactType, class ApproxType>
class LazyWrapExactRV
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType; // Note that we actually return a rvalue-ref!
  using ApproxEvalType= ApproxType;
  using EvalOnce = std::true_type;
  constexpr static bool EvalOnce_v = true;
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

/**
  Both Exact and Approximate Object in separate classes, Const-Ref passing.
*/

template<class ExactType, class ApproxType>
class LazyWrapBothCR
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType;
  using ApproxEvalType = ApproxType;
  using EvalOnce = std::false_type;
  constexpr static bool EvalOnce_v = false;
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

/**
  Both Exact and Approximate Object in separate classes, Delayed-Move
*/


template<class ExactType, class ApproxType>
class LazyWrapBothRV
{
  public:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  using ExactEvalType = ExactType;
  using ApproxEvalType = ApproxType;
  using EvalOnce = std::true_type;
  constexpr static bool EvalOnce_v = true;
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

/**
  Wrap combined class, Const-Ref
*/

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
  using EvalOnce = std::false_type;
  constexpr static bool EvalOnce_v = false;
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

/**
  Wrap combined class, Delayed-Move
*/

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
  using EvalOnce = std::true_type;
  constexpr static bool EvalOnce_v = true;
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

  The class is static (i.e. cannot / should not / need not be instantiated)
  It provides the following public typedefs / static constants:
  static constexpr int nargs : number of arguments.
  typename ExactEvalType : Result Type of Exact Evaluation
  typename ApproxEvalType: Result Type of Approximate Evaluation
  (these should be non-reference, non cv-qualified)
  a function fun_name() returning a string identifying the function. For debug purposes.
  a typedef IsLazyFunction = std::true_type. For traits detection.
  static functions call_exact(args...) and call_approx(args...). These functions should take
  nargs may arguments and return a (possibly reference to cv-qualified) ExactEvalType / ApproxEvalType
  Note that call_exact / call_approx might be templates or overloaded on various argument types.

*/

/**
  Identity function. F(x) = x. Can be used to promote a
  LazyWrap* to a SieveLazyEval<ELP,Approxiomation,Lazy_Identity<...>,...>
  */

#define GAUSS_SIEVE_LAZY_FUN(nargs_,namestring)                             \
  constexpr static unsigned int nargs = (nargs_);                           \
  using IsLazyFunction = std::true_type;                                    \
  static std::string fun_name() {return namestring ;}

template<class ExactType, class ApproxType>
class Lazy_Identity
{
  public:
  GAUSS_SIEVE_LAZY_FUN(1,"Identity function")
  using ExactEvalType = mystd::decay_t<ExactType>;
  using ApproxEvalType= mystd::decay_t<ApproxType>;
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ExactType>::value + 1;
  static_assert(!(std::is_same<ExactEvalType,ApproxEvalType>::value),"");
  template<class Arg, TEMPL_RESTRICT_DECL2(std::is_same<ExactEvalType,mystd::decay_t<Arg>>)>
  inline static Arg && call_exact(Arg &&exact) { return std::forward<Arg>(exact);}
  template<class Arg, TEMPL_RESTRICT_DECL2(std::is_same<ApproxEvalType,typename mystd::decay_t<Arg>>)>
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
  GAUSS_SIEVE_LAZY_FUN(2,"Scalar Product")
  using ExactEvalType  = ExactScalarType;
  using ApproxEvalType = ApproxScalarType;
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ELP>::value + 1;
  template<class LHS, class RHS>
  inline static auto call_exact(LHS &&lhs, RHS &&rhs)
  -> decltype( compute_sc_product_exact( std::declval<LHS &&>(), std::declval<RHS &&>()   )  )
  {
    static_assert(std::is_same<mystd::decay_t<LHS>,ExactVectorType>::value,"LHS wrong type.");
    static_assert(std::is_same<mystd::decay_t<RHS>,ExactVectorType>::value,"RHS wrong type.");
    return compute_sc_product_exact(std::forward<LHS>(lhs),std::forward<RHS>(rhs));
  }
  template<class LHS, class RHS>
  inline static ApproxScalarType call_approx(LHS &&lhs, RHS &&rhs)
  {
    static_assert(std::is_same<mystd::decay_t<LHS>, ApproxVectorType>::value,"LHS wrong type.");
    static_assert(std::is_same<mystd::decay_t<RHS>, ApproxVectorType>::value,"RHS wrong type.");
    return compute_sc_product_approx(std::forward<LHS>(lhs),std::forward<RHS>(rhs));
  }
};

template<class ELP, class Approximation> class Lazy_Norm2
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  GAUSS_SIEVE_LAZY_FUN(1,"Get Norm2")
  constexpr static unsigned int ApproxLevel = ApproxLevelOf<ELP>::value + 1;
  using ExactEvalType  = ExactScalarType;
  using ApproxEvalType = ApproxScalarType;

  inline static auto call_exact(ExactVectorType const &arg)
  -> decltype( std::declval<ExactVectorType const>().get_norm2_exact()  )
  {
    return arg.get_norm2_exact();
  }
  inline static ApproxScalarType call_approx(ApproxVectorType const &arg)
  {
    return arg.get_approx_norm2();
  }
};

}} //end namespaces

#undef CONSTEXPR_IN_NON_DEBUG_TC

#endif
