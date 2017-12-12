#ifndef GAUSS_SIEVE_LAZY_H
#define GAUSS_SIEVE_LAZY_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include "TraitChecks.h"

#include "Lazy_macros.h"

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

/**
  This file also specifies the notion of leveled objects. These are objects that have several
  "levels", intended to be used for approximations. An object with two levels encapsulates a
  level-0 object and a level-1 object. For us, the level-0 object is an exact value (of e.g. a
  lattice point) and the level-1 object an approximation. Higher levels would correspond to coarser
  approximations. Such levels are closely tied to lazy evaluation, because the whole point of lazy
  evaluation is that computations on higher levels are faster than on lower levels and may often
  make the lower-level computations unneccessary.
  In particular, comparison operations on lazy objects first compute everthing on higher levels
  and if that result is false, the whole result is considered false (without even touching the
  lower levels).
*/



// forward declarations:
namespace GaussSieve::LazyEval{     // sub-namespace to inject free functions like abs
template<class LazyFunction, class... Args> class SieveLazyEval;
template<class ExactClass, class ApproximationClass> class ObjectWithApproximation;
namespace Helpers
{
template<unsigned int level, class ExactClass, class Approximation> class ObjectWithApproximationHelper;
}
}

// If DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS is set, we perform some debug output.
// As a consequence, some functions have side-effects (notably, output) and are no longer constexpr.
// This #define is local (i.e. #undef'd at the end of the file).

#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
  #define CONSTEXPR_IN_NON_DEBUG_TC
#else
  #define CONSTEXPR_IN_NON_DEBUG_TC constexpr
#endif

/**
  Leveled Objects:

  A leveled object has to declare the following members / typedefs:

  - using LeveledObject = std::true_type  for trait detection. Use Has_LeveledObject<Obj> to query.

  - static constexpr unsigned int ApproxLevel to denote the maximal level.

  - a template ObjectAtLevel<level> to denote the (std::decay'ed) object encapsulated at the desired
    level. This specifies (up to const- and reference- ness) the return types of the functions below

  - using LeveledObject_Base = std::true_type / false_type.
    Query with Has_LeveledObject_Base<Obj>;
    Denotes whether this is a "base" object. Base objects have access<level> - members and
    are *not* Lazy Objects.
    (Currently, every leveled object is either base or lazy -- this is just to make sure )

  - using LeveledComparison = std::true_type / false_type.
    This controls whether comparison operators work by comparing level-by-level, starting with the
    highest. As soon as any level compares false, the result is false.

  template member functions
  - get_value_at_level<level> : These may return values or (typically const) references to the objects.

  - access<level> : Returns reference / const-reference (depending of const-ness of the this pointer)
                  This is only defined for base objects.
*/


namespace GaussSieve
{
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(IsLazyNode, std::true_type, Has_IsLazyNode);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(LeveledObject, std::true_type, Has_LeveledObject);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(LeveledObject_Base, std::true_type, Has_LeveledObject_Base);
// CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(DelayedDefaultFunctions, std::true_type, Has_DelayedDefaultFunctions);
}

namespace GaussSieve{
//GAUSS_SIEVE_LAZY_UNARY_FUNCTION_LEVEL_DETECTION_FORWARD_DECLARE(get_norm2,get_norm2)
//GAUSS_SIEVE_LAZY_UNARY_FUNCTION_CREATE_LAZY_WRAPPER_FORWARD_DECLARE(get_norm2,get_norm2)
//GAUSS_SIEVE_LAZY_UNARY_FUNCTION_DIRECT_LAZY_FORWARD_DECLARE(call_get_norm2,get_norm2, std::true_type)
}

namespace GaussSieve{
//GAUSS_SIEVE_LAZY_UNARY_FUNCTION_FOR_DELAYED_OBJECTS_FORWARD_DECLARE(abs,abs)
//GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_FOR_DELAYED_OBJECTS_FORWARD_DECLARE(get_norm2)
}
// ---------------

/**
  Example Leveled Object:
  This class defines an interface for objects that store an exact value and an approximation.
  This class itself is only used for testing.
*/

namespace GaussSieve::LazyEval{

template<class ExactClass, class ApproximationClass>
struct ObjectWithApproximation
{
  private:
  // helper classes defined below
  template<unsigned int level> using  Helper = Helpers::ObjectWithApproximationHelper<level,ExactClass, ApproximationClass>;
  template<unsigned int level, class ExactClass2, class ApproximationClass2> friend class Helpers::ObjectWithApproximationHelper;

  // TODO: Rename things and remove these typedefs.
  using ExactType = ExactClass;
  using ApproxType= ApproximationClass;
  public:
  using LeveledComparison = std::true_type;
  using LeveledObject     = std::true_type;
  using LeveledObject_Base= std::true_type;
//  using DelayedDefaultFunctions = std::true_type;
  template<unsigned int level> using ObjectAtLevel = typename Helper<level>::Object;
  ExactType  exact_object;
  ApproxType approx_object;
  static_assert(std::is_same<ExactType, ObjectAtLevel<0>>::value,"");
  static_assert(std::is_same<ApproxType,ObjectAtLevel<1>>::value,"");

  static constexpr unsigned int ApproxLevel = 1;

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

  template<unsigned int level>
  constexpr       ObjectAtLevel<level> const & access() const { return Helper<level>::get(*this); }
  template<unsigned int level>
  CPP14CONSTEXPR  ObjectAtLevel<level>       & access()       { return Helper<level>::get(*this); }
  template<unsigned int level>
  constexpr       ObjectAtLevel<level> const & get_value_at_level() const { return Helper<level>::get(*this); }
  template<unsigned int level>
  CPP14CONSTEXPR  ObjectAtLevel<level>       & get_value_at_level()       { return Helper<level>::get(*this); }

  // TODO: Remove
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
  Helper class for the class above:

  Since C++ does not allow to partially specialize template functions, we cannot directly write a
  get<level>() function to get the object at a certain level.

  Instead, we define a helper class with a get member function, so we can partially specialize the
  class:

  ObjectWithApproximationHelper::Class<level,Exact,Approx>::get(combined_object)
  will get the exact object of combined_object for level==0 and the approximate object for level==1.
*/


namespace Helpers{
template<unsigned int level, class ExactClass, class ApproximationClass> class ObjectWithApproximationHelper
{
  // this should never be instantiated!
  static_assert(level<=1,"wrong usage"); // we only have levels 0 and 1 for ObjectWithApproximation
  friend ObjectWithApproximation<ExactClass,ApproximationClass>;
};
template<class ExactClass, class ApproximationClass>
class ObjectWithApproximationHelper<0,ExactClass,ApproximationClass>
{
  public:
  using CombinedObject = ObjectWithApproximation<ExactClass,ApproximationClass>;
  using Object = ExactClass;
  friend CombinedObject;
  private:
  static inline Object const & get(CombinedObject const &obj) { return obj.exact_object;}
  static inline Object &       get(CombinedObject & obj) { return obj.exact_object; }
};
template<class ExactClass, class ApproximationClass>
class ObjectWithApproximationHelper<1,ExactClass,ApproximationClass>
{
  public:
  using CombinedObject = ObjectWithApproximation<ExactClass,ApproximationClass>;
  using Object = ApproximationClass;
  friend CombinedObject;
  private:
  static inline Object const & get(CombinedObject const &obj) { return obj.approx_object; }
  static inline Object &       get(CombinedObject & obj) { return obj.approx_object; }
};
} // end of Helpers namespace

template<bool Good> struct LevelIsGood{static_assert(Good,"Calling at invalid level");};


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

  Every node of the tree, including SieveLazyEval<LazyFunction, Args> is a leveled object.
  Such objects support an eval<level> function that triggers evaluation at the desired level.
  level 0 corresponds to the "exact" object and higher levels model approximations.

  The maximal level is given by ApproxLevel. Typically ApproxLevel == 1.
  We emphasize that these levels are completely unrelated to the depth of the evaluation tree.

  In order for templates to recognize them. each node has the tag IsLazyNode set to std::true_type
  and IsLazyLeaf declares whether it is a leaf.

  We support (TODO!)
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
  // LazyFunction must have a specified tag defined.
  static_assert(LazyFunction::IsLazyFunction::value,"No Lazy Function");
  static_assert(sizeof...(Args) == LazyFunction::nargs, "Wrong number of arguments");
  // Each child in the expression tree has to be a lazy node.
  static_assert(mystd::conjunction< Has_IsLazyNode<Args>... >::value,"Some argument is wrong.");
  // The template paratmers Args must be std::decay'ed.
  static_assert(mystd::conjunction< std::is_same<Args,mystd::decay_t<Args>>...>::value,"Args are of wrong type");

  public:

  // tags used for various static_assert's and template overload selection:
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::false_type; // not a leaf of the expression tree.
  using LeveledComparison = std::true_type;
  using LeveledObject     = std::true_type;
  using LeveledObject_Base= std::false_type;
//  using DelayedDefaultFunctions = std::true_type;
  static constexpr unsigned int ApproxLevel = LazyFunction::ApproxLevel;
//  static_assert(ApproxLevel >0, "Approximation level is 0.");

// The type encapsulated is encoded in the function. Note that the line below is equivalent to
// template<unsigned int level> using ObjectAtLevel = typename LazyFunction::template EvalType<level>;
// This more complicated way gives better error messages (Notably, aliasing LazyFunction>::EvalType<level>
// might SFINAE-fail, removing functions from overloads. This way ensures we hit static_asserts
  template<unsigned int level> using ObjectAtLevel =
    typename std::conditional<level<=ApproxLevel,typename LazyFunction::template EvalType<level<=ApproxLevel?level:0>,LevelIsGood<level<=ApproxLevel>>::type;


#ifdef DEBUG_SIEVE_ALL_APPROX_AT_SAME_LEVEL
  static_assert(mystd::conjunction<mystd::bool_constant<ApproxLevel == ApproxLevelOf<Args>::value>...>::value,
    "All arguments must have at the same approximation level");
#else
  static_assert(mystd::conjunction<mystd::bool_constant<ApproxLevel <= ApproxLevelOf<Args>::value>...>::value,
    "All arguments must have at least the given approximation level");
#endif

  // EvalOnce means that calling eval_* might actually invalidate the data stored to / refered to.
  // If this is set, we
  // - store a non-const tuple of arguments
  // - optionally (implicitly) disable copying of SieveLazyEval
  // (via recursively disabling the automatic generation of the copy constructors / copy assignment)
  // (Note that disabling automagically generated constructors based on traits is quite non-trivial
  // in C++ and subject to many gotchas... The issue being that if you template these functions, the
  // compiler may no longer recognize them as "special".)
  // - eval may only be called on rvalue.
  // - may only call eval once. More precisely, calling eval<level> forbids any future calls to
  //   eval<level'> with level'>=level.
  // These restrictions might not be enforced.
  // (Note that these restriction could be slightly relaxed, separating the restrictions, but it
  //  seems not worth the hassle.)
  // This happens if the leaves of the trees encode pass by rvalue-semantics.

  // EvalOnce is either std::true_type or std::false_type, EvalOnce_v is either true or false.
  using EvalOnce = mystd::disjunction< typename Args::EvalOnce...>; //OR of Args
  static constexpr bool EvalOnce_v = EvalOnce::value;

  using TreeType = std::tuple<MaybeConst<!(Args::EvalOnce_v),Args>...>; // const-ness depends on EvalOnce.
  MaybeConst<!EvalOnce_v,TreeType> args;

  // This restriction can easily be removed. The issue is just that we want no default-constructor.
  // (at least for nargs > 0)
  SieveLazyEval() = delete; static_assert( sizeof...(Args)>0,"0-ary functions not supported yet." );

// Main Constructor: This takes a list of arguments and constructs a new lazy evaluation object
// that wraps calling function LazyFunction on those arguments.
// Technically, we "just" take any set of arguments and and use them to construct the
// std::tuple args from it.

  template<class... FunArgs,
  // There is a problem with clashing with the copy / move - constructors:
  // SieveLazyEval mave have automatically generated copy / move constructors, depending on whether
  // the std::tuple args is copyable/moveable, which depends on the types of Args...
  // This is as it should be.
  //
  // However, note that the (possibly default) copy constructor has signature
  // SieveLazyEval(SieveLazyEval const &).
  // The variadic template also contains (for FunArgs == SieveLazyEval &) a constructor with
  // signature SieveLazyEval(SieveLazyEval &). Now if we copy from a non-const object, the
  // variadic template is a better match.
  // Of course, this does not work (and the static_assert's will fire).
  // Hence, we deactivate this constructor for that case.
  TEMPL_RESTRICT_DECL
    (!( (sizeof...(Args)==1) &&
        (mystd::conjunction< std::is_same<mystd::decay_t<FunArgs>,SieveLazyEval>...>::value)
      )
    )>
  CONSTEXPR_IN_NON_DEBUG_TC explicit SieveLazyEval(FunArgs &&... fun_args)
    : args(std::forward<FunArgs>(fun_args)...) // the only line that actually does something here...
  {
    // Args and FunArgs must be the same (up to const and reference-ness)
    static_assert(sizeof...(Args) == sizeof...(FunArgs),"wrong number of arguments to constructor.");
    static_assert(mystd::conjunction< std::is_same<Args,mystd::decay_t<FunArgs>>...>::value,"Wrong type of arguments to constructor.");

    // If Args_i::EvalOnce is set, then fun_args_i must be passed via rvalue-reference
    // to the constructor of SieveLazyEval (possibly by explicit std::move).
    // If fun_args_i is passed via rvalue-ref, then (via the way forwarding reference capture works)
    // FunArgs_i itself is a non-refence type.
    // This means that we must have (FunArgs_i is non reference || Args_i::EvalOnce_v == false) for each i
    static_assert(mystd::conjunction< mystd::bool_constant<
        (std::is_reference<FunArgs>::value == false) || (Args::EvalOnce_v == false) >... >::value,
        "Use (possibly explicit) move semantics for Lazy objects that may encapsulate RValues." );
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Lazy function wrapper object for function " << LazyFunction::fun_name() << std::endl;
#endif
  }

  // the functions below have both const and non-const versions, which are identical. This
  // is required to propagate const-ness through the expression tree when evaluating.

  // do_eval is a helper function, whose sole purpose is to introduce a unnamed (dummy) variable,
  // whose template parameters iarg_0, iarg_1,... can be used to un-tuple args.
  // The parameter itself is then removed via inlining.

  template<unsigned int level, std::size_t... iarg>
  [[gnu::always_inline]] inline ObjectAtLevel<level> do_eval(mystd::index_sequence<iarg...>)
  {
    static_assert(sizeof...(iarg) == sizeof...(Args),"This cannot happen.");
#ifdef DEBUG_SIEVE_LAZY_TRACE_EVALS
    std::cout << "Calling function " << LazyFunction::fun_name() << " at level " << level << std::endl;
#endif
    return LazyFunction::template call<level>( std::get<iarg>(args).eval<level>()... );
  }

  template<unsigned int level, std::size_t... iarg>
  [[gnu::always_inline]] inline ObjectAtLevel<level> do_eval(mystd::index_sequence<iarg...>) const
  {
    static_assert(sizeof...(iarg) == sizeof...(Args),"This cannot happen.");
#ifdef DEBUG_SIEVE_LAZY_TRACE_EVALS
    std::cout << "Calling function " << LazyFunction::fun_name() << " at level " << level << std::endl;
#endif
    return LazyFunction::template call<level>( std::get<iarg>(args).eval<level>()... );
  }

// make_index_sequence generates a(n empty) struct MyIndexSeq<0,1,2,...>. This is used as a dummy paramter
// to select the correct version of do_eval_approx and allows to actually un-std::tuple the argument.
// (This is the least roundabout way of doing it and corresponds to the implementation of std::apply)

  template<unsigned int level>
  [[gnu::always_inline]] inline ObjectAtLevel<level> eval()
  {
    static_assert(level <= ApproxLevel, "");
    return do_eval<level>(mystd::make_index_sequence<sizeof...(Args)>{} );
  }

  template<unsigned int level>
  [[gnu::always_inline]] inline ObjectAtLevel<level> eval() const
  {
    static_assert(level <= ApproxLevel, "");
    return do_eval<level>(mystd::make_index_sequence<sizeof...(Args)>{} );
  }

  template<unsigned int level>
  [[gnu::always_inline]] inline ObjectAtLevel<level> get_value_at_level()
  {
    static_assert(level <= ApproxLevel, "");
    return eval<level>();
  }
  template<unsigned int level>
  [[gnu::always_inline]] inline ObjectAtLevel<level> get_value_at_level() const
  {
    static_assert(level <= ApproxLevel, "");
    return eval<level>();
  }

  // force_enable is just to make sure the argument to TEMPL_RESTRICT_DECL is not always false.
  template<unsigned int level, bool force_enable = false, TEMPL_RESTRICT_DECL(force_enable || (level<=ApproxLevel && level>0))>
  inline explicit operator ObjectAtLevel<level>()       { return eval<level>(); }
  template<unsigned int level, bool force_enable = false, TEMPL_RESTRICT_DECL(force_enable || (level<=ApproxLevel && level>0))>
  inline explicit operator ObjectAtLevel<level>() const { return eval<level>(); }

// not explicit (by design)
  inline operator ObjectAtLevel<0>() const  { return eval<0>(); }
  inline operator ObjectAtLevel<0>()        { return eval<0>(); }
  GAUSS_SIEVE_FORWARD_UNARY_MEMBER_FUNCTION_TO_DELAYED(get_norm2,call_get_norm2,SieveLazyEval)
};

/**
  LazyWrap* model the leaves of the expression trees.
  These classes are just wrappers around a reference to a *
  and expose the same interface as SieveLazyEval.

  Note that the functionality may be a bit more limited than that of SieveLazyEval (due to lazyness
  of the coder). As a workaround, you may use SieveLazyEval with an IdentityFunction.

  There are several different wrappers that differ in the calling conventions:

  We support (up to) four different calling conventions:
  - Const-Ref (CR): The wrapper is initialized with an lvalue and stores a const-reference.
                    upon eval_*, the lvalue used to initialize must still be valid.
                    The function used in eval_ must pass this by value or const-ref

  - Delayed-move(RV): This wrapper models delayed (explicit) rvalue-passing.
                      The wrapper is initialized with an lvalue.
                      Upon eval_*, the lvalue used to initialize must still be valid.
                      The function used in eval_ will be called on an rvalue.
                      (This means that the wrapper (lazily) takes ownership of the wrapped object,
                      and usage of eval_* might modify the object.
                      If the evaluation tree contains such an object, EvalOnce is set.)

  Untested, not used so far:

  - Move :  This models immediate rvalue-passing. The wrapper is initialized with an rvalue and
            directly move-stores the wrapped object. The function used in eval_* will be called on
            an rvalue.
            This also sets EvalOnce

  Untested, not used so far:

  - Capture:  This mode stores the values directly (not as reference). The wrapper may be initialized
              with either rvalues or lvalues. The argument will be passed to eval_* by value or const-ref.

  Additionally, we have a wrapper
  - Value:  Stores the value of a NON-LEVELED object itself. Calling eval<level> will return the
            object for any level. Used to encapsulate things like "x<<=2", where x is leveled, but
            2 is not. This is essentially equivalent to Capture except for the leveled-ness.

*/


/**
  Wrap via Const-Ref
*/

// joint code for each wrapper.

#define LAZY_WRAP_ALL \
public: \
static_assert(Has_LeveledObject_Base<CombinedObject>::value,"Can only wrap base leveled objects."); \
static constexpr unsigned int ApproxLevel = maxlevel; \
template<unsigned int level> using ObjectAtLevel = typename CombinedObject::template ObjectAtLevel<level>; \
static_assert(ApproxLevelOf<CombinedObject>::value == maxlevel,""); /* TODO: Relax this */ \
using IsLazyNode = std::true_type; \
using IsLazyLeaf = std::true_type; \
static constexpr bool EvalOnce_v = EvalOnce::value; \
using LeveledObject = std::true_type; \
using LeveledComparison = std::true_type; \
using LeveledObject_Base = std::false_type;

//using DelayedDefaultFunctions = std::true_type


template<class CombinedObject, unsigned int maxlevel = ApproxLevelOf<CombinedObject>::value >
struct LazyWrapCR
{
  using EvalOnce   = std::false_type;
  LAZY_WRAP_ALL;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapCR(CombinedObject const &init_combined):combined_ref(init_combined)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating CR Lazy Wrapper at level " << maxlevel << std::endl; // TODO: Level of CombType
#endif
  }
  LazyWrapCR(CombinedObject const &&) = delete; // No move!
  template<unsigned int level>
  inline ObjectAtLevel<level> const & eval() const
  {
    static_assert(level<=maxlevel,"Cannot evaluate at this level");
    return combined_ref.template access<level>();
  }
  template<unsigned int level>
  inline ObjectAtLevel<level> const & get_value_at_level() const
  {
    static_assert(level<=maxlevel,"Cannote evaluate at this level");
    return combined_ref.template access<level>();
  }
  CombinedObject const & combined_ref;
};

template<class CombinedObject, unsigned int maxlevel = ApproxLevelOf<CombinedObject>::value >
struct LazyWrapRV
{
  public:
  using EvalOnce   = std::true_type;
  LAZY_WRAP_ALL;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapRV(CombinedObject &init_combined):combined_ref(init_combined)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating RV Lazy Wrapper at level " << maxlevel << std::endl; //TODO!
#endif
  }
  LazyWrapRV(CombinedObject &&) = delete;
  template<unsigned int level> inline ObjectAtLevel<level> && eval()
  {
    static_assert(level <= maxlevel,"Cannot evaluate at this level");
    return std::move( combined_ref.template access<level>() );
  }
  template<unsigned int level> inline ObjectAtLevel<level> && get_value_at_level()
  {
    static_assert(level <= maxlevel,"Cannot evaluate at this level");
    return std::move( combined_ref.template access<level>() );
  }
  CombinedObject & combined_ref;
};

// untested:

template<class CombinedObject, unsigned int maxlevel  = ApproxLevelOf<CombinedObject>::value>
struct LazyWrapMove
{
  using EvalOnce   = std::true_type;
  LAZY_WRAP_ALL;
  LazyWrapMove(CombinedObject &&init_combined):combined_store(std::move(init_combined))
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Move Lazy Wrapper at level " << maxlevel << std::endl;
#endif
  }
  LazyWrapMove() = delete;
  LazyWrapMove(LazyWrapMove const &) = delete;
  LazyWrapMove(LazyWrapMove &&) = default;
  template<unsigned int level> inline ObjectAtLevel<level> && eval()
  {
    static_assert(level<= maxlevel,"Cannot evaluate at this level");
    return std::move( combined_store.template access<level>() );
  }
  template<unsigned int level> inline ObjectAtLevel<level> && get_value_at_level()
  {
    static_assert(level<= maxlevel,"Cannot evaluate at this level");
    return std::move( combined_store.template access<level>() );
  }
  CombinedObject combined_store;
};

// untested:

template<class CombinedObject, unsigned int maxlevel = ApproxLevelOf<CombinedObject>::value>
struct LazyWrapCapture
{
  using EvalOnce = std::false_type;
  LAZY_WRAP_ALL;
  CONSTEXPR_IN_NON_DEBUG_TC LazyWrapCapture(CombinedObject const &init_combined)
    :combined_store(init_combined)
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Capture Lazy Wrapper at level " << maxlevel << std::endl;
#endif
  }
  LazyWrapCapture(CombinedObject &&init_combined)
    :combined_store(std::move(init_combined))
  {
#ifdef DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
    std::cout << "Creating Capture Lazy Wrapper at level " << maxlevel << std::endl;
#endif
  }
  template<unsigned int level> inline ObjectAtLevel<level> const & eval() const
  {
    static_assert(level<= maxlevel,"Cannot evaluate at this level");
    return combined_store;
  }
  template<unsigned int level> inline ObjectAtLevel<level> const & get_value_at_level() const
  {
    static_assert(level<= maxlevel,"Cannot evaluate at this level");
    return combined_store;
  }
  CombinedObject const combined_store;
};

// untested:
// this one stores NON-Leveled objects by value. The purpose is to store secondary arguments.
template<class Object, unsigned int maxlevel>
struct LazyWrapValue
{
  using EvalOnce = std::false_type;
  // well, we *could* wrap those, but it's most likely an error. We would rather wrap access<0>().
  static_assert(!(Has_LeveledObject_Base<Object>::value),"Can NOT wrap base leveled objects.");
  static constexpr unsigned int ApproxLevel = maxlevel;
  template<unsigned int level> using ObjectAtLevel = Object;
  using IsLazyNode = std::true_type;
  using IsLazyLeaf = std::true_type;
  static constexpr bool EvalOnce_v = false;
  using LeveledObject = std::true_type;
  using LeveledComparison = std::true_type;
  using LeveledObject_Base = std::false_type;

  // No debug output. retaining constexpr is too important here.
  constexpr LazyWrapValue( Object const & obj) : stored_obj(obj) {}
  LazyWrapValue (Object && obj): stored_obj(std::move(obj)) {}
  Object const stored_obj;
  template<unsigned int level> constexpr inline Object get_value_at_level() const
  {
    static_assert(level<=maxlevel,"");
    return stored_obj;
  }
  template<unsigned int level> constexpr inline Object eval() const
  {
    static_assert(level<=maxlevel,"");
    return stored_obj;
  }
};

#undef LAZY_WRAP_ALL




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

#define GAUSS_SIEVE_LAZY_FUN(nargs_,namestring)                             \
  constexpr static unsigned int nargs = (nargs_);                           \
  using IsLazyFunction = std::true_type;                                    \
  static std::string fun_name() {return namestring ;}

/**
  Identity function. F(x) = x. Can be used to promote a
  LazyWrap* to a SieveLazyEval<ELP,Approxiomation,Lazy_Identity<...>,...>
  */


template<class CombinedType>
class Lazy_Identity
{
  public:
  GAUSS_SIEVE_LAZY_FUN(1,"Identity function")
  static constexpr unsigned int ApproxLevel = ApproxLevelOf<CombinedType>::value;
  template<unsigned int level> using EvalType = mystd::decay_t<typename CombinedType::template ObjectAtLevel<level>>;

  // Note: TEMPL_RESTRICT_DECL2's implicit and short-circuits template instantiation.
  // This means we will not instantiate EvalType<level> for too large level.
  template<unsigned int level, class Arg, TEMPL_RESTRICT_DECL2(
    mystd::bool_constant<level<=ApproxLevel>,
    std::is_same<EvalType<level>, mystd::decay_t<Arg>>
    )>
  inline static Arg call(Arg &&arg) { return std::forward<Arg>(arg); }
};

template<unsigned int maxlevel, template<unsigned int> class LHS, template<unsigned int> class RHS>
struct Lazy_ScalarProduct
{
  constexpr static unsigned int nargs = 2;
  using IsLazyFunction = std::true_type;
  static std::string fun_name() { return "scalar product";}
  static constexpr unsigned int ApproxLevel = maxlevel;
  template<unsigned int level> using EvalType=
  decltype( compute_sc_product( std::declval<LHS<level>>(), std::declval<RHS<level>>())  );

  static_assert(std::is_same<EvalType<0>,mystd::decay_t<EvalType<0>>>::value,"");
  template<unsigned int level, class Arg1, class Arg2>
  [[gnu::always_inline]] inline static EvalType<level> call(Arg1 &&arg1, Arg2 &&arg2)
  {
    static_assert(level <= maxlevel, "Cannot call at this level");
    static_assert(std::is_same<mystd::decay_t<Arg1>,LHS<level> >::value,"Invalid lhs");
    static_assert(std::is_same<mystd::decay_t<Arg2>,RHS<level> >::value,"Invalid rhs");
    return compute_sc_product(std::forward<Arg1>(arg1),std::forward<Arg2>(arg2));
  }
};


} //end namespace GaussSieve::LazyEval


namespace GaussSieve
{
//GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_FOR_DELAYED_OBJECTS(get_norm2)

GAUSS_SIEVE_LAZY_UNARY_FUNCTION_LEVEL_DETECTION    (abs,abs,using std::abs;)
GAUSS_SIEVE_LAZY_UNARY_FUNCTION_CREATE_LAZY_WRAPPER(abs,abs,using std::abs;)
GAUSS_SIEVE_LAZY_UNARY_FUNCTION_DIRECT_LAZY(abs,abs,std::true_type)

GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_LEVEL_DETECTION(get_norm2,get_norm2)
GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_CREATE_LAZY_WRAPPER(get_norm2,get_norm2)
GAUSS_SIEVE_LAZY_UNARY_FUNCTION_DIRECT_LAZY(call_get_norm2,get_norm2,std::true_type)

GAUSS_SIEVE_LAZY_BINARY_OP_TT_LEVEL_DETECTION(+,addition_tt)
GAUSS_SIEVE_LAZY_BINARY_OP_TT_CREATE_LAZY_WRAPPER(+,addition_tt)
GAUSS_SIEVE_LAZY_BINARY_FUNCTION_DIRECT_LAZY(operator +, addition_tt,std::true_type)

GAUSS_SIEVE_LAZY_BINARY_OP_TO_LEVEL_DETECTION(>>,rsh_to)
GAUSS_SIEVE_LAZY_BINARY_OP_TO_CREATE_LAZY_WRAPPER(>>,rsh_to)
GAUSS_SIEVE_LAZY_BINARY_FUNCTION_DIRECT_LAZY_VALUE(operator >>,rsh_to, std::true_type )
GAUSS_SIEVE_LAZY_BINARY_OP_TO_LEVEL_DETECTION(<<,lsh_to)
GAUSS_SIEVE_LAZY_BINARY_OP_TO_CREATE_LAZY_WRAPPER(<<,lsh_to)
GAUSS_SIEVE_LAZY_BINARY_FUNCTION_DIRECT_LAZY_VALUE(operator <<,lsh_to, std::true_type )


//GAUSS_SIEVE_LAZY_UNARY_FUNCTION_FOR_DELAYED_OBJECTS(abs, abs, using std::abs;)
//GAUSS_SIEVE_LAZY_BINARY_OP_FOR_DELAYED_OBJECTS_BOTH(+ ,operator_add_both_delayed)
}

/**
  We compare objects which have a public typedef LeveledComparison set to std::true_type in the
  following way:
  - We assume that there exists a get_value_at_level<level> function that returns something
    for which comparison is meaningful.
  We compare them level by level, starting with the highest level that both objects have.
  (As defined by ApproxLevelOf<Type>)
  If any comparison turns out false, the result is false.
  Note that A < B does not imply B > A in this manner.
  This is useful if comparisons may sometimes turn out wrong. For efficiency, write comparisons such
  that a *false* result occurs more often.

  Note:  If only one side of the comparison has LeveledComparison set, we only use the 0th level.
  Note2: Lattice Points themselves are not supposed to set LeveledComparison to std::true_type
         The reason is that we want LatticePoint1 < LatticePoint2 to compare exactly (by norm), in
         order to ensure that this ordering is actually well-behaved. Such comparison might occur
         for sorting lists. Use LatticePoint1.get_norm2() < LatticePoint2.get_norm2() and set
         LeveledComparison for the return type of get_norm2() if you want to use leveled comparisons
*/



namespace GaussSieve{


// clang-format off

// LHSLeveled / RHSLeveled are shorthands for testing if LHS / RHS have LeveledComparsion set.
// #defines are local (#undef'd after use)

CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(LeveledComparison,std::true_type,Has_LeveledComparison);
#define LHSLeveled (Has_LeveledComparison<mystd::decay_t<LHS>>::value)
#define RHSLeveled (Has_LeveledComparison<mystd::decay_t<RHS>>::value)

/**
  In order to avoid code duplication, we define a macro that is to be instantiated with
  OP in { <, >, <=, >= }. helper_name is just a name of a dummy helper class that needs to be
  different for each OP.
*/

#define DEFINE_DEFAULT_LEVELED_COMPARISON(OP, helper_name)                                        \
                                                                                                  \
/* Case 1: Only one the arguments has LeveledComparison set. We evaluate that at the 0th level    \
           and compare with the other side */                                                     \
template<class LHS, class RHS, TEMPL_RESTRICT_DECL(LHSLeveled && (!RHSLeveled))>                  \
CPP14CONSTEXPR inline bool operator OP (LHS &&lhs, RHS &&rhs)                                     \
{                                                                                                 \
  return std::forward<LHS>(lhs).template get_value_at_level<0>() OP std::forward<RHS>(rhs);       \
}                                                                                                 \
/* Case 2: Like case 1, just the other argument */                                                \
template<class LHS, class RHS, TEMPL_RESTRICT_DECL((!LHSLeveled) && RHSLeveled)>                  \
CPP14CONSTEXPR inline bool operator OP (LHS &&lhs, RHS &&rhs)                                     \
{                                                                                                 \
  return std::forward<LHS>(lhs) OP std::forward<RHS>(rhs).template get_value_at_level<0>();       \
}                                                                                                 \
                                                                                                  \
/* Case 3: Both arguments have LeveledComparison set. In this case, we do a kind-of for loop over \
           the levels, stating with the largest level. Now, observe that the level is a template  \
           argument and hence needs to be known at compile time, so we can't use a for-loop.      \
           Rather, we write it recursively via a compare_at_level<i> function (that may call a    \
           compare_at_level<i-1> function to simulate the for loop. Since C++ does not allow to   \
           partially specialize functions (such as our would-be compare_at_level<i>), we cannot   \
           easily write the (special) base case for the 0th level (lacking C++17 constexpr if)    \
           To work around that, we use a template helper *class* with a compare function, so      \
           we use helper_name<i>::compare rather than compare_at_level<i>.                        \
           The template arguments LHS, RHS are just forwarded through.                            \
  NOTE: Regarding move semantics, we assume that std::move(lhs).access_at_level<i> may only       \
        invalidate lhs at levels >=i. So in particular,                                           \
        std::move(lhs).access_at_level<1>() followed by std::move(lhs).access_at_level<2>()       \
        is valid.                                                                                 \
*/                                                                                                \
namespace ComparisonHelper /* Save the environment! Dont pollute your namespaces! */              \
{                                                                                                 \
template<class LHS, class RHS, unsigned int level> struct helper_name                             \
{                                                                                                 \
  static_assert(level>0,""); /* specialization for level 0 below */                               \
  static_assert(Has_LeveledComparison<mystd::decay_t<LHS>>::value,"");                            \
  static_assert(Has_LeveledComparison<mystd::decay_t<RHS>>::value,"");                            \
  helper_name() = delete;                                                                         \
  [[gnu::always_inline]]static inline bool compare(LHS &&lhs, RHS &&rhs)                          \
  {                                                                                               \
  /* note:  This uses std::forward (i.e. possibly std::move) twice. This is consistent with the   \
            move semantics specification of Leveled objects. See remark above                  */ \
  /* note2: The && short-circuits, so if the comparison at the current level is false, we don't   \
            try the next one.                                                                  */ \
    return (  std::forward<LHS>(lhs).template get_value_at_level<level>()                         \
          OP  std::forward<RHS>(rhs).template get_value_at_level<level>() ) &&                    \
            helper_name<LHS,RHS,level-1>::compare(std::forward<LHS>(lhs),std::forward<RHS>(rhs)); \
  }                                                                                               \
};                                                                                                \
template<class LHS, class RHS> struct helper_name<LHS,RHS,0>                                      \
{                                                                                                 \
  static_assert(Has_LeveledComparison<mystd::decay_t<LHS>>::value,"");                            \
  static_assert(Has_LeveledComparison<mystd::decay_t<RHS>>::value,"");                            \
  helper_name() = delete;                                                                         \
  [[gnu::always_inline]]CPP14CONSTEXPR static inline bool compare(LHS &&lhs, RHS &&rhs)           \
  {                                                                                               \
    return   std::forward<LHS>(lhs).template get_value_at_level<0>()                              \
          OP std::forward<RHS>(rhs).template get_value_at_level<0>();                             \
  }                                                                                               \
};                                                                                                \
} /* end of helper namespace */                                                                   \
                                                                                                  \
template<class LHS, class RHS, TEMPL_RESTRICT_DECL(LHSLeveled && RHSLeveled)>                     \
CPP14CONSTEXPR inline bool operator OP (LHS &&lhs, RHS &&rhs)                                     \
{                                                                                                 \
  return ComparisonHelper::helper_name                                                            \
  <                                                                                               \
    LHS,RHS,                                                                                      \
    /* std::max is not constexpr in C++11..., so we have to compute the max by ?: */              \
    ApproxLevelOf<mystd::decay_t<LHS>>::value >= ApproxLevelOf<mystd::decay_t<RHS>>::value ?      \
    ApproxLevelOf<mystd::decay_t<LHS>>::value : ApproxLevelOf<mystd::decay_t<RHS>>::value         \
  >::compare(std::forward<LHS>(lhs), std::forward<RHS>(rhs));                                     \
}

DEFINE_DEFAULT_LEVELED_COMPARISON(<, helper_less)
DEFINE_DEFAULT_LEVELED_COMPARISON(>, helper_greater)
DEFINE_DEFAULT_LEVELED_COMPARISON(<=,helper_leq)
DEFINE_DEFAULT_LEVELED_COMPARISON(>=,helper_geq)

/** TODO: == and !=, unsure whether to use approximations for this */

#undef DEFINE_DEFAULT_LEVELED_COMPARISON
#undef LHSLeveled
#undef RHSLeveled

// clang-format on


} // end namespace GaussSieve

#undef GAUSS_SIEVE_LAZY_FUN
#undef CONSTEXPR_IN_NON_DEBUG_TC

#endif // of include guard
