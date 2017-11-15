
#define GAUSS_SIEVE_FORWARD_UNARY_MEMBER_FUNCTION_TO_DELAYED(function_name)                        \
template<bool enabled=true>                                                                        \
[[gnu::always_inline]] inline auto function_name()  const &                                        \
-> typename std::enable_if<enabled,decltype(delayed_call_##function_name(*this))>::type            \
{                                                                                                  \
  return delayed_call_##function_name(*this);                                                      \
}                                                                                                  \
                                                                                                   \
template<bool enabled=true>                                                                        \
[[gnu::always_inline]] inline auto function_name()  &&                                             \
-> typename std::enable_if<enabled,decltype(delayed_call_##function_name(std::move(*this)))>::type \
{                                                                                                  \
  return delayed_call_##function_name(std::move(*this));                                           \
}

/***************************************************************************************************
 * Creates a CanCall_function_id struct to detect up to which level function_name(arg) is valid.   *
 *                                                                                                 *
 * Notably, it creates a CanCall_function_id<LeveledObjects> struct into the subnamespace          *
 * LazyEval::function_id_helper_namespace. It contains a public static member function template    *
 * has_function_id_maxlevel<maxlevel>() which returns the maximal level, for which                 *
 * function_name(LeveledObjects<level>) is valid, capped at maxlevel.                              *
 *                                                                                                 *
 * Note: function_id should be a unique (in this context) identifier and is used to compose the    *
 *       names of the helper structures via ##. By default, use function_id == function_name.      *
 *       Use the same function_id in other macros to make them work together.                      *
 *       namespace_injection is used to inject names from other namespaces to affect validity      *
 *       of function_name(LeveledObjects<level>) expression (needed for std::abs and std::swap)    *
 **************************************************************************************************/

#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_LEVEL_DETECTION(function_name,function_id,namespace_injection)\
namespace LazyEval::function_id##_helper_namespace{                                                \
namespace_injection                                                                                \
                                                                                                   \
/** allows to detect whether function_name(arg) is valid at various levels */                      \
template<template<unsigned int> class TestObjectAtLevel>                                           \
struct CanCall_##function_id                                                                       \
{                                                                                                  \
  private:                                                                                         \
  /** has_fun_name<level>(0) returns true if fun_name(TestObjectAtLevel<level>) is valid */        \
  template<unsigned int level> constexpr static bool has_##function_id(...) { return false; }      \
  template<unsigned int level> /* better match, but only if the argument inside declval is valid */\
  constexpr static auto has_##function_id(int) /* argument is unused, but needed */                \
  -> mystd::decay_t<decltype ( static_cast<void>(                                                  \
    function_name(std::declval<TestObjectAtLevel<level> >() )                                      \
  ), bool(true) )>                                                                                 \
  { return true; }                                                                                 \
  static_assert(has_##function_id<0>(1),"");                                                       \
                                                                                                   \
  /** has_fun_name_upto<level>() returns true if function_name(TestObjectAtLevel<level'>)          \
      is valid for all level' <= level */                                                          \
  /* all ?:'s are just to workaround limitations of C++...; */                                     \
  template<unsigned int level> constexpr static bool has_##function_id##_upto()                    \
  {                                                                                                \
    return has_##function_id<level>(1) &&                                                          \
           ( (level==0) ? true : has_##function_id##_upto<( (level>0)?level-1:0)>() );             \
  }                                                                                                \
  /** has_function_name_maxlevel<maxlevel> returns min(maxlevel, maximal level for which           \
   function_name(TestObjectAtLevel<level>)  is valid ) */                                          \
  /* at least the first ?: is actually meaningful. */                                              \
  public:                                                                                          \
  template<unsigned int maxlevel> constexpr static unsigned int has_##function_id##_maxlevel()     \
  {                                                                                                \
    return has_##function_id##_upto<maxlevel>() ?                                                  \
       maxlevel :                                                                                  \
      ( (maxlevel==0) ? 0 : has_##function_id##_maxlevel< ( (maxlevel>0) ? maxlevel-1 : 0)>() );   \
  }                                                                                                \
};                                                                                                 \
} /* end namespace */

/***************************************************************************************************
 * Creates a Lazy_function_name wrapper to be used as a LazyFunction is SieveLazyEval              *                                                                                      *
 *                                                                                                 *
 * Notably, creates a struct Lazy_function_id into the LazyEval::function_id_helper_namespace ns   *
 * Lazy_function_id is templated with a maximal level and an ArgAtLevel class template, where      *
 * maxlevel is the maximal level supported and ArgAtLevel<level> is the argument to function_name  *
 * when calling at the given level.                                                                *
 *                                                                                                 *
 * The class encapsulates a call to function_name by exposing a public static call<level>(arg)     *
 * member function template, which is of the form required by SieveLazyEval.                       *
 **************************************************************************************************/

#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_CREATE_LAZY_WRAPPER(function_name, function_id, namespace_injection)\
namespace LazyEval::function_id##_helper_namespace{                                                \
template<unsigned int maxlevel, template<unsigned int> class ArgAtLevel>                           \
struct Lazy_##function_id                                                                          \
{                                                                                                  \
  constexpr static unsigned int nargs = 1;                                                         \
  using IsLazyFunction = std::true_type;                                                           \
  static std::string fun_name() {return #function_name;}                                           \
  static constexpr unsigned int ApproxLevel = maxlevel;                                            \
  template<unsigned int level> using EvalType = mystd::decay_t<                                    \
    decltype( function_name(std::declval<ArgAtLevel<level>>() ))                                   \
    >;                                                                                             \
                                                                                                   \
  template<unsigned int level, class Arg, TEMPL_RESTRICT_DECL2(                                    \
    std::is_same<ArgAtLevel<level>, mystd::decay_t<Arg>>                                           \
    )>                                                                                             \
  inline static EvalType<level> call(Arg &&arg)                                                    \
  {                                                                                                \
    static_assert(level <= maxlevel,"Cannot call at this level");                                  \
    static_assert(std::is_same<ArgAtLevel<level>,mystd::decay_t<Arg>>::value,"Invalid argument");  \
    return function_name(std::forward<Arg>(arg));                                                  \
  }                                                                                                \
};                                                                                                 \
} /* end LazyEval::function_name_helper_namespace */

/***************************************************************************************************
 * Defines the function function_name(Arg) for Args satisfying a certain predicate by creating     *
 * a SieveLazyEval object. The call is performed directly (without calling convention wrappers     *
 * such as the LazyWrap* wrappers from Lazy.h).                                                    *
 *                                                                                                 *
 * Notably, if the predicate on Args (supplied via the third and further macro argument as comma-  *
 * separated list of types exposing a ::value with an implicit AND between the arguments) holds    *
 * true, function_name(Args) is defined by creating a SieveLazyEval object wrapping the call.      *
 * Note that we include a test for LazyNode==std::true_type automatically, because it would not    *
 * work otherwise. If you want no further conditions, set the third argument to                    *
 * std::true_type                                                                                  *
 *                                                                                                 *
 * The second argument function_id has to be prepared by creating a Lazy_function and a CanCall    *
 * object by the macros above.                                                                     *
 *                                                                                                 *
 * Note: The third argument is eventually fed into a variadic template. In particular, commas in   *
 *       the argument (which are always a problem for variadic macros taking types) are fine.      *
 **************************************************************************************************/

#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_DIRECT_LAZY(function_name, function_id, ...)               \
template<class Arg, TEMPL_RESTRICT_DECL2(Has_IsLazyNode<Arg>,__VA_ARGS__)>                         \
auto function_name(Arg &&arg)                                                                      \
/* decltype(auto) would be REALLY useful */                                                        \
-> LazyEval::SieveLazyEval<                                                                        \
    LazyEval::function_id##_helper_namespace::Lazy_##function_id< /* LazyFun */                    \
      LazyEval::function_id##_helper_namespace::CanCall_##function_id<                             \
          mystd::decay_t<Arg>::template ObjectAtLevel>::template has_##function_id##_maxlevel /* newarglevel */ \
      <                                                                                            \
        ApproxLevelOf<mystd::decay_t<Arg>>::value /* arglevel */                                   \
      >(),                                                                                         \
      mystd::decay_t<Arg>::template ObjectAtLevel>,                                                \
    mystd::decay_t<Arg> >                                                                          \
{                                                                                                  \
  /* approximation level of argument: */                                                           \
  static constexpr unsigned int arglevel = ApproxLevelOf<mystd::decay_t<Arg>>::value;              \
  /* Use CanCall_function_id to detect maximal level for which leveled evaluation is possible */   \
    static constexpr unsigned int newarglevel=                                                     \
    LazyEval::function_id##_helper_namespace::CanCall_##function_id<                               \
      mystd::decay_t<Arg>::template ObjectAtLevel>::template has_##function_id##_maxlevel<arglevel>(); \
  /* Use Lazy_function_id with the correct template arguments, capped at level newarglevel */      \
  using LazyFun = LazyEval::function_id##_helper_namespace::Lazy_##function_id<                    \
      newarglevel, mystd::decay_t<Arg>::template ObjectAtLevel>;                                   \
  /* Return type is a SieveLazyEval object with this function and Arg as argument type */          \
  using RetType = LazyEval::SieveLazyEval<LazyFun, mystd::decay_t<Arg>>;                           \
  return RetType{std::forward<Arg>(arg)};                                                          \
}




/**

  OUTDATED INFO:

  Macros that actually perform the lazy evaluation.
  GAUSS_SIEVE_LAZY_UNARY_FUNCTION_FOR_DELAYED_OBJECTS(fun,fun_id,inject)
  creates a function fun(Arg) in the current namespace for
  Args that have Has_DelayedDefaultFunctions set.
  It is assumed that Arg is a LeveledObject.
  fun(Arg) is constructed by detecting until which level fun(arg_at_level<level>) is defined for its
  individual levels and automagically creates a delayed objects that wraps an evaluation for this
  number of levels.
  fun_id is any unique (in that context) identifier and may usually be chosen identical to fun.
  (It is used to construct names of helper objects via preprocessor ##.
  Only if fun is of the form operator-, operator+ etc. should one choose fun_id different)

  The variant
  GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_FOR_DELAYED_OBJECTS(fun)
  creates a function delayed_call_fun(Arg) in the current namespace that wraps around
  a call to Arg.fun(). It is otherwise identical to the above.
  (Note that the "unary" includes the implicit *this, notationally)
*/


#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_FOR_DELAYED_OBJECTS_FORWARD_DECLARE(function_name, function_id) \
namespace LazyEval::function_id##_helper_namespace \
{                                    \
template<template<unsigned int> class TestObjectAtLevel> struct CanCall_##function_id; \
template<unsigned int maxlevel, template<unsigned int> class ArgAtLevel> struct Lazy_##function_id; \
} \
template<class LazyObject,                                                                         \
         TEMPL_RESTRICT_DECL2(LazyEval::Has_DelayedDefaultFunctions<mystd::decay_t<LazyObject>>)>  \
\
auto function_name(LazyObject &&obj)                                                               \
-> LazyEval::SieveLazyEval<                                                                        \
    LazyEval::function_id##_helper_namespace::Lazy_##function_id< /* LazyFun */                    \
      LazyEval::function_id##_helper_namespace::CanCall_##function_id<                             \
          mystd::decay_t<LazyObject>::template ObjectAtLevel>::template has_##function_id##_maxlevel /* newarglevel */ \
      <                                                                                            \
        ApproxLevelOf<mystd::decay_t<LazyObject>>::value /* arglevel */                            \
      >(),                                                                                         \
      mystd::decay_t<LazyObject>::template ObjectAtLevel>,                                         \
    mystd::decay_t<LazyObject> >;                                                                  \


/*************************
  member function variant
*************************/


#define GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_FOR_DELAYED_OBJECTS(function_name)                  \
namespace LazyEval::function_name##_helper_namespace{                                              \
                                                                                                   \
/** allows to detect whether arg<level>.function_name() is valid at various levels */              \
template<template<unsigned int> class TestObjectAtLevel>                                           \
struct CanCall_##function_name                                                                     \
{                                                                                                  \
  private:                                                                                         \
  /** has_fun_name<level>(0) returns true if TestObjectAtLevel<level>.function_name() is valid */  \
  template<unsigned int level> constexpr static bool has_##function_name(...) { return false; }    \
  template<unsigned int level> /* better match, but only if the argument inside declval is valid */\
  constexpr static auto has_##function_name(int) /* argument is unused, but needed */              \
  -> mystd::decay_t<decltype ( static_cast<void>(                                                  \
     std::declval<TestObjectAtLevel<level> >().function_name()                                     \
  ), bool(true) )>                                                                                 \
  { return true; }                                                                                 \
  /*static_assert(has_##function_name<0>(1),""); */                                                    \
                                                                                                   \
  /** has_fun_name_upto<level>() returns true if function_name(TestObjectAtLevel<level'>)          \
      is valid for all level' <= level */                                                          \
  /* all ?:'s are just to workaround limitations of C++...; */                                     \
  template<unsigned int level> constexpr static bool has_##function_name##_upto()                  \
  {                                                                                                \
    return has_##function_name<level>(1) &&                                                        \
           ( (level==0) ? true : has_##function_name##_upto<( (level>0)?level-1:0)>() );           \
  }                                                                                                \
  /** has_function_name_maxlevel<maxlevel> returns min(maxlevel, maximal level for which           \
   function_name(TestObjectAtLevel<level>)  is valid ) */                                          \
  /* at least the first ?: is actually meaningful. */                                              \
  public:                                                                                          \
  template<unsigned int maxlevel> constexpr static unsigned int has_##function_name##_maxlevel()   \
  {                                                                                                \
    return has_##function_name##_upto<maxlevel>() ?                                                \
       maxlevel :                                                                                  \
      ( (maxlevel==0) ? 0 : has_##function_name##_maxlevel< ( (maxlevel>0) ? maxlevel-1 : 0)>() ); \
  }                                                                                                \
};                                                                                                 \
                                                                                                   \
/** Actual Lazy Function wrapper */                                                                \
template<unsigned int maxlevel, template<unsigned int> class ArgAtLevel>                           \
struct Lazy_##function_name                                                                        \
{                                                                                                  \
  GAUSS_SIEVE_LAZY_FUN(1,#function_name)                                                           \
  static constexpr unsigned int ApproxLevel = maxlevel;                                            \
  template<unsigned int level> using EvalType = mystd::decay_t<                                    \
    decltype( std::declval<ArgAtLevel<level>>().function_name() )                                  \
    >;                                                                                             \
                                                                                                   \
  template<unsigned int level, class Arg>                                                          \
  inline static EvalType<level> call(Arg &&arg)                                                    \
  {                                                                                                \
    static_assert(level<=maxlevel, "Cannot call at this level");                                   \
    static_assert(std::is_same<ArgAtLevel<level>,mystd::decay_t<Arg>>::value,"wrong argument");    \
    return std::forward<Arg>(arg).function_name();                                                 \
  }                                                                                                \
};                                                                                                 \
                                                                                                   \
} /* end namespace LazyEval::function_name_helper_namespace */                                     \
                                                                                                   \
/** actual function: declares delayed_call_function_name(Arg) for all Args for which               \
    Has_DelayedDefaultFunctions is set. It detects the maximal level, for which fun_name is defined\
     for all approximations (static_asserting it is defined at least for level 0) and creates an   \
     object wrapping the  evaluation */                                                            \
template<class LazyObject,                                                                         \
         TEMPL_RESTRICT_DECL2(Has_DelayedDefaultFunctions<mystd::decay_t<LazyObject>>)>  \
auto delayed_call_##function_name(LazyObject &&obj)                                                \
/* decltype(auto) would be REALLY useful */                                                        \
-> LazyEval::SieveLazyEval<                                                                        \
    LazyEval::function_name##_helper_namespace::Lazy_##function_name< /* LazyFun */                \
      LazyEval::function_name##_helper_namespace::CanCall_##function_name<                         \
          mystd::decay_t<LazyObject>::template ObjectAtLevel>::template has_##function_name##_maxlevel /* newarglevel */ \
      <                                                                                            \
        ApproxLevelOf<mystd::decay_t<LazyObject>>::value /* arglevel */                            \
      >(),                                                                                         \
      mystd::decay_t<LazyObject>::template ObjectAtLevel>,                                         \
    mystd::decay_t<LazyObject> >                                                                   \
{                                                                                                  \
  static constexpr unsigned int arglevel = ApproxLevelOf<mystd::decay_t<LazyObject>>::value;       \
  static constexpr unsigned int newarglevel=                                                       \
    LazyEval::function_name##_helper_namespace::CanCall_##function_name<                           \
        mystd::decay_t<LazyObject>::template ObjectAtLevel>::template has_##function_name##_maxlevel<arglevel>(); \
  using LazyFun = LazyEval::function_name##_helper_namespace::Lazy_##function_name<newarglevel, mystd::decay_t<LazyObject>::template ObjectAtLevel>; \
  using RetType = LazyEval::SieveLazyEval<LazyFun, mystd::decay_t<LazyObject>>;                    \
  return RetType{std::forward<LazyObject>(obj)};                                                   \
}


#define GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_FOR_DELAYED_OBJECTS_FORWARD_DECLARE(function_name)  \
namespace LazyEval::function_name##_helper_namespace{                                              \
template<template<unsigned int> class TestObjectAtLevel>  struct CanCall_##function_name;          \
template<unsigned int maxlevel, template<unsigned int> class ArgAtLevel> struct Lazy_##function_name;\
} \
template<class LazyObject,                                                                         \
         TEMPL_RESTRICT_DECL2(Has_DelayedDefaultFunctions<mystd::decay_t<LazyObject>>)>  \
auto delayed_call_##function_name(LazyObject &&obj)                                                \
/* decltype(auto) would be REALLY useful */                                                        \
-> LazyEval::SieveLazyEval<                                                                        \
    LazyEval::function_name##_helper_namespace::Lazy_##function_name< /* LazyFun */                \
      LazyEval::function_name##_helper_namespace::CanCall_##function_name<                         \
          mystd::decay_t<LazyObject>::template ObjectAtLevel>::template has_##function_name##_maxlevel \
      <                                                                                            \
        ApproxLevelOf<mystd::decay_t<LazyObject>>::value /* arglevel */                            \
      >(),                                                                                         \
      mystd::decay_t<LazyObject>::template ObjectAtLevel>,                                         \
    mystd::decay_t<LazyObject> >;

/**
 Binary functions:
*/

#define GAUSS_SIEVE_LAZY_BINARY_OP_FOR_DELAYED_OBJECTS_BOTH(op_name, function_id)    \
namespace LazyEval::function_id##_helper_namespace{                                                \
                                                                                                   \
/* allows to detect whether function_name(arg1,arg2) is valid at various levels */                \
template<template<unsigned int> class Arg1,template<unsigned int> class Arg2>                      \
struct CanCall_##function_id                                                                       \
{                                                                                                  \
  private:                                                                                         \
  /* has_function_id<level>(0) returns true if function_id(Arg1<level>,Arg2<level>) is valid */      \
  template<unsigned int level> constexpr static bool has_##function_id(...) { return false; }      \
  template<unsigned int level> /* better match, but only if the argument inside declval is valid */\
  constexpr static auto has_##function_id(int) /* argument is unused, but needed */                \
  -> mystd::decay_t<decltype ( static_cast<void>(                                                  \
    std::declval<Arg1<level> >() op_name std::declval<Arg2<level> >()                       \
  ), bool(true) )>                                                                                 \
  { return true; }                                                                                 \
  static_assert(has_##function_id<0>(0),"");                                                       \
                                                                                                   \
  /* has_fun_name_upto<level>() returns true if function_name(TestObjectAtLevel<level'>)          \
      is valid for all level' <= level */                                                          \
  /* all ?:'s are just to workaround limitations of C++...; */                                     \
  template<unsigned int level> constexpr static bool has_##function_id##_upto()                    \
  {                                                                                                \
    return has_##function_id<level>(1) &&                                                          \
           ( (level==0) ? true : has_##function_id##_upto<( (level>0)?level-1:0)>() );             \
  }                                                                                                \
  /* has_function_name_maxlevel<maxlevel> returns min(maxlevel, maximal level for which           \
   function_name(TestObjectAtLevel<level>)  is valid ) */                                          \
  /* at least the first ?: is actually meaningful. */                                              \
  public:                                                                                          \
  template<unsigned int maxlevel> constexpr static unsigned int has_##function_id##_maxlevel()     \
  {                                                                                                \
    return has_##function_id##_upto<maxlevel>() ?                                                  \
       maxlevel :                                                                                  \
      ( (maxlevel==0) ? 0 : has_##function_id##_maxlevel< ( (maxlevel>0) ? maxlevel-1 : 0)>() );   \
  }                                                                                                \
};                                                                                                 \
                                                                                                   \
/* Actual Lazy Function wrapper */                                                                \
template<unsigned int maxlevel, template<unsigned int> class Arg1AtLevel, template<unsigned int> class Arg2AtLevel> \
struct Lazy_##function_id                                                                          \
{                                                                                                  \
  GAUSS_SIEVE_LAZY_FUN(2,#op_name)                                                           \
  static constexpr unsigned int ApproxLevel = maxlevel;                                            \
  template<unsigned int level> using EvalType = mystd::decay_t<                                    \
    decltype( std::declval<Arg1AtLevel<level>>() op_name std::declval<Arg2AtLevel<level>>() )      \
    >;                                                                                             \
                                                                                                   \
  template<unsigned int level, class Arg1, class Arg2, TEMPL_RESTRICT_DECL2(                       \
    std::is_same<Arg1AtLevel<level>, mystd::decay_t<Arg1>>,                                        \
    std::is_same<Arg2AtLevel<level>, mystd::decay_t<Arg2>>                                         \
    )>                                                                                             \
  inline static EvalType<level> call(Arg1 &&arg1, Arg2 &&arg2)                                     \
  {                                                                                                \
    static_assert(level<=maxlevel,"cannot evaluate at this level"); \
    return std::forward<Arg1>(arg1) op_name std::forward<Arg2>(arg2);                       \
  }                                                                                                \
};                                                                                                 \
                                                                                                   \
} /* end LazyEval::function_name_helper_namespace */                                               \
                                                                                                   \
/* actual function: declares function_name(Arg1,Arg2) for all Arg1, Arg2s for which               \
    Has_DelayedDefaultFunctions is set for BOTH Arg1, Arg2.                                        \
    It detects the maximal level, for which fun_name is defined for all approximations             \
    (static_asserting it is defined at least for level 0) and creates an object wrapping the       \
    evaluation */                                                                                  \
template<class LazyObject1, class LazyObject2,                                                     \
         TEMPL_RESTRICT_DECL2(Has_DelayedDefaultFunctions<mystd::decay_t<LazyObject1>>,  \
                              Has_DelayedDefaultFunctions<mystd::decay_t<LazyObject2>>)> \
auto operator op_name(LazyObject1 &&obj1, LazyObject2 &&obj2)                                         \
/* decltype(auto) would be REALLY useful. It's just RetType from the function definition below after plugging everything in */ \
-> LazyEval::SieveLazyEval<                                                                        \
    LazyEval::function_id##_helper_namespace::Lazy_##function_id< /* LazyFun */                    \
      LazyEval::function_id##_helper_namespace::CanCall_##function_id<                             \
          mystd::decay_t<LazyObject1>::template ObjectAtLevel, mystd::decay_t<LazyObject2>::template ObjectAtLevel>::template has_##function_id##_maxlevel /* newarglevel */ \
      <                                                                                            \
        mystd::constexpr_min( ApproxLevelOf<mystd::decay_t<LazyObject1>>::value,ApproxLevelOf<mystd::decay_t<LazyObject2>>::value) /* arglevel */                            \
      >(),                                                                                         \
      mystd::decay_t<LazyObject1>::template ObjectAtLevel,                                         \
      mystd::decay_t<LazyObject2>::template ObjectAtLevel>,                                        \
    mystd::decay_t<LazyObject1>,                                                                   \
    mystd::decay_t<LazyObject2> >                                                                  \
{                                                                                                  \
  static constexpr unsigned int arglevel = mystd::constexpr_min(ApproxLevelOf<mystd::decay_t<LazyObject1>>::value,ApproxLevelOf<mystd::decay_t<LazyObject2>>::value);       \
  static constexpr unsigned int newarglevel=                                                       \
    LazyEval::function_id##_helper_namespace::CanCall_##function_id<                               \
        mystd::decay_t<LazyObject1>::template ObjectAtLevel, mystd::decay_t<LazyObject2>::template ObjectAtLevel>::template has_##function_id##_maxlevel<arglevel>(); \
  using LazyFun = LazyEval::function_id##_helper_namespace::Lazy_##function_id<newarglevel, mystd::decay_t<LazyObject1>::template ObjectAtLevel, mystd::decay_t<LazyObject2>::template ObjectAtLevel>; \
  using RetType = LazyEval::SieveLazyEval<LazyFun, mystd::decay_t<LazyObject1>,mystd::decay_t<LazyObject2> >;  \
  return RetType{std::forward<LazyObject1>(obj1),std::forward<LazyObject2>(obj2)};                 \
}

/*****************
 No forward declaration for now...
*****************/

