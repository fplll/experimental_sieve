
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

/**
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


// function_id is typically identical to function_name, except if function_name contains is
// sth. that might cause trouble if concatenated to form a new identifier.
//  e.g. if function_name == operator+, we set function_id as operator_plus
#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_FOR_DELAYED_OBJECTS(function_name, function_id, namespace_injection)    \
namespace LazyEval::function_id##_helper_namespace{                                    \
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
      ( (maxlevel==0) ? 0 : has_##function_id##_maxlevel< ( (maxlevel>0) ? maxlevel-1 : 0)>() ); \
  }                                                                                                \
};                                                                                                 \
                                                                                                   \
/** Actual Lazy Function wrapper */                                                                \
template<unsigned int maxlevel, template<unsigned int> class ArgAtLevel>                           \
struct Lazy_##function_id                                                                          \
{                                                                                                  \
  GAUSS_SIEVE_LAZY_FUN(1,#function_name)                                                           \
  static constexpr unsigned int ApproxLevel = maxlevel;                                            \
  template<unsigned int level> using EvalType = mystd::decay_t<                                    \
    decltype( function_name(std::declval<ArgAtLevel<level>>() ))                                        \
    >;                                                                                             \
                                                                                                   \
  template<unsigned int level, class Arg, TEMPL_RESTRICT_DECL2(                                    \
    mystd::bool_constant<level<=ApproxLevel>,                                                      \
    std::is_same<ArgAtLevel<level>, mystd::decay_t<Arg>>                                           \
    )>                                                                                             \
  inline static EvalType<level> call(Arg &&arg) { return function_name(std::forward<Arg>(arg)); }  \
};                                                                                                 \
                                                                                                   \
} /* end LazyEval::function_name_helper_namespace */                                               \
                                                                                                   \
/** actual function: declares function_name(Arg) for all Args for which Has_DelayedDefaultFunctions\
    is set. It detects the maximal level, for which fun_name is defined for all approximations     \
    (static_asserting it is defined at least for level 0) and creates an object wrapping the       \
    evaluation */                                                                                  \
template<class LazyObject,                                                                         \
         TEMPL_RESTRICT_DECL2(LazyEval::Has_DelayedDefaultFunctions<mystd::decay_t<LazyObject>>)>  \
auto function_name(LazyObject &&obj)                                                               \
/* decltype(auto) would be REALLY useful */                                                        \
-> LazyEval::SieveLazyEval<                                                                        \
    LazyEval::function_id##_helper_namespace::Lazy_##function_id< /* LazyFun */                  \
      LazyEval::function_id##_helper_namespace::CanCall_##function_id<                             \
          mystd::decay_t<LazyObject>::template ObjectAtLevel>::template has_##function_id##_maxlevel /* newarglevel */ \
      <                                                                                            \
        ApproxLevelOf<mystd::decay_t<LazyObject>>::value /* arglevel */                            \
      >(),                                                                                         \
      mystd::decay_t<LazyObject>::template ObjectAtLevel>,                                         \
    mystd::decay_t<LazyObject> >                                                                   \
{                                                                                                  \
  static constexpr unsigned int arglevel = ApproxLevelOf<mystd::decay_t<LazyObject>>::value;       \
  static constexpr unsigned int newarglevel=                                                       \
    LazyEval::function_id##_helper_namespace::CanCall_##function_id<                               \
        mystd::decay_t<LazyObject>::template ObjectAtLevel>::template has_##function_id##_maxlevel<arglevel>(); \
  using LazyFun = LazyEval::function_id##_helper_namespace::Lazy_##function_id<newarglevel, mystd::decay_t<LazyObject>::template ObjectAtLevel>; \
  using RetType = LazyEval::SieveLazyEval<LazyFun, mystd::decay_t<LazyObject>>;                    \
  return RetType{std::forward<LazyObject>(obj)};                                                   \
}

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
         TEMPL_RESTRICT_DECL2(LazyEval::Has_DelayedDefaultFunctions<mystd::decay_t<LazyObject>>)>  \
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
         TEMPL_RESTRICT_DECL2(LazyEval::Has_DelayedDefaultFunctions<mystd::decay_t<LazyObject>>)>  \
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
