/***************************************************************************************************
 * Creates a CanCall_function_id struct to detect up to which level function_name(args) is valid.  *
 *                                                                                                 *
 * Notably, it creates a CanCall_function_id<LeveledObject(s)> struct into the subnamespace        *
 * LazyEval::function_id_helper_namespace. It contains a public static constexpr member function   *
 * template has_function_id_maxlevel<maxlevel>() which returns the maximal level, for which        *
 * function_name() is valid, capped at maxlevel.                                                   *
 *                                                                                                 *
 * For GAUSS_SIEVE_LAZY_UNARY_FUNCTION_LEVEL_DETECTION, this means                                 *
 * validity of function_name(LeveledObject<level>)                                                 *
 * For GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_LEVEL_DETECTION, it means                            *
 * validity of LeveledObject<level>.function_name()                                                *
 * ("unary" always includes the implict *this argument here)                                       *
 * For GAUSS_SIEVE_LAZY_BINARY_OP_**_LEVEL_DETECTION, this means                                   *
 * validity of [Leveled]Object1[<level>] op [Leveled]Object2[<level>], where ** comes in 3 variants*
 * that determine whether Object1/2 are leveled.                                                   *
 *                                                                                                 *
 * Note: function_id should be a unique (in this context) identifier and is used to compose the    *
 *       names of the helper structures via ##. By default, use function_id == function_name.      *
 *       (except for operators)                                                                    *
 *       Use the same function_id in other macros to make them work together.                      *
 *                                                                                                 *
 *       namespace_injection is used to inject names from other namespaces to affect validity      *
 *       of function_name(LeveledObjects<level>) expression (needed for std::abs and std::swap)    *
 **************************************************************************************************/

/***************************************************************************************************
 * Helper macro to unify the unary / binary, member function/non-member function / operator cases  *
 *                                                                                                 *
 * the 4th argument to the macro is the template arguments to the resulting CanCall objects. These *
 * shall include at least one class template, templated by level. Note that these are a variadic,  *
 * because they may include a list of types, hence unshielded commas.                              *
 * expression is the expression whose validity is to be checked. It needs to depend on "level".    *
 * function_id a unique name (used to ## )                                                         *
 * namespace_injection injects into the helper namespace to affect visibility.                     *
 *                                                                                                 *
 * Do not use this macro directly, use the macros below defined in terms of this one instead.      *
 **************************************************************************************************/

#define GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION(expression,function_id,namespace_injection,...)\
namespace LazyEval::function_id##_helper_namespace{                                                \
namespace_injection                                                                                \
                                                                                                   \
/** allows to detect whether function_name(arg) is valid at various levels */                      \
template<__VA_ARGS__>                                                                                \
struct CanCall_##function_id                                                                       \
{                                                                                                  \
  private:                                                                                         \
  /** has_fun_name<level>(0) returns true if expression is valid (expression shall involve )*/        \
  template<unsigned int level> constexpr static bool has_##function_id(...) { return false; }      \
  template<unsigned int level> /* better match, but only if the argument inside declval is valid */\
  constexpr static auto has_##function_id(int) /* argument is unused, but needed */                \
  -> mystd::decay_t<decltype ( static_cast<void>( expression ), bool(true) )>                      \
  { return true; }                                                                                 \
  static_assert(has_##function_id<0>(1)==true,"");                                                 \
                                                                                                   \
  /** has_fun_name_upto<level>() returns true if expression is valid for all level' <= level */    \
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

/**
  For forward declarations.
*/

#define GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION_FORWARD_DECLARE(function_id,...)                  \
namespace LazyEval::function_id##_helper_namespace{                                                \
template<__VA_ARGS__> struct CanCall_##function_id;                                                \
}

/*******************
  Use these macros:
*******************/

/** Detect function_name(arg) */
#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_LEVEL_DETECTION(function_name,function_id,namespace_injection)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION( function_name(std::declval<Arg1<level>>()),function_id,namespace_injection, template<unsigned int> class Arg1)
/** Detect arg.function_name() */
#define GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_LEVEL_DETECTION(function_name,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION( std::declval<Arg1<level>>().function_name(),function_id, ,template<unsigned int> class Arg1 )

/** Detect arg1 op arg2. T-T / T-O / O-T denotes whether the argument is templated by level or
    an oridinary argument (that does not participate in maxlevel-detection) */
#define GAUSS_SIEVE_LAZY_BINARY_OP_TT_LEVEL_DETECTION(op,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION( std::declval<TArg1<level>>() op std::declval<TArg2<level>>(), function_id, ,template<unsigned int>class TArg1, template<unsigned int>class TArg2  )
#define GAUSS_SIEVE_LAZY_BINARY_OP_TO_LEVEL_DETECTION(op,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION( std::declval<TArg1<level>>() op std::declval<Arg2>(), function_id, ,template<unsigned int>class TArg1, class Arg2)
#define GAUSS_SIEVE_LAZY_BINARY_OP_OT_LEVEL_DETECTION(op,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION( std::declval<Arg1>() op std::declval<TArg2<level>>(), function_id, ,class Arg1, template<unsigned int>class TArg2)


/** Forward declarations */
#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_LEVEL_DETECTION_FORWARD_DECLARE(function_name,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION_FORWARD_DECLARE(function_id, template<unsigned int> class Arg1)
#define GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_LEVEL_DETECTION_FORWARD_DECLARE(function_name,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION_FORWARD_DECLARE(function_id, template<unsigned int> class Arg1)

#define GAUSS_SIEVE_LAZY_BINARY_OP_TT_LEVEL_DETECTION_FORWARD_DECLARE(op, function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION_FORWARD_DECLARE(function_id, template<unsigned int> class TArg1, template<unsigned int> class TArg2)
#define GAUSS_SIEVE_LAZY_BINARY_OP_TO_LEVEL_DETECTION_FORWARD_DECLARE(op, function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION_FORWARD_DECLARE(function_id, template<unsigned int> class TArg1, class Arg2)
#define GAUSS_SIEVE_LAZY_BINARY_OP_OT_LEVEL_DETECTION_FORWARD_DELCARE(op, function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_LEVEL_DETECTION_FORWARD_DECLARE(function_id, class Arg1, template<unsigned int> class TArg2)


// Does not work, because the above (iterated) macro expansion is delayed until final use...
// #undef GAUSS_SIEVE_LAZY_UNARY_THINGY_LEVEL_DETECTION

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

/* Helper macro to unify the free function and the member function cases */
#define GAUSS_SIEVE_LAZY_GENERAL_UNARY_CREATE_LAZY_WRAPPER(function_name, function_id, namespace_injection, function_call, function_call_type)\
namespace LazyEval::function_id##_helper_namespace{                                                \
namespace_injection                                                                                \
template<unsigned int maxlevel, template<unsigned int> class ArgAtLevel>                           \
struct Lazy_##function_id                                                                          \
{                                                                                                  \
  constexpr static unsigned int nargs = 1;                                                         \
  using IsLazyFunction = std::true_type;                                                           \
  static std::string fun_name() {return #function_name;}                                           \
  static constexpr unsigned int ApproxLevel = maxlevel;                                            \
  template<unsigned int level> using EvalType = mystd::decay_t< decltype( function_call_type ) >;  \
  /* This is to enforce instantiation of EvalType<0>, checking validity. */                        \
  static_assert(std::is_same<EvalType<0>,mystd::decay_t<EvalType<0>>>::value,"");                  \
                                                                                                   \
  template<unsigned int level, class Arg>                                                          \
  [[gnu::always_inline]] inline static EvalType<level> call(Arg &&arg)                             \
  {                                                                                                \
    static_assert(level <= maxlevel,"Cannot call at this level");                                  \
    static_assert(std::is_same<ArgAtLevel<level>,mystd::decay_t<Arg>>::value,"Invalid argument");  \
    return function_call;                                                                          \
  }                                                                                                \
};                                                                                                 \
} /* end LazyEval::function_id_helper_namespace */

#define GAUSS_SIEVE_LAZY_GENERAL_BINARY_CREATE_LAZY_WRAPPER(function_name, function_id, function_call, function_call_type, extra_static_assertions, namespace_injection,...)\
namespace LazyEval::function_id##_helper_namespace{ \
namespace_injection \
template<unsigned int maxlevel, __VA_ARGS__> \
struct Lazy_##function_id \
{ \
  constexpr static unsigned int nargs = 2; \
  using IsLazyFunction = std::true_type; \
  static std::string fun_name() { return #function_name;} \
  static constexpr unsigned int ApproxLevel = maxlevel; \
  template<unsigned int level> using EvalType = mystd::decay_t< decltype( function_call_type ) >;  \
  /* Will never fire, but causes instantiation of EvalType<0>, which contains static_asserts. */   \
  static_assert(std::is_same<EvalType<0>,mystd::decay_t<EvalType<0>>>::value,""); \
  \
  template<unsigned int level, class Arg1, class Arg2> \
  [[gnu::always_inline]] inline static EvalType<level> call (Arg1 &&arg1, Arg2 &&arg2) \
  { \
    static_assert(level <= maxlevel, "Cannot call at this level"); \
    extra_static_assertions \
    return function_call; \
  } \
};\
} /* End of namespace LazyEval::function_id_helper_namespace*/


/**
  For forward declarations:
*/

#define GAUSS_SIEVE_LAZY_GENERAL_UNARY_CREATE_LAZY_WRAPPER_FORWARD_DECLARE(function_id)             \
namespace LazyEval::function_id##_helper_namespace{                                                \
template<unsigned int maxlevel, template<unsigned int> class ArgAtLevel> struct Lazy_##function_id;\
}


#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_CREATE_LAZY_WRAPPER(function_name,function_id, namespace_injection)\
  GAUSS_SIEVE_LAZY_GENERAL_UNARY_CREATE_LAZY_WRAPPER(function_name,function_id,namespace_injection,function_name(std::forward<Arg>(arg)),function_name(std::declval<ArgAtLevel<level>>()) )
#define GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_CREATE_LAZY_WRAPPER(function_name,function_id) \
  GAUSS_SIEVE_LAZY_GENERAL_UNARY_CREATE_LAZY_WRAPPER(function_name,function_id, ,std::forward<Arg>(arg).function_name(), std::declval<ArgAtLevel<level>>().function_name() )

#define GAUSS_SIEVE_LAZY_BINARY_OP_TT_CREATE_LAZY_WRAPPER(op,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_BINARY_CREATE_LAZY_WRAPPER(op,function_id,\
  std::forward<Arg1>(arg1) op std::forward<Arg2>(arg2),\
  std::declval<Arg1AtLevel<level>>() op std::declval<Arg2AtLevel<level>>(), \
  static_assert(std::is_same<Arg1AtLevel<level>,mystd::decay_t<Arg1>>::value,"Arg1 invalid");  \
  static_assert(std::is_same<Arg2AtLevel<level>,mystd::decay_t<Arg2>>::value,"Arg2 invalid"); ,\
  , \
  template<unsigned int> class Arg1AtLevel, template<unsigned int> class Arg2AtLevel )
#define GAUSS_SIEVE_LAZY_BINARY_OP_TO_CREATE_LAZY_WRAPPER(op,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_BINARY_CREATE_LAZY_WRAPPER(op,function_id,\
  std::forward<Arg1>(arg1) op std::forward<Arg2>(arg2),\
  std::declval<Arg1AtLevel<level>>() op std::declval<Arg2NoLevel>(),\
  static_assert(std::is_same<Arg1AtLevel<level>,mystd::decay_t<Arg1>>::value,"Arg1 invalid");\
  static_assert(std::is_same<Arg2NoLevel,       mystd::decay_t<Arg2>>::value,"Arg2 invalid"); ,\
  , \
  template<unsigned int> class Arg1AtLevel, class Arg2NoLevel)
#define GAUSS_SIEVE_LAZY_BINARY_OP_OT_CREATE_LAZY_WRAPPER(op,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_BINARY_CREATE_LAZY_WRAPPER(op,function_id,\
  std::forward<Arg1>(arg1) op std::forward<Arg2>(arg2),\
  std::declval<Arg1NoLevel>() op std::declval<Arg2AtLevel<level>>(),\
  static_assert(std::is_same<Arg1NoLevel,       mystd::decay_t<Arg1>>::value,"Arg1 invalid");\
  static_assert(std::is_same<Arg2AtLevel<level>,mystd::decay_t<Arg2>>::value,"Arg2 invalid"); ,\
  , \
  class Arg1AtLevel, template<unsigned int> class Arg2NoLevel)




#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_CREATE_LAZY_WRAPPER_FORWARD_DECLARE(function_name,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_UNARY_CREATE_LAZY_WRAPPER_FORWARD_DECLARE(function_id)
#define GAUSS_SIEVE_LAZY_UNARY_MEMBER_FUNCTION_CREATE_LAZY_WRAPPER_FORWARD_DECLARE(function_name,function_id)\
  GAUSS_SIEVE_LAZY_GENERAL_UNARY_CREATE_LAZY_WRAPPER_FORWARD_DELCARE(function_id)


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
template<class Arg, TEMPL_RESTRICT_DECL2(Has_IsLazyNode<mystd::decay_t<Arg>>,__VA_ARGS__)>         \
auto inline function_name(Arg &&arg)                                                               \
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
  Forward declaration:
*/

#define GAUSS_SIEVE_LAZY_UNARY_FUNCTION_DIRECT_LAZY_FORWARD_DECLARE(function_name,function_id,...) \
template<class Arg, TEMPL_RESTRICT_DECL2(Has_IsLazyNode<mystd::decay_t<Arg>>,__VA_ARGS__)>         \
auto inline function_name(Arg &&arg)                                                               \
-> LazyEval::SieveLazyEval<                                                                        \
    LazyEval::function_id##_helper_namespace::Lazy_##function_id< /* LazyFun */                    \
      LazyEval::function_id##_helper_namespace::CanCall_##function_id<                             \
          mystd::decay_t<Arg>::template ObjectAtLevel>::template has_##function_id##_maxlevel /* newarglevel */ \
      <                                                                                            \
        ApproxLevelOf<mystd::decay_t<Arg>>::value /* arglevel */                                   \
      >(),                                                                                         \
      mystd::decay_t<Arg>::template ObjectAtLevel>,                                                \
    mystd::decay_t<Arg> >;



/***************************************************************************************************
Forward unary member function to unary free function. This is just to simplify things.
***************************************************************************************************/

#define GAUSS_SIEVE_FORWARD_UNARY_MEMBER_FUNCTION_TO_DELAYED(member_function_name,free_function_name,class_name) \
/* The template parameter RealClass=class_name serves the purpose that the function is only        \
   instantiated if it is actually used. */                                                         \
template<class RealClass = class_name>                                                             \
[[gnu::always_inline]]  inline auto member_function_name()  const &                                \
-> decltype(free_function_name(std::declval<RealClass>()) )                                        \
{                                                                                                  \
  static_assert(std::is_same<RealClass,class_name>::value,"");                                     \
  return free_function_name(static_cast<RealClass>(*this));                                        \
}                                                                                                  \
                                                                                                   \
template<class RealClass = class_name>                                                             \
[[gnu::always_inline]] inline auto member_function_name()  &&                                      \
-> decltype(free_function_name(std::move(std::declval<RealClass>()   )))                           \
{                                                                                                  \
  static_assert(std::is_same<RealClass,class_name>::value,"");                                     \
  return free_function_name(std::move(static_cast<RealClass>(*this)));                             \
}


/**
 Binary functions:
*/

#define GAUSS_SIEVE_LAZY_BINARY_FUNCTION_DIRECT_LAZY(function_name, function_id,...)\
template<class Arg1, class Arg2, TEMPL_RESTRICT_DECL2(Has_IsLazyNode<mystd::decay_t<Arg1>>,Has_IsLazyNode<mystd::decay_t<Arg2>>,__VA_ARGS__)>         \
auto inline function_name(Arg1 &&arg1, Arg2 &&arg2)                                                \
-> LazyEval::SieveLazyEval<                                                                        \
    LazyEval::function_id##_helper_namespace::Lazy_##function_id< /* LazyFun */                    \
      LazyEval::function_id##_helper_namespace::CanCall_##function_id<                             \
          mystd::decay_t<Arg1>::template ObjectAtLevel,mystd::decay_t<Arg2>::template ObjectAtLevel>::template has_##function_id##_maxlevel /* newarglevel */ \
      <                                                                                            \
        mystd::constexpr_min(ApproxLevelOf<mystd::decay_t<Arg1>>::value,ApproxLevelOf<mystd::decay_t<Arg2>>::value) /* arglevel */                                  \
      >(),                                                                                         \
      mystd::decay_t<Arg1>::template ObjectAtLevel,                                                \
      mystd::decay_t<Arg2>::template ObjectAtLevel>,                                               \
    mystd::decay_t<Arg1>,                                                                          \
    mystd::decay_t<Arg2> >                                                                         \
{                                                                                                  \
  static constexpr unsigned int arglevel = mystd::constexpr_min(ApproxLevelOf<mystd::decay_t<Arg1>>::value,ApproxLevelOf<mystd::decay_t<Arg2>>::value);             \
  static constexpr unsigned int newarglevel =                                                      \
    LazyEval::function_id##_helper_namespace::CanCall_##function_id                                \
      <                                                                                            \
        mystd::decay_t<Arg1>::template ObjectAtLevel,mystd::decay_t<Arg2>::template ObjectAtLevel  \
      >::template has_##function_id##_maxlevel<arglevel>();                                        \
  using LazyFun = LazyEval::function_id##_helper_namespace::Lazy_##function_id                     \
                    <newarglevel,mystd::decay_t<Arg1>::template ObjectAtLevel,mystd::decay_t<Arg2>::template ObjectAtLevel>;\
  using RetType = LazyEval::SieveLazyEval<LazyFun, mystd::decay_t<Arg1>,mystd::decay_t<Arg2>>;     \
  return RetType{std::forward<Arg1>(arg1),std::forward<Arg2>(arg2)};                               \
}

#define GAUSS_SIEVE_LAZY_BINARY_FUNCTION_DIRECT_LAZY_VALUE(function_name,function_id,...)          \
template<class Arg1, class Arg2, TEMPL_RESTRICT_DECL2(Has_IsLazyNode<mystd::decay_t<Arg1>>,mystd::negation<Has_IsLazyNode<mystd::decay_t<Arg2>>>,__VA_ARGS__)>\
decltype(auto) inline function_name(Arg1 &&arg1, Arg2 &&arg2)                                      \
{                                                                                                  \
  static constexpr unsigned int arglevel = ApproxLevelOf<mystd::decay_t<Arg1>>::value;             \
  static constexpr unsigned int newarglevel =                                                      \
    LazyEval::function_id##_helper_namespace::CanCall_##function_id                                \
    <                                                                                              \
      mystd::decay_t<Arg1>::template ObjectAtLevel,                                                \
      mystd::decay_t<Arg2>                                                                         \
    >::template has_##function_id##_maxlevel<arglevel>();                                          \
  using LazyFun = LazyEval::function_id##_helper_namespace::Lazy_##function_id                     \
    <                                                                                              \
      newarglevel,                                                                                 \
      mystd::decay_t<Arg1>::template ObjectAtLevel,                                                \
      mystd::decay_t<Arg2>                                                                         \
    >;                                                                                             \
  using ValueCapture = LazyEval::LazyWrapValue                                                     \
    <                                                                                              \
      mystd::decay_t<Arg2>,                                                                        \
      newarglevel                                                                                  \
    >;                                                                                             \
  using RetType = LazyEval::SieveLazyEval                                                          \
    <                                                                                              \
      LazyFun,                                                                                     \
      mystd::decay_t<Arg1>,                                                                        \
      ValueCapture                                                                                 \
    >;                                                                                             \
  return RetType{std::forward<Arg1>(arg1),ValueCapture{std::forward<Arg2>(arg2)} };                \
}
