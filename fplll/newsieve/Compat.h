// Definitions for C++ version compatibility

// clang-format adjustments finished -- Gotti

/**
  This file defines some compatibility macros / workaround for C++-14 / C++-17 features.
  Notably, it defines the namespace GaussSieve::mystd::, which contains alternative
  implementations for some std:: functions/typedefs/etc.
  If we detect that the C++ implemtation has these, the corresponding mystd:: variants are just
  aliased to the std:: ones.
*/

#ifndef SIEVE_GAUSS_COMPAT_H
#define SIEVE_GAUSS_COMPAT_H

#include <type_traits>
#include <utility>

#if __cpp_lib_experimental_detect >= 201505
#include <experimental/type_traits>
#endif

/**
  CPP14CONSTEXPR is used to declare functions as constexpr, provided the (considerable!) relaxations
  to the constexpr specifier from C++14 are supported.

  CPP17CONSTEXPRIF is used to support C++17 if constexpr. If unsupported, falls back to a normal if.
  Note that the main point of if constexpr(expr){ foo() } is that expr might depend on a template
  argument and if expr is false, foo does not need to even compile (with some caveats).
  This feature would indeed be extremely useful and simplify a lot of code.
  (The workaround is to define a class template X, templated by a bool b. Specialize X for b==true
  with static member function f that implements foo and specialize X for b==false with a member
  function that does nothing. Then we can call X<epxr>::f(arg), passing required arguments arg.
  This is the state-of-the-art in C++11, which is extremely annoying and clutters the code.)
  Unfortunately, to support the fallback to normal if, we cannot use CPP17CONSTEXPRIF that way.
  Instead, CPP17CONSTEXPRIF is supposed to be used for if's that should depend on compile-time
  expressions in order to trigger compiler errors if this is not the case and to document intent.
  Efficiency-wise, there is no benefit if the compile is able to elimate if(false){ ... } code
  (which seems a safe assumption).
*/

// Note that undefined symbols in an #if get replaced by 0, so this works even if
// __cpp_constexpr is not even defined.

// unfortunately, clang-format does not like indendation of (nested) #ifdefs,
// which means it's deactivated 90% of this file.
// clang-format off

#if __cpp_constexpr >= 201304
  #define CPP14CONSTEXPR constexpr
#else
  #define CPP14CONSTEXPR
#endif

#if !defined(__cpp_constexpr)
  #warning Your compiler does not support feature testing for constexpr.
#endif

#if __cpp_if_constexpr >= 201606
  #define CPP17CONSTEXPRIF if constexpr
#else
  #define CPP17CONSTEXPRIF if
#endif

// attributes

#if defined(__has_cpp_attribute)
  #if __has_cpp_attribute(nodiscard)
    // indicates that return value must not be discarded.
    // standardized attribute starting from C++17
    #define NODISCARD [[nodiscard]]
  #else
    #define NODISCARD
  #endif
  #if __has_cpp_attribute(gnu::always_inline)
    // Enforces that the function gets inlined.
    // Note that we only use it if any reasonable compiler *ought* to inline anyway
    // (1-line perfect-forward dispatch functions mimicking constexpr if, for example)
    // The purpose is that we want the compiler to yell at us if it cannot.
    // relevant for some expression-template constructions, where failure to inline can become very
    // expensive.
    #define FORCE_INLINE [[gnu::always_inline]]
  #else
    #define FORCE_INLINE
  #endif
#else
  // Note: per C++ standard, unrecognized attributes are ignored (it may generate warnings, though).
  // This is exactly the behaviour we want.
  // So in case the compiler does not support __has_cpp_attribute, we might also just go ahead and
  // define those as above.
  // Case in point is some clang-versions, which recognize gnu::foo attributes, but do not have
  // __has_cpp_attribute
  #define NODISCARD    [[nodiscard]]
  #define FORCE_INLINE [[gnu::always_inline]]
  #warning "Compiler does support feature testing for attributes."
#endif

// clang-format on

/**
  Poor man's "requires" (C++20 / concepts-lite TS):
  We want a lot of templates to only be instantiatable if the template parameters satify certain
  boolean constraints (e.g. "Is a lattice Point", "Is an integer", "Is an Iterator"...)
  Lacking the "requires" keyword, we use std::enable_if and SFINAE.
  To improve readability and making the intent clearer, we define some macros
  TEMPL_RESTRICT_DECL, TEMPL_RESTRICT_IMPL
  Usage:
  template<class Integer, TEMPL_RESTRICT_DECL((std::is_integral<Integer>::value))> some_declaration;
  template<class Integer, TEMPL_RESTRCIT_IMPL((std::is_integral<Integer>::value))> some_definition;

  This restricts the template to only be instantiatable if the argument to TEMPL_RESTRICT_DECL/IMPL
  is convertible to boolean true at compile time. Use as the last template argument.
  (For variadic templates, you should know what you are doing)

  The variants TEMPL_RESTRICT_DECL2 and TEMPL_RESTRICT_IMPL2 are similar, but take a comma-separated
  list of *types* T_1,T_2,... encapsulating a boolean value instead.
  (i.e. types with static constexpr bool T_i::value defined -- This is compatible with all such
  types and manipulators from std:: ). It restricts to the AND of these types.

  The difference between the DECL and IMPL variants is that you need to use DECL in declarations or
  declarations + definitions, whereas you need to use IMPL in definitions of previously declared
  template member functions (i.e. when declaring a member function inside a class, use DECL. If
  you use a separate file to then define the function, you need to use IMPL -- This is because DECL
  sets default parameters and IMPL does not)

  NOTE: Since these are macros, be aware of unshielded commas in template parameters. Use ((arg))
  for the TEMPL_RESTRICT_DECL/IMPL((arg)) variants. The DECL2/IMPL2 variant automagically work.
  (The arguments are actually chopped up by the preprocessor, but it still works)
*/

/**
  NOTE: If you get an compiler error "no type named "type" in std::enable_if<false,int>,
        then you are using TEMPL_RESTRICT_* wrongly. Due to some choices in the C++ standard
        the condition argument(s) to TEMPL_RESTRIC_* must evaluate to true for at least one possible
        choice of template arguments.

        In the important case of member functions templates of class templates, this means that for
        each instantiations of the class (i.e. class template parameters fixed), there has to be a
        set of template arguments for the function, such that this holds.

        A workaround to use template arguments of the class is as follows:
        Instead of using

        template<class ClassArg>
        class Myclass
        {
          template<TEMPL_RESTRICT_DECL2(Condition<ClassArg>)> foo()
          {
            ...
          }
        };

        you need to use

        template<class ClassArg>
        class Myclass
        {
          template<class dummy = ClassArg, TEMPL_RESTRICT_DECL2(Condition<dummy>)> foo ()
          {
            // you might want to static_assert(std::is_same<dummy,ClassArg>::value, "");
            // because the caller might set the dummy template parameter explitly by calling
            // foo<Bar>()
            // (If foo has other template paramters as well, there could be a mix-up)
            ...
          }
        };
*/

#define TEMPL_RESTRICT_DECL(condition) typename std::enable_if<condition, int>::type = 0
#define TEMPL_RESTRICT_IMPL(condition) typename std::enable_if<condition, int>::type

#define TEMPL_RESTRICT_DECL2(...)                                                                  \
  typename std::enable_if<GaussSieve::mystd::conjunction<__VA_ARGS__>::value, int>::type = 0
#define TEMPL_RESTRICT_IMPL2(...)                                                                  \
  typename std::enable_if<GaussSieve::mystd::conjunction<__VA_ARGS__>::value, int>::type

/**
  Replacements for C++14 (and beyond) standard library features that are missing in C++11.
  See the C++ - documentation for their meaning.
*/

namespace GaussSieve
{
namespace mystd
{
// some often-used shorthands to avoid having to use typename ...
// completely identical to the corresponding std::conditonal_t etc.
// clang-format off
template<bool b, class T, class F> using conditional_t   = typename std::conditional<b, T, F>::type;
template<class T>                  using decay_t         = typename std::decay<T>::type;
template<bool b>                   using bool_constant   = std::integral_constant<bool,b>;
template<bool b, class T = void>   using enable_if_t     = typename std::enable_if<b, T>::type;
template<class T>                  using make_unsigend_t = typename std::make_unsigned<T>::type;
// clang-format o

// std::max is not (and can't be) constexpr until C++14 (where the def. of constexpr was relaxed).
// This version is always constexpr, but does not support custom comparators or initializer lists.
// Hence, this differs from anything in std::
template <class T> constexpr const T &constexpr_max(const T &x1, const T &x2)
{
  return (x1 < x2) ? x2 : x1;
}
template <class T> constexpr const T &constexpr_min(const T &x1, const T &x2)
{
  return (x1 < x2) ? x1 : x2;
}

// clang-format off
#if __cpp_lib_logical_traits >= 201510
  template <class... Bs> using conjunction = std::conjunction<Bs...>;  // AND
  template <class... Bs> using disjunction = std::disjunction<Bs...>;  // OR
  template <class B>     using negation    = std::negation<B>;         // NOT
#else
  // just implement std::conjunction, disjunction and negation myself:
  template <class...> struct conjunction     : std::true_type {};
  template <class B1> struct conjunction<B1> : B1 {};
  template <class B1,class... Bs> struct conjunction<B1,Bs...>
      : conditional_t<static_cast<bool>(B1::value), conjunction<Bs...>, B1 > {};

  template <class...> struct disjunction     : std::false_type {};
  template <class B1> struct disjunction<B1> : B1 {};
  template <class B1, class... Bs> struct disjunction<B1,Bs...>
      : conditional_t<static_cast<bool>(B1::value), B1, disjunction<Bs...>> {};

  template <class B> struct negation : bool_constant<!static_cast<bool>(B::value)> {};
#endif

// clang-format on

// clang-format off

#if __cpp_lib_integer_sequence >= 201304
  // index_sequence<size_t... Ints> and friends are aliases for integer_sequence<size_t, Ints...>
  // (i.e. integer_sequence with type std::size_t). Look up std::integer_sequence if you cannot
  // find std::index_sequence
  template <std::size_t... Ints> using index_sequence      = std::index_sequence<Ints...>;
  template <std::size_t N>       using make_index_sequence = std::make_index_sequence<N>;
  template <class... T>          using index_sequence_for  = std::index_sequence_for<T...>;
#else
  template <std::size_t... Ints> class index_sequence
  {
  public:
    static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
    using value_type = std::size_t;
  };
  namespace IndexSeqHelper
  {
  // encapsulates a integer sequence 0, ..., N-1, RestArgs
  template <std::size_t N, std::size_t... RestArgs>
  struct GenIndexSeq : GenIndexSeq<N - 1, N - 1, RestArgs...>
  {
  };
  template <std::size_t... RestArgs>
  struct GenIndexSeq<0, RestArgs...>
  {
    using type = index_sequence<RestArgs...>;
  };
  }  // end namespace IndexSeqHelper
  template <std::size_t N> using make_index_sequence = typename IndexSeqHelper::GenIndexSeq<N>::type;
  template <class... T>    using index_sequence_for  = make_index_sequence<sizeof...(T)>;
#endif

// clang-format on

// void_t takes any number of arguments and ignores them.
// The point is that arguments have to be valid (to be used in SFINAE contexts).
// This is seemingly trivial, but actually extremely useful.

// clang-format off
#if __cpp_lib_void_t >= 201411
  template <class... Args> using void_t = std::void_t<Args...>;
#else
  // proper definition:
  /*
  template <class... Args> using void_t = void;
  */

  // workaround for compiler bug in GCC4.9
  namespace Void_tHelper
  {
    template <class... > struct make_void { using type = void; };
  }
  template <class... T> using void_t = typename Void_tHelper::make_void<T...>::type;
#endif
// clang-format on

// std::experimental::*_detected_* - features.
// Part of Library fundamentals V2 TS.
// These can be used to query the existence of e.g. typedefs (and more).
// Used for traits checking (More elegant than the old macros)
// Taken from http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/n4502.pdf

// clang-format off
#if __cpp_lib_experimental_detect >= 201505

  using nonesuch = std::experimental::nonesuch;

  template <template<class...> class Op, class... Args>
  using is_detected = std::experimental::is_detected<Op, Args...>;

  template<template<class...> class Op, class... Args>
  using detected_t  = std::experimental::detected_t<Op, Args...>;

  template<class Default, template<class...> class Op, class... Args>
  using detected_or = std::experimental::detected_or<Default, Op, Args...>;

  template<class Default, template<class...> class Op, class... Args>
  using detected_or_t = std::experimental::detected_or_t<Default,Op,Args...>;

#else

  struct nonesuch // indicating "Not detected"
  {
    nonesuch()                       = delete;
    ~nonesuch()                      = delete;
    nonesuch(nonesuch const &)       = delete;
    void operator=(nonesuch const &) = delete;
  };

  namespace detection_impl
  {
  template <class Default, class AlwaysVoid, template<class...> class Op, class... Args>
  struct detector
  {
    using value_t = std::false_type;
    using type    = Default;
  };

  template<class Default, template<class...> class Op, class... Args>
  struct detector<Default, void_t<Op<Args...>>, Op, Args...>
  {
    using value_t = std::true_type;
    using type    = Op<Args...>;
  };
  }  // end namespace detection_impl

  template<template<class...> class Op, class... Args>
  using is_detected   = typename detection_impl::detector<void, void, Op, Args...>::value_t;

  template<template<class...> class Op, class... Args>
  using detected_t    = typename detection_impl::detector<mystd::nonesuch, void, Op, Args...>::type;

  template<class Default, template<class...> class Op, class... Args>
  using detected_or   = typename detection_impl::detector<Default, void, Op, Args...>;

  template<class Default, template<class...> class Op, class... Args>
  using detected_or_t = typename detected_or<Default, Op, Args...>::type;

#endif  // __cpp_lib_experimental_detect >= 201505

// clang-format on

}  // end namespace mystd
}  // end namespace GaussSieve

#endif  // include-guards
