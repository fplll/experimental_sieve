// Definitions for C++ version compatibility

// clang-format off
// clang-format makes readability worse: It does not like indentation of #if's or {} for empty classes.

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
  expressions in order to trigger compiler errors if this is not the case. Efficiency-wise,
  there is no difference if the compile is able to elimate if(false){ ... } code (which seems a safe
  assumption).
*/

#if defined(__cpp_constexpr)
  #if __cpp_constexpr >= 201304
    #define CPP14CONSTEXPR constexpr
  #else
    #define CPP14CONSTEXPR
  #endif
#else
  #define CPP14CONSTEXPR
  #warning "Your compiler does not support feature testing for constexpr."
#endif

#if defined(__if_constexpr)
  #if __if_constexpr
    #define CPP17CONSTEXPRIF if constexpr
  #else
    #define CPP17CONSTEXPRIF if
  #endif
#else
  //#warning "Your compiler does not support feature testing for if constexpr."
  #define CPP17CONSTEXPRIF if
#endif

#if defined(__has_cpp_attribute)
  #if __has_cpp_attribute(nodiscard)
    #define NODISCARD [[nodiscard]]
  #else
    #define NODISCARD
  #endif
  #if __has_cpp_attribute(gnu::always_inline)
    #define FORCE_INLINE [[gnu::always_inline]]
  #else
    #define FORCE_INLINE
  #endif
#else
  // Note: per C++ standard, unrecognized attributes are ignored (it may generate warnings, though)
  // So in case the compiler does not support __has_cpp_attribute, we might also just go ahead and define those as above.
  // Case in point is some clang-versions, which recognize gnu::foo attributes, but do not have __has_cpp_attribute
  #define NODISCARD
  #define FORCE_INLINE
  #warning "Your compiler does support have feature testing for attributes."
#endif

/**
  Replacements for C++14 (and beyond) standard library features that are missing in C++11.
  See the C++ - documentation for their meaning.
*/

namespace GaussSieve
{
namespace mystd
{
  // some often-used shorthands to avoid having to use typename ...
  template<bool B, class T, class F> using conditional_t = typename std::conditional<B,T,F>::type;
  template<class T>                  using decay_t       = typename std::decay<T>::type;

  template<bool b>      using bool_constant = std::integral_constant<bool,b>;

  // std::max is not (and cannote be) constexpr until C++14 (where the def. of constexpr was relaxed).
  // This version is *always* constexpr, but does not support custom comparators or initializer lists.
  template<class T> constexpr const T& constexpr_max(const T &x1, const T &x2)
  { return (x1 < x2) ? x2 : x1; }
  template<class T> constexpr const T& constexpr_min(const T &x1, const T &x2)
  { return (x1 < x2) ? x1 : x2; }

#if __cpp_lib_logical_traits >= 201510
  template<class... Bs> using conjunction   = std::conjunction<Bs...>;           //AND
  template<class... Bs> using disjunction   = std::disjunction<Bs...>;           //OR
  template<class B>     using negation      = std::negation<B>;                  //NOT
#else
  // just implement std::conjunction, disjunction and negation myself:
  template<class...> struct conjunction     : std::true_type{};
  template<class B1> struct conjunction<B1> : B1 {};
  template<class B1,class... Bs> struct conjunction<B1,Bs...>
    : conditional_t< static_cast<bool>(B1::value), conjunction<Bs...>, B1 > {};

  template<class...> struct disjunction     : std::false_type{};
  template<class B1> struct disjunction<B1> : B1 {};
  template<class B1, class... Bs> struct disjunction<B1,Bs...>
    : conditional_t< static_cast<bool>(B1::value), B1, disjunction<Bs...> > {};

  template<class B> struct negation : bool_constant<!static_cast<bool>(B::value)>{};
#endif

#if __cpp_lib_integer_sequence >= 201304
// index_sequence<size_t... Ints> and friends are aliases for integer_sequence<size_t, Ints...>
// (i.e. integer_sequence with type std::size_t). Look up std::integer_sequence if you cannot
// find std::index_sequence
  template<std::size_t... Ints> using index_sequence      = std::index_sequence<Ints...>;
  template<std::size_t N>       using make_index_sequence = std::make_index_sequence<N>;
  template<class... T>          using index_sequence_for  = std::index_sequence_for<T...>;
#else
  template<std::size_t... Ints> class index_sequence
  {
  public:
    static constexpr std::size_t size() noexcept { return sizeof...(Ints); }
    using value_type = std::size_t;
  };
  namespace IndexSeqHelper
  {
  // encapsulates a integer sequence 0,...,N-1, RestArgs
    template<std::size_t N, std::size_t... RestArgs> struct GenIndexSeq
        : GenIndexSeq <N-1,N-1,RestArgs...> {};
    template<std::size_t... RestArgs> struct GenIndexSeq<0,RestArgs...>
    {
      using type = index_sequence<RestArgs...>;
    };
  }
  template<std::size_t N> using make_index_sequence = typename IndexSeqHelper::GenIndexSeq<N>::type;
  template<class... T>    using index_sequence_for  = make_index_sequence<sizeof...(T)>;
#endif

#if __cpp_lib_void_t >= 201411
// takes any number of arguments and ignores them. The point is that arguments have to be valid (to be used in SFINAE contexts)
// seemingly trivial, but actually extremely useful.
  template<class... Args> using void_t = std::void_t<Args...>;
#else
  template<class... Args> using void_t = void;
#endif

// std::experimental::*_detected_* - features.
// Part of Library fundamentals V2 TS.
// These can be used to query the existence of e.g. typedefs (and more).
// Used for traits checking (More elegant than the old macros)
// Taken from http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/n4502.pdf

#if __cpp_lib_experimental_detect >= 201505

using nonesuch = std::experimental::nonesuch;
template<template<class...> class Op, class... Args> using is_detected = std::experimental::is_detected<Op,Args...>;
template<template<class...> class Op, class... Args> using detected_t  = std::experimental::detected_t <Op,Args...>;
template<class Default, template<class...> class Op, class... Args> using detected_or   = std::experimental::detected_or<Default,Op,Args...>;
template<class Default, template<class...> class Op, class... Args> using detected_or_t = std::experimental::detected_or_t<Default,Op,Args...>;

#else

struct nonesuch // indicating "Not detected"
{
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const &) = delete;
  void operator=(nonesuch const &) = delete;
};

namespace detection_impl
{
template<class Default, class AlwaysVoid, template<class...> class Op, class... Args>
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
}

template<template<class...> class Op, class... Args>
using is_detected = typename detection_impl::detector<void, void, Op, Args...>::value_t;

template<template<class...> class Op, class... Args>
using detected_t  = typename detection_impl::detector<mystd::nonesuch, void, Op, Args...>::type;

template<class Default, template<class...> class Op, class... Args>
using detected_or = typename detection_impl::detector<Default, void, Op, Args...>;

template<class Default, template<class...> class Op, class... Args>
using detected_or_t = typename detected_or<Default, Op, Args...>::type;

#endif


}  // end namespace mystd
}  // end namespace GaussSieve

#endif

// clang-format on
