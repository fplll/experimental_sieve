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

/**
  CPP14CONSTEXPR is used to declare functions as constexpr, provided the (considerable!) relaxations
  to the constexpr specifier from C++14 are supported.

  CPP17CONSTEXPRIF is used to support C++17 if constexpr. If unsupported, falls back to a normal if.
  Note that the main point of if constexpr(expr){ foo() } is that expr might depend on a template
  argument and if expr is false, foo does not need to even compile. This feature would indeed be
  extremely useful and simplify a lot of code.
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

#if __cpp_constexpr >= 201304
  #define CPP14CONSTEXPR constexpr
#else
  #define CPP14CONSTEXPR
#endif

#if __if_constexpr
  #define CPP17CONSTEXPRIF if constexpr
#else
  #define CPP17CONSTEXPRIF if
#endif

/**
  Replacements for C++14 standard library features that are missing in C++11.
  See the C++ - documentation for their meaning.

  TODO: The MyNAME variants rather than mystd::name variants is confusing.
        Simplify to only use mystd::
*/

namespace GaussSieve
{
#if __cpp_lib_logical_traits >= 201510
  template<class... Bs> using MyConjunction = std::conjunction<Bs...>;           //AND
  template<class... Bs> using MyDisjunction = std::disjunction<Bs...>;           //OR
  template<class B>     using MyNegation    = std::negation<B>;                  //NOT
  template<class... Bs> using MyNAND        = MyNegation<MyConjunction<Bs...>>;  //NAND
  template<class... Bs> using MyNOR         = MyNegation<MyDisjunction<Bs...>>;  //NOR
#else
// just implement std::conjunction and friends myself:
  template<class...> struct MyConjunction     : std::true_type{};
  template<class B1> struct MyConjunction<B1> : B1 {};
  template<class B1,class... Bs> struct MyConjunction<B1,Bs...>
    : std::conditional<static_cast<bool>(B1::value), MyConjunction<Bs...>,B1>::type {};

  template<class...> struct MyDisjunction     : std::false_type{};
  template<class B1> struct MyDisjunction<B1> : B1 {};
  template<class B1, class... Bs> struct MyDisjunction<B1,Bs...>
    : std::conditional<static_cast<bool>(B1::value), B1, MyDisjunction<Bs...> >::type {};

  template<class B> struct MyNegation : std::integral_constant<bool,!static_cast<bool>(B::value)>{};
  template<class... Bs> using MyNAND = MyNegation<MyConjunction<Bs...>>;
  template<class... Bs> using MyNOR  = MyNegation<MyDisjunction<Bs...>>;
#endif

namespace mystd
{
  template<bool B, class T, class F> using conditional_t = typename std::conditional<B,T,F>::type;
  template<class T>                  using decay_t       = typename std::decay<T>::type;

  template<bool b>      using bool_constant = std::integral_constant<bool,b>;
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
}


#if __cpp_lib_integer_sequence >= 201304
  template<std::size_t... Ints> using MyIndexSeq         = std::index_sequence<Ints...>;
  template<std::size_t N>       using MyMakeIndexSeq     = std::make_index_sequence<N>;
  template<class... T>          using MyIndexSequenceFor = std::index_sequence_for<T...>;
#else
  template<std::size_t... Ints> class MyIndexSeq {}; // not equivalent to the above, lacks size()
  namespace Helpers
  {
  // encapsulates a integer sequence 0,...,N-1, RestArgs
    template<std::size_t N, std::size_t... RestArgs> struct GenIndexSeq
        : GenIndexSeq <N-1,N-1, RestArgs...> {};
    template<std::size_t... RestArgs> struct GenIndexSeq<0,RestArgs...>
    {
      using type = MyIndexSeq<RestArgs...>;
    };
  }
  template<std::size_t N> using MyMakeIndexSeq = typename Helpers::GenIndexSeq<N>::type;
  template<class... T> using MyIndexSequenceFor = MyMakeIndexSeq<sizeof...(T)>;
#endif

namespace mystd
{
//  template<class... Bs> using conjunction   = MyConjunction<Bs...>;    //AND
//  template<class... Bs> using disjunction   = MyDisjunction<Bs...>;    //OR
//  template<class B>     using negation      = MyNegation<B>;           //NOT

  template<bool B, class T, class F> using conditional_t = typename std::conditional<B,T,F>::type;

  // std::max is not constexpr until C++14. This version is always constexpr, but does not support
  // custom comparators or initializer lists.
  template<class T> constexpr const T& constexpr_max(const T& x1, const T& x2)
  { return (x1<x2) ? x2 : x1; }
  template<class T> constexpr const T& constexpr_min(const T& x1, const T& x2)
  { return (x1<x2) ? x1 : x2; }
}  // end namespace mystd
}  // end namespace GaussSieve

#endif

// clang-format on
