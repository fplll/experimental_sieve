// Definitions for C++ version compatibility

#ifndef SIEVE_GAUSS_COMPAT_H
#define SIEVE_GAUSS_COMPAT_H

#include<type_traits>

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

#if __cpp_lib_logical_traits >= 201510
// untested:
namespace GaussSieve
{
  template<class... Bs> using MyConjunction = std::conjunction<Bs...>;    //AND
  template<class... Bs> using MyDisjunction = std::disjunction<Bs...>;    //OR
  template<class B>     using MyNegation    = std::negation<B>;           //NOT
  template<class... Bs> using MyNAND = MyNegation<MyConjunction<Bs...>>;  //NAND
  template<class... Bs> using MyNOR  = MyNegation<MyDisjunction<Bs...>>;  //NOR
}
#else

// just implement std::conjunction and friends myself:
namespace GaussSieve
{
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
}
#endif

#endif
