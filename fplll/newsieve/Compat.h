// Macro definitions for C++ version compatibility

#ifndef SIEVE_GAUSS_COMPAT_H
#define SIEVE_GAUSS_COMPAT_H

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
  #define DECLARECONJUNCTION(ClassName) template<class... B> using ClassName = std::conjunction<B...>;
#else
// define our own version of std::conjunction
  #define DECLARECONJUNCTION(ClassName) \
  template<class...> struct ClassName : std::true_type{}; \
  template<class B1> struct ClassName<B1> : B1 {}; \
  template<class B1,class... Bs> struct ClassName<B1,Bs...> \
    : std::conditional<static_cast<bool>(B1::value), ClassName<Bs...>,B1>::type {};
#endif

#endif
