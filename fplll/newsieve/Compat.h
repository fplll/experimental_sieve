// Macro definitions for C++ version compatibility

#ifndef SIEVE_GAUSS_COMPAT_H
#define SIEVE_GAUSS_COMPAT_H

#if __cpp_constexpr >= 201304
  #define CXX14CONSTEXPR constexpr
#else
  #define CXX14CONSTEXPR
#endif

#if __if_constexpr
  #define CXX17CONSTEXPR constexpr
#else
  #define CXX17CONSTEXPR
#endif

#endif
