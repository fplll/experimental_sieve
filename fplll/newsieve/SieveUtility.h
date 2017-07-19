// This file contains utility (boilerplate) functions and classes only. It should not have
// dependencies on other files within the Gauss Sieve and be header-only.

#ifndef GAUSS_SIEVE_UTILITY_H
#define GAUSS_SIEVE_UTILITY_H

#include <iostream>
#include <istream>
#include <string>
#include <type_traits>

/**
This macro is used to test the presence of a (public) member typedef in a class
Args:   TypeToCheck - typename whose presence to check
        CheckerClassName - Name of the checker class
This macro emits a new template class definition with the name CheckerClassName.
TypeToCheck must not be void.

Usage:
CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TypeToCheck, CheckerClassName);
This creates(!) the template class CheckerClassName.

Then CheckerClassName<SomeSuspiciousClass>::value will be true if
SomeSuspicousClass::TypeToCheck exists, false otherwise

Note that the missing semicolon at the end is intentional.
The user needs to put it to emphasize that this is a declaration.
*/

// clang-format off

#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TypeToCheck, CheckerClassName)                           \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> void                             foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t =                                                                                \
        std::integral_constant<bool, !(std::is_void<decltype(foo<ClassToCheck>(0))>::value)>;      \
    static bool constexpr value = value_t::value;                                                  \
  }

/**
Similar to the above, creates a checker template class that checks wether
TypeToCheck exists and is equal to TypeShouldBe
*/

#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(TypeToCheck, TypeShouldBe, CheckerClassName)      \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> void                             foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t =                                                                                \
        std::integral_constant<bool,                                                               \
                               std::is_same<TypeShouldBe, decltype(foo<ClassToCheck>(0))>::value>; \
    static bool constexpr value = value_t::value;                                                  \
  }

// clang-format on

// class that ignores its argument. Can be used to optimize away unused parameters in function
// templates...
class IgnoreAnyArg
{
public:
  template <class T> constexpr IgnoreAnyArg(T val){};
  constexpr IgnoreAnyArg() = default;
};

// same, but enforces the type of the ignored argument.
template <class T> class IgnoreArg
{
public:
  inline constexpr IgnoreArg(T val){};
  constexpr IgnoreArg() = default;
};

template <int nfixed = -1> class Dimension;

template <> class Dimension<-1>
{
public:
  using IsFixed = std::false_type;

  Dimension(unsigned int const new_dim) : dim(new_dim){};
  Dimension() = default;  // Not sure whether we should allow uninitialized dims here. The issue is that we want the same interface in both cases.
  inline operator unsigned int() const { return dim; };
  unsigned int dim;
};

template <int nfixed> class Dimension
{
public:
  using IsFixed = std::true_type;

  Dimension() = default;
  Dimension(IgnoreArg<unsigned int> new_dim){};  // assert(new_dim==nfixed);}
  //    Dimension(unsigned int){};
  inline constexpr operator unsigned int() const { return nfixed; };
  static constexpr unsigned int dim = nfixed;
};

namespace GaussSieve  // helper functions
{
/**
string_consume(is, str, elim_ws, verbose) reads from stream is.
If the next read on is is not the string str, it returns false,
otherwise it throws away str from is.
elim_ws==true means that string_consume removes whitespace from the stream before and after its
operations.
If verbose is set to true, in the case that the function returns false, we additionally write
diagnostics to cerr.

More precisely, we first remove whitespace (optional), then read and remove len(str) chars from is.
If there is not enough data on the stream, we read and remove all there is.
If we read what was expected, return true, otherwise false.
We also remove whitespace from is afterwards (optional)

This utility function is used to parse dumps.

string_consume assumes that str itself does not start/end with whitespace.
*/

inline bool string_consume(std::istream &is, std::string const &str, bool elim_ws = true,
                           bool verbose = true);  // helper function for dumping/reading
}

inline bool GaussSieve::string_consume(std::istream &is, std::string const &str, bool elim_ws,
                                       bool verbose)
{
  unsigned int len = str.length();
  char *buf        = new char[len + 1];
  buf[len]         = 0;  // for error message.
  if (elim_ws)
  {
    is >> std::ws;
  }
  is.read(buf, len);
  if (is.gcount() != len)
  {
    if (verbose)
    {
      std::cerr << "Failure reading header: Expected to read" << str << std::endl;
      std::cerr << "Read only " << is.gcount() << "bytes. String read was" << buf << std::endl;
    }
    return false;
  }
  if (elim_ws)
  {
    is >> std::ws;
  }
  if (str.compare(0, len, buf, len) != 0)
  {
    if (verbose)
    {
      std::cerr << "Failure reading header: Expected to read" << str << std::endl;
      std::cerr << "Read instead:" << buf << std::endl;
    }
    return false;
  }
  return true;
}

#endif  // GAUSS_SIEVE_UTILITY_H
