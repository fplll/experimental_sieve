// This file contains utility (boilerplate) functions and classes only. It should not have
// dependencies on other files within the Gauss Sieve and be header-only.

#ifndef GAUSS_SIEVE_UTILITY_H
#define GAUSS_SIEVE_UTILITY_H

#include "DefaultIncludes.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"
#include "gmp.h"
#include "gmpxx.h"
// #include "TraitChecks.h"  -- No longer needed

namespace GaussSieve
{

// forward declarations
template <class T> class IgnoreArg;
class IgnoreAnyArg;
template <class T> class IsZNRClass;
template <class Integer> struct ConvertMaybeMPZ;
template <class T> class UnZNR;

/*****************************************************************************
is_a_power_of_two(n) checks whether n is a power of 2.
The constexpr version may be slow; it is designed for use in static asserts
(This is why it's a recursive one-line function, to be C++11 - constexpr)
******************************************************************************/
template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
constexpr bool is_a_power_of_two_constexpr(Integer const n)
{
  // one-line function to be C++11 - constexpr. Slow, but only used in static_asserts anyway.
  return (n > 0) && ((n == 1) || ((n % 2 == 0) && is_a_power_of_two_constexpr(n / 2)));
}
template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>, std::is_unsigned<Integer>)>
CPP14CONSTEXPR bool is_a_power_of_two(Integer const n)
{
  // essentially performs a popcount and checks whether it's 1.
  // clang-format off
  return std::bitset< std::numeric_limits<Integer>::digits >{n}.count() == 1;
  // clang-format on
}

// retrieves first type of a sequence of types. Useful in variadic templates.
// may be unused.

// clang-format off
template <class Arg1, class...>
class GetFirstArg
{
  using type = Arg1;
};
template <class Arg1, class... Args> using GetFirstArg_t = GetFirstArg<Arg1, Args...>;
// clang-format on

// class that ignores its argument. Can be used to optimize away unused parameters in function
// Can probably go away...
class IgnoreAnyArg
{
public:
  template <class T> constexpr IgnoreAnyArg(T) noexcept {};
  constexpr IgnoreAnyArg() = default;
};

// same, but enforces the type of the ignored argument.
template <class T> class IgnoreArg
{
public:
  inline constexpr IgnoreArg(T) noexcept {};
  constexpr IgnoreArg() = default;
};

// MaybeFixed<-1> encapsulates a nonnegative runtime integer
// MaybeFixed<nfixed> encapsulates a compile-time known integer nfixed.
// This is used as a template parameter to propagate a fixed dimension through
// a lot of our classes. Some data structures make use of a known dimension.

template <int nfixed = -1, class UIntClass = unsigned int> class MaybeFixed;

template <class UIntClass>  // specialization for nfixed == -1
class MaybeFixed<-1, UIntClass>
{
  static_assert(std::is_unsigned<UIntClass>::value, "MaybeFixed only works with unsigned types.");

public:
  static constexpr bool IsFixed = false;
  using IsFixed_t               = std::false_type;
  using type                    = UIntClass;
  constexpr MaybeFixed(UIntClass const new_val) noexcept : value(new_val) {}
  MaybeFixed() = default;  // Not sure whether we should allow uninitialized dims here. The issue is
                           // that we want the same interface in both cases.
  inline operator UIntClass() const { return value; }
  inline constexpr UIntClass const &get_num() { return value; }
  UIntClass value;
};

template <int nfixed, class UIntClass> class MaybeFixed
{
  static_assert(std::is_unsigned<UIntClass>::value, "MaybeFixed only works with unsigned types.");
  static_assert(nfixed >= 0, "nfixed negative and / or wrong specialization used.");

public:
  using type                    = UIntClass;
  using IsFixed_t               = std::true_type;
  static constexpr bool IsFixed = true;

  constexpr MaybeFixed() = default;
  //#ifdef DEBUG_SIEVE_LP_MATCHDIM
  //  constexpr  MaybeFixed(UIntClass const new_value) { assert(new_value == nfixed); }
  //#else
  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  constexpr MaybeFixed(Integer const)
  {
  }
  // #endif
  FORCE_INLINE inline constexpr operator UIntClass() const { return nfixed; }
  static inline constexpr UIntClass get_num() { return nfixed; }
  static constexpr UIntClass value = nfixed;
};

// Type normalization:
// We turn any Trait class T with T::value==false into a standard std::false_type
// and T::value==true into a standard std::true_type
// Note that such non-standard classes T encapsulating a bool constexpr may appear due to
// processing such types, depending on post-processing of types.
// Note: Currently used, but not really needed due to changes to trait getters.
template <class T> using NormalizeTrait = mystd::bool_constant<static_cast<bool>(T::value)>;

// MaybeConst<true,  T> = T const;
// MaybeConst<false, T> = T;
template <bool IsConst, class T>
using MaybeConst = mystd::conditional_t<IsConst, typename std::add_const<T>::type, T>;

// Z_NR - detection and modification...

/**
  IsZNRClass detects whether a class T is of the form Z_NR<ET> and allows to obtain ET.
  IsZZMatClass detects whether a class T is of the form ZZ_Mat<ET> and allows to obtain ET.

  PossiblyRemoveZNR  changes Z_NR<ET> into ET, leaving non-Z_NR<...> class untouched.
  PossiblyMpztToMpzClass changes mpz_t into mpz_class, leaving other classes untouched.
  PossiblyMpzClassToMpzt is the other way round.
  AddZNR turns ET into Z_NR<ET>, changing mpz_class to mpz_t as needed.
  AddZNR checks that ET is one of double, long or mpz_*.

  Except for IsZNRClass / IsZZMatClass, these are used as Result = Transform<Input>.
  no typename or ::type is needed.
*/

/**
  Detects whether a class T is of the form Z_NR<ET> and allows to obtain ET.
*/
template <class T> class IsZNRClass
{
public:
  using type                  = std::false_type;
  static bool constexpr value = false;
  constexpr operator bool() const { return false; }
};

template <class T> class IsZNRClass<fplll::Z_NR<T>>
{
public:
  using type                  = std::true_type;
  static bool constexpr value = true;
  constexpr operator bool() const { return true; }
  using GetUnderlyingType = T;
};

/**
  Detects whether T is of the form T = ZZ_mat<ET> and allows to recover ET.
*/
template <class T> struct IsZZMatClass
{
  using type                  = std::false_type;
  static bool constexpr value = false;
  operator bool() { return value; }
};

template <class ET> struct IsZZMatClass<fplll::ZZ_mat<ET>>
{
  using type                  = std::true_type;
  static bool constexpr value = true;
  operator bool() { return true; }
  using GetET = ET;
};

namespace MpzConversionHelper  // implementation details
{
template <class T> struct RemoveZNRHelper
{
  using type = T;
};
template <class T> struct RemoveZNRHelper<fplll::Z_NR<T>>
{
  using type = T;
};
template <class T> struct VerifyValidZNRArg
{
  static_assert(std::is_same<T, double>::value || std::is_same<T, mpz_t>::value ||
                    std::is_same<T, long>::value,
                "Invalid Argument to Z_NR");
  using type = T;
};
}  // end namespace MpzConversionHelper

// clang-format off
template<class T> using PossiblyRemoveZNR      = typename MpzConversionHelper::RemoveZNRHelper<T>::type;
template<class T> using PossiblyMpzClassToMpzt = mystd::conditional_t< std::is_same<T,mpz_class>::value, mpz_t, T>;
template<class T> using PossiblyMpztToMpzClass = mystd::conditional_t< std::is_same<T,mpz_t>::value, mpz_class, T>;
template<class T> using AddZNR = fplll::Z_NR<  typename MpzConversionHelper::VerifyValidZNRArg< PossiblyMpzClassToMpzt<T> >::type  >;
// clang-format on

/**
  Conversion to double that works with both mpz_class and anything convertible to double.
*/

// template<class Source> double convert_to_double(Source const & source);

template <class Source> double convert_to_double(Source const &source)
{
  static_assert(!std::is_same<Source, mpz_class>::value, "Source is mpz_class");
  return static_cast<double>(source);
}

double convert_to_double(mpz_class const &source) { return source.get_d(); }

/**
  Conversion to an integer type or mpz_class that works with both mpz_class and non-mpz types
*/

// helper class (wrapped in class to enable partial specializations)
// Used to define ConvertMaybeMPZ<Target>::do_convert<Source>()

namespace ConversionHelpers  // namespace for implementation details
{
// Integer is the convertion TARGET. There is specialization for mpz_class below
template <class Integer> struct ConvertMaybeMPZ
{
  static_assert(std::is_integral<Integer>::value, "Use only for integral classes.");

  // this restriction comes from mpz_* :
  static_assert(std::numeric_limits<Integer>::digits <= std::numeric_limits<long>::digits,
                "Converter does not work properly for types larger than long.");

  // clang-format off
  template <class Source,
  TEMPL_RESTRICT_DECL2(mystd::negation<std::is_same<mystd::decay_t<Source>,mpz_t>>,
                       mystd::negation<std::is_same<mystd::decay_t<Source>,mpz_class>>)>
  static Integer do_convert(Source &&source)
  {
    return static_cast<Integer>(std::forward<Source>(source));
  }
  // clang-format on

  static Integer do_convert(mpz_class const &source)
  {
    if (std::numeric_limits<Integer>::is_signed)
    {
      return static_cast<Integer>(source.get_si());
    }
    else
    {
      return static_cast<Integer>(source.get_ui());
    }
  }

  static Integer do_convert(mpz_t const &source)
  {
    if (std::numeric_limits<Integer>::is_signed)
    {
      return static_cast<Integer>(mpz_get_si(source));
    }
    else
    {
      return static_cast<Integer>(mpz_get_ui(source));
    }
  }
};

// specialization when the target type is mpz_class:
template <> struct ConvertMaybeMPZ<mpz_class>
{
  template <class Source> static mpz_class do_convert(Source const &source)
  {
    return static_cast<mpz_class>(source);
  }
  static mpz_class do_convert(mpz_class const &source) { return source; }
  static mpz_class do_convert(mpz_t const &source) { return static_cast<mpz_class>(source); }
};
}  // end namespace ConversionHelpers

template <class Target, class Source> Target convert_to_inttype(Source &&source)
{
  return ConversionHelpers::ConvertMaybeMPZ<Target>::do_convert(std::forward<Source>(source));
}

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

inline bool string_consume(std::istream &is, std::string const &str, bool elim_ws, bool verbose)
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

}  // end namespace

#endif  // GAUSS_SIEVE_UTILITY_H
