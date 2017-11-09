// This file contains utility (boilerplate) functions and classes only. It should not have
// dependencies on other files within the Gauss Sieve and be header-only.

#ifndef GAUSS_SIEVE_UTILITY_H
#define GAUSS_SIEVE_UTILITY_H

#include "DefaultIncludes.h"
#include <istream>
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"
#include "gmpxx.h"
#include "gmp.h"


namespace GaussSieve
{
template <class T> class IgnoreArg;
class IgnoreAnyArg;
template<class T> class IsZNRClass;
template<class Integer> struct ConvertMaybeMPZ;
template<class T> class UnZNR;
template<class Integer> struct ConvertMaybeMPZ;
}

/**
This macro is used to test the presence of a (public) member typedef in a class
Args:   TypeToCheck - typename whose presence to check
        CheckerClassName - Name of the checker class
This macro emits a new template class definition with the name CheckerClassName.
TypeToCheck must not be void (or another incomplete type)

Usage:
CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TypeToCheck, CheckerClassName);
This creates(!) the template class CheckerClassName.

Then CheckerClassName<SomeSuspiciousClass>::value will be true if
SomeSuspicousClass::TypeToCheck exists, false otherwise

The missing semicolon at the end of the macro is intentional.
The user needs to put it to emphasize that this is a declaration.
*/

// clang-format off

#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TypeToCheck, CheckerClassName)                           \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> static void                      foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t =                                                                                \
        mystd::bool_constant< !(std::is_void<decltype(foo<ClassToCheck>(0))>::value)>;      \
    static bool constexpr value = value_t::value;                                                  \
    constexpr operator bool() const { return value; };                                             \
  }

// clang-format on

/**
Similar to the above, creates a checker template class that checks wether
TypeToCheck exists and is equal to TypeShouldBe
*/

// clang-format off

#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(TypeToCheck, TypeShouldBe, CheckerClassName)      \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> static void                      foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t =                                                                                \
        mystd::bool_constant<                                                               \
                               std::is_same<TypeShouldBe, decltype(foo<ClassToCheck>(0))>::value>; \
    static bool constexpr value = value_t::value;                                                  \
    constexpr operator bool() const { return value; };                                             \
  }

// clang-format on

/**
  Checks whether TypeToCheck exists in TraitClass<ClassToCheck>.
*/

// clang-format off

#define CREATE_TRAIT_CHECK_CLASS(TraitClass, TypeToCheck, CheckerClassName)                    \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> static void                      foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t = mystd::bool_constant<                                                   \
        !(std::is_void<decltype(foo<TraitClass<ClassToCheck>>(0))>::value)>;                       \
    static bool constexpr value = value_t::value;                                                  \
    constexpr operator bool() const { return value; };                                             \
  }

// clang-format on

/**
  Checks whether TraitClass<ClassToCheck>::TypeToCheck exists and equals TypeShouldBe.
*/

// clang-format off

#define CREATE_TRAIT_EQUALS_CHECK(TraitClass, TypeToCheck, TypeShouldBe, CheckerClassName)         \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> static void                      foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t =                                                                                \
        mystd::bool_constant<                                                               \
                  std::is_same<TypeShouldBe, decltype(foo<TraitClass<ClassToCheck>>(0))>::value>;  \
    static bool constexpr value = value_t::value;                                                  \
    constexpr operator bool() const { return value; };                                             \
  }

// clang-format on

/**
  This is used to obtain traits from a trait class, with default settings.
  Notably CheckerClassName<T>::type is equal to
    TraitClass<T>::TypeToCheck if this exists,
    DefaultType otherwise.
*/

#define MAKE_TRAIT_GETTER(TraitClass, TypeToCheck, DefaultType, CheckerClassName)                  \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class...> static DefaultType foo(...);                                               \
                                                                                                   \
  public:                                                                                          \
    using type = decltype(foo<TraitClass<ClassToCheck>>(0));                                       \
  }


/**
  This is to improve readability of declarations.
  Usage:
  template<class Integer, TEMPL_RESTRICT_DECL((std::is_integral<Integer>::value))>
*/

/**
  NOTE: If you get an compiler error "no type named "type" in std::enable_if<false,int>,
        then you are using TEMPL_RESTRICT_* wrongly. Due to some choices in the C++ standard
        (to ensure backwards compatibility with certain internal workings of compilers),
        the condition argument(s) to TEMPL_RESTRIC_* must evaluate to true for at least one possible
        choice of template arguments. In case of template member functions of template classes,
        this means that for each instantiations of the class (i.e. class template parameters fixed),
        there has to be a set of template arguments for the function, such that this holds.
*/


#define TEMPL_RESTRICT_DECL(condition) typename std::enable_if<static_cast<bool>(condition),int>::type = 0
#define TEMPL_RESTRICT_IMPL(condition) typename std::enable_if<static_cast<bool>(condition),int>::type

// same as above, but accepts a type rather than a value. This allows for easier syntax.
// Limitation: If the argument contains a comma in its arguments (e.g. template<A,B>), it does not
// work

#define TEMPL_RESTRICT_DECL2(...) typename std::enable_if<static_cast<bool>(GaussSieve::mystd::conjunction<__VA_ARGS__>::value),int>::type = 0
#define TEMPL_RESTRICT_IMPL2(...) typename std::enable_if<static_cast<bool>(GaussSieve::mystd::conjunction<__VA_ARGS__>::value),int>::type

namespace GaussSieve
{

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

template <int nfixed = -1, class UIntClass = unsigned int> class MaybeFixed;

template <class UIntClass> class MaybeFixed<-1,UIntClass>
{
static_assert(std::is_unsigned<UIntClass>::value,"MaybeFixed only works with unsigned types.");
public:
  static constexpr bool IsFixed = false;
  using IsFixed_t  = std::false_type;
  using type = UIntClass;
  constexpr MaybeFixed(UIntClass const new_val) : value(new_val){};
  MaybeFixed() = default;  // Not sure whether we should allow uninitialized dims here. The issue is
                          // that we want the same interface in both cases.
  inline operator UIntClass() const { return value; };
  UIntClass value;
};

template <int nfixed, class UIntClass> class MaybeFixed
{
static_assert(std::is_unsigned<UIntClass>::value,"MaybeFixed only works with unsigned types.");
static_assert(nfixed>=0,"nfixed negative and / or wrong specialization used.");
public:
  using type = UIntClass;
  static constexpr bool IsFixed = true;

  using IsFixed_t  = std::true_type;
  constexpr MaybeFixed()         = default;
//#ifdef DEBUG_SIEVE_LP_MATCHDIM
//  constexpr  MaybeFixed(UIntClass const new_value) { assert(new_value == nfixed); }
//#else
  template<class Integer, typename std::enable_if<std::is_integral<Integer>::value,int>::type =0>
  constexpr MaybeFixed(Integer const){};
//#endif
  inline constexpr operator UIntClass() const { return nfixed; };
  static constexpr UIntClass value = nfixed;
};

// Type normalization:
// We turn any Trait class T with T::value==false into a standard std::false_type
// and T::value==true into a standard std::true_type
// Note that such non-standard classes T encapsulating a bool constexpr may appear due to
// processing such types. In particular, there are issues with the TRAIT_CHECK* macros:
// These internally use a std::is_same< > and the checker classes generated may not be standard.
// (This may be unnecessariy after recent changes to the checker classes)
// These classes help to remedy that:

template<bool> struct TypeNormalize_Helper;
template<> struct TypeNormalize_Helper<true >{using type=std::true_type; };
template<> struct TypeNormalize_Helper<false>{using type=std::false_type;};
template<class T> using NormalizeTrait = typename TypeNormalize_Helper<static_cast<bool>(T::value)>::type;

// MaybeConst<true,  T> = T const;
// MaybeConst<false, T> = T;
template<bool IsConst, class T>
using MaybeConst = typename std::conditional<IsConst,typename std::add_const<T>::type,T>::type;

// Z_NR - detection and modification...

/**
  Detects whether a class T is of the form Z_NR<ET> and allows to obtain ET.
*/

template<class T> class IsZNRClass
{
public:
  using type = std::false_type;
  static bool constexpr value = false;
  constexpr operator bool() const { return false; };
};

template<class T> class IsZNRClass<fplll::Z_NR<T>>
{
public:
  using type = std::true_type;
  static bool constexpr value = true;
  constexpr operator bool() const { return true; };
  using GetUnderlyingType = T;
};

/**
  Turns a Z_NR<long> into a long,
          Z_NR<double> into a double and
          Z_NR<mpz_t> into a mpz_class.
          The latter conversion is to get an arithmetic type that supports +, *, etc.
          If you want mpz_t, just use typename IsZNRClass<...>::GetUnderlyingType.

          Other types are unchanged.

          Usage: using NewType = typename UnZNR<OldType>::type
*/



template<class T> class UnZNR
{
  public: using type = T;
};

template<> class UnZNR<fplll::Z_NR<long>>
{
  public: using type = long;
};

template<> class UnZNR<fplll::Z_NR<double>>
{
  public: using type = double;
};

template<> class UnZNR<fplll::Z_NR<mpz_t>>
{
  public: using type = mpz_class;
};

/**
  Turns T into Z_NR<T>
  mpz_class is turned into Z_NR<mpz_t>
*/

template<class T> class AddZNR
{
  static_assert(std::is_same<T,long>::value || std::is_same<T,double>::value,
                "Unsupported type for ZNR");
  public: using type = fplll::Z_NR<T>;
};

template<> class AddZNR<mpz_t>
{
  public: using type = fplll::Z_NR<mpz_t>;
};

template<> class AddZNR<mpz_class>
{
  public: using type = fplll::Z_NR<mpz_t>;
};

/**
  Turns mpz_class into mpz_t
  Other classes are unchanged.
*/

template<class T> class FixMpz_classToMpz_t
{
  public: using type = T;
};

template<> class FixMpz_classToMpz_t<mpz_class>
{
  public: using type = mpz_t;
};

/**
  Detects whether T is of the form T = ZZ_mat<ET>
  and allows to recover ET.
*/

template<class T> class IsZZMatClass
{
  public:
  using type = std::false_type;
  static bool constexpr value = false;
  operator bool() {return value;}
};

template<class ET> class IsZZMatClass<fplll::ZZ_mat<ET>>
{
  public:
  using type = std::true_type;
  static bool constexpr value = true;
  operator bool() {return true;}
  using GetET = ET;
};

/**
  Conversion to double
*/

//template<class Source> double convert_to_double(Source const & source);

template<class Source>
double convert_to_double(Source const & source)
{
  static_assert(!std::is_same<Source,mpz_class>::value, "Source is mpz_class");
  return static_cast<double>(source);
}

double convert_to_double(mpz_class const & source)
{
  return source.get_d();
}

// ConvertMaybeMPZ<Integer>::convert_to_inttype(source)
// is a static_cast<Integer>(source) that also works for mpz_class
template<class Integer> struct ConvertMaybeMPZ;

template<class Integer> struct ConvertMaybeMPZ
{
  static_assert(std::is_integral<Integer>::value, "Use only for integral classes.");
  static_assert(std::numeric_limits<Integer>::digits <= std::numeric_limits<long>::digits,"Converter does not work properly for types larger than long.");
  // and in fact, there is issues with signed / unsigned.

  template<class Source>
  static Integer convert_to_inttype(Source const & source)
  {
    static_assert(!std::is_same<mystd::decay_t<Source>,mpz_class>::value, "Source is mpz_class");
    return static_cast<Integer>(source);
  }

  static Integer convert_to_inttype(mpz_class const & source)
  {
    if(std::numeric_limits<Integer>::is_signed)
    {
      return static_cast<Integer>(source.get_si() );
    }
    else
    {
      return static_cast<Integer>(source.get_ui() );
    }
  }

  static Integer convert_to_inttype(mpz_t const & source)
  {
    if(std::numeric_limits<Integer>::is_signed)
    {
      return static_cast<Integer>(mpz_get_si(source) );
    }
    else
    {
      return static_cast<Integer>(mpz_get_ui(source) );
    }
  }
};

template<>
struct ConvertMaybeMPZ<mpz_class>
{
  template<class Source>
  static mpz_class convert_to_inttype(Source const & source)
  {
    return static_cast<mpz_class>(source);
  }
  static mpz_class convert_to_inttype(mpz_class const & source) {return source;}
  static mpz_class convert_to_inttype(mpz_t const & source) {return static_cast<mpz_class>(source);}
};



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
