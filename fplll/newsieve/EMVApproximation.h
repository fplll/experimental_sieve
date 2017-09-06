#ifndef EMV_APPROXIMATION_H
#define EMV_APPROXIMATION_H

#include "DebugAll.h"
#include "Typedefs.h"
#include "SieveUtility.h"
#include <cstdint>
#include <type_traits>
#include "gmpxx.h"
#include <limits>
#include <cmath>
#include <iostream>

/** This defines a lattice point approximation, where the approximation consists of
a (shared) exponenent and a vector of >= 16-bit mantissas.
*/

namespace GaussSieve{

template<class DimensionType> class EMVApproximation;
class EMVApproximationTraits;
class EMVScalarProduct;

/**
  global constants:
*/
class EMVApproximationTraits
{
  public:
  // Note: This is only about the mantissas; there is an exponent in addition to those data.
  using ApproxEntryType = int_least16_t; // we *do* care about space.
  using ApproxNorm2Type = int_fast32_t;
};

//template<class DimensionType>
class EMVScalarProduct
{
  public:
  using MantissaType = typename EMVApproximationTraits::ApproxNorm2Type;
  static_assert(std::numeric_limits<MantissaType>::is_specialized,"Wrong ApproxNorm2Type");

  int exponent;
  MantissaType mantissa;
  inline double get_double() const;

  // default constructor
  explicit EMVScalarProduct(int const new_exponent, MantissaType const new_mantissa):
    exponent(new_exponent), mantissa(new_mantissa) {};

  // construct from integral or floating type
  template<class Integer, typename std::enable_if<std::is_integral<Integer>::value, int>::type=0>
  explicit EMVScalarProduct(Integer const source_arithmetic);

  template<class FloatType, typename std::enable_if<std::is_floating_point<FloatType>::value, int>::type=0>
  explicit EMVScalarProduct(FloatType source_float);

  explicit EMVScalarProduct(mpz_class const &source_mpz);
};


template<class DimensionType>
class EMVApproximation
{
  public:
  using ApproxEntryType = typename EMVApproximationTraits::ApproxEntryType;
  using ApproxNorm2Type = typename EMVApproximationTraits::ApproxNorm2Type;

  using AuxData = DimensionType; // No need for a traits class.

  template<class LatticePoint> // TODO : enable_if to select Lattice Points only.
  explicit EMVApproximation(LatticePoint const &exact_point);
  explicit EMVApproximation(DimensionType const &dim); //
  EMVApproximation() = delete;
  EMVApproximation(EMVApproximation const &old) = delete;
  EMVApproximation(EMVApproximation && old) = default;
  EMVApproximation& operator=(EMVApproximation const &other) = delete;
  EMVApproximation& operator=(EMVApproximation && other) = default;
  ~EMVApproximation();

  private:


};


/**
  Implementations:
*/

// EMVScalarProduct:

template<class Integer, typename std::enable_if<std::is_integral<Integer>::value, int>::type>
EMVScalarProduct::EMVScalarProduct(Integer const source_integer)
{
  static_assert(std::numeric_limits<Integer>::radix == 2,"Should never happen");
  static_assert(std::numeric_limits<MantissaType>::radix == 2,"Should never happen");
  constexpr int input_digits = std::numeric_limits<Integer>::digits;
  constexpr int mantissa_digits =std::numeric_limits<MantissaType>::digits;

  if (mantissa_digits >=input_digits)
  {
    exponent = 0;
    mantissa = source_integer;
  }
  else
  {
    // inefficient, but works:
    bool const sign = (source_integer < 0);
    typename std::make_unsigned<Integer>::type tmp = sign ? -source_integer : source_integer;
//    assert(tmp >=0);
    exponent = 0;
    while(tmp > (std::numeric_limits<MantissaType>::max() ))
    {
      ++exponent;
      tmp >>= 1;
    }
    mantissa = sign ? -tmp : tmp;
  }
}

template<class FloatType, typename std::enable_if<std::is_floating_point<FloatType>::value, int>::type>
EMVScalarProduct::EMVScalarProduct(FloatType source_float)
{

  // This is probably slow...

  static_assert(std::numeric_limits<MantissaType>::radix == 2,"Should never happen");
  constexpr int mantissa_digits =std::numeric_limits<MantissaType>::digits;

  source_float = std::frexp(source_float, &exponent); // This means that 2^exponent * source_float stores the correct value and abs(source_float) is in [1/2,1)

  source_float = std::ldexp(source_float, mantissa_digits);
  exponent-=mantissa_digits;
  // 2^exponent * mantissa is still correct.
  mantissa = std::trunc(source_float); // Note : Rounding towards zero is required to prevent overflow in corner cases.
}

EMVScalarProduct::EMVScalarProduct(mpz_class const & source_mpz)
{
  static_assert(std::numeric_limits<MantissaType>::radix == 2,"Should never happen");
  constexpr int mantissa_digits =std::numeric_limits<MantissaType>::digits;


  // This is equivalent to source_float  = std::frexp(source_float, & exponent) above:
  long tmp_exponent;
  double source_float = mpz_get_d_2exp(&tmp_exponent, source_mpz.get_mpz_t() );
  exponent = tmp_exponent; // convert to int.

  exponent-=mantissa_digits;
  mantissa = std::trunc( std::ldexp(source_float, mantissa_digits) );
}

inline bool operator< (EMVScalarProduct const & lhs, EMVScalarProduct const & rhs)
{
// We compare 2^lhs.exponent * mantissa < 2^rhs.exponent * mantissa
// The following works, but might need improvement:

  if(lhs.exponent > rhs.exponent)
  {
    return lhs.mantissa < (rhs.mantissa >> (lhs.exponent - rhs.exponent));
  }
  else
  {
    return (lhs.mantissa >> (rhs.exponent - lhs.exponent)) < rhs.mantissa;
  }
}

// inefficient, but should work:

template<class T,
typename std::enable_if < !std::is_same<typename std::decay<T>::type,EMVScalarProduct>::value, int>::type =0
>
inline bool operator< (EMVScalarProduct const & lhs, T && rhs)
{
  return lhs <  static_cast<EMVScalarProduct>(rhs);
}

template<class T,
typename std::enable_if < !std::is_same<typename std::decay<T>::type,EMVScalarProduct>::value, int>::type =0
>
inline bool operator< (T && lhs, EMVScalarProduct const & rhs)
{
  return static_cast<EMVScalarProduct>(lhs) < rhs;
}


template<class T,
typename std::enable_if < !std::is_same<typename std::decay<T>::type,EMVScalarProduct>::value, int>::type =0
>
inline bool operator> (EMVScalarProduct const & lhs, T && rhs)
{
  return std::forward<T>(rhs) < lhs;
}

template<class T,
typename std::enable_if < !std::is_same<typename std::decay<T>::type,EMVScalarProduct>::value, int>::type =0
>
inline bool operator> (T && lhs, EMVScalarProduct const & rhs)
{
  return rhs < std::forward<T>(lhs);
}

inline bool operator> (EMVScalarProduct const & lhs, EMVScalarProduct const rhs)
{
  return rhs < lhs;
}


inline std::ostream & operator<<(std::ostream &os, EMVScalarProduct const &approximated_number)
{
  os << approximated_number.mantissa << "x2^" << approximated_number.exponent;
  return os;
}




} // end namespace

#endif
