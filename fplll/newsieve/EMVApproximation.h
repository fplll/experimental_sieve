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
#include <vector>

/** This defines a lattice point approximation, where the approximation consists of
a (shared) exponenent and a vector of >= 16-bit mantissas.
*/

namespace GaussSieve{

template<int nfixed> class EMVApproximation;
class EMVApproximationTraits;
class EMVScalar;

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

class EMVScalar
{
  public:
  using MantissaType = typename EMVApproximationTraits::ApproxNorm2Type;
  static_assert(std::numeric_limits<MantissaType>::is_specialized,"Wrong ApproxNorm2Type");

  int exponent;
  MantissaType mantissa;
  inline double get_double() const { return std::ldexp(mantissa,exponent); };

  // default constructor
  explicit EMVScalar(int const new_exponent, MantissaType const new_mantissa):
    exponent(new_exponent), mantissa(new_mantissa) {};

  // construct from integral or floating type
  template<class Integer, TEMPL_RESTRICT_DECL((std::is_integral<Integer>::value))>
  explicit EMVScalar(Integer const source_arithmetic);

  template<class FloatType, TEMPL_RESTRICT_DECL((std::is_floating_point<FloatType>::value))>
  explicit EMVScalar(FloatType source_float);

  explicit EMVScalar(mpz_class const &source_mpz);


  // get_exponent(source) Returns an exponent, such that source = 2^exponent * x,
  // where 1/2 - eps <= |x| < 1.
  // (eps > 0 is a small error term that only occurs due to numerical issues in the current
  // for Integral source.)
  template<class Integer, TEMPL_RESTRICT_DECL((std::is_integral<Integer>::value))>
  static signed int get_exponent(Integer const source_int);
  template<class FloatType, TEMPL_RESTRICT_DECL((std::is_floating_point<FloatType>::value))>
  static signed int get_exponent(FloatType const source_float);
  static signed int get_exponent(mpz_class const &source_mpz);

  // divides by 2^exp. In the integer and mpz_class version, we assume that exp is positive
  template<class Integer, TEMPL_RESTRICT_DECL((std::is_integral<Integer>::value))>
  static Integer divide_by_power_of_2(Integer const source_int, unsigned int exponent);
  template<class FloatType, TEMPL_RESTRICT_DECL((std::is_floating_point<FloatType>::value))>
  static FloatType divide_by_power_of_2(FloatType const source_float, int exponent);
  static mpz_class divide_by_power_of_2(mpz_class const &source_mpz, unsigned int exponent);

};

#define FOR_FIXED_DIM template <int X = nfixed, typename std::enable_if<X >= 0, int>::type = 0>
#define FOR_VARIABLE_DIM template <int X = nfixed, typename std::enable_if<X == -1, int>::type = 0>

template<int nfixed>
class EMVApproximation
{
  friend
  StaticInitializer<EMVApproximation<nfixed>>;
  public:

  using AuxData = MaybeFixed<nfixed>; // No need for a traits class.
  using ScalarProductType = EMVScalar;

  using ApproxEntryType = typename EMVApproximationTraits::ApproxEntryType;
  using ApproxNorm2Type = typename EMVApproximationTraits::ApproxNorm2Type;

  using Container = typename std::conditional<nfixed >= 0,
                      std::array<ApproxEntryType, nfixed >=0 ? nfixed:0>,  // if nfixed >= 0
                      std::vector<ApproxEntryType> >                       // if nfixed <  0
                      ::type;


  template<class LatticePoint> // TODO : enable_if to select Lattice Points only (having [])
  explicit EMVApproximation(LatticePoint const &exact_point);

  //explicit EMVApproximation(AuxData const &dim); //
  EMVApproximation() = delete;
  EMVApproximation(EMVApproximation const &old) = delete;
  EMVApproximation(EMVApproximation && old) = default;
  EMVApproximation& operator=(EMVApproximation const &other) = delete;
  EMVApproximation& operator=(EMVApproximation && other) = default;
  ~EMVApproximation() {};

  ApproxEntryType &operator[](uint_fast16_t idx) { return data[idx]; };
  ApproxEntryType const &operator[](uint_fast16_t idx) const { return data[idx]; };

  FOR_FIXED_DIM
  static constexpr MaybeFixed<nfixed> get_dim() { return MaybeFixed<nfixed>(nfixed); }

  FOR_VARIABLE_DIM
  static MaybeFixed<-1> get_dim() { return dim; }

  FOR_FIXED_DIM
  void reserve_size() {}

  FOR_VARIABLE_DIM
  void reserve_size() {data.reserve(static_cast<typename Container::size_type>(dim));}

  public:

  signed int exponent; // shared exponent
  Container data; // array or vector of 16-bit ints
  private:
  static MaybeFixed<nfixed> dim; // dimension of data
  static unsigned int max_bits; // maximum number of bits inside each entry.
                                // This depends on the bitlenght of both ApproxEntryType
                                // and ApproxNorm2Type, as well as on dim.
#ifdef DEBUG_SIEVE_LP_INIT
  static bool class_initialized;
#endif // DEBUG_SIEVE_LP_INIT

};



// initialization of static data of template class.
template<int nfixed> MaybeFixed<nfixed> EMVApproximation<nfixed>::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);
template<int nfixed> unsigned int EMVApproximation<nfixed>::max_bits  = 0;
#ifdef DEBUG_SIEVE_LP_INIT
template<int nfixed> bool EMVApproximation<nfixed>::class_initialized = false;
#endif

// Static Initializer:
template<int nfixed> class StaticInitializer<EMVApproximation<nfixed> >
{
  public:
  StaticInitializer(MaybeFixed<nfixed> const new_dim)
  {
    if(user_counter>0)
    {
      assert((new_dim == EMVApproximation<nfixed>::dim));
      // TODO: Throw exception!
    }
    else
    {
      EMVApproximation<nfixed>::dim = new_dim;

      static_assert (std::numeric_limits<typename EMVApproximationTraits::ApproxEntryType>::is_specialized, "bad ApproxType");
      static_assert (std::numeric_limits<typename EMVApproximationTraits::ApproxNorm2Type>::is_specialized, "bad ApproxTypeNorm2");

      unsigned int constexpr entry_digits = std::numeric_limits<typename EMVApproximationTraits::ApproxEntryType>::digits;
      unsigned int constexpr result_digits= std::numeric_limits<typename EMVApproximationTraits::ApproxNorm2Type>::digits;

      // We need to ensure that (2^max_bits - 1)^2 * dim <= 2^result_digits
      // and max_bits <= entry_digits
      // in order to avoid overflows.
      unsigned int const extra_bits  = std::floor(std::log2(static_cast<int>(new_dim)));
      assert(result_digits - extra_bits>=0);
      EMVApproximation<nfixed>::max_bits= std::min(entry_digits, static_cast<unsigned int>(( result_digits - extra_bits) /2) ) ;
    }
    ++user_counter;
#ifdef DEBUG_SIEVE_LP_INIT
    EMVApproximation<nfixed>::class_initialized = true;
#endif
  }

  ~StaticInitializer()
  {
    assert(user_counter > 0);
    --user_counter;
#ifdef DEBUG_SIEVE_LP_INIT
    EMVApproximation<nfixed>::class_initialized = (user_counter > 0);
#endif
  }

  static bool is_initialized(){ return user_counter > 0; }; // Does an object exist?

  private:
  static unsigned int user_counter; // counts number of StaticInitializer instances.
};

template<int nfixed>
unsigned int StaticInitializer<EMVApproximation<nfixed>>::user_counter = 0;



/*******************
 * Implementations *
 *******************/

/**
  EMVScalar Implementation
*/

// static helper functions.

template<class Integer, TEMPL_RESTRICT_IMPL((std::is_integral<Integer>::value))>
signed int EMVScalar::get_exponent(Integer const source_int)
{
  static_assert(std::numeric_limits<Integer>::radix == 2,"Should never happen");
  signed int ret;
  std::frexp(source_int,&ret); // Note: This is inefficient and has rounding issues.
  // Alas, there is no easy, fast and portable way to obtain the (signed) msb in c++ (c++17, even!)
  // If this is time-critical, consider __builtin's (Clang and gcc offer them for this!)

  // The issue with rounding is that if source_int is a 64 bit long slightly smaller than a large power of 2,
  // then the result is too large.
  return ret;
}

template<class FloatType, TEMPL_RESTRICT_IMPL((std::is_floating_point<FloatType>::value))>
signed int EMVScalar::get_exponent(FloatType const source_float)
{
  signed int ret;
  std::frexp(source_float, &ret);
  return ret;
}

signed int EMVScalar::get_exponent(mpz_class const &source_mpz)
{
  long ret; // mpz_get_d_2exp requires long
  mpz_get_d_2exp(&ret, source_mpz.get_mpz_t() );
  return ret;
}

template<class Integer, TEMPL_RESTRICT_IMPL((std::is_integral<Integer>::value))>
Integer EMVScalar::divide_by_power_of_2(Integer const source_int, unsigned int exponent)
{
  return source_int / (1<<exponent); // Note: We don't right-shift source_int, because it might be
                                     // negative
}
template<class FloatType, TEMPL_RESTRICT_IMPL((std::is_floating_point<FloatType>::value))>
FloatType EMVScalar::divide_by_power_of_2(FloatType const source_float, int exponent)
{
  return std::ldexp(source_float,-exponent);
}

mpz_class EMVScalar::divide_by_power_of_2(mpz_class const &source_mpz, unsigned int exponent)
{
  mpz_class ret;
  mpz_tdiv_q_2exp(ret.get_mpz_t(), source_mpz.get_mpz_t(),exponent );
  return ret;
}


// constructors of EMVScalar

template<class Integer, typename std::enable_if<std::is_integral<Integer>::value, int>::type>
EMVScalar::EMVScalar(Integer const source_integer)
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
EMVScalar::EMVScalar(FloatType source_float)
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

EMVScalar::EMVScalar(mpz_class const & source_mpz)
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

// comparison operators

inline bool operator< (EMVScalar const & lhs, EMVScalar const & rhs)
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

template<class T,
typename std::enable_if < !std::is_same<typename std::decay<T>::type,EMVScalar>::value, int>::type =0
>
inline bool operator< (EMVScalar const & lhs, T && rhs)
{
  return lhs <  static_cast<EMVScalar>(rhs);
}

template<class T,
typename std::enable_if < !std::is_same<typename std::decay<T>::type,EMVScalar>::value, int>::type =0
>
inline bool operator< (T && lhs, EMVScalar const & rhs)
{
  return static_cast<EMVScalar>(lhs) < rhs;
}


template<class T,
typename std::enable_if < !std::is_same<typename std::decay<T>::type,EMVScalar>::value, int>::type =0
>
inline bool operator> (EMVScalar const & lhs, T && rhs)
{
  return std::forward<T>(rhs) < lhs;
}

template<class T,
typename std::enable_if < !std::is_same<typename std::decay<T>::type,EMVScalar>::value, int>::type =0
>
inline bool operator> (T && lhs, EMVScalar const & rhs)
{
  return rhs < std::forward<T>(lhs);
}

inline bool operator> (EMVScalar const & lhs, EMVScalar const rhs)
{
  return rhs < lhs;
}


// output of EMVScalars

inline std::ostream & operator<<(std::ostream &os, EMVScalar const &approximated_number)
{
  os << approximated_number.mantissa << "x2^" << approximated_number.exponent;
  return os;
}


/**
  EMVApproximation
*/

// Constructor:

template<int nfixed>
template<class LatticePoint>
EMVApproximation<nfixed>::EMVApproximation(LatticePoint const &exact_point):
exponent(0),
data()
{
  reserve_size();
  auto const dimension = get_dim();
#ifdef DEBUG_SIEVE_LP_MATCHDIM
  assert(dimension == exact_point.get_dim() );
  assert(dimension == exact_point.get_vec_size() );
#endif
#ifdef DEBUG_SIEVE_LP_INIT
  assert(class_initialized);
#endif
  using CooType = typename GetCooType<LatticePoint>::type;
  using std::abs;

  CooType max_entry = 0;
  for(uint_fast16_t i=0;i<dimension;++i)
  {
    // The static_cast<CooType> is there to disable lazy evaluation (via expression templates) done
    // in mpz_class. The issue is that this confuses template argument deduction for std::max.
    max_entry = std::max(max_entry, static_cast<CooType>(abs(exact_point[i])) );
  }

  if(max_entry ==0 )
  {
    exponent = 0;
    for(uint_fast16_t i=0; i<dimension; ++i)
    {
      data[i] = 0;
    }
    std::cerr << "Warning: approximating all-zero vector." << std::endl; // TODO : Remove
  }
  else // max_entry !=0
  {
    exponent = EMVScalar::get_exponent(max_entry) - max_bits;

    // This avoids approximating small numbers by big ones, a la 1 = 2^-15 * 32728,
    // which is correct, but hurts readability

    // Note that std::numeric_limits<mpz_class> is specialized by gmpxx, so this case is caught.
    if(std::numeric_limits<CooType>::is_integer) // constexpr if
    {
      exponent = std::max(0,exponent);
    }

    for(uint_fast16_t i=0;i<dimension;++i)
    {
      data[i] = ConvertMaybeMPZ<ApproxEntryType>::convert_to_inttype(
            EMVScalar::divide_by_power_of_2(exact_point[i],exponent) );
    }
  }
}

// actual scalar product

template<int nfixed>
inline auto approximate_scalar_product(EMVApproximation<nfixed> const &lhs, EMVApproximation<nfixed> const &rhs)
-> EMVScalar
{
  static_assert(std::is_same<EMVScalar, typename EMVApproximation<nfixed>::ScalarProductType>::value,"");
  using ProductType = typename EMVApproximationTraits::ApproxNorm2Type;

  ProductType scp =0;
  auto const dim = rhs.get_dim();
  for(int i=0;i<dim;++i)
  {
    scp += lhs[i] * rhs[i];
  }

  EMVScalar result(lhs.exponent + rhs.exponent, scp);
  return result;
}




// output

template<int nfixed>
inline std::ostream & operator<<(std::ostream &os, EMVApproximation<nfixed> const &approximated_vector)
{
  os << "Approximation :2^" << approximated_vector.exponent <<"x [";
  auto const dim = approximated_vector.get_dim();
  for(uint_fast16_t i = 0;i<dim;++i)
  {
    os << approximated_vector.data[i] << " ";
  }
  os << "]" << std::endl;
  return os;
}


} // end namespace

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM


#endif
