#warning This file is not supposed to be used at the moment

#ifndef EMV_APPROXIMATION_H
#define EMV_APPROXIMATION_H

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"
#include "SieveUtility.h"
#include "gmpxx.h"

/**
  This defines a lattice point approximation, where the approximation consists of
  a (shared) exponenent and a vector of >= 16-bit mantissas, i.e. we approximate by
  2^exponent * (vector of 16-bit ints)
  The idea is that in order to test whether a scalar product between two vectors is larger than
  a certain threshold (that may depend on the lengths of the vectors), it is enough to consider
  the most significant bits of the invidual vectors.
  The advantage is that scalar product computations of 16-bit ints are much faster than of, say,
  mpz_class. Also, 16-bit scalar product computations can make use of vectorized instructions.
  (provided the architecture has vectorized instructions for ints rather than floats, which is
  somewhat uncommon at the time of writing, but exists).
  Essentially, such an approximation makes the length of the input data almost irrevalant and allows
  to optimize for ONE given data size.

  An alternative would be to simply approximate by a vector of doubles. Note that the difference
  here is that the exponent is shared and we use fewer bits.

  TODO: Test an alternative with floats. The main downside is that this could be the main
        contribution to total memory consumption. Ideally, we would like half-precision floats...
        (which exists either as slow software implementations or compiler && architecture-dependent
        C++ extensions (e.g. __fp16 with GCC on ARM using the -mfp16-format compiler option)
*/

namespace GaussSieve
{

template <int nfixed> class EMVApproximation;
class EMVApproximationTraits;
class EMVScalar;

/**
  global constants:
*/
class EMVApproximationTraits
{
public:
  // Note: ApproxEntryType and ApproxNorm2Type are only about the mantissas;
  // there is an exponent in addition to those data.
  // Also note that we use slightly less than 16 bits (depending on the dimension) to avoid over-
  // flows when computing scalar products.
  using ApproxEntryType = int_least16_t;  // we *do* care about space.
  using ApproxNorm2Type = int_fast32_t;
};

/************************
Approximation of a scalar:

This stores an exponent and a mantissa. Of course, we could just use floating point types here,
but this is more compatible with the rest of the code. Converts to a double, still.

Usage: EMVScalar can be constructed from any integral or floating type or from mpz_class.
       It currently supports only limited arithmetic.
       Can convert to double / float.
*************************/

class EMVScalar
{
public:
  using HasDefaultStaticInitializer = std::true_type;
  using MantissaType                = typename EMVApproximationTraits::ApproxNorm2Type;
  static_assert(std::numeric_limits<MantissaType>::is_specialized, "Wrong ApproxNorm2Type");

  int exponent;
  MantissaType mantissa;
  // clang-format off
  inline double get_double() const { return std::ldexp(mantissa, exponent); }
  explicit operator double()       { return std::ldexp(mantissa, exponent); }
  explicit operator float()        { return std::ldexp(mantissa, exponent); }
  // clang-format on

  // constructors
  constexpr explicit EMVScalar(int const new_exponent, MantissaType const new_mantissa)
      : exponent(new_exponent), mantissa(new_mantissa){};
  // default copy constructors (they are automatic anyway)
  // constexpr EMVScalar(EMVScalar const &) = default;
  // EMVScalar(EMVScalar &&) = default;

  // construct from integral or floating type
  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  explicit EMVScalar(Integer const source_arithmetic);

  template <class FloatType, TEMPL_RESTRICT_DECL2(std::is_floating_point<FloatType>)>
  explicit EMVScalar(FloatType source_float);

  explicit EMVScalar(mpz_class const &source_mpz);

  constexpr EMVScalar operator-() const & { return EMVScalar{exponent, -mantissa}; }
  EMVScalar operator-() &&
  {
    mantissa = -mantissa;
    return *this;
  }

  // TODO: Return value of operators does not adhere usual conventions.
  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline void operator>>=(Integer const &shift)
  {
    exponent -= shift;
  }

  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline void operator<<=(Integer const &shift)
  {
    exponent += shift;
  }

  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline EMVScalar operator>>(Integer const &shift) const &
  {
    EMVScalar retval{*this};
    retval >>= shift;
    return retval;
  }

  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline EMVScalar operator>>(Integer const &shift) &&
  {
    (*this) >>= shift;
    return *this;
  }
  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline EMVScalar operator<<(Integer const &shift) const &
  {
    EMVScalar retval{*this};
    retval <<= shift;
    return retval;
  }
  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline EMVScalar operator<<(Integer const &shift) &&
  {
    (*this) <<= shift;
    return *this;
  }

  inline friend EMVScalar abs(EMVScalar const &arg)
  {
    using std::abs;
    return EMVScalar{arg.exponent, abs(arg.mantissa)};
  }
  inline friend EMVScalar abs(EMVScalar &&arg)
  {
    using std::abs;
    arg.mantissa = abs(arg.mantissa);
    return arg;
  }

  // helper functions: included as static functions tied to the class:

  // get_exponent(source) Returns an exponent, such that source = 2^exponent * x,
  // where 1/2 - eps <= |x| < 1.
  // (eps > 0 is a small error term that only occurs due to numerical issues in the current
  // for Integral source.)
  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  static signed int get_exponent(Integer const source_int);
  template <class FloatType, TEMPL_RESTRICT_DECL2(std::is_floating_point<FloatType>)>
  static signed int get_exponent(FloatType const source_float);
  static signed int get_exponent(mpz_class const &source_mpz);

  // divides by 2^exp. In the integer and mpz_class version, we assume that exp is positive
  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  constexpr static Integer divide_by_power_of_2(Integer const source_int, unsigned int exponent);
  template <class FloatType, TEMPL_RESTRICT_DECL2(std::is_floating_point<FloatType>)>
  static FloatType divide_by_power_of_2(FloatType const source_float, int exponent);
  static mpz_class divide_by_power_of_2(mpz_class const &source_mpz, unsigned int exponent);
};

// clang-format off
#define FOR_FIXED_DIM    template<int dummy_nfixed = nfixed, TEMPL_RESTRICT_DECL(dummy_nfixed >= 0)>
#define FOR_VARIABLE_DIM template<int dummy_nfixed = nfixed, TEMPL_RESTRICT_DECL(dummy_nfixed ==-1)>
// clang-format on

/************************
Approximation to a vector
************************/

// clang-format off
template<int nfixed>
class EMVApproximation
{
  friend StaticInitializer<EMVApproximation<nfixed>>;
public:
  using ScalarProductType = EMVScalar;

private:
  using ApproxEntryType = typename EMVApproximationTraits::ApproxEntryType;
  using ApproxNorm2Type = typename EMVApproximationTraits::ApproxNorm2Type;

  using Container = mystd::conditional_t<  (nfixed >= 0),
                        std::array<ApproxEntryType, (nfixed >= 0 ? nfixed : 0)>,  // if nfixed >= 0
                        std::vector<ApproxEntryType>  >;                          // if nfixed <  0
  // clang-format on
public:
  template <class LatticePoint, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatticePoint>)>
  explicit EMVApproximation(LatticePoint const &exact_point);

  // clang-format off
  EMVApproximation()                         		= delete;
  EMVApproximation(EMVApproximation const &) 		= delete;
  EMVApproximation(EMVApproximation &&)      		= default;
  EMVApproximation &operator=(EMVApproximation const &) = delete;
  EMVApproximation &operator=(EMVApproximation &&) 	= default;
  ~EMVApproximation(){};
  // clang-format on

  ApproxEntryType &operator[](uint_fast16_t idx) { return data[idx]; };
  ApproxEntryType const &operator[](uint_fast16_t idx) const { return data[idx]; };

  FOR_FIXED_DIM
  static constexpr MaybeFixed<nfixed> get_dim() { return MaybeFixed<nfixed>(nfixed); }

  FOR_VARIABLE_DIM
  static MaybeFixed<-1> get_dim() { return dim; }

  FOR_FIXED_DIM
  inline void reserve_size() {}

  FOR_VARIABLE_DIM
  inline void reserve_size() { data.reserve(static_cast<typename Container::size_type>(dim)); }
  // clang-format off
public:  // unsure about data members
  signed int exponent;  // shared exponent
  Container data;       // array or vector of 16-bit ints
  // Note: Might want to avoid storing exponent twice.
  EMVScalar get_norm2() const { return approx_norm2; }
  // clang-format on

private:
  EMVScalar approx_norm2;
  inline static bool is_class_initialized()  // checks static initialization
  {
    return StaticInitializer<EMVApproximation<nfixed>>::is_initialized();
  }
  static MaybeFixed<nfixed> dim;  // dimension of data
  static unsigned int max_bits;   // maximum number of bits inside each entry.
                                  // This depends on the bitlenght of both ApproxEntryType
                                  // and ApproxNorm2Type, as well as on dim.
};

// initialization of static data of template class.
// clang-format off
template<int nfixed> MaybeFixed<nfixed> EMVApproximation<nfixed>::dim      = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);
template<int nfixed> unsigned int       EMVApproximation<nfixed>::max_bits = 0;
// clang-format on

// Static Initializer:
template <int nfixed>
class StaticInitializer<EMVApproximation<nfixed>> final
    : public DefaultStaticInitializer<EMVApproximation<nfixed>>
{
public:
  using Parent = DefaultStaticInitializer<EMVApproximation<nfixed>>;

  // clang-format off
  template <class InitArg, TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<InitArg>)>
  StaticInitializer(InitArg const &init_arg)
      : StaticInitializer(init_arg.dim) {}
  // clang-format on

  StaticInitializer(MaybeFixed<nfixed> const new_dim)
  {
    // note: user_count increased by one in initializer list prior to this constructor body.
    assert(Parent::user_count > 0);
    if (Parent::user_count > 1)
    {
      assert((new_dim == EMVApproximation<nfixed>::dim));  // TODO: Throw exception!
    }
    else
    {
      EMVApproximation<nfixed>::dim = new_dim;
      // clang-format off
      static_assert(std::numeric_limits<typename EMVApproximationTraits::ApproxEntryType>::is_specialized, "bad ApproxType");
      static_assert(std::numeric_limits<typename EMVApproximationTraits::ApproxNorm2Type>::is_specialized, "bad ApproxTypeNorm2");

      unsigned int constexpr entry_digits  = std::numeric_limits<typename EMVApproximationTraits::ApproxEntryType>::digits;
      unsigned int constexpr result_digits = std::numeric_limits<typename EMVApproximationTraits::ApproxNorm2Type>::digits;
      // clang-format on

      // We need to ensure that (2^max_bits - 1)^2 * dim <= 2^result_digits
      // and max_bits <= entry_digits
      // in order to avoid overflows.
      // TODO: Use const-expr-compatible computation (Note std::min is not in C++11)
      unsigned int const extra_bits = std::floor(std::log2(static_cast<int>(new_dim)));
      assert(result_digits - extra_bits >= 0);
      EMVApproximation<nfixed>::max_bits =
          std::min(entry_digits, static_cast<unsigned int>((result_digits - extra_bits) / 2));
    }
  }

  ~StaticInitializer() {}

  inline static bool is_initialized() { return Parent::user_count > 0; };  // Does an object exist?
};

/*******************
 * Implementations *
 *******************/

/**
  EMVScalar Implementation
*/

// static helper functions.

template <class Integer, TEMPL_RESTRICT_IMPL2(std::is_integral<Integer>)>
signed int EMVScalar::get_exponent(Integer const source_int)
{
  static_assert(std::numeric_limits<Integer>::radix == 2, "Should never happen");
  signed int ret;
  std::frexp(source_int, &ret);  // Note: This is inefficient and has rounding issues.
  // Alas, there is no easy, fast and portable way to obtain the (signed) msb in c++ (c++17, even!)
  // If this is time-critical, consider __builtin's (Clang and gcc offer them for this!)

  // The issue with rounding is that if source_int is a 64 bit long slightly smaller than a large
  // power of 2, then the result is too large.
  return ret;
}

template <class FloatType, TEMPL_RESTRICT_IMPL2(std::is_floating_point<FloatType>)>
signed int EMVScalar::get_exponent(FloatType const source_float)
{
  signed int ret;
  std::frexp(source_float, &ret);
  return ret;
}

signed int EMVScalar::get_exponent(mpz_class const &source_mpz)
{
  long ret;  // mpz_get_d_2exp requires long
  mpz_get_d_2exp(&ret, source_mpz.get_mpz_t());
  return ret;
}

template <class Integer, TEMPL_RESTRICT_IMPL2(std::is_integral<Integer>)>
constexpr Integer EMVScalar::divide_by_power_of_2(Integer const source_int, unsigned int exponent)
{
  return source_int / (1 << exponent);  // Note: We don't right-shift source_int, because it
                                        // might be negative
}
template <class FloatType, TEMPL_RESTRICT_IMPL2(std::is_floating_point<FloatType>)>
FloatType EMVScalar::divide_by_power_of_2(FloatType const source_float, int exponent)
{
  return std::ldexp(source_float, -exponent);
}

mpz_class EMVScalar::divide_by_power_of_2(mpz_class const &source_mpz, unsigned int exponent)
{
  mpz_class ret;
  mpz_tdiv_q_2exp(ret.get_mpz_t(), source_mpz.get_mpz_t(), exponent);
  return ret;
}

// constructors of EMVScalar

template <class Integer, TEMPL_RESTRICT_IMPL2(std::is_integral<Integer>)>
EMVScalar::EMVScalar(Integer const source_integer)
{
  static_assert(std::numeric_limits<Integer>::radix == 2, "Should never happen");
  static_assert(std::numeric_limits<MantissaType>::radix == 2, "Should never happen");
  constexpr int input_digits    = std::numeric_limits<Integer>::digits;
  constexpr int mantissa_digits = std::numeric_limits<MantissaType>::digits;

  if (mantissa_digits >= input_digits)
  {
    exponent = 0;
    mantissa = source_integer;
  }
  else
  {
    // clang-format off
    // inefficient, but works:
    bool const                               sign = (source_integer < 0);
    typename mystd::make_unsigned_t<Integer> tmp  = sign ? -source_integer : source_integer;
    // clang-format on
    exponent = 0;
    while (tmp > (std::numeric_limits<MantissaType>::max()))
    {
      ++exponent;
      tmp >>= 1;
    }
    mantissa = sign ? -tmp : tmp;
  }
}

template <class FloatType, TEMPL_RESTRICT_IMPL2(std::is_floating_point<FloatType>)>
EMVScalar::EMVScalar(FloatType source_float)
{
  // This is probably slow...

  static_assert(std::numeric_limits<MantissaType>::radix == 2, "Should never happen");
  constexpr int mantissa_digits = std::numeric_limits<MantissaType>::digits;

  // The following means that 2^exponent * source_float stores the correct value and
  // abs(source_float) is in the range [1/2,1)
  source_float = std::frexp(source_float, &exponent);

  source_float = std::ldexp(source_float, mantissa_digits);
  exponent -= mantissa_digits;
  // 2^exponent * mantissa is still correct.
  mantissa = std::trunc(source_float);  // Note : Rounding towards zero is required to prevent
                                        // overflow in corner cases.
}

EMVScalar::EMVScalar(mpz_class const &source_mpz)
{
  static_assert(std::numeric_limits<MantissaType>::radix == 2, "Should never happen");
  constexpr int mantissa_digits = std::numeric_limits<MantissaType>::digits;

  // This is equivalent to source_float  = std::frexp(source_float, & exponent) above:
  long tmp_exponent;
  double source_float = mpz_get_d_2exp(&tmp_exponent, source_mpz.get_mpz_t());
  exponent            = tmp_exponent;  // convert to int.

  exponent -= mantissa_digits;
  mantissa = std::trunc(std::ldexp(source_float, mantissa_digits));
}

/*************************************
Non-member functions:
*************************************/

// comparison operators

inline bool operator<(EMVScalar const &lhs, EMVScalar const &rhs)
{
  // We compare 2^lhs.exponent * mantissa < 2^rhs.exponent * mantissa
  // The following works, but might need improvement:

  if (lhs.exponent > rhs.exponent)
  {
    return lhs.mantissa < (rhs.mantissa >> (lhs.exponent - rhs.exponent));
  }
  else
  {
    return (lhs.mantissa >> (rhs.exponent - lhs.exponent)) < rhs.mantissa;
  }
}

inline bool operator<=(EMVScalar const &lhs, EMVScalar const &rhs)
{
  // We compare 2^lhs.exponent * mantissa < 2^rhs.exponent * mantissa
  // The following works, but might need improvement:

  if (lhs.exponent > rhs.exponent)
  {
    return lhs.mantissa <= (rhs.mantissa >> (lhs.exponent - rhs.exponent));
  }
  else
  {
    return (lhs.mantissa >> (rhs.exponent - lhs.exponent)) <= rhs.mantissa;
  }
}

// other cases forward to the two cases above, possibly after conversion:

template <class T, TEMPL_RESTRICT_DECL(!(std::is_same<mystd::decay_t<T>, EMVScalar>::value))>
inline bool operator<(EMVScalar const &lhs, T &&rhs)
{
  return lhs < static_cast<EMVScalar>(std::forward<T>(rhs));
}

template <class T, TEMPL_RESTRICT_DECL(!(std::is_same<mystd::decay_t<T>, EMVScalar>::value))>
inline bool operator<=(EMVScalar const &lhs, T &&rhs)
{
  return lhs <= static_cast<EMVScalar>(std::forward<T>(rhs));
}

template <class T, TEMPL_RESTRICT_DECL(!(std::is_same<mystd::decay_t<T>, EMVScalar>::value))>
inline bool operator<(T &&lhs, EMVScalar const &rhs)
{
  return static_cast<EMVScalar>(std::forward<T>(lhs)) < rhs;
}

template <class T, TEMPL_RESTRICT_DECL(!(std::is_same<mystd::decay_t<T>, EMVScalar>::value))>
inline bool operator<=(T &&lhs, EMVScalar const &rhs)
{
  return static_cast<EMVScalar>(std::forward<T>(lhs)) <= rhs;
}

template <class T, TEMPL_RESTRICT_DECL(!(std::is_same<mystd::decay_t<T>, EMVScalar>::value))>
inline bool operator>(EMVScalar const &lhs, T &&rhs)
{
  return std::forward<T>(rhs) < lhs;
}

template <class T, TEMPL_RESTRICT_DECL(!(std::is_same<mystd::decay_t<T>, EMVScalar>::value))>
inline bool operator>(T &&lhs, EMVScalar const &rhs)
{
  return rhs < std::forward<T>(lhs);
}

inline bool operator>(EMVScalar const &lhs, EMVScalar const rhs) { return rhs < lhs; }

template <class T, TEMPL_RESTRICT_DECL(!(std::is_same<mystd::decay_t<T>, EMVScalar>::value))>
inline bool operator>=(EMVScalar const &lhs, T &&rhs)
{
  return std::forward<T>(rhs) <= lhs;
}

template <class T, TEMPL_RESTRICT_DECL(!(std::is_same<mystd::decay_t<T>, EMVScalar>::value))>
inline bool operator>=(T &&lhs, EMVScalar const &rhs)
{
  return rhs <= std::forward<T>(lhs);
}

// clang-format off
inline bool operator>=(EMVScalar const &lhs, EMVScalar const rhs)
{
  return rhs <= lhs;
}
// clang-format on

// output of EMVScalars

inline std::ostream &operator<<(std::ostream &os, EMVScalar const &approximated_number)
{
  os << approximated_number.mantissa << "x2^" << approximated_number.exponent;
  return os;
}

/**
  EMVApproximation implementation
*/

// Constructor:

// clang-format off
template <int nfixed>
template <class LatticePoint, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatticePoint>)>
EMVApproximation<nfixed>::EMVApproximation(LatticePoint const &exact_point)
    : exponent(0),
      data(),
      approx_norm2(static_cast<EMVScalar>(exact_point.get_norm2()))
{
  // TODO: static_assert() Traits for LatticePoint?
  // clang-format on
  reserve_size();
  auto const dimension = get_dim();
#ifdef DEBUG_SIEVE_LP_MATCHDIM
  assert(dimension == exact_point.get_dim());
//  assert(dimension == exact_point.get_internal_rep_size() );
#endif
#ifdef DEBUG_SIEVE_LP_INIT
  assert(is_class_initialized());
#endif
  using AbsoluteCooType = Get_AbsoluteCooType<LatticePoint>;
  using std::abs;

  AbsoluteCooType max_entry = 0;
  for (uint_fast16_t i = 0; i < dimension; ++i)
  {
    // clang-format off
    // The static_cast<AbsoluteCooType> is there to disable lazy evaluation (via expression templates)
    // done inside mpz_class. The issue is that this confuses template argument deduction for std::max.
    max_entry = std::max(max_entry, static_cast<AbsoluteCooType>(abs(exact_point.get_absolute_coo(i))));
    // clang-format on
  }

  if (max_entry == 0)
  {
    exponent = 0;
    for (uint_fast16_t i = 0; i < dimension; ++i)
    {
      data[i] = 0;
    }
    // std::cerr << "Warning: approximating all-zero vector." << std::endl; // TODO : Remove
  }
  else  // max_entry !=0
  {
    exponent = EMVScalar::get_exponent(max_entry) - max_bits;

    // This avoids approximating small numbers by big ones, a la 1 = 2^-15 * 32728,
    // which is correct, but hurts readability

    // Note that std::numeric_limits<mpz_class> is specialized by gmpxx, so this case is caught.
    CPP17CONSTEXPRIF (std::numeric_limits<AbsoluteCooType>::is_integer)
    {
      exponent = std::max(0, exponent);
    }

    for (uint_fast16_t i = 0; i < dimension; ++i)
    {
      data[i] = ConvertMaybeMPZ<ApproxEntryType>::convert_to_inttype(
          EMVScalar::divide_by_power_of_2(exact_point.get_absolute_coo(i), exponent));
    }
  }
}

// actual scalar product

template <int nfixed>
inline auto compute_sc_product(EMVApproximation<nfixed> const &lhs,
                               EMVApproximation<nfixed> const &rhs) -> EMVScalar
{
  static_assert(
      std::is_same<EMVScalar, typename EMVApproximation<nfixed>::ScalarProductType>::value, "");

  using ProductType = typename EMVApproximationTraits::ApproxNorm2Type;

  ProductType scp = 0;
  auto const dim  = rhs.get_dim();
  for (uint_fast16_t i = 0; i < dim; ++i)
  {
    scp += lhs[i] * rhs[i];
  }
  return EMVScalar{lhs.exponent + rhs.exponent, scp};
}

// output

template <int nfixed>
inline std::ostream &operator<<(std::ostream &os,
                                EMVApproximation<nfixed> const &approximated_vector)
{
  os << "Approximation :2^" << approximated_vector.exponent << " x [";
  auto const dim = approximated_vector.get_dim();
  for (uint_fast16_t i = 0; i < dim; ++i)
  {
    os << approximated_vector.data[i] << " ";
  }
  os << "]";  // std::endl by caller;
  return os;
}

}  // end namespace GaussSieve

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM

#endif

// clang-format on
