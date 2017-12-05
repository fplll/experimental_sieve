#ifndef GAUSS_SIEVE_BITAPPROX_NEW_H
#define GAUSS_SIEVE_BITAPPROX_NEW_H

#include <bitset>
#include <boost/dynamic_bitset.hpp> // maybe remove at some point.
#include <random>
#include <math.h> //for log2, sqrt

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"

// clang-format off

namespace GaussSieve
{
namespace SimHashNew
{

/*****************************************************************************
is_a_power_of_two(n) checks whether n is a power of 2.
Slow, but designed for use in static asserts
(This is why it's a recursive one-line function, to be C++11 - constexpr)
******************************************************************************/

template<class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
constexpr bool is_a_power_of_two(Integer const n)
{
  // one-line function to be C++11 - constexpr. Slow, but only used in static_asserts anyway.
  return (n > 0) && ( (n == 1) || ( (n%2 == 0) && is_a_power_of_two(n/2) ) );
}

template<class SieveTraits, bool MT> class CoordinateSelection;  // forward declare
template<class SieveTraits, bool MT> using SimHashBlock = std::bitset<SieveTraits::sim_hash_len>;


/**
 PMatrix creates and stores a random permutation (matrix) in dim dimensions.
 Can be apply()'d to a vector.
 */
class PMatrix
{
private:
  std::vector<uint_fast16_t> permutation;
public:
  PMatrix() = default;
  PMatrix(unsigned int dim, std::mt19937 &rng);  // initialize with a random dim-dimensional permutation.
  // apply the stored permutation to the given vector. Changes the argument.
  template<class T> inline void apply(std::vector<T> &vec) const;
  inline void print(std::ostream &os = std::cout) const;
};

/**
  Creates and stores a random diagonal matrix with +/-1 on the diagonal.
  Can be apply()'ed to a vector.
*/

class DMatrix
{
private:
  std::vector<uint_fast8_t> diagonal;  // 0 means no sign-flip, 1 means sign-flip:
public:
  DMatrix() = default;
  DMatrix(unsigned int const dim, std::mt19937 &rng);
  template<class T> inline void apply(std::vector<T> &vec) const;
  inline void print(std::ostream &os = std::cout) const;
};

/********************************************************************
SimHash::CoordinateSelection<SieveTraits,IsMultithreaded>
is a class template holding the data related to bitapproximation.

Notably, it determines the random directions we are using for SimHash.
This is essentially encoded in a set of orthogonal matrices (up to scaling), each of which is
(for efficiency reason) defined and stored as a sequence of special matrices (permutations,
diagonal matrices and the Walsh-Hadamard matrix) that can be applied to a vector fast.

Note: We may template this by the user lattice point class as well or
store these data inside some other class which uses the data.
Initializiation is performed by StaticInitializer<CoordinateSelection>.
This is called by ExactLatticePoint at the moment, but this is subject
to change.

Note2: Template arguments are subject to change.

Note3: Currently, we only use CoordinateSelection<void,false>
       This will change once control over bitapproximations is moved into the main list.

TODO: Consider writing a wrapper around either a std::bitset or a std::dynamic_bitset
      that actually has these as static data and initialize it via the main sieve.

TODO: Update documentation to reflect refactoring due to inclusion into main list
********************************************************************/

template<class SieveTraits, bool MT>
class CoordinateSelection
{
public:
  static unsigned int constexpr num_of_transforms = 2;  // affects the "quality" vs. speed tradeoff
  static unsigned int constexpr &sim_hash_num = SieveTraits::sim_hash_num;
  static unsigned int constexpr &sim_hash_len = SieveTraits::sim_hash_len;
  CoordinateSelection() = delete;  // cannot be instantiated.
  CoordinateSelection(typename SieveTraits::DimensionType const dim, unsigned int random_seed);
  CoordinateSelection(typename SieveTraits::DimensionType const dim)
      :CoordinateSelection(dim, std::random_device{}()) {}

  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  inline auto transform_and_bitapprox(LatP const &point)
      -> std::array< SimHashBlock<SieveTraits,MT>, sim_hash_num >;

  //template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  //inline static auto transform_and_bitapprox_simple(LatP const &point)
  //    -> std::array< std::bitset<sim_hash_len>, num_of_levels >;

  /*
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  inline static auto transform_and_bitapprox_2nd_layer(LatP const &point)
      -> std::array< SimHashBlock<SieveTraits,MT>, num_of_levels >;
  */
private:
  template<class T>
  auto inline fast_partial_walsh_hadamard(std::vector<T> input)
    -> std::vector<T>;
  unsigned int number_of_blocks;
  unsigned int fast_walsh_hadamard_len;
  std::vector< std::array<PMatrix,num_of_transforms> > pmatrices;
  std::vector< std::array<DMatrix,num_of_transforms> > dmatrices;
  // unsigned int ambient_dimension;  // TODO: Remove this in favor of a static scratch space

  //static std::array<RMatrix,SimHash::num_of_levels> rmatrices;
};

/**
  This class stores the result of computing a scalar product of bitwise
  approximations.
  Essentially, just wraps around an int.
  Note: We might want to wrap an approximate scalar product of t as value = #bits - t. -- Gotti
*/

class BitApproxScalarProduct
{
public:
  using BitApproxScalarProduct_WrappedType = uint_fast32_t;


  BitApproxScalarProduct(BitApproxScalarProduct const &) = delete;  // why not copy ints? -- Gotti
  BitApproxScalarProduct(BitApproxScalarProduct &&)      = default;

  explicit constexpr BitApproxScalarProduct(BitApproxScalarProduct_WrappedType const rhs):value(rhs) {}
  explicit operator  BitApproxScalarProduct_WrappedType() { return value; }

  BitApproxScalarProduct &operator=(BitApproxScalarProduct const &other) = delete;  // Why? -- Gotti
  BitApproxScalarProduct &operator=(BitApproxScalarProduct &&other) = default;

  // TODO: operator >=, <=

  inline bool operator<=(BitApproxScalarProduct_WrappedType &&rhs) const
  {
    return (this->value <= rhs);
  }

  BitApproxScalarProduct_WrappedType value;
};

template <class SieveTraits, bool MT>
[[gnu::always_inline]] inline constexpr BitApproxScalarProduct compute_simhash_scalar_product_block(
    SimHashBlock<SieveTraits,MT> const &lhs, SimHashBlock<SieveTraits,MT> const &rhs)
{
  return BitApproxScalarProduct {static_cast<uint_fast32_t>(SieveTraits::sim_hash_len - (lhs ^ rhs).count()) };
}

}  // end namespace (GaussSieve::)SimHash
}  // end namespace GaussSieve

#include "BitApproximationNew_impl.h"

#endif // include guards

// clang-format on
