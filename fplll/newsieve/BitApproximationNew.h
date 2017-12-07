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

// forward declarations:
template<unsigned int sim_hash_len_arg, unsigned int sim_hash_num_arg, bool MT, class DimensionType>
class CoordinateSelection;
template<class CooSelection> class ObtainSimHashBlock;

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

template<unsigned int sim_hash_len_arg, unsigned int sim_hash_num_arg, bool MT, class DimensionType> // Nfixed?
class CoordinateSelection
{
public:
  static unsigned int constexpr num_of_transforms = 2;  // affects the "quality" vs. speed tradeoff
  static unsigned int constexpr sim_hash_num = sim_hash_num_arg;
  static unsigned int constexpr sim_hash_len = sim_hash_len_arg;
  static_assert(is_a_power_of_two(sim_hash_len),"");
  using SimHashBlock  = std::bitset<sim_hash_len>;
  using SimHashes     = std::array<SimHashBlock,sim_hash_num>;
  using Get_DimensionType = DimensionType;

  static_assert(MT==false,""); using SimHashesForList = SimHashes; // TODO!

  CoordinateSelection() = default;
  CoordinateSelection(DimensionType const dim, unsigned int random_seed);
  CoordinateSelection(DimensionType const dim)
      :CoordinateSelection(dim, std::random_device{}()) {}

  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  inline SimHashes compute_all_bitapproximations(LatP const &point) const;

  template<class LHS, class RHS, class LowerThresholds, class UpperThresholds>
  FORCE_INLINE static inline bool check_simhash_scalar_product(LHS const &lhs, RHS const &rhs, LowerThresholds const &lb, UpperThresholds const &ub)
  {
    uint_fast32_t approx_scprod = 0;
    for (unsigned int level = 0; level < sim_hash_num; ++level)
    {
      approx_scprod += (ObtainSimHashBlock<CoordinateSelection>::get(lhs) ^ ObtainSimHashBlock<CoordinateSelection>::get(rhs)).count();
      if (approx_scprod >= ub[level] || approx_scprod <= lb[level])
      {
        continue;
      }
      else
      {
        return false;
      }
    }
    return true;
  }

private:
  template<class T>
  auto inline fast_partial_walsh_hadamard(std::vector<T> input) const
    -> std::vector<T>;

  unsigned int number_of_blocks;  // TODO: rename
  unsigned int fast_walsh_hadamard_len;
  std::vector< std::array<PMatrix,num_of_transforms> > pmatrices;
  std::vector< std::array<DMatrix,num_of_transforms> > dmatrices;

  //static std::array<RMatrix,SimHash::num_of_levels> rmatrices;
};

/**
  This class stores the result of computing a scalar product of bitwise
  approximations.
  Essentially, just wraps around an int.
  Note: We might want to wrap an approximate scalar product of t as value = #bits - t. -- Gotti
*/

/*
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
*/

/*
template <class SieveTraits, bool MT>
[[gnu::always_inline]] inline constexpr BitApproxScalarProduct compute_simhash_scalar_product_block(
    SimHashBlock<SieveTraits,MT> const &lhs, SimHashBlock<SieveTraits,MT> const &rhs)
{
  return BitApproxScalarProduct {static_cast<uint_fast32_t>(SieveTraits::sim_hash_len - (lhs ^ rhs).count()) };
}
*/


}  // end namespace (GaussSieve::)SimHash
}  // end namespace GaussSieve

namespace GaussSieve
{

template<class CooSelection> class ObtainSimHashBlock
{
  template<class Arg>
  FORCE_INLINE static typename CooSelection::SimHashBlock const & get(Arg const &arg, unsigned int const level)
  {
    return arg.access_bitapproximation(level);
  }
  FORCE_INLINE static typename CooSelection::SimHashBlock const & get(typename CooSelection::SimHashes const &arg, unsigned int const level)
  {
    return arg[level];
  }
};


template<class LatticePoint, class CooSelect>
struct LPWithBitapprox
{
public:
  static_assert(IsALatticePoint<LatticePoint>::value, "Wrong template argument");

  using SimHashes = typename CooSelect::SimHashes;
  using SimHashBlock = typename CooSelect::SimHashBlock;

  LatticePoint latp;
  SimHashes sim_hashes;
  explicit LPWithBitapprox(LatticePoint &&new_latp, SimHashes const &new_sim_hashes) noexcept
      : latp(std::move(new_latp)), sim_hashes(new_sim_hashes) {}
  explicit LPWithBitapprox(LatticePoint &&new_latp, SimHashes &&new_sim_hashes) noexcept
      : latp(std::move(new_latp)), sim_hashes(std::move(new_sim_hashes)) {}
  explicit LPWithBitapprox(LatticePoint &&new_latp, CooSelect const &coo_select) noexcept
      : latp(std::move(new_latp)), sim_hashes(coo_select.compute_all_bitapproximations(latp)) {}
  constexpr operator LatticePoint() const & { return latp; }
            operator LatticePoint() &&      { return std::move(latp); }
  constexpr operator SimHashes()            { return sim_hashes; }
  void update_sim_hashes(CooSelect const &coo_select) noexcept
  {
    sim_hashes = coo_select.compute_all_bitapproximations(latp);
  }
  SimHashBlock const & access_bitapproximation(unsigned int level) const
  {
    return sim_hashes[level];
  }
};


}  // end namespace GaussSieve

#include "BitApproximationNew_impl.h"

#endif // include guards

// clang-format on
