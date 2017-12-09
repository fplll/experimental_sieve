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
template<std::size_t sim_hash_len_arg, std::size_t sim_hash_num_arg, bool MT, class DimensionType>
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

// Note args should be size_t, because that is what std::bitset and std::array expect
// (Otherwise, certain templates might not work)
template<std::size_t sim_hash_len_arg, std::size_t sim_hash_num_arg, bool MT, class DimensionType> // Nfixed?
class CoordinateSelection
{
public:
  static unsigned int constexpr num_of_transforms = 2;  // affects the "quality" vs. speed tradeoff
  static std::size_t constexpr sim_hash_num = sim_hash_num_arg;
  static std::size_t constexpr sim_hash_len = sim_hash_len_arg;
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

  // LHS, RHS: Either SimHashBlock or an iterator or a lattice point that contains a bitapproximation
  // TODO: static_assert those conditions.

private:
  template<class T>
  [[deprecated]] // use faster_partial_walsh_hadamard, which performs inplace operation.
  auto inline fast_partial_walsh_hadamard(std::vector<T> input) const
    -> std::vector<T>;

  // performs WH-Trafo (well, almost: up to a fixed permutation / sign-flips) in
  // the first fast_walsh_hadamard_len coordinates. Modifies the argument (hence inplace).
  // This version is an order of magnitude faster than the above.
  template<class T>
  void inline faster_almost_partial_walsh_hadamard_inplace(std::vector<T> & input) const;

  unsigned int number_of_blocks;  // TODO: rename
  unsigned int fast_walsh_hadamard_len;
  unsigned int fast_walsh_hadamard_loglen;
  std::vector< std::array<PMatrix,num_of_transforms> > pmatrices;
  std::vector< std::array<DMatrix,num_of_transforms> > dmatrices;
  //static std::array<RMatrix,SimHash::num_of_levels> rmatrices;
};

}  // end namespace (GaussSieve::)SimHash
}  // end namespace GaussSieve


#include "BitApproximationNew_impl.h"

#endif // include guards

// clang-format on
