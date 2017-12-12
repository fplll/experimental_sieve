#ifndef GAUSS_SIEVE_BITAPPROX_NEW_H
#define GAUSS_SIEVE_BITAPPROX_NEW_H


#include "DefaultIncludes.h"
#include "GlobalStaticData.h"

// clang-format off

namespace GaussSieve
{
namespace SimHashNew
{

// forward declarations:
template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
class BlockOrthogonalSimHash;
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

// Note args should be size_t, because that is what std::bitset and std::array expect
// (Otherwise, certain templates might not work)
// Note: The _arg suffix is used to be able to mirror the template argument in a local typedef.
template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg> // Nfixed?
class BlockOrthogonalSimHash
{
public:
  using IsCooSelection = std::true_type;
  static_assert(is_a_power_of_two_constexpr(sim_hash_len),"");
  static unsigned int constexpr num_of_transforms = 2;  // affects the "quality" vs. speed tradeoff
  using SimHashBlock  = std::bitset<sim_hash_len>;
  using SimHashes     = std::array<SimHashBlock,sim_hash_num>;
  using DimensionType = DimensionType_arg;
  static std::size_t constexpr get_sim_hash_num() { return sim_hash_num; }
  static std::size_t constexpr get_sim_hash_len() { return sim_hash_len; }

  BlockOrthogonalSimHash() = default;
  BlockOrthogonalSimHash(DimensionType const dim, unsigned int random_seed);
  BlockOrthogonalSimHash(DimensionType const dim)
      :BlockOrthogonalSimHash(dim, std::random_device{}()) {}

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

  unsigned int number_of_orthogonal_blocks;  // TODO: rename
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
