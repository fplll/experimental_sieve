#ifndef GAUSS_SIEVE_BITAPPROX_NEW_H
#define GAUSS_SIEVE_BITAPPROX_NEW_H

// TODO: Rename file

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"

/**
  This files defines the class BlockOrthogonalSimHash, which is a "CoordinateSelection" in the sense
  explained in SimHash.h
  Roughly, such a CoordinateSelection stores the parameters for and allows the application of a hash
  function H that maps lattice points to sequences of bits (stored as a vector of blocks of bits).

  The hash function family (approximately) defined by BlockOrthogonalSimHash is given as follows:
  (NOTE:  This just defines a hash function family; the actual way we compute things differs greatly
          from that description for efficiency reasons, see below)
  Choose n * m unit vectors v_i, organized in m blocks of othonormal vectors, i.e. within each block
  the vectors are pairwise orthogonal. We choose n equal to the dimension we operate in, so each
  block is a basis. Then the i-th bit of H(p) is given by the sign of <p,v_i>.
  The v_i's constitute the parameters of the hash function.

  Note that the output is given by a vector of blocks of bits, but this does not correspond to the
  decomposition of the v_i's into orthonormal bases; rather, we truncate the last block of v_i's.
  The selection of v_i's defines a concrete function within the family.

  Note: We want to choose the v_i's "evenly spread" out on the sphere, we want to actually find such
        a good selection of v_i's quickly and we want to compute the <p,v_i>'s quickly. The choices
        taken are a comprise between these goals. In particular, block-wise orthogonal v_i's are
        better than uniformly random v_i's and allows the efficiency-improvement explained below.

  There are various ways to uniformly sample orthonormal bases.
  For efficiency reasons, we do the following (implicitly): We take the standard basis and apply a
  sufficiently random uniform orthogonal transformation A.
  To get such a "sufficiently random" uniform transformation, we define A as
  A^(-1) = P_0 *D_0 * WH * P_1 * D_1 * WH * P_2 * D_2 * WH * ... * D_{num_of_transforms},
  where P_i are permutation matrices, D_i are diagonal matrices with +/-1's on the diagonal and
  WH is the Walsh-Hadamard transform. Note that, due to their special shape, applying such a matrix
  on a vector is much faster than a generic matrix-vector multiplication.
  Note: For efficiency reasons, we only perform WH in a subset of the dimensions that is a power of
        two; the primary purpose of the permuations is to ensure all coos get eventually affected.

  Now, instead of computing <p, v_i> = <p, A(e_i)>, where e_i are the standard basis, we compute the
  equivalent <A^(-1)p, e_i> for each i. This means we do not store the v_i's. Rather, we store the
  sequence of P_i, D_i's. To compute a given block (wrt to the orthogonal decomposition) of the hash
  function output H(p), we compute A(p) by applying P_i's, D_i's and WH on p and then just look at
  the sign bits of A(p). More precisely, to account for differing decompositions, we consider a
  concatenation p||p||p||... of sufficient length, and apply such an (independently random)
  A-transform on each block. The resuting output is then decomposed into blocks according to the
  decomposition required for the output.

  Note that this is actually *faster* than a uniformly random selection of v_i's.
  While asymptotically, the computation of such hashes does not matter, the time spent is actually
  a significant part of the running time (dominating, if done naively!) for small dims <= 60...
  (Then again, SimHashes do not give any asymptotic benefit either in the way we use them)

  The number num_of_transforms parameterizes a tradeoff between the randomness of the output and
  computational complexity. We choose num_of_transforms small.

  TODO: The mixing of coordinates might be improved by a non-random permutation.
*/

namespace GaussSieve
{

// forward declarations:
template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
class BlockOrthogonalSimHash;
template<class CooSelection> class ObtainSimHashBlock;

/**
 PMatrix creates a random permutation matrix.
 Stored as a permutation of the sequence of coordinates, not as a matrix.
 Can be apply()'d to a vector.
 */
class PMatrix
{
private:
  std::vector<uint_fast16_t> permutation;  // result[i] = input[permutation[i]]
public:
  PMatrix() = default;
  // creates a random dim-dimensional permutation using rng as a source of randomness.
  PMatrix(unsigned int dim, std::mt19937 &rng);
  // apply the stored permutation to the given vector. Changes the argument.
  template<class T> inline void apply(std::vector<T> &vec) const;
  // prints the permutation (as a sequence of permuation[i], not as a matrix)
  inline void print(std::ostream &os = std::cout) const;
};

/**
  Creates a random diagonal matrix with +/-1 on the diagonal.
  We only store the diagonal, of course.
  Can be apply()'ed to a vector.
*/

class DMatrix
{
private:
  // 0 means no sign-flip, 1 means sign-flip,
  // i.e. the actual matrix entry is Matrix[i,i] = (-1)^{diagonal[i])
  std::vector<uint_fast8_t> diagonal;
public:
  DMatrix() = default;
  // creates a random such DMatrix, using rng as randomness source.
  DMatrix(unsigned int const dim, std::mt19937 &rng);
  template<class T> inline void apply(std::vector<T> &vec) const;
  // prints diagonal (as a sequence of 0's and 1's, not as +/-1's)
  inline void print(std::ostream &os = std::cout) const;
};

/**
  Actual BlockOrthogonalSimHash: Stores parameters for a sim_hash function from the family of
  block-orthogonal matrices and allows applying it to a point.
  template parameters sim_hash_len and sim_hash_num determine the requested number of output bits
  as sim_hash_len * sim_hash_len, output in blocks of sim_hash_len bits.
  MT denotes whether we are multi-threaded. This parameter is ignored, computation of a hash is
  always thread-safe anyway. (Parameter initialization is not, but does not have to be)
  DimensionType is a type capable of holding the dimension of the input vectors.
  Can be nfixed<dim> to hard-wire the dimension (this might give potential improvements for the
  computation of the WH trafo due to better loop unrolling)

  Note: The time spent on computing the hashes is visible for the 2-sieve in dims up to at least 60,
        so we care about efficiency.
*/

// Note args should be size_t, because that is what std::bitset and std::array expect
// (Otherwise, certain templates might not work)
// Note: The _arg suffix is used to be able to mirror the template argument in a local typedef.
template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
class BlockOrthogonalSimHash
{
public:
  using IsCooSelection = std::true_type;  // Satisfies the CoordinateSelection concept
  static unsigned int constexpr num_of_transforms = 2;  // affects the "quality" vs. speed tradeoff

  // required to tell users the type of outputs, required by the interface specified in
  // SimHash.h
  using SimHashBlock  = std::bitset<sim_hash_len>;
  using SimHashes     = std::array<SimHashBlock,sim_hash_num>;
  using DimensionType = DimensionType_arg;  // convertible to unsigned int
  static std::size_t constexpr get_sim_hash_num() { return sim_hash_num; }
  static std::size_t constexpr get_sim_hash_len() { return sim_hash_len; }

  BlockOrthogonalSimHash() = default;
  // dim is the ambient dimension of the points on which we call compute_all_bitapproximations
  // randomly samples parameters for a hash functions fit on points of that dimension, using
  // random_seed as (optional) input randomness. If the user does not provide it, use random_device
  BlockOrthogonalSimHash(DimensionType const dim, unsigned int random_seed);
  BlockOrthogonalSimHash(DimensionType const dim)
      :BlockOrthogonalSimHash(dim, std::random_device{}()) {}

  // computes the SimHashes of a LatticePoint point
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  inline SimHashes compute_all_bitapproximations(LatP const &point) const;

private:
  template<class T>
  [[deprecated]] // use faster_partial_walsh_hadamard instead, which performs inplace operation.
  auto inline fast_partial_walsh_hadamard(std::vector<T> input) const
    -> std::vector<T>;

  // performs WH-Trafo (well, almost: up to a fixed permutation / sign-flips) in
  // the first fast_walsh_hadamard_len coordinates. Modifies the argument (hence inplace).
  // This version is an order of magnitude faster than the above.
  // Note: Allowing a fixed permutation / sign-flips allows some micro-optimizations.
  template<class T>
  void inline faster_almost_partial_walsh_hadamard_inplace(std::vector<T> & input) const;

  // to compute a bitapproximation on a vector p, we duplicate p number_of_orthogonal_blocks many
  // times and operate on each block individually.
  unsigned int number_of_orthogonal_blocks;
  // length in which we perform WH trafos. Currently, the largest power of 2 <= dim.
  unsigned int fast_walsh_hadamard_len;
  unsigned int fast_walsh_hadamard_loglen;  // log_2 of the above.
  // number_of_orthogonal_blocks x num_of_transforms may pmatrices / dmatrices.
  std::vector< std::array<PMatrix,num_of_transforms> > pmatrices;  // permutations
  std::vector< std::array<DMatrix,num_of_transforms> > dmatrices;  // diagonal matrices
};

}  // end namespace GaussSieve


#include "BlockOrthogonalSimHash_impl.h"

#endif // include guards
