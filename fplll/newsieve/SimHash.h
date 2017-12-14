#ifndef GLOBAL_BITAPPROX_DATA_H
#define GLOBAL_BITAPPROX_DATA_H

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"

/******
  This file defines "general" aspect of SimHashes (i.e. bit-vectors that are used as a "hash" of a
  point to approximately tell whether a vector is a good candiate for a reduction)

  - This file specifies the CoordinateSelection concept, which encapsulates the (usually randomized)
    *parameters* of the concrete SimHash in use.

  - This file defines the GlobalBitApproxData template, which holds the actual (global, static)
    parameters we are using. Note that those parameters (should) act essentially as a Singleton.

  - This file defines some functions to operate on SimHashes.

  Notably, such a SimHash can be applied to a lattice point p, resulting in a sequence of bits H(p).
  For 2 lattice points p, q on the sphere (Note that H is invariant under scaling), we have that
  the number of coinciding bits of their hashes is a rough measure for their distance, i.e.
  H(p) XOR H(q) is correlated with the (angular) distance between p,q.
  We use that to determine whether points are promising candidates for a reduction.

  Note: The difference between this usage (SimHash) and locality sensitive hashing / filtering
        techniques is that for the latter, we would organize our storage according to the value of
        the sim_hash, in order to directly obtain all points with a given (set of restrictions on
        the) sim_hash without looking at other points. This requires more elaborate data structure
        with a larger overhead. By constrast, for SimHashes we store the hashes along the points
        and *test* each candidate.

******/

// TODO: Unify naming "bit_approximation", "bitapproximation", "bitapprox", "sim_hashes", both in
// documentation and naming of classes / functions
// Note : These are the same, really!

/**
  CONCEPT: CoordinateSelection

  For our use of SimHashes, we need to define a class that collects the (usually randomized)
  parameters of the hash functions(s).

  Such a class CoordinateSelection is called a Coordinate Selection (CooSelection).
  For a typical SimHash, we choose #output_bits many directions v_i and set the i-th bit of H(i) as
  the sign of <v_i, p>. The parameters of the hash functions then correspond to a selection of v_i's
  This is why we chose that name.
  Note that a given class that satisfies this concept might choose the v_i's in a structured, maybe
  implicit, way and does not need to store such v_i's explicitly (in fact, our choices don't).

  A CoordinateSelection is a class (concept) that has the following interface:
  - It has a public member typedef IsCooSelection = std::true_type
    (Check with IsACoordinateSelection<T>)
  - It has public member typedefs SimHashBlock and SimHashes
    SimHashBlock needs to be a std::bitset<block_len> or compatible,
    SimHashes    needs to be a random access container of SimHashBlocks
    (Note: "or compatible" means we might support e.g. boost/dynamic_bitset
           random access container means access via sim_hash[index],
           i.e. a std::array, std::vector or C-style array of SimHashBlocks)

    Note: We usually process a given SimHash block by block (in particular, when using a SimHash to
          detect whether a lattice vector is promising, we may want to terminate early and not even
          look at all bits of the SimHash to save time. This is why we consider a block structure,
          where the block_len is to be optimized for the target architecture, typically 64 bits.
  - public member functions get_sim_hash_len(), get_sim_hash_num() that return the number of bits
    per block (len) respectively the number of blocks (num).
    (Note: These are typically constexpr static.)
  - It has a public member typedef DimensionType,
    which is Default-constructible, copy-assignable and convertible to / from integral types.
  - CoordinateSelection is default-constructible and copy-assignable.
    (Note: it need not be usable after default-construction, this requirement is purely syntactic)
  - It has constructors CoordinateSelection(DimensionType const dim, unsigned int seed)
                    and CoordinateSelection(DimensionType const dim)
    If seed it set, all functions (including those constructors) should be deterministic.
    (Note: DimensionType typically is either some MaybeFixed<nfixed> or some unsigned int class)
  - it has a public member template function
    SimHashes compute_all_bitapproximations(LatticePoint const &p) const
    After a CoordinateSelection object was created by one of the two mentioned constructors above,
    this deterministically computes all sim_hashes of the point p. LatticePoint may be any lattice
    point class (see LatticePointConcept.h for what is a "LatticePoint" in this context).
*/

/**
  Usage inside the GaussSieve:

  GlobalBitApproxData<CoordinateSelection> holds a static object of type CoordinateSelection.

  The SieveTraits class that is passed as a parameter to most modules of the GaussSieve selects a
  concrete CoordinateSelection class that is used. Classes then can use
  GlobalBitApproxData<concrete CoordinateSelection>::coo_selection to access the global static.
  This global static is initialized using the CoordinateSelection(dim, seed) - constructor by
  the main sieve.

  An object of type StaticInitializer<GlobalBitApproxData<CoordinateSelection>> defined below is
  used for this. Note that this ensure that you do not (re-)initialize CoordinateSelection with
  different paramters as long as the old StaticInitializer<GlobalBitApproxData<CoordinateSelection>>
  is alive. In particular, to change the coordinate selection, you have to destroy all objects of
  the GaussSieve class beforehand.
  (See GlobalStaticData.h for details on StaticInitializer's)
*/

/**
  Currently, we support and use the BlockwiseOrthogonalSimHash class as the selected
  CoordinateSelection, defined in BlockOrthogonalSimHash.h
*/

namespace GaussSieve
{

// To ensure validity of template arguments, in order to give meaningful compiler errors.
template <class T> using Predicate_IsCooSelection = mystd::enable_if_t<  std::is_same<typename T::IsCooSelection, std::true_type>::value  >;
template <class T> using IsACoordinateSelection   = mystd::is_detected<Predicate_IsCooSelection,T>;

/**
  This is a class that encapsulates (as a static member) a singleton of type CooSelection.
  Various (otherwise loosely related) have to use the SAME hash function, hence we want a global.
*/
template<class CooSelection>
class GlobalBitApproxData
{
  static_assert(IsACoordinateSelection<CooSelection>::value, "Invalid template argument");
  friend StaticInitializer<GlobalBitApproxData>;

private:
  // we store the data used to construct coo_selection during its initialization.
  // This is done to detect re-initialization attempts with different dim / seed while we still
  // use the current coo_selection.
  static typename CooSelection::DimensionType dim_used;
  static unsigned int random_seed_used;

public:
  static CooSelection coo_selection;  // the actual global object.
};

// Note: Recall that static member variables of class *templates* are automatically inline
//       (pre C++-17, this is the only way to get inline variables)

// compile-time initialization of coo_selection. This will be overwritten at runtime before its
// first use by the StaticInitializer. Note that our StaticInitializer keep a "WasInitialized" bit
// around, so we do not need one here.
template<class CooSelection> CooSelection                         GlobalBitApproxData<CooSelection>::coo_selection;
template<class CooSelection> typename CooSelection::DimensionType GlobalBitApproxData<CooSelection>::dim_used;
template<class CooSelection> unsigned int                         GlobalBitApproxData<CooSelection>::random_seed_used;

/**
  Static Initializer for run-time initialization (RAII-style wrapper).
  An object of this type is created by the GaussSieve very early during its initialization process.
  (cf. GlobalStaticData.h)
 */

template<class CooSelection> class StaticInitializer< GlobalBitApproxData<CooSelection> > final
    : public DefaultStaticInitializer< GlobalBitApproxData<CooSelection> >
{
  static_assert(IsACoordinateSelection<CooSelection>::value, "Invalid template argument");
  using Parent        = DefaultStaticInitializer< GlobalBitApproxData<CooSelection> >;
  using DimensionType = typename CooSelection::DimensionType;

public:
  // See GlobalStaticData.h for the meaning of IsArgForStaticInitializer.
  // Essentially, we allow certain types that have a .dim member.
  template<class T, TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const &initializer, unsigned int random_seed)
      : StaticInitializer(static_cast<DimensionType>(initializer.dim), random_seed) {}

  StaticInitializer(DimensionType const new_dim, unsigned int random_seed)
  {
    assert(Parent::user_count > 0);  // Parent::user_count includes this object itself as well.
    if (Parent::user_count > 1)
    {
      // ensure we not reinitialize as long as there are previous users around.
      // TODO: Throw exception rather than assert!
      assert(new_dim == GlobalBitApproxData<CooSelection>::dim_used);
      assert(random_seed == GlobalBitApproxData<CooSelection>::random_seed_used);
    }
    else  // we are the first / only user of GlobalBitApproxData<CooSelection>
    {
      GlobalBitApproxData<CooSelection>::dim_used         = new_dim;
      GlobalBitApproxData<CooSelection>::random_seed_used = random_seed;
      GlobalBitApproxData<CooSelection>::coo_selection    = CooSelection{new_dim, random_seed};
    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing Global data for Bitapproximations."
                                      << " Counter is" << Parent::user_count)
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing Global data for Bitapproximations."
                                      << " Counter is " << Parent::user_count)
  }
};

// TODO: Currently, we only support std::arrays of std::bitset's, not different containers.

/**
 Outputs a SimHash (i.e. an array of bitsets)
*/
template<std::size_t sim_hash_num, std::size_t sim_hash_len>
std::ostream & operator<<(std::ostream                                               &os,
                          std::array< std::bitset<sim_hash_len>, sim_hash_num> const &sim_hashes)
{
  for (unsigned int i = 0; i < sim_hash_num; ++i)
  {
    os << sim_hashes[i] << " ";
  }
  return os;
}

/**
  takes a SimHash and returns a copy of it, with all bits flipped.
*/
template<std::size_t sim_hash_num, std::size_t sim_hash_len>
auto flip_all_bits(std::array< std::bitset<sim_hash_len>, sim_hash_num > const &sim_hashes)
    -> std::array< std::bitset<sim_hash_len>, sim_hash_num >
{
  auto sim_hash_ret = sim_hashes;
  for (size_t i = 0; i < sim_hash_num; ++i)
  {
    sim_hash_ret[i] = ~sim_hashes[i];
  }
  return sim_hash_ret;
}

/**
  Helper for the function below.
  ObtainSimHashBlock<CooSelection>get(arg, level)
  is equivalent to arg[level] if arg *is* a SimHash, but also works for certain classes that
  instead *contain* a SimHash, if they provide an access_bitapproximation(unsigned int) - member.
  (Some Lattice point classes and our list iterators do that.) Used to unify code.
*/
namespace Helpers
{
template<class CooSelection> struct ObtainSimHashBlock
{
  static_assert(IsACoordinateSelection<CooSelection>::value, "Wrong template argument");

  template<class Arg>
  FORCE_INLINE static auto get(Arg const &arg, unsigned int const level)
      -> typename CooSelection::SimHashBlock const &
  {
    return arg.access_bitapproximation(level);
  }

  FORCE_INLINE static auto get(typename CooSelection::SimHashes const &arg, unsigned int const level)
      -> typename CooSelection::SimHashBlock const &
  {
    return arg[level];
  }
};
}  // end namespace Helpers

/**
  Checks whether popcount(lhs XOR rhs) is >= ub or <= lb (cumulatively) for each level.
  LHS / RHS may be either a SimHash or something having an access_bitapproximation function
  (i.e. currently either a SimHash, a list iterator or a lattice point with BitApproximation)

  Cumulative means the following:

  We first check lhs[0] XOR rhs[0] against ub[0], lb[0].
  On the next level, we check (lhs[0] XOR rhs[0]) + (lhs[1] XOR rhs[1]) against ub[1] and lb[1] etc.
  We return true iff for each level : (we are either >= ub or <= lb).

  Note: a "scalar product" corresponds rougly to 1 - 2*(lhs XOR rhs)/len, in the sense that
        for points on the unit sphere, the expression 1-2*(lhs XOR rhs)/len is close to 1 if the
        points are close, close to -1 if the points are near-antipodes and close to 0 if the points
        are near-orthogonal.
        In particular, lhs XOR rhs itself becomes *smaller* if the scalar product gets *larger*, so
        the lower thresholds on the XOR means an upper threshold on the scalar product.


*/
template<class CoordinateSelection, class LHS, class RHS,
	 class LowerThresholds, class UpperThresholds>
FORCE_INLINE static inline bool CPP14CONSTEXPR check_simhash_scalar_product(
    LHS const &lhs, RHS const &rhs, LowerThresholds const &lb, UpperThresholds const &ub)
{
  uint_fast16_t approx_scprod = 0;  // holds accumulated XOR - value
  for (unsigned int level = 0;
       level < GlobalBitApproxData<CoordinateSelection>::coo_selection.get_sim_hash_num(); ++level)
  {
    approx_scprod += ( Helpers::ObtainSimHashBlock<CoordinateSelection>::get(lhs,level)
                      ^Helpers::ObtainSimHashBlock<CoordinateSelection>::get(rhs,level) ).count();
    if (approx_scprod >= ub[level] || approx_scprod <= lb[level])
    {
      continue;  // We are outside of the bounds, so we go to the next iteration of for loop
    }
    else  // otherwise, we are inside the bounds and already know we return false.
    {
      return false;
    }
  }
  return true;
}

/*
  same as above but also returns (in case of true) which of the if-conditions was satisfied
  I.e., tell whether the two points are close or far apart
  If (approx_scprod <= lb[level]), bool is_close is set to true
 */
template<class CoordinateSelection, class LHS, class RHS,
         class LowerThresholds, class UpperThresholds>
FORCE_INLINE static inline bool CPP14CONSTEXPR check_simhash_scalar_product_ext(
                                                        LHS const &lhs, RHS const &rhs, LowerThresholds const &lb, UpperThresholds const &ub, bool &is_close)
{
  uint_fast16_t approx_scprod = 0;  // holds accumulated XOR - value
  for (unsigned int level = 0;
       level < GlobalBitApproxData<CoordinateSelection>::coo_selection.get_sim_hash_num(); ++level)
  {
    approx_scprod += ( Helpers::ObtainSimHashBlock<CoordinateSelection>::get(lhs,level)
                      ^Helpers::ObtainSimHashBlock<CoordinateSelection>::get(rhs,level) ).count();
    if (approx_scprod >= ub[level])
    {
      is_close = false;
      continue;  // We are outside of the bounds, so we go to the next iteration of for loop
    }
    else if (approx_scprod <= lb[level])
    {
      is_close = true;
      continue;
    }
    else  // otherwise, we are inside the bounds and already know we return false.
    {
      return false;
    }
  }
  return true;
}

}  // end namespace GaussSieve

#endif
