#ifndef GLOBAL_BITAPPROX_DATA_H
#define GLOBAL_BITAPPROX_DATA_H

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"
#include <bitset>
#include <array>

namespace GaussSieve
{

// This is a class that encapsulates (as a static member) a singleton
// of type CooSelection.
template<class CooSelection>
class GlobalBitApproxData
{
  friend StaticInitializer<GlobalBitApproxData>;
private:
  static typename CooSelection::Get_DimensionType dim_used;
  static unsigned int random_seed_used;
public:
  static CooSelection coo_selection;
};

// compile-time initialization of coo_selection;
template<class CooSelection>
CooSelection GlobalBitApproxData<CooSelection>::coo_selection;
template<class CooSelection>
typename CooSelection::Get_DimensionType GlobalBitApproxData<CooSelection>::dim_used;
template<class CooSelection>
unsigned int GlobalBitApproxData<CooSelection>::random_seed_used;

/**
 Outputs a SimHash (i.e. an array of bitsets)
*/
template<std::size_t sim_hash_num, std::size_t sim_hash_len>
std::ostream& operator<<(std::ostream &os, std::array< std::bitset<sim_hash_len>, sim_hash_num> const &sim_hashes)
{
  for(unsigned int i=0; i < sim_hash_num; ++i)
  {
    os << sim_hashes[i] << " ";
  }
  return os;
}

template<std::size_t sim_hash_num, std::size_t sim_hash_len>
auto flip_all_bits( std::array< std::bitset<sim_hash_len>, sim_hash_num > const &sim_hashes)
    -> std::array< std::bitset<sim_hash_len>, sim_hash_num >
{
  auto sim_hash_ret = sim_hashes;
  for(size_t i = 0; i < sim_hash_num; ++i )
  {
    sim_hash_ret[i] = ~sim_hashes[i];
  }
  return sim_hash_ret;
}

/**
  Helper for the function below
*/
// TODO: Encapsulation
template<class CooSelection> struct ObtainSimHashBlock
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


/**
  Checks whether popcount(lhs XOR rhs) is > ub or < lb for each level.
  lhs / rhs may be either a SimHash or somthing having an access_bitapproximation function
  (i.e. either a SimHash, a list iterator or a lattice point with Bitapproximation)
*/
template<class CoordinateSelection, class LHS, class RHS, class LowerThresholds, class UpperThresholds>
FORCE_INLINE static inline bool CPP14CONSTEXPR check_simhash_scalar_product(LHS const &lhs, RHS const &rhs, LowerThresholds const &lb, UpperThresholds const &ub)
{
  uint_fast32_t approx_scprod = 0;
  for (unsigned int level = 0; level < CoordinateSelection::sim_hash_num; ++level)
  {
    approx_scprod += (   ObtainSimHashBlock<CoordinateSelection>::get(lhs,level)
                       ^ ObtainSimHashBlock<CoordinateSelection>::get(rhs,level)  ).count();
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



// Static Initializer for run-time initialization:
template<class CooSelection> class StaticInitializer< GlobalBitApproxData<CooSelection> >
  : public DefaultStaticInitializer< GlobalBitApproxData<CooSelection> >
{
  using Parent = DefaultStaticInitializer< GlobalBitApproxData<CooSelection> >;
  using DimensionType = typename CooSelection::Get_DimensionType;
  public:

  template<class T,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const & initializer, unsigned int random_seed) : StaticInitializer(initializer.dim,random_seed) {}

  StaticInitializer(DimensionType const new_dim, unsigned int random_seed)
  {
    assert(Parent::user_count > 0);
    if(Parent::user_count>1)
    {
      assert(new_dim == GlobalBitApproxData<CooSelection>::dim_used);
      assert(random_seed == GlobalBitApproxData<CooSelection>::random_seed_used);
      // TODO: Throw exception!
    }
    else
    {
      GlobalBitApproxData<CooSelection>::dim_used = new_dim;
      GlobalBitApproxData<CooSelection>::random_seed_used = random_seed;
      GlobalBitApproxData<CooSelection>::coo_selection = CooSelection{new_dim,random_seed};
    }
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing Global data for Bitapproximations." << " Counter is" << Parent::user_count )
  }
  ~StaticInitializer()
  {
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing Global data for Bitapproximations." << " Counter is " << Parent::user_count )
  }

//    StaticInitializer<SimHash::CoordinateSelection<void,false>> init_bitapprox;
//  GaussSieve::StaticInitializer<DMatrix> init_D_matrices;
//  GaussSieve::StaticInitializer<PMatrix> init_P_matrices;
  //GaussSieve::StaticInitializer<RelevantCoordinates> init_relevant_coo_matrix;
};



}  // end namespace GaussSieve

#endif
