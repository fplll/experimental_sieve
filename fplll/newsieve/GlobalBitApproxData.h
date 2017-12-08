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

template<std::size_t sim_hash_num, std::size_t sim_hash_len>
std::ostream& operator<<(std::ostream &os, std::array< std::bitset<sim_hash_len>, sim_hash_num> const &sim_hashes)
{
  for(unsigned int i=0; i < sim_hash_num; ++i)
  {
    os << sim_hashes[i] << " ";
  }
  return os;
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
