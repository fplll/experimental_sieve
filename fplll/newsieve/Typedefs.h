// This header defines (global) traits. We might want to encapsulte these in some trait classes that
// get forwarded, eventually.
// This header should not depend on anything within the sieve.

// Consider renaming the file to avoid clashes.

#ifndef GAUSS_SIEVE_TYPEDEFS_H
#define GAUSS_SIEVE_TYPEDEFS_H

#include "DefaultIncludes.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"
#include "gmpxx.h"
#include "ExactLatticePoint.h"
#include "PlainLatticePoint.h"
#include "FilteredPoint.h"
#include "FilteredPoint2.h"
#include "SieveUtility.h"
#include "GlobalStaticData.h"
//#include "HashedLatticePoint.h"
//#include "ApproximatedPoint.h"
#include "EMVApproximation.h"
#include "BitApproximationNew.h"
#include "GlobalBitApproxData.h"
#include "PointWithBitapprox.h"


namespace GaussSieve
{
// unfortunately, trigonometric functions to compute pi don't have constexpr variants on all
// compilers we want to support, so we just define pi directly
long double constexpr pi_long = 3.14159265358979323846264338327950288419716939937510L;
double constexpr pi_double    = 3.14159265358979323846264338327950288419716939937510;
long double constexpr pi      = 3.14159265358979323846264338327950288419716939937510L;

double constexpr list_size_k2 = 0.2075187494;
double constexpr list_size_k3 = 0.1887218757;
//double constexpr list_size_k3 = 0.195;
double constexpr list_size_k4 = 0.1723692862;




// forward-declarations:
template <class ET, int nfixed> class MyLatticePoint;
template <class ET, int nfixed> class PlainLatticePoint;
template <class ET, int nfixed> class ExactLatticePoint;
template <class ET, int nfixed> class HashedLatticePoint;
template <class ELP, class Approximation> class VectorWithApproximation;
template <class ELP, class CooSelection>  class AddBitApproximationToLP;

// Note: ET does *not* include Z_NR<...> here

template
<
  class ET, bool MT, int nfixed,
  class InputBT = typename fplll::ZZ_mat< typename FixMpz_classToMpz_t<ET>::type>
>
class DefaultSieveTraits
{
  public:

  //static constexpr int max_bucket_size = 250;

  using IsSieveTraitsClass = std::true_type;

  using GaussSampler_ReturnType = ExactLatticePoint<ET,nfixed>;
  using GaussList_StoredPoint   = ExactLatticePoint<ET,nfixed>;
  using GaussList_ReturnType    = ExactLatticePoint<ET,nfixed>;
  using GaussQueue_ReturnType   = GaussSampler_ReturnType;
  using GaussQueue_DataType     = GaussQueue_ReturnType;

//#endif

#ifdef USE_LSH
  using GaussSampler_ReturnType = HashedLatticePoint<ET,nfixed>;
  using GaussList_StoredPoint   = HashedLatticePoint<ET,nfixed>;
  using GaussList_ReturnType    = HashedLatticePoint<ET,nfixed>;
  using FastAccess_Point        = HashedLatticePoint<ET,nfixed>;

  //--------HYPERPLANE LSH SPECIFIC----------
  static constexpr unsigned short number_of_hash_tables = HashedLatticePoint<ET,nfixed>::number_of_hash_tables;
  static constexpr int number_of_hash_functions = 11;
#endif

  static std::size_t constexpr sim_hash_len = 64;  // number of bits per simhash block
  static std::size_t constexpr sim_hash_num = 2;   // number of simhash blocks/levels per vector
  // -> Total number of bits is given by sim_hash_len * sim_hash_num

  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_2sieve_lb = {{64-5,128-8}};
  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_2sieve_ub = {{64+5,128+8}};

//  constexpr std::array<unsigned int, sim_hash_num> threshold_lvls_3sieve = {{0}};

  using DimensionType           = MaybeFixed<nfixed>;

  using SimHashGlobalDataType   = SimHashNew::CoordinateSelection<sim_hash_len,sim_hash_num,MT,DimensionType>;
  using SimHashBlock            = typename SimHashGlobalDataType::SimHashBlock;
  using SimHashes               = typename SimHashGlobalDataType::SimHashes;


  using EntryType               = ET;

  //for the class FilteredPoint, the template paremeters are: <Entry type, if_dim_is_fixed, scalar_prod. type>
  //using FlilteredPointType      = FilteredPoint<ET, nfixed, EntryType>;
  using FlilteredPointType      = FilteredPointPointer<ET, nfixed, EntryType>;
  using FilteredListType        = std::vector<FlilteredPointType>;
  using GlobalStaticDataInitializer = StaticInitializerArg<DimensionType>;

  using FastAccess_Point        = AddBitApproximationToLP< ExactLatticePoint<ET,nfixed>, SimHashGlobalDataType >;


  // note that if ET = mpz_class, then ZNREntryType::underlying_data_type = mpz_t,
  // otherwise ET == ZNREntryType::underlying_data_type
//  using InputBasisType = typename fplll::ZZ_mat< typename ZNREntryType::underlying_data_type>;
//  using InputBasisType          = fplll::ZZ_mat< mpz_t>;
  using InputBasisType          = InputBT;
  using PlainPoint              = PlainLatticePoint<ET,nfixed>;
  static int constexpr get_nfixed = nfixed; // TODO: Remove and forward DimensionType throughout...

};
  template< class ET, bool MT, int nfixed, class InputBT>
  constexpr std::array<unsigned int, 2> DefaultSieveTraits<ET,MT,nfixed, InputBT>::threshold_lvls_2sieve_lb;
  template< class ET, bool MT, int nfixed, class InputBT>
  constexpr std::array<unsigned int, 2> DefaultSieveTraits<ET,MT,nfixed, InputBT>::threshold_lvls_2sieve_ub;

//  template< class ET, bool MT, int nfixed, class InputBT>
//  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_2sieve_ub; = {{64+5,128+8}};


// clang-format on
};

#endif
