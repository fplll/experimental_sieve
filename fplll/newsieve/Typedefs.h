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
#include "SieveUtility.h"
#include "GlobalStaticData.h"
// #include "HashedLatticePoint.h"
// #include "ApproximatedPoint.h"
// #include "EMVApproximation.h"  // currently unused
#include "BlockOrthogonalSimHash.h"
#include "SimHash.h"
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
// double constexpr list_size_k3 = 0.195;
double constexpr list_size_k4 = 0.1723692862;




// forward-declarations:
template<class ET, int nfixed> class MyLatticePoint;
template<class ET, int nfixed> class PlainLatticePoint;
template<class ET, int nfixed> class ExactLatticePoint;
template<class ET, int nfixed> class HashedLatticePoint;
template<class ELP, class Approximation> class VectorWithApproximation;
template<class ELP, class CooSelection>  class AddBitApproximationToLP;

// Note: ET does *not* include Z_NR<...> here

template
<
  class ET, bool MT, int nfixed,
  class InputBT = typename fplll::ZZ_mat< PossiblyMpzClassToMpzt<ET> >
>
class DefaultSieveTraits
{
  public:

  using IsSieveTraitsClass = std::true_type;

  using DimensionType           = MaybeFixed<nfixed>;

#if defined(OVERRIDE_SIM_HASH_LEN)
  static std::size_t constexpr sim_hash_len = OVERRIDE_SIM_HASH_LEN;  // number of bits per simhash block
#else
  static std::size_t constexpr sim_hash_len = 64;  // number of bits per simhash block
#endif
#if defined(OVERRIDE_SIM_HASH_NUM)
  static std::size_t constexpr sim_hash_num = OVERRIDE_SIM_HASH_NUM;   // number of simhash blocks/levels per vector
#else
  static std::size_t constexpr sim_hash_num = 2;   // number of simhash blocks/levels per vector
#endif
  // -> Total number of bits is given by sim_hash_len * sim_hash_num

  using CoordinateSelectionUsed = BlockOrthogonalSimHash<sim_hash_len, sim_hash_num, MT, DimensionType>;
  using SimHashBlock            = typename CoordinateSelectionUsed::SimHashBlock;
  using SimHashes               = typename CoordinateSelectionUsed::SimHashes;

  using LengthType               = ET;

  using GlobalStaticDataInitializer = StaticInitializerArg<DimensionType>;

  using FastAccess_Point        = AddBitApproximationToLP< ExactLatticePoint<ET,nfixed>, CoordinateSelectionUsed >;


  using GaussSampler_ReturnType = ExactLatticePoint<ET,nfixed>;
  using GaussList_StoredPoint   = ExactLatticePoint<ET,nfixed>;
//  using GaussList_ReturnType    = ExactLatticePoint<ET,nfixed>;
  using GaussList_ReturnType    = FastAccess_Point;
  using GaussQueue_ReturnType   = GaussSampler_ReturnType;
  using GaussQueue_DataType     = GaussQueue_ReturnType;

//#ifdef USE_LSH
//  using GaussSampler_ReturnType = HashedLatticePoint<ET,nfixed>;
//  using GaussList_StoredPoint   = HashedLatticePoint<ET,nfixed>;
//  using GaussList_ReturnType    = HashedLatticePoint<ET,nfixed>;
//  using FastAccess_Point        = HashedLatticePoint<ET,nfixed>;
//
//  //--------HYPERPLANE LSH SPECIFIC----------
//  static constexpr unsigned short number_of_hash_tables = HashedLatticePoint<ET,nfixed>::number_of_hash_tables;
//  static constexpr int number_of_hash_functions = 11;
//#endif


  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_2sieve_lb = {{32-7, 64-12}};
  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_2sieve_ub = {{32+7, 64+12}};

  // for 3-sieve: outer loop
  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_3sieve_lb_out = {{32-5, 64-11}};
  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_3sieve_ub_out = {{32+5, 64+11}};

  // for 3-sieve: outer loop
  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_3sieve_lb_inn = {{32-3, 64-6}};
  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_3sieve_ub_inn = {{32+3, 64+6}};

  // for 3-sieve: FilteredList is implemented as vector
  // we reserve filtered_list_size_max inside for its length

  constexpr static int filtered_list_size_max = 500;

  // for 3-sieve exact check, squared, normalized
  //constexpr static double x1x2_target = .1111;
  //constexpr static double x2x3_target = .1111;
  constexpr static double x1x2_target = .0911;
  constexpr static double x2x3_target = .0911;

  /* for printing routines */
  // print after every print_step_*sieve iterations
  constexpr static int print_step_2sieve  = 2000;
  constexpr static int print_step_3sieve  = 1000;



  // note that if ET = mpz_class, then ZNREntryType::underlying_data_type = mpz_t,
  // otherwise ET == ZNREntryType::underlying_data_type
//  using InputBasisType = typename fplll::ZZ_mat< typename ZNREntryType::underlying_data_type>;
//  using InputBasisType          = fplll::ZZ_mat< mpz_t>;
  using InputBasisType          = InputBT;
  using PlainPoint              = PlainLatticePoint<ET,nfixed>;
  static int constexpr get_nfixed = nfixed; // TODO: Remove and forward DimensionType throughout...

};
  template< class ET, bool MT, int nfixed, class InputBT>
  constexpr std::array<unsigned int, DefaultSieveTraits<ET,MT,nfixed,InputBT>::sim_hash_num> DefaultSieveTraits<ET,MT,nfixed,InputBT>::threshold_lvls_2sieve_lb;
  template< class ET, bool MT, int nfixed, class InputBT>
  constexpr std::array<unsigned int, DefaultSieveTraits<ET,MT,nfixed,InputBT>::sim_hash_num> DefaultSieveTraits<ET,MT,nfixed,InputBT>::threshold_lvls_2sieve_ub;
  template< class ET, bool MT, int nfixed, class InputBT>
  constexpr std::array<unsigned int, DefaultSieveTraits<ET,MT,nfixed,InputBT>::sim_hash_num> DefaultSieveTraits<ET,MT,nfixed,InputBT>::threshold_lvls_3sieve_lb_out;
  template< class ET, bool MT, int nfixed, class InputBT>
  constexpr std::array<unsigned int, DefaultSieveTraits<ET,MT,nfixed,InputBT>::sim_hash_num> DefaultSieveTraits<ET,MT,nfixed,InputBT>::threshold_lvls_3sieve_ub_out;

  template< class ET, bool MT, int nfixed, class InputBT>
  constexpr std::array<unsigned int, DefaultSieveTraits<ET,MT,nfixed,InputBT>::sim_hash_num> DefaultSieveTraits<ET,MT,nfixed,InputBT>::threshold_lvls_3sieve_lb_inn;
  template< class ET, bool MT, int nfixed, class InputBT>
  constexpr std::array<unsigned int, DefaultSieveTraits<ET,MT,nfixed,InputBT>::sim_hash_num> DefaultSieveTraits<ET,MT,nfixed,InputBT>::threshold_lvls_3sieve_ub_inn;

//  template< class ET, bool MT, int nfixed, class InputBT>
//  constexpr static std::array<unsigned int, sim_hash_num> threshold_lvls_2sieve_ub; = {{64+5,128+8}};


// clang-format on
};

#endif
