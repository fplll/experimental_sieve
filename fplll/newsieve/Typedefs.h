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


namespace GaussSieve
{
// unfortunately, trigonometric functions to compute pi don't have constexpr variants on all
// compilers we want to support, so we just define pi directly
long double constexpr pi_long = 3.14159265358979323846264338327950288419716939937510L;
double constexpr pi_double    = 3.14159265358979323846264338327950288419716939937510;
long double constexpr pi      = 3.14159265358979323846264338327950288419716939937510L;



// forward-declarations:
template <class ET, int nfixed> class MyLatticePoint;
template <class ET, int nfixed> class PlainLatticePoint;
template <class ET, int nfixed> class ExactLatticePoint;

class JustSomeExampleSieveTraitsThatDoNotWork
{
  public:
  using IsSieveTraitsClass = std::true_type;   // set to true to indicate purpose of trait class.
  using GaussSampler_ReturnType = void;        // What the Sampler has to return

  using GaussList_ReturnType = void;           // What the list has to return.
  using GaussList_StoredPoint= void;           // What is actually stored inside the list.

  using GaussQueue_ReturnType= void;           // What pop()-ing the queue returns.
  using GaussQueue_DataType  = void;           // Storage inside queue.
  using FastAccess_Point     = void;           // For internal use
  using DimensionType        = unsigned int;   // Type to store the dimension.

  using InputBasisType       = void;           // Input basis type that is used.

  // Indicates whether we need to remove Z_NR from input.
  // While we don't use Z_NR classes internally (much), the input may be.
  using UnZNRInput           = std::true_type;

// Note that MT is intentionally missing here. MT is its own template argument.
};

// Note: ET does *not* include Z_NR<...> here

template
<
  class ET, bool MT, int nfixed,
  class InputBT = typename fplll::ZZ_mat< typename FixMpz_classToMpz_t<ET>::type>
>
class DefaultSieveTraits
{
  public:
  using IsSieveTraitsClass = std::true_type;
  using GaussSampler_ReturnType = ExactLatticePoint<ET,nfixed>;
  using GaussList_StoredPoint   = ExactLatticePoint<ET,nfixed>;
  using GaussList_ReturnType    = ExactLatticePoint<ET,nfixed>;
  using GaussQueue_ReturnType   = GaussSampler_ReturnType;
  using GaussQueue_DataType     = GaussQueue_ReturnType;
  using FastAccess_Point        = ExactLatticePoint<ET,nfixed>;
  using DimensionType           = MaybeFixed<nfixed>;
  using EntryType               = ET;
  using ZNREntryType            = typename AddZNR<ET>::type; // should be unused

  //for the class FilteredPoint, the template paremeters are: <Entry type, if_dim_is_fixed, scalar_prod. type>
  //using FlilteredPointType      = FilteredPoint<ET, nfixed, EntryType>;
  using FlilteredPointType      = FilteredPointPointer<ET, nfixed, EntryType>;
  using FilteredListType        = std::vector<FlilteredPointType>;
  using GlobalStaticDataInitializer = StaticDataInitializer<DimensionType>;

  // note that if ET = mpz_class, then ZNREntryType::underlying_data_type = mpz_t,
  // otherwise ET == ZNREntryType::underlying_data_type
//  using InputBasisType = typename fplll::ZZ_mat< typename ZNREntryType::underlying_data_type>;
//  using InputBasisType          = fplll::ZZ_mat< mpz_t>;
  using InputBasisType          = InputBT;
  using PlainPoint              = PlainLatticePoint<ET,nfixed>;
  static int constexpr get_nfixed = nfixed; // TODO: Remove and forward DimensionType throughout...


  //--------HYPERPLANE LSH SPECIFIC----------
  static constexpr unsigned short number_of_hash_tables = 12;
  static constexpr int number_of_hash_functions = 11;
  static constexpr int max_bucket_size = 250;
};


// clang-format on
};

#endif
