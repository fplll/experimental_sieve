// clang-format off

#ifndef POINT_LIST_WITH_BITAPPROX_H
#define POINT_LIST_WITH_BITAPPROX_H

#include "DefaultIncludes.h"

#include <mutex>
#include <atomic>
#include <forward_list>
#include <queue>
#include <stack>
#include "Typedefs.h"
#include <list>
#include "SieveUtility.h"
#include "BitApproximation.h"

namespace GaussSieve{

//forward declarations
template<class SieveTraits, bool MT> class GaussListWithBitApprox;
template<class SieveTraits, bool MT> class GaussIteratorBitApprox;


// MT implementation:


template<class SieveTraits>
class MTNodeWithBitapprox
{

};

template<class SieveTraits>
class GaussListWithBitApprox<SieveTraits, true>
{
public:
  using StoredPoint = typename SieveTraits::GaussList_StoredPoint;
  using ReturnType  = typename SieveTraits::GaussList_ReturnType;
  using Iterator    = GaussIteratorBitApprox<SieveTraits,true>;
  friend Iterator;
};




}  // end namespace GaussSieve

#endif

//clang-format on
