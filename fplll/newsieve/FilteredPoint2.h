// clang-format off

//
//  FilteredPoint.h
//
//
//  Created by Elena on 06/03/17.
//
//

#ifndef FILTERED_POINT2_H
#define FILTERED_POINT2_H

#include "DefaultIncludes.h"
#include "DebugAll.h"
#include "SieveUtility.h"
#include "ExactLatticePoint.h"
//#include "HashedLatticePoint.h"
#include "GaussListBitapprox.h"
#include "BitApproximationNew.h"
#include "Typedefs.h"

namespace GaussSieve{

//template <class ET, int nfixed> class ExactLatticePoint;
//template <class ET, int nfixed> class HashedLatticePoint;

//template <class SieveTraits, class ET> class FilteredPoint;
template<class SieveTraits, bool MT> class FilteredPoint;

// Template parameters are:
//  ET: entry type
//  nfixed: indicates whether the dimension is fixed or not
//  SC: scalar-product type

template<class SieveTraits>
class FilteredPoint<SieveTraits, false>
//template <class SieveTraits, class ET>
//class FilteredPoint
{
public:

  using LengthType    = typename SieveTraits::LengthType;
  using StoredData    = STNode<SieveTraits>;

  FilteredPoint()                       = delete;
  FilteredPoint(const FilteredPoint &)  = delete;
  FilteredPoint(FilteredPoint &&)       = default;

  explicit constexpr FilteredPoint(StoredData const * pt, LengthType const & sc_pr) noexcept
      : point(pt), sc_prod(sc_pr) {}

  FilteredPoint& operator=(FilteredPoint const &) = delete;
  FilteredPoint& operator=(FilteredPoint &&     ) = default;

  ~FilteredPoint()  = default; // Note: Pointer is NOT owning.

  inline StoredData const& get_point()   const {return *point;}
  inline LengthType        get_sc_prod() const {return sc_prod;}
  // inline StoredData& get_ptr_toexact() const {return ptr_to_exact;}
  // inline bool get_sign() const {return minus;}

private:
  // StoredPoint* ptr_to_exact;

  StoredData const * point;
  LengthType sc_prod;
};

}

#endif
