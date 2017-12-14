
//
//  FilteredPoint.h
//
//
//  Created by Elena on 06/03/17.
//
//

#ifndef FILTERED_POINT2_H
#define FILTERED_POINT2_H

#include "DebugAll.h"
#include "DefaultIncludes.h"
#include "ExactLatticePoint.h"
#include "GaussListBitapprox.h"
#include "SieveUtility.h"
#include "Typedefs.h"

namespace GaussSieve
{

/**
  Old version, unused atm.
*/

// Might not work due to changes to NodeST internals.
/*
template<class SieveTraits, bool MT> class FilteredPoint;
template<class SieveTraits>
class FilteredPoint<SieveTraits, false>
{
public:

  using LengthType = typename SieveTraits::LengthType;
  using StoredData = STNode<SieveTraits>;

  FilteredPoint()                        = delete;
  FilteredPoint(FilteredPoint const  &)  = delete;
  FilteredPoint(FilteredPoint       &&)  = default;

  explicit constexpr FilteredPoint(StoredData const * pt, LengthType const & sc_pr) noexcept
      : point(pt), sc_prod(sc_pr) {}

  FilteredPoint& operator=(FilteredPoint const  &) = delete;
  FilteredPoint& operator=(FilteredPoint       &&) = default;

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
*/

/**
  New implementation:
*/

template <class SieveTraits, bool MT> class GaussListWithBitApprox;
template <class SieveTraits, bool MT> class GaussIteratorBitApprox;
template <class SieveTraits, bool MT> struct FilteredPoint2;

// clang-format off
template<class SieveTraits>
struct FilteredPoint2<SieveTraits,false>
{
  using LengthType   = typename SieveTraits::LengthType;
  using StoredPoint  = typename SieveTraits::GaussList_StoredPoint;  // or pointer
  using SimHashBlock = typename SieveTraits::SimHashBlock;
  using SimHashes    = typename SieveTraits::SimHashes;
  // clang-format on

  SimHashes sim_hashes;  // stores sim_hashes to *ptr_to_exact if sign_flip == false
                         // otherwise, stores bit-negated sim_hashes to *ptr_to_exact
  bool sign_flip;
  typename SieveTraits::GaussList_StoredPoint const *ptr_to_exact;  // non-owning pointer
  LengthType cond;  // stores -||x||^2 - 2 * <p, +/-x>
                    // where  x==*ptr_to_exact and p is the point with respect
                    // to which the filtered list is computed.
                    // This expression is optimized for the algorithm.
                    // Note that this will always be <=0.
  FilteredPoint2()                       = delete;
  FilteredPoint2(FilteredPoint2 const &) = delete;
  FilteredPoint2(FilteredPoint2 &&)      = default;

  // clang-format off
  explicit constexpr FilteredPoint2(GaussIteratorBitApprox<SieveTraits, false> const &list_iterator,
                                    bool const flip,
                                    LengthType const &precompute) noexcept
      : sim_hashes(flip ? flip_all_bits(list_iterator.get_all_bitapproximations())
                        :               list_iterator.get_all_bitapproximations() ),
        sign_flip(flip),
        ptr_to_exact(static_cast<StoredPoint const *>(list_iterator)),
        cond(precompute)
  {
  }
  // clang-format on
};

}  // end namespace GaussSieve

#endif
