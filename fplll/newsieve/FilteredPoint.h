//
//  FilteredPoint.h
//
//
//  Created by Elena on 06/03/17.
//

/**
  In the k-sieve for k>=3, we need to create sub-lists of the main list insider the iteration.
  Notably, in the sieve iteration for a given candidate p=x_0, we need (for k=3) we iterate over
  the main list and find candidates for x_1, x_2, given by some criteria depending only on p.
  We then try (for p fixed) x_1, x_2 - pairs.
  These sublists are stored as a list of filtered points, which are are class defined in this file.
  Note that there is the option of storing some precomputed values, such as the scalar product with
  p, to avoid recomputation / aid in determining good triples more quickly.
  This class is subject to experimentation what data is worth storing.
*/

#ifndef FILTERED_POINT2_H
#define FILTERED_POINT2_H

#include "DebugAll.h"
#include "DefaultIncludes.h"
#include "ExactLatticePoint.h"
//#include "GaussListBitapprox.h"
#include "GaussVectorBitApprox.h"
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
  using LengthType            = typename SieveTraits::LengthType;
  // we store a pointer to those:
  using GaussList_StoredPoint = typename SieveTraits::GaussList_StoredPoint;
  using SimHashBlock          = typename SieveTraits::SimHashBlock;
  using SimHashes             = typename SieveTraits::SimHashes;
  // clang-format on

  SimHashes sim_hashes;  // stores sim_hashes to *ptr_to_exact if sign_flip == false
                         // otherwise, stores bit-negated sim_hashes to *ptr_to_exact
  bool sign_flip;
  GaussList_StoredPoint const *ptr_to_exact;  // non-owning pointer

  /*
    cond stores -||x||^2 - 2 * <p, +/-x>, where  x==*ptr_to_exact and p is the point with respect to
    which the filtered list is computed. This expression is optimized for the algorithm and known at
    the call site.
    Note: The name cond is because this quantity determines whether a two-reduction is possible. In
    particular, cond will always be <= 0. (Otherwise, we perform a 2-reduction directly and do not
    need to use FilteredPoint)
  */
  LengthType cond;
  FilteredPoint2()                       = delete;
  FilteredPoint2(FilteredPoint2 const &) = delete;
  FilteredPoint2(FilteredPoint2 &&)      = default;

  // clang-format off
  // we take an iterator to the main lists as an argument rather than a point itself, because
  // the iterator actually holds more data (like sim_hashes)
  // flip determines whether the stored point implicitly has a - sign.
  // by flipping the sign inside the list of filtered points, we do not need to perform a case
  // distinction when iterating over pairs from the filtered list.
  // precompute equals -||x||^2 -/+ 2<p,x>, which we record for later use.
  explicit constexpr FilteredPoint2(GaussIteratorBitApprox<SieveTraits, false> const &list_iterator,
                                    bool const flip,
                                    LengthType const &precompute) noexcept
      : sim_hashes(flip ? flip_all_bits(list_iterator.get_all_bitapproximations())
                        :               list_iterator.get_all_bitapproximations() ),
        sign_flip(flip),
        ptr_to_exact(static_cast<GaussList_StoredPoint const *>(list_iterator)),
        cond(precompute)
  {
  }
  // clang-format on
};

}  // end namespace GaussSieve

#endif  // include guards
