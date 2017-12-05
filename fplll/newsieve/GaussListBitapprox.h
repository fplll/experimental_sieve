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
#include "BitApproximationNew.h"

namespace GaussSieve{

//forward declarations
template<class SieveTraits, bool MT> class GaussListWithBitApprox;
template<class SieveTraits, bool MT> class GaussIteratorBitApprox;

using SimHashApproxNorm2 = double;

// MT implementation:


// internal element stored in the list, single-threaded cases
template<class SieveTraits>
struct STNode
{
  std::array< SimHashNew::SimHashBlock<SieveTraits,false>, SieveTraits::sim_hash_num > bit_approximations;  // or C-Style?
  SimHashApproxNorm2 approx_norm2;  // consider making it a float
  typename SieveTraits::GaussList_StoredPoint* ptr_to_exact;
};


// single-threaded variant
template<class SieveTraits>
class GaussListWithBitApprox<SieveTraits, false>
{
public:
  using StoredPoint  = typename SieveTraits::GaussList_StoredPoint;
  using ReturnType   = typename SieveTraits::GaussList_ReturnType;
  using Iterator     = GaussIteratorBitApprox<SieveTraits,true>;
  using SimHashBlock = SimHashNew::SimHashBlock<SieveTraits,false>;
  friend Iterator;
  using UnderlyingContainer = std::list<STNode< SieveTraits> >;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;
  GaussListWithBitApprox()                                          = delete;
  GaussListWithBitApprox(GaussListWithBitApprox const &)            = delete;
  GaussListWithBitApprox(GaussListWithBitApprox &&)                 = delete;
  GaussListWithBitApprox& operator= (GaussListWithBitApprox const&) = delete;
  GaussListWithBitApprox& operator= (GaussListWithBitApprox &&)     = delete;
  explicit GaussListWithBitApprox(GlobalStaticDataInitializer const &static_data)
      :sim_hash_data(static_data.dim), init_stored_point(static_data),
       init_return_type(static_data), actual_list() {}

  Iterator cbegin() const;
  Iterator cend() const;

  //  Note: LatticePoint && is a forwarding/universal reference, not a rvalue reference.
  //        The function static_asserts that it is called on a rvalue.
  //        This means you might have to use insert_before(pos, std::move(new_lp));
  template<class LatticePoint, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatticePoint>>)>
  Iterator insert_before(Iterator const &pos, LatticePoint && new_point);

  Iterator erase(Iterator pos);
  void sort();

  static SimHashApproxNorm2 get_approx_norm2(Iterator const &);
  static std::array<SimHashBlock,SieveTraits::sim_hash_num> get_all_bitapproximations(Iterator const &it);
  // static std::array<SimHashBlock,SieveTraits::sim_hash_num> get_all_bitapproximations_ref(Iterator const &it);
  static SimHashBlock const & get_bitapproximation(Iterator const &it, unsigned int level);
  static StoredPoint const & get_point_r(Iterator const &it);
  static ReturnType pop_point(Iterator &it);

  private:
  SimHashNew::CoordinateSelection<SieveTraits,false> const sim_hash_data;
  StaticInitializer<StoredPoint> const init_stored_point;
  StaticInitializer<ReturnType>  const init_return_type;
  UnderlyingContainer actual_list;
};

//  Note: This is *NOT* an iterator in the sense of the C++ concept.
template<class SieveTraits>
class GaussIteratorBitApprox<SieveTraits, false>
{
private:
  using ListType = GaussListWithBitApprox<SieveTraits,false>;
  using UnderlyingIterator  = typename ListType::UnderlyingContainer::iterator;
  using CUnderlyingIterator = typename ListType::UnderlyingContainer::const_iterator;
  UnderlyingIterator it;
  friend GaussListWithBitApprox<SieveTraits, false>;
public:
  GaussIteratorBitApprox() = delete;
  GaussIteratorBitApprox(GaussIteratorBitApprox const & )           = default;
  GaussIteratorBitApprox(GaussIteratorBitApprox &&)                 = default;
  GaussIteratorBitApprox& operator=(GaussIteratorBitApprox const &) = default;
  GaussIteratorBitApprox& operator=(GaussIteratorBitApprox &&)      = default;
  explicit GaussIteratorBitApprox(UnderlyingIterator const &new_it) : it(new_it) {}
  bool operator==(GaussIteratorBitApprox const &other) const { return it==other.it; }
  bool operator!=(GaussIteratorBitApprox const &other) const { return it!=other.it; }
  GaussIteratorBitApprox& operator++()    { ++it; return *this; }  // prefix version
  GaussIteratorBitApprox  operator++(int) { return it++; }         // postfix version
  // is_end
};


}  // end namespace GaussSieve

#include "GaussListBitapprox_impl.h"

#endif

//clang-format on
