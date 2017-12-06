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
// The !is_reference<Arg> makes this template only valid for actual rvalues.
  template<class Arg,
      TEMPL_RESTRICT_DECL2(IsALatticePoint< mystd::decay_t<Arg> >, mystd::negation<std::is_reference<Arg>>)>
  explicit constexpr STNode(Arg &&arg, SimHashNew::CoordinateSelection<SieveTraits,false> const &coo_selection) noexcept
      : bit_approximations(coo_selection.compute_all_bitapproximations(arg)),
        approx_norm2(convert_to_double(arg.get_norm2())),
        ptr_to_exact(new typename SieveTraits::GaussList_StoredPoint(std::move(arg)) ) {}
  explicit constexpr STNode(typename SieveTraits::GaussList_StoredPoint * const point_ptr, SimHashNew::CoordinateSelection<SieveTraits,false> const &coo_selection) noexcept
      : bit_approximations(coo_selection.compute_all_bitapproximations(*point_ptr)),
        approx_norm2(convert_to_double(point_ptr->get_norm2())),
        ptr_to_exact(point_ptr) {}

  ~STNode() noexcept { delete ptr_to_exact; }

  std::array< SimHashNew::SimHashBlock<SieveTraits,false>, SieveTraits::sim_hash_num > bit_approximations;  // or C-Style?
  SimHashApproxNorm2 approx_norm2;  // consider making it a float
  typename SieveTraits::GaussList_StoredPoint* ptr_to_exact; // owning pointer
  bool operator<(STNode const &other) const {return *ptr_to_exact < *(other.ptr_to_exact); }
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
  explicit GaussListWithBitApprox(GlobalStaticDataInitializer const &static_data) noexcept
      :sim_hash_data(static_data.dim), init_stored_point(static_data),
       init_return_type(static_data), actual_list() {}

  CPP14CONSTEXPR Iterator cbegin() const noexcept { return actual_list.cbegin(); }
  CPP14CONSTEXPR Iterator cend() const noexcept   { return actual_list.cend(); }

  //  Note: LatticePoint && is a forwarding/universal reference, not a rvalue reference.
  //        The function static_asserts that it is called on a rvalue, to be sure.
  //        This means you might have to use insert_before(pos, std::move(new_lp));
  template<class LatticePoint, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatticePoint>>)>
  Iterator insert_before(Iterator const &pos, LatticePoint && new_point)
  {
    static_assert(std::is_reference<LatticePoint>::value == false,"Must call on rvalues");
    return actual_list.emplace(pos.it, std::move(new_point), sim_hash_data);
  }

  Iterator insert_before(Iterator const &pos, StoredPoint * const point_ptr)
  {
    return actual_list.emplace(pos.it, point_ptr, sim_hash_data);
  }

  ReturnType true_pop_point(Iterator &pos)
  {
    ReturnType retval = static_cast<ReturnType>(*pos);
    pos.it->ptr_to_exact = nullptr;
    pos = actual_list.erase(pos);
    return retval;
  }

  Iterator erase(Iterator pos) { return actual_list.erase(pos.it); }
  void sort() { actual_list.sort(); }

  private:
  SimHashNew::CoordinateSelection<SieveTraits,false> const sim_hash_data;
  StaticInitializer<StoredPoint> const init_stored_point;
  StaticInitializer<ReturnType>  const init_return_type;
  UnderlyingContainer actual_list;
};

//  Note: This may or may not be an iterator in the sense of the C++ concept.
//        The issue is that we might deprecate dereferencing...
//  TOOD: If we allow dereferencing, set iterator_traits for compatibility with std::
template<class SieveTraits>
class GaussIteratorBitApprox<SieveTraits, false>
{
private:
  using ListType = GaussListWithBitApprox<SieveTraits,false>;
  using StoredPoint  = typename SieveTraits::GaussList_StoredPoint;
  using ReturnType   = typename SieveTraits::GaussList_ReturnType;
  using UnderlyingIterator  = typename ListType::UnderlyingContainer::iterator;
  using CUnderlyingIterator = typename ListType::UnderlyingContainer::const_iterator;
  using SimHashBlock = SimHashNew::SimHashBlock<SieveTraits,false>;

  CUnderlyingIterator it;

  friend GaussListWithBitApprox<SieveTraits, false>;
public:
  GaussIteratorBitApprox() = delete;
  GaussIteratorBitApprox(GaussIteratorBitApprox const &)            = default;
  GaussIteratorBitApprox(GaussIteratorBitApprox &&)                 = default;
  GaussIteratorBitApprox& operator=(GaussIteratorBitApprox const &) = default;
  GaussIteratorBitApprox& operator=(GaussIteratorBitApprox &&)      = default;
  explicit constexpr GaussIteratorBitApprox(UnderlyingIterator const &new_it) : it(new_it) {}
  bool operator==(GaussIteratorBitApprox const &other) const { return it==other.it; }
  bool operator!=(GaussIteratorBitApprox const &other) const { return it!=other.it; }
  GaussIteratorBitApprox& operator++()    { ++it; return *this; }  // prefix version
  GaussIteratorBitApprox  operator++(int) { return it++; }         // postfix version
  // is_end

  SimHashApproxNorm2 get_approx_norm2() const;
  auto get_all_bitapproximations() const
      -> std::array< SimHashBlock, SieveTraits::sim_hash_num >
  { return it->bit_approximations; }

  // the distinction between these two is not very significant in the single-threaded case
  // in the multi-threaded cases, there is a difference, because access* returns a
  // const-reference to atomic<SimHashBlock>, whereas get_ returns a (new, non-atomic copy of)
  // SimHashBlock by value.
  // Note that atomic<SimHashBlock> is convertible to SimHashBlock, but that conversion uses
  // std::memory_order_seq_cst, which may be slower than needed.

  // TODO: Wrap atomic<SimHashBlock> into another type that is convertible to SimHashBlock via
  // a load with a custom memory order

  SimHashBlock const & access_bitapproximation(unsigned int level) const
  { return it->bitapproximations[level]; }
  //SimHashBlock         get_bitapproximation(unsigned int level) const;
  StoredPoint const &  operator*() const    { return *(it->ptr_to_exact); }
  StoredPoint const *  operator->() const   { return it->ptr_to_exact;    }
  explicit operator StoredPoint* ()         { return it->ptr_to_exact;    }
};


}  // end namespace GaussSieve

#include "GaussListBitapprox_impl.h"

#endif

//clang-format on
