#ifndef POINT_LIST_WITH_BITAPPROX_H
#define POINT_LIST_WITH_BITAPPROX_H

/**
  New version of GaussList (i.e. the main list of the lattice points we keep during the algorithm)
  Think of this as similar to a std::list of lattice points.
  Note, however that we directly integrate and expose SimHashes into the interface of the list.
  (rather than keeping them as members of the list elements)
  This is done for two main reasons:
  1.) During an iteration, most of the time we only look at the sim_hashes.
      In particular, it makes sense to consider storing a list of sim_hashes separately
      to optimize for looking at the sim_hashes: this gives better memory access patters.
  2.) For the multi-threaded implementation, we want to keep the sim_hashes directly inside the node
      The reason is that the thread synchronization requirements for sim_hashes are different from
      the requirements for the point: If we read a sim_hash that was only partially updated / messed
      up by concurrent threads, the resulting sim_hash is still a valid object (any sequence of bits
      is, after all). It will at worst make us give a false positive / false negative in prediction.
      Consequently, we can work with more relaxed atomics. Since memory access is actually a
      bottleneck for the 2-sieve, this matters a lot.

  We also try to maintain the same interface for the single- and multi-threaded versions.
  Since this is a bottleneck for the algorithm, we allow somewhat weird interface definitions (at
  least to some extent) for the sake of efficiency.
  NOTE: The interface is NOT stable at all.
*/

#include "DefaultIncludes.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include "SimHash.h"
#include "Typedefs.h"

namespace GaussSieve
{

// forward declarations
template <class SieveTraits, bool MT> class GaussListWithBitApprox;
template <class SieveTraits, bool MT> class GaussIteratorBitApprox;

// we directly store an approximation (which is actually exact in most cases) to norm2 in the list
// nodes. The class storing this should have a fixed size in memory (no dynamic allocation!)
// and NOT depend on the input data types.
// Currently, this is not used.
using SimHashApproxNorm2 = double;

/**
  STNode is what we currently internally store inside the list for the single-threaded case.
  It is an implementation detail of the list class below and not supposed to be exposed.
  NOTE: There are several choices here that matter for efficiency, which need some profiling and
        experimentation. This is subject to change.
        (E.g. we could store 2 separate lists instead and make our iterators to the main list
        internally consist of a pair of iterators.)

  STNode currently stores an approximation to a point (sim_hashes and an approximation to the norm2)
  and an (owning) pointer to the exact lattice point.

  The class used for the latter is GaussList_StoredPoint (given by SieveTraits), and is
  supposed NOT to have a bitapproximation.
*/

namespace GaussListDetails
{

// clang-format off
template <class SieveTraits>
class STNode
{
  friend GaussListWithBitApprox<SieveTraits, false>;
  friend GaussIteratorBitApprox<SieveTraits, false>;
  // retrieve typedefs to avoid having to write long names.

private:  // shorthands
  using CoordinateSelectionUsed = typename SieveTraits::CoordinateSelectionUsed;
  using SimHashes               = typename CoordinateSelectionUsed::SimHashes;
  using GlobalSimHashClass      = GlobalBitApproxData<CoordinateSelectionUsed>;
  // This would store 2 independently maintained sim_hashes.
  static_assert(Has_BitApprox<typename SieveTraits::GaussList_StoredPoint>::value == false, "");

public:
  // we should never need to copy. If we do, we want the compiler to tell us.
  STNode(STNode const  &)          = delete;
  STNode(STNode       &&) noexcept = default;
  // clang-format on

  // We construct a STNode from a lattice point.
  // We differentiate whether the input lattice point already stores a sim_hash
  // (if no, we compute a sim_hash, if yes, we take it from the point).
  // Note that the initializer list also (tries to) convert the point, so if the argument
  // type differs from GaussList_StoredPoint, it needs to be convertible
  // (AddBitApproximationToPoint (cf. PointWithApproximation.h) ensures the latter)
  // The !is_reference<Arg> in the restrictions makes this template only valid for actual rvalues.

  // Variant for Arg's without SimHashes
  template <class Arg, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<Arg>>,
                                            mystd::negation<std::is_reference<Arg>>,
                                            mystd::negation<Has_BitApprox<mystd::decay_t<Arg>>>)>
  explicit STNode(Arg &&arg) noexcept  // Note: new might actually throw.
      : bit_approximations(GlobalSimHashClass::coo_selection.compute_all_bitapproximations(arg)),
        approx_norm2(convert_to_double(arg.get_norm2())),
        ptr_to_exact(new typename SieveTraits::GaussList_StoredPoint(std::move(arg)))
  {
  }

  // Variant for Arg's with SimHashes
  template <class Arg, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<Arg>>,
                                            mystd::negation<std::is_reference<Arg>>,
                                            Has_BitApprox<mystd::decay_t<Arg>>)>
  explicit STNode(Arg &&arg) noexcept  // Note: new might actually throw...
      : bit_approximations(std::move(arg).take_bitapproximations()),
        approx_norm2(convert_to_double(arg.get_norm2())),
        ptr_to_exact(new typename SieveTraits::GaussList_StoredPoint(std::move(arg)))
  {
  }

  // currently unused: Takes a pointer to a Lattice Point (of the correct type)
  //                   This version has the advantage of not reallocating memory.
  //                   NOTE: THIS TAKES OWNERSHIP OF THE POINTER
  explicit constexpr STNode(typename SieveTraits::GaussList_StoredPoint *&&point_ptr) noexcept
      : bit_approximations(
            GlobalSimHashClass::coo_selection.compute_all_bitapproximations(*point_ptr)),
        approx_norm2(convert_to_double(point_ptr->get_norm2())),
        ptr_to_exact(std::move(point_ptr))
  {
  }

  ~STNode() noexcept { delete ptr_to_exact; }

private:
  SimHashes bit_approximations;
  SimHashApproxNorm2 approx_norm2;  // currently a double, consider making it a float
  typename SieveTraits::GaussList_StoredPoint *ptr_to_exact;  // owning pointer

public:
  // we keep the list sorted according to length. We define operator< to use sort() from std::list.
  // (Note that comparions for the lattice points *ptr_to_exact are by norm2)
  bool operator<(STNode const &other) const { return *ptr_to_exact < *(other.ptr_to_exact); }
};

}  // end namespace GaussListDetails

/**
  Single-threaded variant of our list class.
  Wraps around a std::list of STNodes.
  Note that it behaves like a forward_list (i.e. iterators can only be incremented) with the
  exception that we have insert_before rather than insert_after. Furthermore, we only have
  const-iterators (i.e. we can not modify the points in the list directly, we need to
  erase-modify-reinsert instead)
  Note that large parts of the functionality are defined as members of the iterator below.
*/
// clang-format off
template<class SieveTraits>
class GaussListWithBitApprox<SieveTraits, false>
{
private:
  using StoredPoint             = typename SieveTraits::GaussList_StoredPoint;
  using ReturnType              = typename SieveTraits::GaussList_ReturnType;
  using CoordinateSelectionUsed = typename SieveTraits::CoordinateSelectionUsed;

  using SimHashBlock            = typename CoordinateSelectionUsed::SimHashBlock;
  using SimHashes               = typename CoordinateSelectionUsed::SimHashes;

  using Iterator                = GaussIteratorBitApprox<SieveTraits, false>;  // custom iterator class below
  // clang-format on

  friend Iterator;
  // The class is essentially just a wrapper around UnderlyingContainer
  using UnderlyingContainer         = std::list<GaussListDetails::STNode<SieveTraits>>;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;

public:
  // clang-format off
  GaussListWithBitApprox()                                          = delete;
  GaussListWithBitApprox(GaussListWithBitApprox const &)            = delete;
  GaussListWithBitApprox(GaussListWithBitApprox &&)                 = delete;
  GaussListWithBitApprox &operator=(GaussListWithBitApprox const &) = delete;
  GaussListWithBitApprox &operator=(GaussListWithBitApprox &&)      = delete;
  // clang-format on

  // Sole constructor. The argument is used to initialize the static data of the used lattice point
  // classes.
  // NOTE / TODO: Making this noexcept means we have to catch reinitializations in the caller.
  explicit GaussListWithBitApprox(GlobalStaticDataInitializer const &static_data) noexcept
      : init_stored_point(static_data),
        init_return_type(static_data),
        actual_list()
  {
  }

  // behaves like cbegin, cend from STL containers, i.e. gives const-iterator to begin/end.
  // clang-format off
  CPP14CONSTEXPR Iterator cbegin() const noexcept { return actual_list.cbegin(); }
  CPP14CONSTEXPR Iterator cend()   const noexcept { return actual_list.cend();   }
  // clang-format on

  // insert_before(pos, point) inserts the point just before pos.
  // the return value is an iterator to the newly inserted point.
  //  Note: LatticePoint && is a forwarding/universal reference, NOT a rvalue reference. Still, the
  //        function is only supposed to be called on rvalues, so we static_asserts that.
  //        This means you might have to use insert_before(pos, std::move(new_lp));
  template <class LatticePoint, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatticePoint>>)>
  Iterator insert_before(Iterator const &pos, LatticePoint &&new_point)
  {
    static_assert(std::is_reference<LatticePoint>::value == false, "Must call on rvalues");
    return actual_list.emplace(pos.it, std::move(new_point));
  }

  // Inserts the lattice point pointed to by a pointer. The list takes ownership of the pointee.
  // TODO: Consider using unique_ptr
  // Untested
  Iterator insert_before(Iterator const &pos, StoredPoint *&&point_ptr)
  {
    return actual_list.emplace(pos.it, std::move(point_ptr));
  }

  // removes the element at position pos from the list and returns (and converts) it.
  // "increments" the iterator passed as argument: it now point to the next element.

  // The implementation differs, depending on whether ReturnType includes SimHashes
  // (in which case we do not recompute them)
  template <class dummy = ReturnType, TEMPL_RESTRICT_DECL2(mystd::negation<Has_BitApprox<dummy>>)>
  auto true_pop_point(Iterator &pos) -> ReturnType;

  template <class dummy = ReturnType, TEMPL_RESTRICT_DECL2(Has_BitApprox<dummy>)>
  auto true_pop_point(Iterator &pos) -> ReturnType;

  // erase, sort, size, empty follow std::list's semantics
  // Note that sort and size may not be available / differ for the multithreaded case.
  Iterator erase(Iterator pos) { return actual_list.erase(pos.it); }
  void sort() { actual_list.sort(); }
  typename UnderlyingContainer::size_type size() const noexcept { return actual_list.size(); }
  NODISCARD bool empty() const noexcept { return actual_list.empty(); }

private:
  // clang-format off
  StaticInitializer<StoredPoint> const init_stored_point;
  StaticInitializer<ReturnType>  const init_return_type;
  UnderlyingContainer actual_list;
  // clang-format on
};

/**
  Iterator for the list above. Note that some essential extra functionality is implemented as member
  functions of the iterators. Notably, we access the sim_hashes by calling functions on the
  iterators directly (without dereferencing)
  Alternatives would be exposing STNode to the user or dereferencing returning a proxy object
  with get_sim_hashes functions etc. Unfortunately, such proxy objects would probably return / hold
  references, which does not work very well with "auto" and get tricky with std::atomic.
  ( The issue is that std::atomic<foo> converts to foo, but this does not use relaxed memory_order.
    This would mean that we either have to return by value or references to atomics, exposing the
    atomics to the user).
  Since we want to keep the interfaces for single- and multi-threaded similar, we chose this weird
  interface.
  Note that this iterator is a const_iterator in the C++ sense.

  Note: This may or may not be an iterator in the sense of the C++ concept.
        The issue is that we might deprecate dereferencing...
  TODO: If we allow dereferencing, set iterator_traits for compatibility with std:: / range-based
        loops etc.
*/

// clang-format off
template<class SieveTraits>
class GaussIteratorBitApprox<SieveTraits, false>
{
friend GaussListWithBitApprox<SieveTraits, false>;

private:
  using ListType                = GaussListWithBitApprox<SieveTraits, false>;

  using StoredPoint             = typename SieveTraits::GaussList_StoredPoint;
  using ReturnType              = typename SieveTraits::GaussList_ReturnType;
  using CoordinateSelectionUsed = typename SieveTraits::CoordinateSelectionUsed;

  using SimHashBlock            = typename CoordinateSelectionUsed::SimHashBlock;
  using SimHashes               = typename CoordinateSelectionUsed::SimHashes;

  // (UnderlyingIterator is only used internally (by the Iterator class in conversions)
  using UnderlyingIterator      = typename ListType::UnderlyingContainer::iterator;
  using CUnderlyingIterator     = typename ListType::UnderlyingContainer::const_iterator;
  // clang-format on

private:
  CUnderlyingIterator it;

public:
  // clang-format off
  GaussIteratorBitApprox()                                          = delete;
  GaussIteratorBitApprox(GaussIteratorBitApprox const &)            = default;
  GaussIteratorBitApprox(GaussIteratorBitApprox &&)                 = default;
  GaussIteratorBitApprox &operator=(GaussIteratorBitApprox const &) = default;
  GaussIteratorBitApprox &operator=(GaussIteratorBitApprox &&)      = default;
  // clang-format on

  // we can convert from a "plain" iterator to the underlying list. Only used internally or by
  // the list class.
private:
  // clang-format off
  constexpr GaussIteratorBitApprox( UnderlyingIterator const &new_it) noexcept : it(new_it) {}
  constexpr GaussIteratorBitApprox(CUnderlyingIterator const &new_it) noexcept : it(new_it) {}
  // clang-format on

public:
  // comparison, needed for for-loops (i.e. compare against list.cend() )
  bool operator==(GaussIteratorBitApprox const &other) const { return it == other.it; }
  bool operator!=(GaussIteratorBitApprox const &other) const { return it != other.it; }

  // increment operators. We have NO decrement, these are forward iterators.
  // (The latter is for interface compatibility with multi-threaded, where that restriction makes
  // everything both faster and easier.)
  // clang-format off
  GaussIteratorBitApprox &operator++()    { ++it; return *this; }  // prefix version
  GaussIteratorBitApprox  operator++(int) { return it++; }         // postfix version
  // clang-format on

  // TODO: Consider adding bool is_end() member function rather than comparing
  //       with cend() (in multithreading, cend() might mutate during an iteration,
  //       which can lead to unexpected issue.)

  // obtains an approximation to norm2 of the pointee.
  SimHashApproxNorm2 get_approx_norm2() const { return it->approx_norm2; }

  // obtains all bitapproximations of the pointee.
  SimHashes get_all_bitapproximations() const { return it->bit_approximations; }

  // gives access to the level'th block the bitapproximations.
  // In the MT case, we might return a reference to atomic.
  // NOTE: the check_simhash_scalar_product function works with *any* type that exposes an
  // access_bitapproximation(unsigned int) function (i.e. some lattice points and this class)
  // this way, we can call these directly on the iterator.
  SimHashBlock const &access_bitapproximation(unsigned int level) const
  {
    return it->bit_approximations[level];
  }

  // derefencing gives us the point stored in the STNodes (i.e. without sim_hashes)
  StoredPoint const &operator*() const { return *(it->ptr_to_exact); }

  // Note that for overloads of the -> operator, the return type is a class, whose -> operator is
  // recursively called (recursion stops at plain pointers).
  StoredPoint const *operator->() const { return it->ptr_to_exact; }

  // Iterator is explicitly convertible to pointer to const-Lattice Point.
  explicit operator StoredPoint const *() const { return it->ptr_to_exact; }
};

/**
  turn_maybe_iterator_to_point(PointOrIterator &&arg) will dereference iterators and do nothing on
  lattice points. This (weird) function is / was used to unify some code which can either take an
  iterator or a point.
*/

namespace Helpers
{
// helper classes to work around the lack of C++17 constexpr if / unwillingness to resort to arcane
// rules how overloads / template argument deduction / *partial* specializations mix...

// general template is for lattice points : get does nothing
// clang-format off
template <class PointOrIterator>
struct ConvertIteratorToPoint_Helper
{
  static_assert(IsALatticePoint<PointOrIterator>::value == true, "");
  using RetType = PointOrIterator;
  FORCE_INLINE static inline constexpr PointOrIterator const &get(PointOrIterator const &arg)
  {
    return arg;
  }
};

// partial specialization is for iterators : get dereferences
// Note: clang actually warns about compliance with FORCE_INLINE here.
// (It might depend on the template arguments; needs further investigation)
template <class SieveTraits, bool MT>
struct ConvertIteratorToPoint_Helper<GaussIteratorBitApprox<SieveTraits, MT>>
{
  using RetType = decltype(   *std::declval< GaussIteratorBitApprox<SieveTraits,MT> >()   );
  FORCE_INLINE static inline RetType const &get(GaussIteratorBitApprox<SieveTraits, MT> const &arg)
  {
    return *arg;
  }
};
}  // end Helper namespace
// clang-format on

// clang-format off
template<class PointOrIterator>
FORCE_INLINE auto CPP14CONSTEXPR turn_maybe_iterator_to_point(PointOrIterator &&arg)
    -> typename Helpers::ConvertIteratorToPoint_Helper<mystd::decay_t<PointOrIterator>>::RetType const &
{
  return Helpers::ConvertIteratorToPoint_Helper<mystd::decay_t<PointOrIterator>>::get(
      std::forward<PointOrIterator>(arg));
}
// clang-format on

}  // end namespace GaussSieve

#include "GaussListBitapprox_impl.h"

#endif  // include guards
