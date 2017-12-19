#ifndef GAUSS_VECTOR_BITAPPROX_H
#define GAUSS_VECTOR_BITAPPROX_H

#include "DefaultIncludes.h"
#include "GaussListBitapprox.h"

namespace GaussSieve
{

template<class SieveTraits>
class GaussVectorWithBitApprox<SieveTraits, false>
{
private:
  using StoredPoint             = typename SieveTraits::GaussList_StoredPoint;
  using ReturnType              = typename SieveTraits::GaussList_ReturnType;
  using CoordinateSelectionUsed = typename SieveTraits::CoordinateSelectionUsed;

  using SimHashBlock            = typename CoordinateSelectionUsed::SimHashBlock;
  using SimHashes               = typename CoordinateSelectionUsed::SimHashes;

  using Iterator                = GaussIteratorBitApproxForVector<SieveTraits, false>;  // custom iterator class below
  // clang-format on

  friend Iterator;
  // The class is essentially just a wrapper around UnderlyingContainer
  using UnderlyingContainer         = std::vector<GaussListDetails::STNode<SieveTraits>>;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;

public:
  // clang-format off
  GaussVectorWithBitApprox()                                          = delete;
  GaussVectorWithBitApprox(GaussVectorWithBitApprox const &)            = delete;
  GaussVectorWithBitApprox(GaussVectorWithBitApprox &&)                 = delete;
  GaussVectorWithBitApprox &operator=(GaussVectorWithBitApprox const &) = delete;
  GaussVectorWithBitApprox &operator=(GaussVectorWithBitApprox &&)      = delete;
  // clang-format on

  // Sole constructor. The argument is used to initialize the static data of the used lattice point
  // classes.
  // NOTE / TODO: Making this noexcept means we have to catch reinitializations in the caller.
  explicit GaussVectorWithBitApprox(GlobalStaticDataInitializer const &static_data) noexcept
      : init_stored_point(static_data),
        init_return_type(static_data),
        actual_vector()
  {
  }

  // behaves like cbegin, cend from STL containers, i.e. gives const-iterator to begin/end.
  // clang-format off
  CPP14CONSTEXPR Iterator begin() const noexcept { return actual_vector.begin(); }
  CPP14CONSTEXPR Iterator end()   const noexcept { return actual_vector.end();   }
  // clang-format on

  // insert_before(pos, point) inserts the point just before pos.
  // the return value is an iterator to the newly inserted point.
  //  Note: LatticePoint && is a forwarding/universal reference, NOT a rvalue reference. Still, the
  //        function is only supposed to be called on rvalues, so we static_asserts that.
  //        This means you might have to use insert_before(pos, std::move(new_lp));
  /*
  template <class LatticePoint, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatticePoint>>)>
  Iterator insert_before(Iterator const &pos, LatticePoint &&new_point)
  {
    static_assert(std::is_reference<LatticePoint>::value == false, "Must call on rvalues");
    return actual_list.emplace(pos.it, std::move(new_point));
  }
  */

  template <class LatticePoint, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatticePoint>>)>
  void emplace_back(LatticePoint &&new_point)
  {
    static_assert(std::is_reference<LatticePoint>::value == false, "Must call on rvalues");
    actual_vector.emplace_back(std::move(new_point));
  }

  // Inserts the lattice point pointed to by a pointer. The list takes ownership of the pointee.
  // TODO: Consider using unique_ptr
  // Untested
  /*
  Iterator insert_before(Iterator const &pos, StoredPoint *&&point_ptr)
  {
    return actual_list.emplace(pos.it, std::move(point_ptr));
  }
  */

  // removes the element at position pos from the list and returns (and converts) it.
  // "increments" the iterator passed as argument: it now point to the next element.

  // The implementation differs, depending on whether ReturnType includes SimHashes
  // (in which case we do not recompute them)
  /*
  template <class dummy = ReturnType, TEMPL_RESTRICT_DECL2(mystd::negation<Has_BitApprox<dummy>>)>
  auto true_pop_point(Iterator &pos) -> ReturnType;

  template <class dummy = ReturnType, TEMPL_RESTRICT_DECL2(Has_BitApprox<dummy>)>
  auto true_pop_point(Iterator &pos) -> ReturnType;
  */

  // erase, sort, size, empty follow std::list's semantics
  // Note that sort and size may not be available / differ for the multithreaded case.
  /*
  Iterator erase(Iterator pos) { return actual_list.erase(pos.it); }
  */
  void sort() { actual_vector.sort(); }
  typename UnderlyingContainer::size_type size() const noexcept { return actual_vector.size(); }
  NODISCARD bool empty() const noexcept { return actual_vector.empty(); }

private:
  // clang-format off
  StaticInitializer<StoredPoint> const init_stored_point;
  StaticInitializer<ReturnType>  const init_return_type;
  UnderlyingContainer actual_vector;
  // clang-format on
};


// clang-format off
template<class SieveTraits>
class GaussIteratorBitApproxForVector<SieveTraits, false>
{
friend GaussVectorWithBitApprox<SieveTraits, false>;

private:
  using ListType                = GaussVectorWithBitApprox<SieveTraits, false>;

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
  UnderlyingIterator it;

public:
  // clang-format off
  GaussIteratorBitApproxForVector()                                          = delete;
  GaussIteratorBitApproxForVector(GaussIteratorBitApproxForVector const &)            = default;
  GaussIteratorBitApproxForVector(GaussIteratorBitApproxForVector &&)                 = default;
  GaussIteratorBitApproxForVector &operator=(GaussIteratorBitApproxForVector const &) = default;
  GaussIteratorBitApproxForVector &operator=(GaussIteratorBitApproxForVector &&)      = default;
  // clang-format on

  // we can convert from a "plain" iterator to the underlying list. Only used internally or by
  // the list class.
private:
  // clang-format off
  constexpr GaussIteratorBitApproxForVector( UnderlyingIterator const &new_it) noexcept : it(new_it) {}
 // constexpr GaussIteratorBitApprox(CUnderlyingIterator const &new_it) noexcept : it(new_it) {}
  // clang-format on

public:
  // comparison, needed for for-loops (i.e. compare against list.cend() )
  bool operator==(GaussIteratorBitApproxForVector const &other) const { return it == other.it; }
  bool operator!=(GaussIteratorBitApproxForVector const &other) const { return it != other.it; }

  // increment operators. We have NO decrement, these are forward iterators.
  // (The latter is for interface compatibility with multi-threaded, where that restriction makes
  // everything both faster and easier.)
  // clang-format off
  GaussIteratorBitApproxForVector &operator++()    { ++it; return *this; }  // prefix version
  GaussIteratorBitApproxForVector  operator++(int) { return it++; }         // postfix version
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
  StoredPoint &operator*() { return *(it->ptr_to_exact); }

  // Note that for overloads of the -> operator, the return type is a class, whose -> operator is
  // recursively called (recursion stops at plain pointers).
  StoredPoint *operator->() { return it->ptr_to_exact; }

  // Iterator is explicitly convertible to pointer to const-Lattice Point.
  explicit operator StoredPoint *() const { return it->ptr_to_exact; }
};



}  // end namespace GaussSieve


#endif // include guard
