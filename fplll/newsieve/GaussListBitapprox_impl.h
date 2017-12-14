#ifndef POINT_LIST_WITH_BITAPPROX_IMPL_H
#define POINT_LIST_WITH_BITAPPROX_IMPL_H

// clang-format off
#ifndef POINT_LIST_WITH_BITAPPROX_H
  #error wrong usage
#endif
// clang-format on

namespace GaussSieve
{

// Removes the lattice point indicated by the iterator pos from the list and returns the lattice
// point.
// Modifies its argument: pos is incremented (so it keeps pointing to something valid)
// variant for lattice points that do not include a SimHash
template <class SieveTraits>
template <class dummy, TEMPL_RESTRICT_IMPL2(mystd::negation<Has_BitApprox<dummy>>)>
auto GaussListWithBitApprox<SieveTraits, false>::true_pop_point(Iterator &pos) -> ReturnType
{
  // retrieve the value before we erase
  ReturnType retval = static_cast<ReturnType>(std::move(*(pos.it->ptr_to_exact)));
  // erases the point and increments the iterator
  pos = erase(pos);
  return retval;  // Uses named-return-value-optimization to avoid copying.
}

// variant for lattice point classes that include a SimHash
template <class SieveTraits>
template <class dummy, TEMPL_RESTRICT_IMPL2(Has_BitApprox<dummy>)>
auto GaussListWithBitApprox<SieveTraits, false>::true_pop_point(Iterator &pos) -> ReturnType
{
  // clang-format off
  // SimHashArgTag is a tag used to select the correct constructor.
  ReturnType retval = ReturnType{ SimHashArgTag{},
                                  std::move(*(pos.it->ptr_to_exact)),
                                  std::move(pos.it->bit_approximations) };
  // clang-format on
  pos = erase(pos);
  return retval;
}

}  // end namespace GaussSieve

#endif  // include guards
