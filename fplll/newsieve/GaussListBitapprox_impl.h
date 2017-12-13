#ifndef POINT_LIST_WITH_BITAPPROX_IMPL_H
#define POINT_LIST_WITH_BITAPPROX_IMPL_H

#ifndef POINT_LIST_WITH_BITAPPROX_H
  #error wrong usage
#endif

namespace GaussSieve
{

template<class SieveTraits>
template<class dummy, TEMPL_RESTRICT_IMPL2(mystd::negation<Has_BitApprox<dummy>>)>
auto GaussListWithBitApprox<SieveTraits,false>::true_pop_point(Iterator &pos) -> ReturnType
{
  ReturnType retval = static_cast<ReturnType>(std::move(*(pos.it->ptr_to_exact))) ;
  pos = erase(pos);
  return retval;  // (mandatory) named-return-value-optimization to avoid copying.
}

template<class SieveTraits>
template<class dummy, TEMPL_RESTRICT_IMPL2(Has_BitApprox<dummy>)>
auto GaussListWithBitApprox<SieveTraits,false>::true_pop_point(Iterator &pos) -> ReturnType
{
  ReturnType retval = ReturnType{ SimHashArgTag{},
                                  std::move(*(pos.it->ptr_to_exact)),
                                  std::move(pos.it->bit_approximations) };
  pos = erase(pos);
  return retval;
}

}  // end namespace GaussSieve

#endif // include guards
