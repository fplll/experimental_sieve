#ifndef DEFAULT_TERM_CONDS_IMPL_H
#define DEFAULT_TERM_CONDS_IMPL_H

#include "DefaultIncludes.h"
#include "DefaultTermConds.h"

namespace GaussSieve
{

template <class SieveTraits, bool MT>
inline int LengthTerminationCondition<SieveTraits, MT>::check(Sieve<SieveTraits, MT> *const sieve)
{
  return (sieve->get_best_length2() <= target_length) ? 1 : 0;
}

template <class SieveTraits, bool MT>
inline int
MinkowskiTerminationCondition<SieveTraits, MT>::check(Sieve<SieveTraits, MT> *const sieve)
{
  return (sieve->get_best_length2() <= target_length) ? 1 : 0;
}

// clang-format off
template <class SieveTraits, bool MT>
inline void MinkowskiTerminationCondition<SieveTraits, MT>::init(Sieve<SieveTraits, MT> *const sieve)
{
  target_length = ConvertMaybeMPZ<LengthType>::convert_to_inttype( sieve->get_basis().get_minkowski_bound() );
}
// clang-format on

}  // namespace GaussSieve

#endif
