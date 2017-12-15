#ifndef DEFAULT_TERM_CONDS_IMPL_H
#define DEFAULT_TERM_CONDS_IMPL_H

#include "DefaultIncludes.h"
#include "DefaultTermConds.h"

namespace GaussSieve
{

template <class SieveTraits, bool MT>
inline int LengthTerminationCondition<SieveTraits, MT>::check()
{
  return (this->sieveptr->get_best_length2() <= target_norm2) ? 1 : 0;
}

template <class SieveTraits, bool MT>
inline int LengthTerminationCondition<SieveTraits, MT>::check_vec(
    typename SieveTraits::TermCond_QueryType const &lattice_point)
{
  return (lattice_point.get_norm2() <= target_norm2) ? 1 : 0;
}

template <class SieveTraits, bool MT>
inline int MinkowskiTerminationCondition<SieveTraits, MT>::check()
{
  return (sieveptr->get_best_length2() <= target_norm2) ? 1 : 0;
}

template <class SieveTraits, bool MT>
inline int MinkowskiTerminationCondition<SieveTraits, MT>::check_vec(
    typename SieveTraits::TermCond_QueryType const &lattice_point)
{
  return (lattice_point.get_norm2() <= target_norm2) ? 1 : 0;
}

template <class SieveTraits, bool MT>
inline void MinkowskiTerminationCondition<SieveTraits, MT>::custom_init()
{
  // clang-format off
  target_norm2 = ConvertMaybeMPZ<LengthType>::convert_to_inttype( sieveptr->get_basis().get_minkowski_bound() );
  // clang-format on
}

}  // namespace GaussSieve

#endif
