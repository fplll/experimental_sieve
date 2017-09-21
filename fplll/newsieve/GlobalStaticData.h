#ifndef GLOBAL_STATIC_DATA_H
#define GLOBAL_STATIC_DATA_H

#include "DefaultIncludes.h"

// This struct holds all data that is used to initialize the static data associated to our various classes.
// Note that it need not be a singleton itself

namespace GaussSieve{

CREATE_MEMBER_TYPEDEF_CHECK_CLASS(StaticDataInitializerTag,IsStaticDataInitializer);

template<class DimensionType>
struct StaticDataInitializer
{
  using StaticDataInitializerTag = std::true_type;
  DimensionType const dim;
  unsigned int const dim_int;

  constexpr StaticDataInitializer(DimensionType const &new_dim) : dim(dim), dim_int(dim) {}
};




} // namespace

#endif
