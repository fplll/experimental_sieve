#ifndef APPROXIMATED_POINT_H
#define APPROXIMATED_POINT_H

#include "DebugAll.h"
#include <type_traits>
#include "SieveUtility.h"

namespace GaussSieve{

template<class InnerLatticePoint, class Approximation> class PointWithApproximation;

template<class InnerLatticePoint, class Approximation>
class LatticePointTraits< PointWithApproximation <InnerLatticePoint, Approximation> >
{
public:
  using AuxDataType             = MaybeFixed<nfixed>;
  using ScalarProductReturnType = ET;
  using CoordinateVector        = std::true_type;
  using CoordinateAccess        = std::true_type;
  using AbsoluteCoos            = std::true_type;
  using CoordinateType          = ET;
  using CheapNorm2              = std::true_type;
  using CheapNegate             = std::true_type;
};



} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
