#ifndef APPROXIMATED_POINT_H
#define APPROXIMATED_POINT_H

#include "DebugAll.h"
#include <type_traits>
#include "SieveUtility.h"
#include "LatticePointConcept.h"

namespace GaussSieve{

template<class ELP, class Approximation> class PointWithApproximation; //ELP = exact lattice point
template<class ELP, class Approximation> class ScalarProductWithApproximation;

template<class ELP, class Approximation>
class LatticePointTraits< PointWithApproximation <ELP, Approximation> >
{
static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
public:
// forwarding traits from ELP
  using CoordinateVector        = typename IsCooVector<ELP>::value_t;
  using CoordinateAccess        = typename HasCoos<ELP>::value_t;
  using AbsoluteCoos            = typename CoosAreAbsolute<ELP>::value_t;
  using CheapNorm2              = typename IsNorm2Cheap<ELP>::value_t;
  using CheapNegate             = typename IsNegateCheap<ELP>::value_t;
  using HasApproximations       = std::true_type;
//  using AuxDataType             = MaybeFixed<nfixed>;
//  using ScalarProductReturnType = ET;
//  using CoordinateType          = ET;

};

template<class ELP, class Approximation>
class ScalarProductWithApproximation
{
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");

  using ExactValueType = typename GetScPType<ELP>::type;
// TODO:
//  using ApproxType     = TODO;
  public:

  private:



};



template<class ELP, class Approximation>
class PointWithApproximation: public GeneralLatticePoint<PointWithApproximation<ELP,Approximation>>;
{
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  public:
  using LatticePointTag = std::true_type;
  using ExactCoos = typename GetCooType<ELP>::type; // may be void

  PointWithApproximation(PointWithApproximation const &old) = delete;
  PointWithApproximation(PointWithApproximation && old) = default;
  PointWithApproximation & operator= (PointWithApproximation const & other) = delete;
  PointWithApproximation & operator= (PointWithApproximation && other) = default;

  template<class ELP2 = ELP, typename std::enable_if<HasCoos<ELP2>::value,int>::type = 0>
  ExactCoos &operator[](uint_fast16_t idx) { return ; };

  template<class ELP2 = ELP, typename std::enable_if<HasCoos<ELP2>::value,int>::type = 0>
  ExactCoos const &operator[](uint_fast16_t idx) const { return exact_point[idx]; };

  static std::string class_name() { return ELP::class_name() + "with approximation"; }; // TODO:class_name for Approximation
  private:

  ELP exact_point;
  Approximation approx;

};



} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
