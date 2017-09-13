#ifndef APPROXIMATED_POINT_H
#define APPROXIMATED_POINT_H

#include "DebugAll.h"
#include <type_traits>
#include "SieveUtility.h"
#include "LatticePointConcept.h"

namespace GaussSieve{

template<class ELP, class Approximation> class PointWithApproximation; //ELP = exact lattice point
//template<class ELP, class Approximation> class ScalarProductWithApproximation;
template<class ELP, class Approximation> class DelayedScProduct;

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
  using CoordinateType          = typename GetCooType<ELP>::type;
  using HasDelayedScProd        = std::true_type;

  using AuxDataType             = typename GetAuxDataType<ELP>::type; // for now. This whole AuxDataType needs to be redesigned.
  using ScalarProductReturnType = typename GetScPType<ELP>::type; // Requires thinking about the concept class.
};

template<class ELP, class Approximation>
struct ScalarProductWithApproximation
{
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  using ExactScProdType   = typename GetScPType<ELP>::type;
  using ApproxScProdType  = typename Approximation::ScalarProductType;

  explicit ScalarProductWithApproximation(ExactScProdType const &exact, ApproxScProdType const &approx)
    :exact_sc_product(exact), approx_sc_product(approx) {};

//  ScalarProductWithApproximation(ScalarProductWithApproximation const &) = default;
//  ScalarProductWithApproximation(ScalarProductWithApproximation &&) = default;
//  ScalarProductWithApproximation& operator=(ScalarProductWithApproximation const &) = default;
//  ScalarProductWithApproximation& operator=(ScalarProductWithApproximation &&) = default;


  ExactScProdType  exact_sc_product;
  ApproxScProdType approx_sc_product;

// TODO:
//  using ApproxType     = TODO;

};



template<class ELP, class Approximation>
class PointWithApproximation: public GeneralLatticePoint<PointWithApproximation<ELP,Approximation>>
{
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  static_assert(std::is_same<typename GetAuxDataType<ELP>::type, typename Approximation::AuxDataType>::value,"AuxDataType must be the same");
  public:
  using LatticePointTag = std::true_type;
  using ExactCoos = typename GetCooType<ELP>::type; // may be void
//  using AuxDataType = typename GetAuxDataType<ELP>::type;
//  using ScalarProductReturnType = typename GetScPType<ELP>::type;
  using typename GeneralLatticePoint<PointWithApproximation<ELP,Approximation>>::ScalarProductReturnType;
  using typename GeneralLatticePoint<PointWithApproximation<ELP,Approximation>>::AuxDataType;


  PointWithApproximation(PointWithApproximation const &old) = delete;
  PointWithApproximation(PointWithApproximation && old) = default;
  PointWithApproximation & operator= (PointWithApproximation const & other) = delete;
  PointWithApproximation & operator= (PointWithApproximation && other) = default;
  explicit PointWithApproximation(ELP && new_exact_point)
    : exact_point(std::move(new_exact_point)), approx(exact_point) {};

  static std::string class_name() { return ELP::class_name() + "with approximation"; }; // TODO:class_name for Approximation

  // operators<,>,<=, >= : No overload for now. Defaults to exact comparison.

  ExactCoos &operator[](uint_fast16_t idx) { return exact_point[idx]; };
  ExactCoos const &operator[](uint_fast16_t idx) const { return exact_point[idx]; };

  template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value)>
  PointWithApproximation& operator+=(LatP2 const &x2) { exact_point+=x2; recompute_approx(); return *this;};


  auto get_dim() -> decltype( std::declval<ELP>().get_dim() ) { return exact_point.get_dim(); }

  void sanitize() { exact_point.sanitize(); approx=static_cast<Approximation>(exact_point); }


  private:

  ELP exact_point;
  Approximation approx;

};



} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
