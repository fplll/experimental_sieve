#ifndef APPROXIMATED_POINT_H
#define APPROXIMATED_POINT_H

#include "DebugAll.h"
#include <type_traits>
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include <iostream>
#include <utility>
#include <tuple>
#include "Lazy.h"

namespace GaussSieve{

template<class ELP, class Approximation> class PointWithApproximation; //ELP = exact lattice point


// Alias : Scalar + approx is identified with a constant function.

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
  using HasDelayedScalarProduct = std::true_type;

  using AuxDataType             = typename GetAuxDataType<ELP>::type; // for now. This whole AuxDataType needs to be redesigned.
  using ScalarProductReturnType = typename GetScPType<ELP>::type; // Requires thinking about the concept class.
};

/*

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

  operator ELP() const { return exact_sc_product;}


  ExactScProdType const  exact_sc_product;
  ApproxScProdType const approx_sc_product;

// TODO:
//  using ApproxType     = TODO;

};

*/



template<class ELP, class Approximation>
class PointWithApproximation: public GeneralLatticePoint<PointWithApproximation<ELP,Approximation>>
{
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  static_assert(std::is_same<typename GetAuxDataType<ELP>::type, typename Approximation::AuxDataType>::value,"AuxDataType must be the same");
//  friend DelayedScProduct<ELP,Approximation>;
  public:
  using LatticePointTag = std::true_type;
  using ExactCoos = typename GetCooType<ELP>::type; // may be void
//  using AuxDataType = typename GetAuxDataType<ELP>::type;
//  using ScalarProductReturnType = typename GetScPType<ELP>::type;
  using typename GeneralLatticePoint<PointWithApproximation<ELP,Approximation>>::ScalarProductReturnType;
  using typename GeneralLatticePoint<PointWithApproximation<ELP,Approximation>>::AuxDataType;

  using ExactScalarProductType    = typename GetScPType<ELP>::type;
  using ApproxScalarProductType   = typename Approximation::ScalarProductType;
  using CombinedScalarProductType = ScalarWithApproximation<ELP,Approximation>;
  using DelayedScalarProductType  = DelayedScalarProduct<ELP,Approximation>;



  PointWithApproximation(PointWithApproximation const &old) = delete;
  PointWithApproximation(PointWithApproximation && old) = default;
  PointWithApproximation & operator= (PointWithApproximation const & other) = delete;
  PointWithApproximation & operator= (PointWithApproximation && other) = default;
  explicit PointWithApproximation(ELP && new_exact_point)
    : exact_point(std::move(new_exact_point)), approx(exact_point) {};

  static std::string class_name() { return ELP::class_name() + "with approximation"; } // TODO:class_name for Approximation

  // operators<,>,<=, >= : No overload for now. Defaults to exact comparison.

  ExactCoos &operator[](uint_fast16_t idx) { return exact_point[idx]; }
  ExactCoos const &operator[](uint_fast16_t idx) const { return exact_point[idx]; }

  template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value)>
  PointWithApproximation& operator+=(LatP2 const &x2) { exact_point+=x2; recompute_approx(); return *this; }

  template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value)>
  PointWithApproximation& operator-=(LatP2 const &x2) { exact_point-=x2; recompute_approx(); return *this; }

  template<class Integer>
  PointWithApproximation& operator*=(Integer const &x2) { exact_point*=x2; recompute_approx(); return *this; }

  PointWithApproximation operator-()&& { return static_cast<PointWithApproximation>(-std::move(exact_point)); }

  template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value)>
  bool operator==(LatP2 const &x2) const { return exact_point== x2;};
  bool operator==(PointWithApproximation const &x2) const
  {
    if(approx!=x2.approx)
    {
      return false;
    }
    else
    {
      return exact_point ==x2.exact_point;
    }
  }

  auto get_dim() const -> decltype( std::declval<ELP>().get_dim() ) { return exact_point.get_dim(); }
  auto get_vec_size() const -> decltype( std::declval<ELP>().get_vec_size() ) { return exact_point.get_vec_size(); }

  std::ostream& write_to_stream(std::ostream &os, bool const include_norm2=true) const
  {
    return os << exact_point << approx;
  }

  //TODO: read_from_stream

  void fill_with_zero() { exact_point.fill_with_zero(); recompute_approx(); }
  void make_negative()  { exact_point.make_negative(); recompute_approx(); }
  bool is_zero() { return exact_point.is_zero(); }
  PointWithApproximation make_copy() const { return PointWithApproximation(exact_point.make_copy()); }

  void sanitize() { exact_point.sanitize(); recompute_approx(); }
  void sanitize(ExactScalarProductType const &norm2) { exact_point.sanitize(norm2); recompute_approx(); }
  void recompute_approx() { approx = static_cast<Approximation>(exact_point); }

  CombinedScalarProductType get_norm2() const { return exact_point.get_norm2_exact(); }
  ExactScalarProductType get_norm2_exact() const {return exact_point.get_norm2_exact(); }

  DelayedScalarProductType do_compute_sc_product(PointWithApproximation const &x2) { };
  ExactScalarProductType   do_compute_sc_product_exact(PointWithApproximation const &x2)
  {
    return exact_point.do_compute_sc_product_exact(x2);
  };

  private:

  ELP exact_point;
  Approximation approx;
};





} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
