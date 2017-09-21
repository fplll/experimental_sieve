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

template<class ELP, class Approximation> class VectorWithApproximation; //ELP = exact lattice point
template<class ELP, class Approximation> class ScalarWithApproximation;

template<class ELP, class Approximation>
struct ScalarWithApproximation
{
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation)
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");

  explicit ScalarWithApproximation(ExactScalarType const &exact, ApproxScalarType const &approx)
    :exact_sc_product(exact), approx_sc_product(approx) {};
  explicit ScalarWithApproximation(ExactScProdt)
    :exact_sc_product(exact), approx_sc_product(static_cast<ApproxScalarType>(exact)){};


//  ScalarWithApproximation(ScalarWithApproximation const &) = default;
//  ScalarWithApproximation(ScalarWithApproximation &&) = default;
//  ScalarWithApproximation& operator=(ScalarWithApproximation const &) = default;
//  ScalarWithApproximation& operator=(ScalarWithApproximation &&) = default;

  operator ExactVectorType() const { return exact_sc_product;}
  explicit operator ApproxVectorType() const { return approx_sc_product; }

  ExactScalarType const  exact_sc_product;
  ApproxScalarType const approx_sc_product;

// TODO:
//  using ApproxType     = TODO;

};


template<class ELP, class Approximation>
class LatticePointTraits< VectorWithApproximation <ELP, Approximation> >
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
//  using ScalarProductStorageType = ET;
  using CoordinateType          = typename GetCooType<ELP>::type;
  using HasDelayedScalarProduct = std::true_type;

  using AuxDataType             = typename GetAuxDataType<ELP>::type; // for now. This whole AuxDataType needs to be redesigned.
  using ScalarProductStorageType = typename GetScalarProductStorageType<ELP>::type; // Requires thinking about the concept class.
};


template<class ELP, class Approximation>
class VectorWithApproximation: public GeneralLatticePoint<VectorWithApproximation<ELP,Approximation>>
{
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  static_assert(std::is_same<typename GetAuxDataType<ELP>::type, typename Approximation::AuxDataType>::value,"AuxDataType must be the same");
//  friend DelayedScProduct<ELP,Approximation>;
  public:
  using LatticePointTag = std::true_type;
  using ExactCoos = typename GetCooType<ELP>::type; // may be void
//  using AuxDataType = typename GetAuxDataType<ELP>::type;
//  using ScalarProductStorageType = typename GetScalarProductStorageType<ELP>::type;
  using typename GeneralLatticePoint<VectorWithApproximation<ELP,Approximation>>::ScalarProductStorageType;
  using typename GeneralLatticePoint<VectorWithApproximation<ELP,Approximation>>::AuxDataType;

  using ExactScalarProductType    = typename GetScalarProductStorageType<ELP>::type;
  using ApproxScalarProductType   = typename Approximation::ScalarProductType;
  using CombinedScalarProductType = ScalarWithApproximation<ELP,Approximation>;
  using DelayedScalarProductType  = DelayedScalarProduct<ELP,Approximation>;



  VectorWithApproximation(VectorWithApproximation const &old) = delete;
  VectorWithApproximation(VectorWithApproximation && old) = default;
  VectorWithApproximation & operator= (VectorWithApproximation const & other) = delete;
  VectorWithApproximation & operator= (VectorWithApproximation && other) = default;
  explicit VectorWithApproximation(ELP && new_exact_point)
    : exact_point(std::move(new_exact_point)), approx(exact_point) {};

  static std::string class_name() { return ELP::class_name() + "with approximation"; } // TODO:class_name for Approximation

  // operators<,>,<=, >= : No overload for now. Defaults to exact comparison.

  ExactCoos &operator[](uint_fast16_t idx) { return exact_point[idx]; }
  ExactCoos const &operator[](uint_fast16_t idx) const { return exact_point[idx]; }

  template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value)>
  VectorWithApproximation& operator+=(LatP2 const &x2) { exact_point+=x2; recompute_approx(); return *this; }

  template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value)>
  VectorWithApproximation& operator-=(LatP2 const &x2) { exact_point-=x2; recompute_approx(); return *this; }

  template<class Integer>
  VectorWithApproximation& operator*=(Integer const &x2) { exact_point*=x2; recompute_approx(); return *this; }

  VectorWithApproximation operator-()&& { return static_cast<VectorWithApproximation>(-std::move(exact_point)); }

  template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value)>
  bool operator==(LatP2 const &x2) const { return exact_point== x2;};
  bool operator==(VectorWithApproximation const &x2) const
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
  VectorWithApproximation make_copy() const { return VectorWithApproximation(exact_point.make_copy()); }

  void sanitize() { exact_point.sanitize(); recompute_approx(); }
  void sanitize(ExactScalarProductType const &norm2) { exact_point.sanitize(norm2); recompute_approx(); }
  void recompute_approx() { approx = static_cast<Approximation>(exact_point); }

  CombinedScalarProductType get_norm2() const { return exact_point.get_norm2_exact(); }
  ExactScalarProductType get_norm2_exact() const {return exact_point.get_norm2_exact(); }

  DelayedScalarProductType do_compute_sc_product(VectorWithApproximation const &x2) { };
  ExactScalarProductType   do_compute_sc_product_exact(VectorWithApproximation const &x2)
  {
    return exact_point.do_compute_sc_product_exact(x2);
  };

  private:

  ELP exact_point;
  Approximation approx;
};





} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
