#ifndef APPROXIMATED_POINT_H
#define APPROXIMATED_POINT_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include <tuple>
#include "Lazy.h"
#include "GlobalStaticData.h"

namespace GaussSieve{

template<class ELP, class Approximation> class VectorWithApproximation; //ELP = exact lattice point
template<class ELP, class Approximation> class ScalarWithApproximation;

/**
  This class stores a scalar with an approximation to that scalar.
  Currently, this class only allows constant objects. We have no arithmetic.
*/
template<class ELP, class Approximation>
struct ScalarWithApproximation
{
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");

  constexpr explicit ScalarWithApproximation(ExactScalarType const &exact, ApproxScalarType const &approx)
    :exact_sc_product(exact), approx_sc_product(approx) {};
  constexpr explicit ScalarWithApproximation(ExactScalarType const &exact)
    :exact_sc_product(exact), approx_sc_product(static_cast<ApproxScalarType>(exact)){};

  constexpr operator ExactVectorType() const { return exact_sc_product;}
  constexpr explicit operator ApproxVectorType() const { return approx_sc_product; }

  // const-ness restriction is for debug purposes, mostly.
  ExactScalarType const  exact_sc_product;
  ApproxScalarType const approx_sc_product;
  template<class LazyFunction>
  constexpr explicit ScalarWithApproximation(LazyEval::SieveLazyEval<LazyFunction> const &lazy_fun )
    :exact_sc_product(lazy_fun.eval_exact()),approx_sc_product(lazy_fun.eval_approx())
    {
      static_assert(LazyEval::SieveLazyEval<LazyFunction>::scalar_or_vector == LazyEval::ScalarOrVector::scalar_type,"Trying to assign a vector to a scalar");
    }
  constexpr explicit ScalarWithApproximation(LazyEval::LazyWrapExactScalar<ELP,Approximation> const &wrapper)
    :exact_sc_product(wrapper.eval_exact()),approx_sc_product(wrapper.eval_approx()) {}
  constexpr explicit ScalarWithApproximation(LazyEval::LazyWrapExactAndApproxScalar<ELP,Approximation> const &wrapper)
    :exact_sc_product(wrapper.eval_exact()),approx_sc_product(wrapper.eval_approx()) {}
};

/**
  Initializes Static Data for the combination (by forwarding to the individual components)
*/
template<class ELP, class Approximation>
class StaticInitializer<ScalarWithApproximation<ELP,Approximation>>
final : public DefaultStaticInitializer<ScalarWithApproximation<ELP,Approximation>>
{
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  StaticInitializer<ExactScalarType> const init_exact_scalar;
  StaticInitializer<ApproxScalarType> const init_approx_scalar;

  template<class X,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<typename std::decay<X>::type>)>
  explicit StaticInitializer(X &&init_arg) :
    init_exact_scalar(std::forward<X>(init_arg)), init_approx_scalar(std::forward<X>(init_arg)){}
};

/**
  VectorWithApproximation<ELP,Approximation> is a lattice point class that is
  made from combining ELP with the approximation.
*/

/**
  The lattice point traits for VectorWithApproximation<ELP,Approximation> are deduced from
  the traits of ELP.
*/

template<class ELP, class Approximation>
class LatticePointTraits< VectorWithApproximation <ELP, Approximation> >
{
static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
public:
// forwarding traits from ELP
  using ScalarProductStorageType = typename GetScalarProductStorageType<ELP>::type;
  using ScalarProductStorageType_Full  = void; //TODO!
  using CoordinateType          = typename GetCooType<ELP>::type;
  using AbsoluteCoos            = typename GetAbsoluteCooType<ELP>::type;
  using RepCooType              = typename GetRepCooType<ELP>::type;
  using ExposesCoos             = NormalizeTrait<DoesExposeCoos<ELP>>;
  using ExposesInternalRep      = NormalizeTrait<HasInternalRep<ELP>>;
  using InternalRepVector       = NormalizeTrait<IsRepLinear<ELP>>;
  using InternalRep_RW          = NormalizeTrait<IsRepRW<ELP>>;
  using InternalRepByCoos       = NormalizeTrait<HasRepByCoos<ELP>>;
  using InternalRepIsAbsolute   = NormalizeTrait<HasAbsoluteRep<ELP>>;
  using CheapNorm2              = NormalizeTrait<IsNorm2Cheap<ELP>>;
  using CheapNegate             = NormalizeTrait<IsNegateCheap<ELP>>;

  using HasApproximations       = std::true_type;
  using HasDelayedScalarProduct = std::true_type;
};

// clang-format off
template<class ELP, class Approximation>
using DelayedScalarProduct = LazyEval::SieveLazyEval<
    LazyEval::Lazy_ScalarProduct // encapsulated function
  < ELP,Approximation,
    LazyEval::LazyWrapExactAndApproxVector<ELP,Approximation>, // type of LHS function arg
    LazyEval::LazyWrapExactAndApproxVector<ELP,Approximation>  // type of RHS function arg
  > >;
// clang-format on

template<class ELP, class Approximation>
class VectorWithApproximation: public GeneralLatticePoint<VectorWithApproximation<ELP,Approximation>>
{
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  public:
  using LatticePointTag = std::true_type;
  using ExactCoos = typename GetCooType<ELP>::type; // may be void
  using RepCooType = typename GetRepCooType<ELP>::type;
  using AbsoluteCooType = typename GetAbsoluteCooType<ELP>::type;
//  using ScalarProductStorageType = typename GetScalarProductStorageType<ELP>::type;
  using typename GeneralLatticePoint<VectorWithApproximation<ELP,Approximation>>::ScalarProductStorageType;

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

  // operators<,>,<=, >= : No overloads. Defaults to exact comparison.

  // forward [] to exact class
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(DoesExposeCoos<T>)>
  ExactCoos &operator[](Arg &&arg) { return exact_point[std::forward<Arg>(arg)]; }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(DoesExposeCoos<T>)>
  ExactCoos const &operator[](Arg &&arg) const { return exact_point[std::forward<Arg>(arg)]; }

  // +=, -=, *= and unary- just forward to the exact class and recompute the approximation.
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  VectorWithApproximation& operator+=(LatP2 const &x2) { exact_point+=x2; recompute_approx(); return *this; }

  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  VectorWithApproximation& operator-=(LatP2 const &x2) { exact_point-=x2; recompute_approx(); return *this; }

  template<class Integer>
  VectorWithApproximation& operator*=(Integer &&x2) { exact_point*=std::forward<Integer>(x2); recompute_approx(); return *this; }

  // TODO: Inefficient
  VectorWithApproximation operator-()&&
  {
    return static_cast<VectorWithApproximation>(-std::move(exact_point));
  }

  // equality comparison.
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  bool operator==(LatP2 const &x2) const { return exact_point == x2;};
  bool operator==(VectorWithApproximation const &x2) const
  {
    if(approx!=x2.approx)
    {
      return false;
    }
    else
    {
      return exact_point == x2.exact_point;
    }
  }

  // forward internal_rep to exact point if it exists.
  template<class T=ELP, TEMPL_RESTRICT_DECL2(HasInternalRep<T>)>
  auto get_internal_rep_size() const -> decltype( std::declval<ELP>().get_internal_rep_size() ) { return exact_point.get_internal_rep_size(); }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(HasInternalRep<T>)>
  RepCooType const & get_internal_rep(Arg &&arg) const { return exact_point.get_internal_rep(std::forward<Arg>(arg)); }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(HasInternalRep<T>, IsRepRW<T>)>
  RepCooType & get_internal_rep(Arg &&arg) {return exact_point.get_internal_rep(std::forward<Arg>(arg));}

  // forward absolute coos to exact point
  template<class Arg>
  AbsoluteCooType get_absolute_coo(Arg &&arg) const { return exact_point.get_absolute_coo(std::forward<Arg>(arg));  }

  // forward get_dim
  auto get_dim() const -> decltype( std::declval<ELP>().get_dim() ) { return exact_point.get_dim(); }


  inline std::ostream& write_lp_to_stream(std::ostream &os, bool const include_norm2=true, bool const include_approx =true) const
  {
    exact_point.write_lp_to_stream(os, include_norm2,include_approx);
    if(include_approx)
    {
      os << approx;
    }
    return os;
  }

  template<class T=ELP, TEMPL_RESTRICT_DECL2(HasInternalRep<ELP>)>
  inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const { return exact_point.write_lp_rep_to_stream(os); }


  //TODO: read_from_stream

  void fill_with_zero() { exact_point.fill_with_zero(); recompute_approx(); }
  void make_negative()  { exact_point.make_negative(); recompute_approx(); } // todo : may optimize
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
