#ifndef APPROXIMATED_POINT_H
#define APPROXIMATED_POINT_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include <tuple>
#include "Lazy.h"
#include "GlobalStaticData.h"
#include "PlainLatticePoint.h"

namespace GaussSieve{

template<class ELP, class Approximation> class VectorWithApproximation; //ELP = exact lattice point
template<class ELP, class Approximation> class ScalarWithApproximation;

/**
  This class stores a scalar with an approximation to that scalar.
  Currently, this class only allows constant objects. We have no arithmetic.
*/
template<class ELP, class Approximation>
struct ScalarWithApproximation
//  :public LazyEval::ObjectWithApproximation<typename Get_ScalarProductStorageType<ELP>::type,typename Approximation::ScalarProductType>
{
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  public:
  using ExactType  = ExactScalarType;
  using ApproxType = ApproxScalarType;
  static constexpr unsigned int ApproxLevel = ApproxLevelOf<ELP>::value + 1;
  static_assert(!(std::is_same<ExactType,ApproxType>::value),"Can not approximate by itself currently");

  private:
  ExactType  exact_scalar;
  ApproxType approx_scalar;
  public:
  constexpr      explicit ScalarWithApproximation(ExactType const& exact,ApproxType const &approx)
    :exact_scalar(exact),approx_scalar(approx){}
  CPP14CONSTEXPR explicit ScalarWithApproximation(ExactType     && exact,ApproxType const &approx)
    :exact_scalar(std::move(exact)),approx_scalar(approx){}
  CPP14CONSTEXPR explicit ScalarWithApproximation(ExactType const& exact,ApproxType      &&approx)
    :exact_scalar(exact),approx_scalar(std::move(approx)){}
  CPP14CONSTEXPR explicit ScalarWithApproximation(ExactType     && exact,ApproxType      &&approx)
    :exact_scalar(std::move(exact)),approx_scalar(std::move(approx)){}
  constexpr      explicit ScalarWithApproximation(ExactType const &exact)
    :exact_scalar(exact), approx_scalar(exact) {}
  CPP14CONSTEXPR explicit ScalarWithApproximation(ExactType && exact)
    :exact_scalar(std::move(exact)), approx_scalar(exact_scalar) {} // Note: Ordering matters here

  // default constructors / move assignment should be done automatically.

  constexpr      explicit operator ExactType()  const & { return exact_scalar;}
  CPP14CONSTEXPR explicit operator ExactType()  &&      { return std::move(exact_scalar);}
  constexpr      explicit operator ApproxType() const & { return approx_scalar;}
  CPP14CONSTEXPR explicit operator ApproxType() &&      { return std::move(approx_scalar);}

  // In case of unknown conversion, evaluate
  template<class T>
  CPP14CONSTEXPR explicit operator T() const & { return static_cast<T>(exact_scalar); }
  template<class T>
  CPP14CONSTEXPR explicit operator T() &&      { return static_cast<T>(std::move(exact_scalar));}

  // Accessing the exact / approximate types. Used by the Lazy Evaluation class
  // (notably, the leaves of the expression trees).
  constexpr      ExactType  const & access_exact()  const { return exact_scalar; }
  CPP14CONSTEXPR ExactType        & access_exact()        { return exact_scalar; }
  constexpr      ApproxType const & access_approx() const { return approx_scalar; }
  CPP14CONSTEXPR ApproxType       & access_approx()       { return approx_scalar; }
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
  using Trait_ScalarProductStorageType = typename Get_ScalarProductStorageType<ELP>::type;
  using Trait_ScalarProductStorageType_Full  = ScalarWithApproximation<ELP,Approximation>;
  using Trait_CoordinateType          = typename Get_CoordinateType<ELP>::type;
  using Trait_AbsoluteCoos            = typename Get_AbsoluteCooType<ELP>::type;
  using Trait_RepCooType              = typename Get_RepCooType<ELP>::type;
  using Trait_ExposesCoos             = NormalizeTrait<Has_ExposesCoos<ELP>>;
  using Trait_Coos_RW                 = NormalizeTrait<Has_Coos_RW<ELP>>;
  using Trait_ExposesInternalRep      = NormalizeTrait<Has_ExposesInternalRep<ELP>>;
  using Trait_InternalRepLinear       = NormalizeTrait<Has_InternalRepLinear<ELP>>;
  using Trait_InternalRep_RW          = NormalizeTrait<Has_InternalRep_RW<ELP>>;
  using Trait_InternalRepByCoos       = NormalizeTrait<Has_InternalRepByCoos<ELP>>;
  using Trait_InternalRepIsAbsolute   = NormalizeTrait<Has_InternalRepIsAbsolute<ELP>>;
  using Trait_CheapNorm2              = NormalizeTrait<Has_CheapNorm2<ELP>>;
  using Trait_CheapNegate             = NormalizeTrait<Has_CheapNegate<ELP>>;
  using Trait_BitApprox               = NormalizeTrait<Has_BitApprox<ELP>>;

  using Trait_Approximations       = std::true_type;
  using Trait_DelayedScalarProduct = std::true_type;
  static constexpr unsigned int ApproxLevel = ApproxLevelOf<ELP>::value + 1;
  using Trait_ApproxLevel          = std::integral_constant<unsigned int, ApproxLevel>;
};

template<class ELP, class Approximation, class Function, class... Args>
class DelayedScalar
  : LazyEval::SieveLazyEval<Function,Args...>
{
  using Parent = LazyEval::SieveLazyEval<Function,Args...>;
  using ExactEvalType = typename Parent::ExactEvalType;
  using ApproxEvalType = typename Parent::ApproxEvalType;
  public:
  // forward all arguments to Parent constructor
  template<class... ConstructorArgs>
  DelayedScalar(ConstructorArgs &&... constructor_args)
    : Parent(std::forward<ConstructorArgs>(constructor_args)... ) {}
  operator ExactEvalType() {return this->eval_exact(); }
};

template<class ELP, class Approximation, class Function, class... Args>
class DelayedVector
  : LazyEval::SieveLazyEval<Function,Args...>
{
  using Parent = LazyEval::SieveLazyEval<Function, Args...>;
};

template<class ELP, class Approximation>
using Get_DelayedScalarProductTypeCR = DelayedScalar
<
  ELP, Approximation,
  LazyEval::Lazy_ScalarProduct<ELP,Approximation>, //Function
  LazyEval::LazyWrapCombinedCR<VectorWithApproximation<ELP,Approximation>>, // LHS
  LazyEval::LazyWrapCombinedCR<VectorWithApproximation<ELP,Approximation>>  // RHS
>;

template<class ELP, class Approximation>
using Get_DelayedNorm2TypeCR = DelayedScalar
<
  ELP, Approximation,
  LazyEval::Lazy_Norm2<ELP,Approximation>, //Function
  LazyEval::LazyWrapCombinedCR<VectorWithApproximation<ELP,Approximation>>  // Argument
>;

// clang-format on

template<class ELP, class Approximation>
class VectorWithApproximation
: public GeneralLatticePoint<VectorWithApproximation<ELP,Approximation>>
{
  static_assert(Has_CheapNorm2<ELP>::value, "ELP does not store its norm.");
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");

  private:
  ELP exact_point;
  Approximation approx;


  public:
  using LatticePointTag = std::true_type;
  using Myself = VectorWithApproximation<ELP,Approximation>;

  using ExactCoos = typename Get_CoordinateType<ELP>::type; // may be void
  using RepCooType = typename Get_RepCooType<ELP>::type;
//  using AbsoluteCooType = typename Get_AbsoluteCooType<ELP>::type;
  using typename GeneralLatticePoint<VectorWithApproximation<ELP,Approximation>>::ScalarProductStorageType;
//
  using ExactScalarProductType    = typename Get_ScalarProductStorageType<ELP>::type;
//  using ApproxScalarProductType   = typename Approximation::ScalarProductType;
  using CombinedScalarProductType = ScalarWithApproximation<ELP,Approximation>;
//  using DelayedScalarProductType  = Get_DelayedScalarProductType<ELP,Approximation>;
//  // Think about this:
//  using DelayedNorm2Type          = Get_DelayedNorm2Type<ELP,Approximation>;
  static constexpr unsigned int ApproxLevel = ApproxLevelOf<ELP>::value + 1;
//
  VectorWithApproximation(VectorWithApproximation const &old) = delete;
  VectorWithApproximation(VectorWithApproximation && old) = default;
  VectorWithApproximation & operator= (VectorWithApproximation const & other) = delete;
  VectorWithApproximation & operator= (VectorWithApproximation && other) = default;
  explicit VectorWithApproximation(ELP && new_exact_point)
    : exact_point(std::move(new_exact_point)), approx(exact_point) {};
  VectorWithApproximation(ELP const & new_exact_point) = delete;

  // construct with precomputed approximation:
  template<class Arg, TEMPL_RESTRICT_DECL2(std::is_same<Approximation, mystd::decay_t<Arg> >)>
  explicit VectorWithApproximation(ELP && new_exact_point, Arg && new_approx)
    : exact_point(std::move(new_exact_point)), approx(std::forward<Arg>(new_approx)) {}

  template<class ET, int nfixed>
  VectorWithApproximation(PlainLatticePoint<ET,nfixed> const &) = delete;
  template<class ET, int nfixed>
  VectorWithApproximation(PlainLatticePoint<ET,nfixed> && plain_point)
    : VectorWithApproximation(static_cast<ELP>(std::move(plain_point))){}


// Implement ObjectWithApproximation's interface as far as meaningful
  using ExactType  = ELP;
  using ApproxType = Approximation;

                 explicit operator ExactType() const & = delete; // would copy point
  CPP14CONSTEXPR explicit operator ExactType()  &&      { return std::move(exact_point);}
  constexpr      explicit operator ApproxType() const & { return approx;}
  CPP14CONSTEXPR explicit operator ApproxType() &&      { return std::move(approx);}

  constexpr      ExactType  const & access_exact()  const { return exact_point; }
  CPP14CONSTEXPR ExactType        & access_exact()        { return exact_point; }
  constexpr      ApproxType const & access_approx() const { return approx; }
  CPP14CONSTEXPR ApproxType       & access_approx()       { return approx; }

  static std::string class_name() { return ELP::class_name() + "with approximation"; } // TODO:class_name for Approximation
//
//  // operators<,>,<=, >= : No overloads. Defaults to exact comparison.
//
//  // forward [] to exact class
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesCoos<T>)>
  ExactCoos &operator[](Arg &&arg) { return exact_point[std::forward<Arg>(arg)]; }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesCoos<T>)>
  ExactCoos const &operator[](Arg &&arg) const { return exact_point[std::forward<Arg>(arg)]; }
//
//  // +=, -=, *= and unary- just forward to the exact class and recompute the approximation.
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  VectorWithApproximation& operator+=(LatP2 const &x2) { exact_point+=x2; recompute_approx(); return *this; }
//
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  VectorWithApproximation& operator-=(LatP2 const &x2) { exact_point-=x2; recompute_approx(); return *this; }
//
  template<class Multiplier>
  VectorWithApproximation& operator*=(Multiplier &&x2) { exact_point*=std::forward<Multiplier>(x2); recompute_approx(); return *this; }
//
//  // TODO: Inefficient
  VectorWithApproximation operator-()&&
  {
    return static_cast<VectorWithApproximation>(-std::move(exact_point));
  }
//
//  // equality comparison.
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
  template<class T=ELP, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>)>
  auto get_internal_rep_size() const -> decltype( std::declval<ELP>().get_internal_rep_size() ) { return exact_point.get_internal_rep_size(); }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>)>
  RepCooType const & get_internal_rep(Arg &&arg) const { return exact_point.get_internal_rep(std::forward<Arg>(arg)); }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>, Has_InternalRep_RW<T>)>
  RepCooType & get_internal_rep(Arg &&arg) {return exact_point.get_internal_rep(std::forward<Arg>(arg));}

  // forward absolute coos to exact point
  template<class Arg> auto get_absolute_coo(Arg &&arg) const
    -> decltype( std::declval<ELP>().get_absolute_coo( std::declval<Arg &&>() ))
  { return exact_point.get_absolute_coo(std::forward<Arg>(arg));  }

  // forward get_dim
  auto get_dim() const -> decltype( std::declval<ELP>().get_dim() ) { return exact_point.get_dim(); }

  // write_lp_to_stream outputs the point, then the approximation
  inline std::ostream& write_lp_to_stream(std::ostream &os, bool const include_norm2=true, bool const include_approx =true) const
  {
    exact_point.write_lp_to_stream(os, include_norm2,include_approx);
    if(include_approx)
    {
      os << approx;
    }
    return os;
  }

  // write_lp_rep_to_stream only outputs the exact point. The approximation does not contribute
  // to the internal representation.
  template<class T=ELP, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<ELP>)>
  inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const { return exact_point.write_lp_rep_to_stream(os); }

//
//  //TODO: read_from_stream
//


  void fill_with_zero() { exact_point.fill_with_zero(); recompute_approx(); } // may optimize
  void make_negative()  { exact_point.make_negative(); recompute_approx(); } // may optimize
  bool is_zero() const { return exact_point.is_zero(); }

  // TODO: Copy approximation
  VectorWithApproximation make_copy() const { return VectorWithApproximation(exact_point.make_copy()); }

// TODO: More efficient sanitize routines that take prior knowledge into account.
  void sanitize() { exact_point.sanitize(); recompute_approx(); }
  void sanitize(ExactScalarProductType const &norm2) { exact_point.sanitize(norm2); recompute_approx(); }
  void recompute_approx() { approx = static_cast<Approximation>(exact_point); }

  // bypass 1 level of approximations
  inline auto get_norm2_exact() const -> decltype( std::declval<ELP>().get_norm2() )
  { return exact_point.get_norm2(); }

  // bypass all levels of approximations
  inline ExactScalarProductType get_norm2_exact_recursive() const { return exact_point.get_norm2_exact_recursive(); }

  // return both
  inline CombinedScalarProductType get_norm2_full() const
  { return CombinedScalarProductType{exact_point.get_norm2_exact(), approx.get_approx_norm2()}; }

  // Delayed evaluation(!)
  inline Get_DelayedNorm2TypeCR<ELP,Approximation> get_norm2() const & // Note: The & is important
  {
    return Get_DelayedNorm2TypeCR<ELP,Approximation> { LazyEval::LazyWrapCombinedCR<Myself> {*this } };
  }

  // scalar products:

  inline auto do_compute_sc_product_exact(Myself const &x2) const
  -> decltype( std::declval<ELP>().do_compute_sc_product(std::declval<ELP>() ) )
  {
    return exact_point.do_compute_sc_product(x2.exact_point);
  }

  inline ExactScalarProductType do_compute_sc_product_exact_recursive(Myself const &x2) const
  {
    return exact_point.do_compute_sc_product_exact_recursive(x2.exact_point);
  }

  inline CombinedScalarProductType do_compute_sc_product_full(Myself const &x2) const
  {
    return CombinedScalarProductType
    { // constructor args
      exact_point.do_compute_sc_product_exact(x2.exact_point),
      compute_sc_product_approx(approx,x2.approx)
    };
  }

  inline Get_DelayedScalarProductTypeCR<ELP,Approximation> do_compute_sc_product(Myself const &x2) const &
  {
    return Get_DelayedScalarProductTypeCR<ELP,Approximation>
    { // constructor args
      LazyEval::LazyWrapCombinedCR<Myself> {*this},
      LazyEval::LazyWrapCombinedCR<Myself> {x2}
    };
  }

  inline Get_DelayedScalarProductTypeCR<ELP,Approximation> do_compute_sc_product(Myself const &&) const & = delete;

  inline auto get_bitapprox_norm2() const -> decltype( std::declval<ELP>().get_bitapprox_norm2() )
  {
    return exact_point.get_bitapprox_norm2();
  }

  inline int do_compute_sc_product_bitapprox(Myself const & x2) const
  {
    return exact_point.do_compute_sc_product_bitapprox(x2.exact_point);
  }

  inline int do_compute_sc_product_bitapprox_2nd_order(Myself const & x2) const
  {
    return exact_point.do_compute_sc_product_bitapprox_2nd_order(x2.exact_point);
  }

};

/**
  Initializes Static Data for the combination vector (by forwarding to the individual components)
  Note: Initializing scalars is the job of the initializers of the vectors.
*/
template<class ELP, class Approximation>
class StaticInitializer<VectorWithApproximation<ELP,Approximation>>
  final : public DefaultStaticInitializer<VectorWithApproximation<ELP,Approximation>>
{
  public:
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  StaticInitializer<ExactVectorType>  const init_exact_vector;
  StaticInitializer<ApproxVectorType> const init_approx_vector;

  template<class X,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<mystd::decay_t<X>>)>
  explicit StaticInitializer(X &&init_arg) :
    init_exact_vector(std::forward<X>(init_arg)), init_approx_vector(std::forward<X>(init_arg)){}
};

} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
