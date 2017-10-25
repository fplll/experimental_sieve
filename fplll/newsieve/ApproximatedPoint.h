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
  :public LazyEval::ObjectWithApproximation<typename Get_ScalarProductStorageType<ELP>::type,typename Approximation::ScalarProductType>
{
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
//
//  constexpr explicit ScalarWithApproximation(ExactScalarType const &exact, ApproxScalarType const &approx)
//    :exact_sc_product(exact), approx_sc_product(approx) {};
//  constexpr explicit ScalarWithApproximation(ExactScalarType const &exact)
//    :exact_sc_product(exact), approx_sc_product(static_cast<ApproxScalarType>(exact)){};
//
//  constexpr operator ExactScalarType() const { return exact_sc_product;}
//  constexpr explicit operator ApproxScalarType() const { return approx_sc_product; }
//
//  // const-ness restriction is for debug purposes, mostly.
//  ExactScalarType const  exact_sc_product;
//  ApproxScalarType const approx_sc_product;
//
//  // initialize from the result of a lazy computation.
//  template<class LazyFunction>
//  constexpr explicit ScalarWithApproximation(LazyEval::SieveLazyEval<LazyFunction> const &lazy_fun )
//    :exact_sc_product(lazy_fun.eval_exact()),approx_sc_product(lazy_fun.eval_approx())
//    {
//      static_assert(LazyEval::SieveLazyEval<LazyFunction>::scalar_or_vector == LazyEval::ScalarOrVector::scalar_type,"Trying to assign a vector to a scalar");
//    }
//
//  // initialize from wrappers
//  constexpr explicit ScalarWithApproximation(LazyEval::LazyWrapExactScalar<ELP,Approximation> const &wrapper)
//    :exact_sc_product(wrapper.eval_exact()),approx_sc_product(wrapper.eval_approx()) {}
//  constexpr explicit ScalarWithApproximation(LazyEval::LazyWrapExactAndApproxScalar<ELP,Approximation> const &wrapper)
//    :exact_sc_product(wrapper.eval_exact()),approx_sc_product(wrapper.eval_approx()) {}
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

  using Trait_Approximations       = std::true_type;
  using Trait_DelayedScalarProduct = std::true_type;
};

template<class ELP, class Approximation, class Function, class... Args>
class DelayedScalar
  : LazyEval::SieveLazyEval<Function, Args...>
{

};

template<class ELP, class Approximation, class Function, class... Args>
class DelayedVector
  : LazyEval::SieveLazyEval<Function, Args...>
{

};

template<class ELP, class Approximation>
using Get_DelayedScalarProductType = DelayedScalar
<
  ELP, Approximation,
  LazyEval::Lazy_ScalarProduct<ELP,Approximation>, //Function
  LazyEval::LazyWrapCombinedCR<VectorWithApproximation<ELP,Approximation>>, // LHS
  LazyEval::LazyWrapCombinedCR<VectorWithApproximation<ELP,Approximation>>  // RHS
>;

template<class ELP, class Approximation>
using Get_DelayedNorm2Type = DelayedScalar
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
  using ExactCoos = typename Get_CoordinateType<ELP>::type; // may be void
  using RepCooType = typename Get_RepCooType<ELP>::type;
  using AbsoluteCooType = typename Get_AbsoluteCooType<ELP>::type;
//  using ScalarProductStorageType = typename Get_ScalarProductStorageType<ELP>::type;
  using typename GeneralLatticePoint<VectorWithApproximation<ELP,Approximation>>::ScalarProductStorageType;
//
  using ExactScalarProductType    = typename Get_ScalarProductStorageType<ELP>::type;
  using ApproxScalarProductType   = typename Approximation::ScalarProductType;
  using CombinedScalarProductType = ScalarWithApproximation<ELP,Approximation>;
  using DelayedScalarProductType  = Get_DelayedScalarProductType<ELP,Approximation>;
//  // Think about this:
  using DelayedNorm2Type          = Get_DelayedNorm2Type<ELP,Approximation>;
//
  VectorWithApproximation(VectorWithApproximation const &old) = delete;
  VectorWithApproximation(VectorWithApproximation && old) = default;
  VectorWithApproximation & operator= (VectorWithApproximation const & other) = delete;
  VectorWithApproximation & operator= (VectorWithApproximation && other) = default;
  explicit VectorWithApproximation(ELP && new_exact_point)
    : exact_point(std::move(new_exact_point)), approx(exact_point) {};

  // construct with precomputed approximation:
  template<class Arg, TEMPL_RESTRICT_DECL2(std::is_same<Approximation, typename std::decay<Arg>::type>)>
  explicit VectorWithApproximation(ELP && new_exact_point, Arg && new_approx)
    : exact_point(std::move(new_exact_point)), approx(std::forward<Arg>(new_approx)) {}

// Implement ObjectWithApproximation's interface as far as meaningful
  using ExactType  = ELP;
  using ApproxType = Approximation;

  //constexpr      explicit operator ExactType()  const & { return exact_point;}
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
//
//  // forward internal_rep to exact point if it exists.
  template<class T=ELP, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>)>
  auto get_internal_rep_size() const -> decltype( std::declval<ELP>().get_internal_rep_size() ) { return exact_point.get_internal_rep_size(); }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>)>
  RepCooType const & get_internal_rep(Arg &&arg) const { return exact_point.get_internal_rep(std::forward<Arg>(arg)); }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>, Has_InternalRep_RW<T>)>
  RepCooType & get_internal_rep(Arg &&arg) {return exact_point.get_internal_rep(std::forward<Arg>(arg));}
//
//  // forward absolute coos to exact point
  template<class Arg>
  AbsoluteCooType get_absolute_coo(Arg &&arg) const { return exact_point.get_absolute_coo(std::forward<Arg>(arg));  }
//
//  // forward get_dim
  auto get_dim() const -> decltype( std::declval<ELP>().get_dim() ) { return exact_point.get_dim(); }
//
//
  inline std::ostream& write_lp_to_stream(std::ostream &os, bool const include_norm2=true, bool const include_approx =true) const
  {
    exact_point.write_lp_to_stream(os, include_norm2,include_approx);
    if(include_approx)
    {
      os << approx;
    }
    return os;
  }
//
  template<class T=ELP, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<ELP>)>
  inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const { return exact_point.write_lp_rep_to_stream(os); }
//
//
//  //TODO: read_from_stream
//
  void fill_with_zero() { exact_point.fill_with_zero(); recompute_approx(); }
  void make_negative()  { exact_point.make_negative(); recompute_approx(); } // todo : may optimize
  bool is_zero() { return exact_point.is_zero(); }

  // TODO: Copy approximation
  VectorWithApproximation make_copy() const { return VectorWithApproximation(exact_point.make_copy()); }

// TODO: More efficient sanitize routines that take prior knowledge into account.
  void sanitize() { exact_point.sanitize(); recompute_approx(); }
  void sanitize(ExactScalarProductType const &norm2) { exact_point.sanitize(norm2); recompute_approx(); }
  void recompute_approx() { approx = static_cast<Approximation>(exact_point); }
//
//  inline DelayedNorm2Type get_norm2() const
//  {
////    return DelayedNorm2Type ( exact_point.get_norm2(), approx.get_approx_norm2()  );
//    return DelayedNorm2Type(  );
//  }
//
//  ExactScalarProductType get_norm2_exact() const {return exact_point.get_norm2_exact(); }
//
//  CombinedScalarProductType get_norm2_full() const
//  {
//    return CombinedScalarProductType(exact_point.get_norm2_exact(),approx.get_approx_norm2()  );
//  }
//
//  DelayedScalarProductType do_compute_sc_product(VectorWithApproximation const &x2) { };
//  ExactScalarProductType   do_compute_sc_product_exact(VectorWithApproximation const &x2)
//  {
//    return exact_point.do_compute_sc_product_exact(x2);
//  };
//
};

/**
  Initializes Static Data for the combination vector (by forwarding to the individual components)
  Note: Initializing scalars is the job of the initializers of the vectors.
*/
template<class ELP, class Approximation>
class StaticInitializer<VectorWithApproximation<ELP,Approximation>>
  final : public DefaultStaticInitializer<VectorWithApproximation<ELP,Approximation>>
{
  BRING_TYPES_INTO_SCOPE_Lazy_GetTypes(ELP,Approximation);
  StaticInitializer<ExactVectorType>  const init_exact_vector;
  StaticInitializer<ApproxVectorType> const init_approx_vector;

  template<class X,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<typename std::decay<X>::type>)>
  explicit StaticInitializer(X &&init_arg) :
    init_exact_vector(std::forward<X>(init_arg)), init_approx_vector(std::forward<X>(init_arg)){}
};



} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
