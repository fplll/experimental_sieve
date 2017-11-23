#ifndef APPROXIMATED_POINT_H
#define APPROXIMATED_POINT_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include <tuple>
#include "Lazy.h"
#include "GlobalStaticData.h"
#include "PlainLatticePoint.h"


#define T_IS_ELP static_assert(std::is_same<T,ELP>::value,"Wrong template argument")

// forward declarations:

namespace GaussSieve{

template<class Scalar> struct MakeLeveledScalar;
template<class OldLevels, class NextLevel> struct AddApproximation;
template<class ELP> class MakeLeveledLatticePoint;
template<class ELP, class ApproxLP> class AddApproximationToLatP;
namespace Helpers
{
template<bool> struct AddApproximation_Helper;
}
}
namespace GaussSieve{

/**
  This turns a non-leveled (scalar) object into a leveled object (cf. Lazy.h for what this means).
  There is no functionality preserved besides storage, but we have implicit conversion from/to
  the non-leveled object.
*/

template<class Scalar>
struct MakeLeveledScalar
{
  static_assert(Has_LeveledObject<Scalar>::value==false,"");
  static_assert(ApproxLevelOf<Scalar>::value==0,"");
  using LeveledComparison = std::true_type; // Unsure about this one.
  using LeveledObject = std::true_type;
  using LeveledObject_Base = std::true_type;
  static constexpr unsigned int ApproxLevel = 0;
  template<unsigned int level> using ObjectAtLevel = mystd::conditional_t<level==0, Scalar, void>;
  Scalar scalar;
  constexpr      MakeLeveledScalar(Scalar const &sc) : scalar(sc) {}
  CPP14CONSTEXPR MakeLeveledScalar(Scalar &&sc) : scalar(std::move(sc)){}
  // default move, copy etc.
  constexpr operator Scalar() const & {return scalar;}
  constexpr operator Scalar()      && {return std::move(scalar);}

  template<unsigned int level>
  constexpr Scalar const & access() const { static_assert(level==0,""); return scalar;}
  template<unsigned int level>
            Scalar       & access()       { static_assert(level==0,""); return scalar;}
  // TODO : RValue variants

  template<unsigned int level>
  constexpr Scalar const & get_value_at_level() const { static_assert(level==0,""); return scalar;}
  template<unsigned int level>
            Scalar       & get_value_at_level()       { static_assert(level==0,""); return scalar;}
};

/**
  This adds another level to a leveled (scalar) object.
  Notably, OldLevels has to be a leveled object with n>=0 levels.
  This then results in a new leveled objects with n+1 levels, with NextLevel at the highest level.
*/

template<class OldLevels, class NextLevel>
struct AddApproximation
{
  static_assert(Has_LeveledObject<OldLevels>::value==true,"");
  static_assert(Has_LeveledObject<NextLevel>::value==false,"");
  using LeveledComparison = std::true_type;
  using LeveledObject = std::true_type;
  using LeveledObject_Base = std::true_type;
  static constexpr unsigned int ApproxLevel = ApproxLevelOf<OldLevels>::value + 1;
  static_assert(std::is_same<typename OldLevels::template ObjectAtLevel<ApproxLevel-1>,NextLevel >::value==false,"");
  template<unsigned int level> using ObjectAtLevel = mystd::conditional_t
    < level==ApproxLevel, NextLevel, typename OldLevels::template ObjectAtLevel<level>  >;
  using BaseScalar = ObjectAtLevel<0>;
  template<unsigned int level> using Helper = Helpers::AddApproximation_Helper<level==ApproxLevel>;
  template<bool> friend class Helpers::AddApproximation_Helper;

  private:
  OldLevels old_levels;
  NextLevel next_level;
  public:

  constexpr      explicit AddApproximation(OldLevels const &olds, NextLevel const &news)
    :old_levels(olds), next_level(news) {}
  CPP14CONSTEXPR explicit AddApproximation(OldLevels      &&olds, NextLevel const &news)
    :old_levels(std::move(olds)),next_level(news) {}
  CPP14CONSTEXPR explicit AddApproximation(OldLevels const &olds, NextLevel      &&news)
    :old_levels(olds), next_level(std::move(news)) {}
  CPP14CONSTEXPR explicit AddApproximation(OldLevels      &&olds, NextLevel      &&news)
    :old_levels(std::move(olds)),next_level(std::move(news)){}
  constexpr      explicit AddApproximation(OldLevels const &olds)
    :old_levels(olds), next_level(static_cast<BaseScalar>(olds)){}
  CPP14CONSTEXPR explicit AddApproximation(OldLevels      &&olds)
    :old_levels(std::move(olds)), next_level(static_cast<BaseScalar>(old_levels)) {}
  constexpr               AddApproximation(BaseScalar const &old_base)
    :old_levels(old_base), next_level(old_base) {}
  CPP14CONSTEXPR          AddApproximation(BaseScalar      &&old_base)
    :old_levels(std::move(old_base)),next_level(static_cast<BaseScalar>(old_levels)) {}

  constexpr      operator BaseScalar() const & { return static_cast<BaseScalar>(old_levels); }
  CPP14CONSTEXPR operator BaseScalar()      && { return static_cast<BaseScalar>(std::move(old_levels));}

  template<unsigned int level>
  constexpr ObjectAtLevel<level> const & access() const
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }
  template<unsigned int level>
  CPP14CONSTEXPR ObjectAtLevel<level> &  access()
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }
  template<unsigned int level>
  constexpr ObjectAtLevel<level> const & get_value_at_level() const
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }
  template<unsigned int level>
  CPP14CONSTEXPR ObjectAtLevel<level> &  get_value_at_level()
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }
  // TODO: Rvalue variants.
};

namespace Helpers{
template<bool AtMaxLevel> struct AddApproximation_Helper;

template<>
struct AddApproximation_Helper<true>
{
  template<unsigned int level, class Argument>
  static inline constexpr auto get(Argument const &ob)
  -> decltype(std::declval<Argument>().next_level) const &
  {
    static_assert(level==Argument::ApproxLevel,"");
    return ob.next_level;
  }
  template<unsigned int level, class Argument>
  static inline auto get(Argument & ob)
  -> decltype(std::declval<Argument>().next_level) &
  {
    static_assert(level==Argument::ApproxLevel,"");
    return ob.next_level;
  }
};

template<>
struct AddApproximation_Helper<false>
{
  template<unsigned int level, class Argument>
  static inline constexpr auto get(Argument const &ob)
  -> typename decltype(std::declval<Argument>().old_levels)::template ObjectAtLevel<level> const &
  {
    static_assert(level < Argument::ApproxLevel,"");
    return ob.old_levels.template access<level>();
  }
  template<unsigned int level, class Argument>
  static inline auto get(Argument &ob)
  -> typename decltype(std::declval<Argument>().old_levels)::template ObjectAtLevel<level> &
  {
    static_assert(level < Argument::ApproxLevel,"");
    return ob.old_levels.template access<level>();
  }
};
} // end namespace Helpers

/******************
 This turns a lattice vector into a leveled lattice point.
*******************/


// Forward traits:
template<class ELP>
class LatticePointTraits< MakeLeveledLatticePoint<ELP> >
{
static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
public:
// forwarding traits from ELP
  using Trait_ScalarProductStorageType = Get_ScalarProductStorageType<ELP>;
  using Trait_ScalarProductStorageType_Full  = MakeLeveledScalar<Get_ScalarProductStorageType<ELP>>;
  using Trait_CoordinateType          = Get_CoordinateType<ELP>;
  using Trait_AbsoluteCoos            = Get_AbsoluteCooType<ELP>;
  using Trait_RepCooType              = Get_RepCooType<ELP>;
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

  static_assert(Has_Leveled<ELP>::value == false,"");
  static_assert(ApproxLevelOf<ELP>::value == 0,"");
  using Trait_Leveled = std::true_type;
  // using Trait_DelayedScalarProduct = std::true_type;
  using Trait_ApproxLevel          = std::integral_constant<unsigned int, 0>;
};

template<class ELP>
class MakeLeveledLatticePoint
: public GeneralLatticePoint<MakeLeveledLatticePoint<ELP>>
{
  public:
  using LatticePointTag         = std::true_type;
  using Myself = MakeLeveledLatticePoint<ELP>;

  using PlainCooType  = Get_CoordinateType<ELP>; // may be void
  using RepCooType    = Get_RepCooType<ELP>;
  using ScalarProductStorageType = Get_ScalarProductStorageType<ELP>;
  using ScalarProductStorageType_Full = Get_ScalarProductStorageType_Full<Myself>;

  using LeveledComparison = std::false_type; // Unsure about this one...
  static_assert(Has_LeveledObject<ELP>::value == false,"");
  using LeveledObject = std::true_type;
  using LeveledObject_Base = std::true_type;
  static_assert(ApproxLevelOf<ELP>::value==0,"");
  static constexpr unsigned int ApproxLevel = 0;
  template<unsigned int level> using ObjectAtLevel = mystd::conditional_t<level==0, ELP, void>;
  private:
  ELP elp;
  public:
  constexpr       MakeLeveledLatticePoint (ELP const &v) = delete;
  CPP14CONSTEXPR  MakeLeveledLatticePoint (ELP      &&v) : elp(std::move(v)) {}
  constexpr       operator ELP() const & {return elp;}
  CPP14CONSTEXPR  operator ELP()      && {return std::move(elp);}

  template<unsigned int level>
  constexpr ELP const & access() const { static_assert(level==0,""); return elp;}
  template<unsigned int level>
            ELP       & access()       { static_assert(level==0,""); return elp;}
  // TODO: Rvalue variants?

  template<unsigned int level>
  constexpr ELP const & get_value_at_level() const { static_assert(level==0,""); return elp;}
  template<unsigned int level>
            ELP       & get_value_at_level()       { static_assert(level==0,""); return elp;}
// forward functionality of ELP. Note that the this is *not* captured by the conversion operators,
// because MakeLeveledLatticePoint<ELP> is derived from GeneralLatticePoint, and those defaults will have
// precendence over any potential conversions.

  static std::string class_name() { return ELP::class_name() + " L"; }

  template<class T=ELP, class Arg>
  inline PlainCooType &      operator[](Arg &&arg)
  {
    T_IS_ELP; static_assert(Has_ExposesCoos<T>::value && Has_Coos_RW<T>::value,"");
    return elp[std::forward<Arg>(arg)];
  }
  template<class T=ELP, class Arg>
  inline PlainCooType const& operator[](Arg &&arg) const
  {
    T_IS_ELP; static_assert(Has_ExposesCoos<T>::value,"");
    return elp[std::forward<Arg>(arg)];
  }

// operators<,>,<=, >= : No overloads. Defaults is correct.
// forward +=,*=,-=,unary-
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline Myself& operator+=(LatP2 &&x2) { elp+=std::forward<LatP2>(x2); return *this; }
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline Myself& operator-=(LatP2 &&x2) { elp-=std::forward<LatP2>(x2); return *this; }
  template<class Multiplier>
  inline Myself& operator*=(Multiplier &&x2) { elp*=std::forward<Multiplier>(x2); return *this; }

  inline Myself operator-() &&   { return static_cast<Myself>(-std::move(elp)); }
  inline bool operator==(Myself const &x2) const { return elp == x2.elp; }
  inline bool operator==(ELP const &x2) const { return elp=x2; }

  // forward get_internal_rep_size, get_internal_rep
  template<class T=ELP>
  CPP14CONSTEXPR inline auto get_internal_rep_size() const -> decltype( std::declval<T>().get_internal_rep_size() )
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return elp.get_internal_rep_size();
  }
  template<class T=ELP, class Arg>
  inline RepCooType const & get_internal_rep(Arg &&arg) const
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return elp.get_internal_rep(std::forward<Arg>(arg));
  }
  template<class T=ELP, class Arg>
  inline RepCooType & get_internal_rep(Arg &&arg)
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value && Has_InternalRep_RW<T>::value,"");
    return elp.get_internal_rep(std::forward<Arg>(arg));
  }

  // forward get_absolute_coo
  template<class Arg> auto inline get_absolute_coo(Arg &&arg) const
    -> decltype( std::declval<ELP>().get_absolute_coo( std::declval<Arg &&>() ))
  { return elp.get_absolute_coo(std::forward<Arg>(arg));  }

  // forward get_dim
  CPP14CONSTEXPR auto inline get_dim() const -> decltype( std::declval<ELP>().get_dim() ) { return elp.get_dim(); }

  // forward write_lp_to_stream
  inline std::ostream& write_lp_to_stream(std::ostream &os, bool const include_norm2=true, bool const include_approx =true) const
  { return elp.write_lp_to_stream(os, include_norm2,include_approx); }

  // forward write_lp_rep_to_stream
  template<class T=ELP>
  inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return elp.write_lp_rep_to_stream(os);
  }

//
//  //TODO: read_from_stream
//

  void fill_with_zero() { elp.fill_with_zero(); }
  void make_negative()  { elp.make_negative(); }
  bool is_zero() const { return elp.is_zero(); }

  MakeLeveledLatticePoint make_copy() const & { return static_cast<MakeLeveledLatticePoint>(elp.make_copy()); }

// TODO: More overloads?
  void sanitize() { elp.sanitize(); }
  void sanitize( ScalarProductStorageType const &norm2) { elp.sanitize(norm2); }
  void sanitize( ScalarProductStorageType_Full const &norm2) { elp.sanitize(static_cast<ScalarProductStorageType>(norm2)); }

// forward get_norm2()
  [[deprecated]] // Note: This function should never be called anyway. Once we add an approximation, we return a delayed object anyway.
                 // TODO: Maybe return a delayed object already?
  inline ScalarProductStorageType_Full get_norm2() const { return static_cast<ScalarProductStorageType_Full>(elp.get_norm2());}

  template<unsigned int level>
  inline auto get_norm2_at_level() const -> decltype (std::declval<ELP>().get_norm2() )
  {
    static_assert(level==0,"");
    return elp.get_norm2();
  }

  inline ScalarProductStorageType_Full get_norm2_full() const { return static_cast<ScalarProductStorageType_Full>(elp.get_norm2()); }

  // forward scalar products:
  template<unsigned int level>
  inline auto do_compute_sc_product_at_level(Myself const &x2) const
    -> decltype( std::declval<ELP>().do_compute_sc_product(std::declval<ELP>() ) )
  {
    static_assert(level==0,"");
    return elp.do_compute_sc_product(x2.elp);
  }

  // Lazy???
  [[deprecated]] inline ScalarProductStorageType_Full do_compute_sc_product(Myself const &x2) const
  { return static_cast<ScalarProductStorageType_Full>(elp.do_compute_sc_product(x2.elp)); }

  inline ScalarProductStorageType_Full do_compute_sc_product_full(Myself const &x2) const
  { return static_cast<ScalarProductStorageType_Full>(elp.do_compute_sc_product(x2.elp)); }

/*
TODO: Sanitize interface
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
*/
};

// static initializer
template<class ELP>
class StaticInitializer<MakeLeveledLatticePoint<ELP>>
final : public DefaultStaticInitializer<MakeLeveledLatticePoint<ELP>>
{
  StaticInitializer<ELP>  const init_elp;
  public:
  template<class X>
  explicit StaticInitializer(X &&init_arg) : init_elp(std::forward<X>(init_arg))
  {
    static_assert(IsArgForStaticInitializer<mystd::decay_t<X>>::value,"");
  }
};

/*********************************************************************************************
This adds new levels to a lattice point, retaining the functionality of being a lattice point.
**********************************************************************************************/

// Forward traits. We add the ScalarProduct of the Approximation to the Full Scalar product type.
template<class ELP, class ApproxLP>
class LatticePointTraits< AddApproximationToLatP<ELP,ApproxLP> >
{
  static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
  static_assert(Has_Leveled<ELP>::value == true,"");
private:
  using ApproxScalar = typename ApproxLP::ScalarProductType;
public:
// forwarding traits from ELP
  using Trait_ScalarProductStorageType = Get_ScalarProductStorageType<ELP>;
  using Trait_ScalarProductStorageType_Full  = AddApproximation<Get_ScalarProductStorageType_Full<ELP>,ApproxScalar>;
  using Trait_CoordinateType          = Get_CoordinateType<ELP>;
  using Trait_AbsoluteCoos            = Get_AbsoluteCooType<ELP>;
  using Trait_RepCooType              = Get_RepCooType<ELP>;
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
  using Trait_Leveled                 = std::true_type; // we static_assert above it's true for ELP.
  using Trait_ApproxLevel          = std::integral_constant<unsigned int, ApproxLevelOf<ELP>::value +1>;
};

template<class ELP, class ApproxLP> class AddApproximationToLatP
: public GeneralLatticePoint<AddApproximationToLatP<ELP,ApproxLP>>
{
  public:
  using LatticePointTag         = std::true_type;
  using Myself                  = AddApproximationToLatP<ELP,ApproxLP>;
  static_assert(Has_LeveledObject<ELP>::value == true,"");
  using LeveledComparison = std::false_type;
  using LeveledObject     = std::true_type;
  using LeveledObject_Base= std::true_type;
  static constexpr unsigned int ApproxLevel = ApproxLevelOf<ELP>::value + 1;
  static_assert(std::is_same<typename ELP::template ObjectAtLevel<ApproxLevel-1>,ApproxLP>::value == false,"");
  template<unsigned int level> using ObjectAtLevel = mystd::conditional_t
    <level == ApproxLevel,ApproxLP, typename ELP::template ObjectAtLevel<level> >;
  template<unsigned int level> using Helper = Helpers::AddApproximation_Helper<level==ApproxLevel>;
  template<bool> friend class Helpers::AddApproximation_Helper;
  using BaseVector = ObjectAtLevel<0>;

  using PlainCooType  = Get_CoordinateType<ELP>; // may be void
  using RepCooType    = Get_RepCooType<ELP>;
  using ScalarProductStorageType = Get_ScalarProductStorageType<ELP>;
  using ScalarProductStorageType_Full = Get_ScalarProductStorageType_Full<Myself>;


  private:
  ELP old_levels;
  ApproxLP next_level;
  public:

  // leveled object functionality, with deleted copy constructor at level 0 (which implies others
  // need to be deleted as well).
                 explicit AddApproximationToLatP(ELP const &olds, ApproxLP const &news) = delete;
  CPP14CONSTEXPR explicit AddApproximationToLatP(ELP      &&olds, ApproxLP const &news)
    :old_levels(std::move(olds)),next_level(news) {}
  CPP14CONSTEXPR explicit AddApproximationToLatP(ELP const &olds, ApproxLP      &&news) = delete;
  CPP14CONSTEXPR explicit AddApproximationToLatP(ELP      &&olds, ApproxLP      &&news)
    :old_levels(std::move(olds)),next_level(std::move(news)){}
                 explicit AddApproximationToLatP(ELP const &olds) = delete;
  CPP14CONSTEXPR explicit AddApproximationToLatP(ELP      &&olds)
    :old_levels(std::move(olds)), next_level(static_cast<BaseVector>(old_levels)) {}
                          AddApproximationToLatP(BaseVector const &old_base) = delete;
  CPP14CONSTEXPR          AddApproximationToLatP(BaseVector      &&old_base)
    :old_levels(std::move(old_base)),next_level(static_cast<BaseVector>(old_levels)) {}

  constexpr      operator BaseVector() const & { return static_cast<BaseVector>(old_levels); }
  CPP14CONSTEXPR operator BaseVector()      && { return static_cast<BaseVector>(std::move(old_levels));}

  template<unsigned int level>
  constexpr ObjectAtLevel<level> const & access() const
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }
  template<unsigned int level>
  CPP14CONSTEXPR ObjectAtLevel<level> &  access()
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }
  template<unsigned int level>
  constexpr ObjectAtLevel<level> const & get_value_at_level() const
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }
  template<unsigned int level>
  CPP14CONSTEXPR ObjectAtLevel<level> &  get_value_at_level()
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }
// TODO: RValue variants?
  private:
  void recompute_approx() { next_level = static_cast<ApproxLP>(old_levels.access<0>() ); }
  public:

  static std::string class_name() { return ELP::class_name() + " with " << ApproxLevel << " approximations"; }

  template<class T=ELP, class Arg>
  inline PlainCooType &      operator[](Arg &&arg)
  {
    static_assert(std::is_same<T,ELP>::value,"");
    static_assert(Has_ExposesCoos<T>::value && Has_Coos_RW<T>::value,"");
    return old_levels[std::forward<Arg>(arg)];
  }
  template<class T=ELP, class Arg>
  inline PlainCooType const& operator[](Arg &&arg) const
  {
    static_assert(std::is_same<T,ELP>::value,"");
    static_assert(Has_ExposesCoos<T>::value,"");
    return old_levels[std::forward<Arg>(arg)];
  }

// operators<,>,<=, >= : No overloads. Defaults is correct.
// forward +=,*=,-=,unary-

// TODO: Consider making those Lazy.
// Note that these are not level-aware wrt x2 and thus might not work at all.

  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  [[deprecated]] inline Myself& operator+=(LatP2 &&x2)
  {
    old_levels+=std::forward<LatP2>(x2); recompute_approx();
    return *this;
  }
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  [[deprecated]]inline Myself& operator-=(LatP2 &&x2)
  {
    old_levels-=std::forward<LatP2>(x2); recompute_approx();
    return *this;
  }
  template<class Multiplier>
  [[deprecated]] inline Myself& operator*=(Multiplier &&x2)
  {
    old_levels*=std::forward<Multiplier>(x2); recompute_approx();
    return *this;
  }

  [[deprecated]] inline Myself operator-() &&   { return static_cast<Myself>(-std::move(old_levels)); }
  inline bool operator==(Myself const &x2) const
  {
    return (next_level == x2.next_level) && (old_levels == x2.old_levels) ;
  }
  inline bool operator==(ELP const &x2) const { return old_levels==x2; }

  // forward get_internal_rep_size, get_internal_rep
  template<class T=ELP>
  CPP14CONSTEXPR inline auto get_internal_rep_size() const -> decltype( std::declval<T>().get_internal_rep_size() )
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return old_levels.get_internal_rep_size();
  }
  template<class T=ELP, class Arg>
  inline RepCooType const & get_internal_rep(Arg &&arg) const
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return old_levels.get_internal_rep(std::forward<Arg>(arg));
  }
  template<class T=ELP, class Arg>
  inline RepCooType & get_internal_rep(Arg &&arg)
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value && Has_InternalRep_RW<T>::value,"");
    return old_levels.get_internal_rep(std::forward<Arg>(arg));
  }

  // forward get_absolute_coo
  template<class Arg> auto inline get_absolute_coo(Arg &&arg) const
    -> decltype( std::declval<ELP>().get_absolute_coo( std::declval<Arg &&>() ))
  { return old_levels.get_absolute_coo(std::forward<Arg>(arg));  }

  // forward get_dim
  CPP14CONSTEXPR auto inline get_dim() const -> decltype( std::declval<ELP>().get_dim() ) { return old_levels.get_dim(); }

  // forward write_lp_to_stream
  inline std::ostream& write_lp_to_stream(std::ostream &os, bool const include_norm2=true, bool const include_approx =true) const
  { return old_levels.write_lp_to_stream(os, include_norm2,include_approx); }

  // forward write_lp_rep_to_stream
  template<class T=ELP>
  inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return old_levels.write_lp_rep_to_stream(os);
  }

//
//  //TODO: read_from_stream
//

  [[deprecated]] void fill_with_zero() { old_levels.fill_with_zero();recompute_approx(); }
  [[deprecated]] void make_negative()  { old_levels.make_negative();recompute_approx(); }
  bool is_zero() const { return old_levels.is_zero(); }

  Myself make_copy() const & { return static_cast<Myself>(old_levels.access<0>().make_copy() ); }
  Myself make_copy() && =delete;

// TODO: More overloads?
  void sanitize() { old_levels.sanitize(); recompute_approx(); }
  void sanitize( ScalarProductStorageType const &norm2) { old_levels.sanitize(norm2); recompute_approx(); }
  void sanitize( ScalarProductStorageType_Full const &norm2) { old_levels.sanitize(static_cast<ScalarProductStorageType>(norm2)); recompute_approx(); }

/*
// forward get_norm2()
  [[deprecated]] // Note: This function should never be called anyway. Once we add an approximation, we return a delayed object anyway.
                 // TODO: Maybe return a delayed object already?
  inline ScalarProductStorageType_Full get_norm2() const { return static_cast<ScalarProductStorageType_Full>(elp.get_norm2());}

  template<unsigned int level>
  inline auto get_norm2_at_level() const -> decltype (std::declval<ELP>().get_norm2() )
  {
    static_assert(level==0,"");
    return elp.get_norm2();
  }

  inline ScalarProductStorageType_Full get_norm2_full() const { return static_cast<ScalarProductStorageType_Full>(elp.get_norm2()); }

  // forward scalar products:
  template<unsigned int level>
  inline auto do_compute_sc_product_at_level(Myself const &x2) const
    -> decltype( std::declval<ELP>().do_compute_sc_product(std::declval<ELP>() ) )
  {
    static_assert(level==0,"");
    return elp.do_compute_sc_product(x2.elp);
  }

  // Lazy???
  [[deprecated]] inline ScalarProductStorageType_Full do_compute_sc_product(Myself const &x2) const
  { return static_cast<ScalarProductStorageType_Full>(elp.do_compute_sc_product(x2.elp)); }

  inline ScalarProductStorageType_Full do_compute_sc_product_full(Myself const &x2) const
  { return static_cast<ScalarProductStorageType_Full>(elp.do_compute_sc_product(x2.elp)); }

*/

};

} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
