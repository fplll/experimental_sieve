#ifndef APPROXIMATED_POINT_H
#define APPROXIMATED_POINT_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "LatticePointConcept.h"
#include <tuple>
#include "Lazy.h"
#include "GlobalStaticData.h"
#include "PlainLatticePoint.h"

// forward declarations:

namespace GaussSieve{

template<class Scalar> struct MakeLeveledScalar;
template<class OldLevels, class NextLevel> struct AddApproximation;
template<class ELP> class MakeLeveledVector;
namespace Helpers
{
template<class,class,bool> struct AddApproximation_Helper;
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
  This adds another level to a leveled object.
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
  template<unsigned int level> using Helper = Helpers::AddApproximation_Helper<OldLevels,NextLevel,level==ApproxLevel>;
  template<class, class, bool> friend class Helpers::AddApproximation_Helper;

  OldLevels old_levels;
  NextLevel next_level;

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
  // TODO: Rvalue variants.
};

namespace Helpers{
template<class OldLevels, class NextLevel, bool AtMaxLevel> struct AddApproximation_Helper;

template<class OldLevels, class NextLevel>
struct AddApproximation_Helper<OldLevels, NextLevel, true>
{
  using Argument = AddApproximation<OldLevels,NextLevel>;
  template<unsigned int level>
  static inline constexpr NextLevel const & get(Argument const &ob)
  {
    static_assert(level==Argument::ApproxLevel,"");
    return ob.approx_scalar;
  }
  template<unsigned int level>
  static inline NextLevel       & get(Argument & ob)
  {
    static_assert(level==Argument::ApproxLevel,"");
    return ob.approx_scalar;
  }
};

template<class OldLevels,class NextLevel>
struct AddApproximation_Helper<OldLevels,NextLevel,false>
{
  using Argument = AddApproximation<OldLevels,NextLevel>;
  template<unsigned int level>
  static inline constexpr typename OldLevels::template ObjectAtLevel<level> const & get(Argument const &ob)
  {
    static_assert(level < Argument::ApproxLevel,"");
    return ob.old_scalar.template access<level>();
  }
  template<unsigned int level>
  static inline typename OldLevels::template ObjectAtLevel<level> &       get(Argument       &ob)
  {
    static_assert(level < Argument::ApproxLevel,"");
    return ob.old_scalar.template access<level>();
  }
};
} // end namespace Helpers

// Forward traits:
template<class ELP>
class LatticePointTraits< MakeLeveledVector<ELP> >
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
class MakeLeveledVector
: public GeneralLatticePoint<MakeLeveledVector<ELP>>
{
  public:
  using LatticePointTag         = std::true_type;
  using Myself = MakeLeveledVector<ELP>;

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
  constexpr       MakeLeveledVector (ELP const &v) : elp(v)            {}
  CPP14CONSTEXPR  MakeLeveledVector (ELP      &&v) : elp(std::move(v)) {}
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
// because MakeLeveledVector<ELP> is derived from GeneralLatticePoint, and those defaults will have
// precendence over any potential conversions.

  static std::string class_name() { return ELP::class_name() + " L"; }

  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesCoos<T>)>
  PlainCooType &      operator[](Arg &&arg)       { return elp[std::forward<Arg>(arg)]; }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesCoos<T>)>
  PlainCooType const& operator[](Arg &&arg) const { return elp[std::forward<Arg>(arg)]; }

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
  template<class T=ELP, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>)>
  CPP14CONSTEXPR inline auto get_internal_rep_size() const -> decltype( std::declval<ELP>().get_internal_rep_size() ) { return elp.get_internal_rep_size(); }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>)>
  inline RepCooType const & get_internal_rep(Arg &&arg) const { return elp.get_internal_rep(std::forward<Arg>(arg)); }
  template<class T=ELP, class Arg, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<T>, Has_InternalRep_RW<T>)>
  inline RepCooType & get_internal_rep(Arg &&arg) {return elp.get_internal_rep(std::forward<Arg>(arg));}

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
  template<class T=ELP, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<ELP>)>
  inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const { return elp.write_lp_rep_to_stream(os); }

//
//  //TODO: read_from_stream
//

  void fill_with_zero() { elp.fill_with_zero(); }
  void make_negative()  { elp.make_negative(); }
  bool is_zero() const { return elp.is_zero(); }

  MakeLeveledVector make_copy() const & { return static_cast<MakeLeveledVector>(elp.make_copy()); }

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
class StaticInitializer<MakeLeveledVector<ELP>>
final : public DefaultStaticInitializer<MakeLeveledVector<ELP>>
{
  StaticInitializer<ELP>  const init_elp;
  public:
  template<class X,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<mystd::decay_t<X>>)>
  explicit StaticInitializer(X &&init_arg) : init_elp(std::forward<X>(init_arg)){}
};




} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
