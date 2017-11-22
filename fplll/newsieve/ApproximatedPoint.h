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

template<class ELP>
class LatticePointTraits< MakeLeveledVector<ELP> >
{
static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
public:
// forwarding traits from ELP
  using Trait_ScalarProductStorageType = Get_ScalarProductStorageType<ELP>;
  using Trait_ScalarProductStorageType_Full  = Get_ScalarProductStorageType_Full<ELP>;
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



} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
