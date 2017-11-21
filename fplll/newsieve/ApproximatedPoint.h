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
template<class OldScalar, class ApproxScalar> struct AddApproximation;
namespace Helpers
{
template<class,class,bool> struct AddApproximation_Helper;
}
}
/**
  This class stores a scalar with an approximation to that scalar.
  Currently, this class only allows constant objects. We have no arithmetic.
*/
namespace GaussSieve{

template<class Scalar>
struct MakeLeveledScalar
{
  static_assert(Has_LeveledObject<Scalar>::value==false,"");
  static_assert(ApproxLevelOf<Scalar>::value==0,"");
  using LeveledComparison = std::true_type;
  using LeveledObject = std::true_type;
  using LeveledObject_Base = std::true_type;
  static constexpr unsigned int ApproxLevel = 0;
  template<unsigned int level> using ObjectAtLevel = mystd::conditional_t<level==0, Scalar, void>;
  Scalar scalar;
  constexpr      MakeLeveledScalar(Scalar const &sc) : scalar(sc) {}
  CPP14CONSTEXPR MakeLeveledScalar(Scalar &&sc) : scalar(std::move(sc)){}
  // default move, copy etc.
  template<unsigned int level>
  constexpr operator Scalar() const & {return scalar;}
  constexpr operator Scalar()      && {return std::move(scalar);}

  template<unsigned int level>
  constexpr Scalar const & access() const { static_assert(level==0,""); return scalar;}
  template<unsigned int level>
            Scalar       & access()       { static_assert(level==0,""); return scalar;}

  template<unsigned int level>
  constexpr Scalar const & get_value_at_level() const { static_assert(level==0,""); return scalar;}
  template<unsigned int level>
            Scalar       & get_value_at_level()       { static_assert(level==0,""); return scalar;}
};

// OldScalar is already a leveled object which may contain approximations
template<class OldScalar, class ApproxScalar>
struct AddApproximation
{
  static_assert(Has_LeveledObject<OldScalar>::value==true,"");
  static_assert(Has_LeveledObject<ApproxScalar>::value==false,"");
  using LeveledComparison = std::true_type;
  using LeveledObject = std::true_type;
  using LeveledObject_Base = std::true_type;
  static constexpr unsigned int ApproxLevel = ApproxLevelOf<OldScalar>::value + 1;
  static_assert(std::is_same<typename OldScalar::template ObjectAtLevel<ApproxLevel-1>,ApproxScalar >::value==false,"");
  template<unsigned int level> using ObjectAtLevel = mystd::conditional_t
    < level==ApproxLevel, ApproxScalar, typename OldScalar::template ObjectAtLevel<level>  >;
  using BaseScalar = ObjectAtLevel<0>;
  template<unsigned int level> using Helper = Helpers::AddApproximation_Helper<OldScalar,ApproxScalar,level==ApproxLevel>;
  template<class, class, bool> friend class Helpers::AddApproximation_Helper;

  OldScalar old_scalar;
  ApproxScalar approx_scalar;

  constexpr      explicit AddApproximation(OldScalar const &exact, ApproxScalar const &approx)
    :old_scalar(exact), approx_scalar(approx) {}
  CPP14CONSTEXPR explicit AddApproximation(OldScalar      &&exact, ApproxScalar const &approx)
    :old_scalar(std::move(exact)),approx_scalar(approx) {}
  CPP14CONSTEXPR explicit AddApproximation(OldScalar const &exact, ApproxScalar      &&approx)
    :old_scalar(exact), approx_scalar(std::move(approx)) {}
  CPP14CONSTEXPR explicit AddApproximation(OldScalar      &&exact, ApproxScalar      &&approx)
    :old_scalar(std::move(exact)),approx_scalar(std::move(approx)){}
  constexpr      explicit AddApproximation(OldScalar const &exact)
    :old_scalar(exact), approx_scalar(static_cast<BaseScalar>(exact)){}
  CPP14CONSTEXPR explicit AddApproximation(OldScalar      &&exact)
    :old_scalar(std::move(exact)), approx_scalar(static_cast<BaseScalar>(old_scalar)) {}
  constexpr               AddApproximation(BaseScalar const &exact_base)
    :old_scalar(exact_base), approx_scalar(exact_base) {}
  CPP14CONSTEXPR          AddApproximation(BaseScalar      &&exact_base)
    :old_scalar(std::move(exact_base)),approx_scalar(static_cast<BaseScalar>(old_scalar)) {}

  constexpr      operator BaseScalar() const & { return static_cast<BaseScalar>(old_scalar); }
  CPP14CONSTEXPR operator BaseScalar()      && { return static_cast<BaseScalar>(std::move(old_scalar));}

  template<unsigned int level>
  constexpr ObjectAtLevel<level> const & access()
  {
    static_assert(level<=ApproxLevel,"");
    return Helper<level>::template get<level>(*this);
  }

};

namespace Helpers{
template<class OldScalar, class ApproxScalar, bool AtMaxLevel> struct AddApproximation_Helper;

template<class OldScalar, class ApproxScalar>
struct AddApproximation_Helper<OldScalar, ApproxScalar, true>
{
  using Argument = AddApproximation<OldScalar,ApproxScalar>;
  template<unsigned int level>
  static inline constexpr ApproxScalar const & get(Argument const &ob)
  {
    static_assert(level==Argument::ApproxLevel,"");
    return ob.approx_scalar;
  }
  template<unsigned int level>
  static inline ApproxScalar       & get(Argument & ob)
  {
    static_assert(level==Argument::ApproxLevel,"");
    return ob.approx_scalar;
  }
};

template<class OldScalar,class ApproxScalar>
struct AddApproximation_Helper<OldScalar,ApproxScalar,false>
{
  using Argument = AddApproximation<OldScalar,ApproxScalar>;
  template<unsigned int level>
  static inline constexpr typename OldScalar::template ObjectAtLevel<level> const & get(Argument const &ob)
  {
    static_assert(level < Argument::ApproxLevel,"");
    return ob.old_scalar.template access<level>();
  }
  template<unsigned int level>
  static inline typename OldScalar::template ObjectAtLevel<level> &       get(Argument       &ob)
  {
    static_assert(level < Argument::ApproxLevel,"");
    return ob.old_scalar.template access<level>();
  }
};
} // end namespace Helpers



} // end namespace GaussSieve

#endif // APPROXIMATED_POINT_H
