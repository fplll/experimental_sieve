#ifndef POINT_WITH_BITAPPROX_H
#define POINT_WITH_BITAPPROX_H

#include "DefaultIncludes.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include "SimHash.h"

/**
  This file declares and defines a class template
  AddBitApproximationToLP<ELP, CooSelection>
  that is used to add / combine a bitapproximation to a given lattice point.

  Notably, if ELP is a Lattice Point class (-> LatticePointConcept.h)
  and CooSelection is a Coordinate Selection ( for SimHashes) (cf. -> SimHash.h)

  then AddBitApproximationToLP<ELP, CooSelection> is another lattice point class
  that includes a sim_hash (whose paramteres are controlled by CooSelection).
  This resulting class is convertible to / from ELP.
*/

// We use class T = ELP dummy parameters for TEMPL_RESTRICT_*
// (see description of TEMPL_RESTRICT_DECL in Compat.h for the reason why)
// This macro is undef'd at the end of the file.
#define T_IS_ELP static_assert(std::is_same<T, ELP>::value, "Wrong template argument")

namespace GaussSieve
{

// clang-format off
class SimHashArgTag {}; // empty class, purely used to disambiguate function signatures.
// clang-format on

template <class ELP, class CooSelection> class AddBitApproximationToLP;

// Forward traits:
template <class ELP, class CooSelection>
struct LatticePointTraits<AddBitApproximationToLP<ELP, CooSelection>>
{
  static_assert(IsALatticePoint<ELP>::value, "ELP is no lattice point");
  static_assert(IsACoordinateSelection<CooSelection>::value, "Wrong 2nd template arg.");

public:
  // forwarding traits from ELP, except for BitApprox
  // Note: NormalizeTraits transforms anything with a bool ::value to a usual
  // std::true_type / std::false_type. Probably unneccessary, but avoids surprises.
  using Trait_ScalarProductStorageType      = Get_ScalarProductStorageType<ELP>;
  using Trait_ScalarProductStorageType_Full = Get_ScalarProductStorageType_Full<ELP>;
  using Trait_CoordinateType                = Get_CoordinateType<ELP>;
  using Trait_AbsoluteCoos                  = Get_AbsoluteCooType<ELP>;
  using Trait_RepCooType                    = Get_RepCooType<ELP>;
  using Trait_ExposesCoos                   = NormalizeTrait<Has_ExposesCoos<ELP>>;
  using Trait_Coos_RW                       = NormalizeTrait<Has_Coos_RW<ELP>>;
  using Trait_ExposesInternalRep            = NormalizeTrait<Has_ExposesInternalRep<ELP>>;
  using Trait_InternalRepLinear             = NormalizeTrait<Has_InternalRepLinear<ELP>>;
  using Trait_InternalRep_RW                = NormalizeTrait<Has_InternalRep_RW<ELP>>;
  using Trait_InternalRepByCoos             = NormalizeTrait<Has_InternalRepByCoos<ELP>>;
  using Trait_InternalRepIsAbsolute         = NormalizeTrait<Has_InternalRepIsAbsolute<ELP>>;
  using Trait_CheapNorm2                    = NormalizeTrait<Has_CheapNorm2<ELP>>;
  using Trait_CheapNegate                   = NormalizeTrait<Has_CheapNegate<ELP>>;
  using Trait_Leveled                       = NormalizeTrait<Has_Leveled<ELP>>;
  using Trait_ApproxLevel                   = Get_ApproxLevel<ELP>;

  // Adding two Bitapproximations would not work as expected:
  static_assert(Has_BitApprox<ELP>::value == false, "Trying to add 2 bitapproximations");
  // clang-format off
  using Trait_BitApprox                     = std::true_type;
  // clang-format on
};

/**
  This creates a new lattice point with a sim_hash from an old one.
  Essentially, it has an ELP member and a sim_hash member and just forwards the whole interface
  (as specified by the lattice point concept) to elp.
  The only difference to elp are the functions relatied to sim_hashes.
*/

template <class ELP, class CooSelection>
class AddBitApproximationToLP final
    : public GeneralLatticePoint<AddBitApproximationToLP<ELP, CooSelection>>
{
  static_assert(IsALatticePoint<ELP>::value, "");
  static_assert(IsACoordinateSelection<CooSelection>::value, "");

public:
  // forward from ELP / CooSelection
  using LatticePointTag = std::true_type;
  using Myself          = AddBitApproximationToLP<ELP, CooSelection>;
  using SimHashBlock    = typename CooSelection::SimHashBlock;
  using SimHashes       = typename CooSelection::SimHashes;

  using PlainCooType                  = Get_CoordinateType<ELP>;  // may be void
  using RepCooType                    = Get_RepCooType<ELP>;
  using ScalarProductStorageType      = Get_ScalarProductStorageType<ELP>;
  using ScalarProductStorageType_Full = Get_ScalarProductStorageType_Full<Myself>;

private:
  // sim_hashes are a sim_hash of the point elp.
  // Note that we do NOT automatically maintain sim_hashes.
  // We rely on the user to call update_bitapprox()
  // This is arguably ugly, but update_bitapprox() matters for the runtime.
  // (The issue is with multiple operations like x+=y; x*=z in a row, where we do not need to
  // update the sim_hash for the intermediates. While we mostly avoid such intermediates in the
  // first place, we do not want to rely on that)
  ELP elp;
  SimHashes sim_hashes;

public:
  // clang-format off

  // move-construct from ELP
  AddBitApproximationToLP(ELP const &v) = delete;
  CPP14CONSTEXPR AddBitApproximationToLP(ELP &&v) noexcept
      : elp(std::move(v)),
        sim_hashes(GlobalBitApproxData<CooSelection>::coo_selection.compute_all_bitapproximations(*this))
  {
  }
  // clang-format on

  // allow conversion back to ELP.
  // clang-format off
  constexpr      operator ELP() const & { return elp; }
  CPP14CONSTEXPR operator ELP()      && { return std::move(elp); }
  // clang-format on

  /*
    Untested, probably does not work anyway.
    template<class Arg1, class... Args>
    CPP14CONSTEXPR AddBitApproximationToLP(Arg1 &&arg1, Args &&... args)
        : elp(std::forward<Arg1>(arg1), std::forward<Args>(args)...),
          sim_hashes(GlobalBitApproxData<CooSelection>::coo_selection.compute_all_bitapproximations(*this))
    {}
  */

  // a constructor that takes an already known sim_hash, to avoid recomputation.
  // Note: There is some serious danger here to accidentally invoke the perferct forward above.
  //       SimHashArgTag's sole purpose is disambiguation and causing easier-to-understand errors.
  CPP14CONSTEXPR AddBitApproximationToLP(SimHashArgTag, ELP &&v, SimHashes const &new_sim_hashes)
      : elp(std::move(v)), sim_hashes(new_sim_hashes)
  {
  }

  CPP14CONSTEXPR AddBitApproximationToLP(SimHashArgTag, ELP &&v, SimHashes &&new_sim_hashes)
      : elp(std::move(v)), sim_hashes(std::move(new_sim_hashes))
  {
  }

  // update the stored sim_hashes. Note that it is the users job to call this function!
  void update_bitapprox()
  {
    sim_hashes =
        GlobalBitApproxData<CooSelection>::coo_selection.compute_all_bitapproximations(*this);
  }

  // gives access to a specific block of the sim_hash
  SimHashBlock const &access_bitapproximation(unsigned int level) const
  {
    return sim_hashes[level];
  }

  // allows to retrieve all sim_hashes. Note that this can only be used on temporaries
  // or with std::move(). The sim_hash stored in the point is invalidated, but elp remains valid.
  SimHashes take_bitapproximations() && { return std::move(sim_hashes); }

  // forward functionality of ELP. Note that the this is *not* captured by the conversion operators,
  // because MakeLeveledLatticePoint<ELP> is derived from GeneralLatticePoint, and those defaults
  // will have precendence over any potential conversions.

  // A LOT OF BORING BOILERPLATE CODE FOLLOWS...

  static std::string class_name() { return ELP::class_name() + " with Bitapproximations"; }

  // clang-format puts the template introduction on a line of its own or not depending
  // on available space... Unfortunately, this makes analogous usages in consecutive lines
  // inconsistent. We keep it consistent within this file.
  // clang-format off
  template<class T = ELP, class Arg>
  inline PlainCooType &operator[](Arg &&arg)
  {
    T_IS_ELP;
    static_assert(Has_ExposesCoos<T>::value && Has_Coos_RW<T>::value, "");
    return elp[std::forward<Arg>(arg)];
  }
  template<class T = ELP, class Arg>
  inline PlainCooType const &operator[](Arg &&arg) const
  {
    T_IS_ELP;
    static_assert(Has_ExposesCoos<T>::value, "");
    return elp[std::forward<Arg>(arg)];
  }
  // clang-format on

  // operators<,>,<=, >= : No overloads. Defaults is correct.
  // forward +=,*=,-=,unary-
  // clang-format off
  template <class LatP2>
  inline Myself &operator+=(LatP2 &&x2) { elp += std::forward<LatP2>(x2); return *this; }
  template <class LatP2>
  inline Myself &operator-=(LatP2 &&x2) { elp -= std::forward<LatP2>(x2); return *this; }
  template <class Multiplier>
  inline Myself &operator*=(Multiplier &&x2) { elp *= std::forward<Multiplier>(x2); return *this; }
  // clang-format on

  // clang-format off
  inline Myself operator-() &&  { elp.make_negative(); return std::move(*this); }
  inline bool operator==(Myself const &x2) const { return elp == x2.elp; }
  inline bool operator==(ELP    const &x2) const { return elp == x2;     }
  // clang-format on

  // forward get_internal_rep_size, get_internal_rep
  template <class T = ELP>
  CPP14CONSTEXPR inline auto get_internal_rep_size() const
      -> decltype(std::declval<T>().get_internal_rep_size())
  {
    T_IS_ELP;
    static_assert(Has_ExposesInternalRep<T>::value, "");
    return elp.get_internal_rep_size();
  }

  // clang-format off
  template <class T = ELP, class Arg>
  inline RepCooType const &get_internal_rep(Arg &&arg) const
  {
    T_IS_ELP;
    static_assert(Has_ExposesInternalRep<T>::value, "");
    return elp.get_internal_rep(std::forward<Arg>(arg));
  }
  // clang-formant on

  // clang-format off
  template <class T = ELP, class Arg>
  inline RepCooType &get_internal_rep(Arg &&arg)
  {
    T_IS_ELP;
    static_assert(Has_ExposesInternalRep<T>::value && Has_InternalRep_RW<T>::value, "");
    return elp.get_internal_rep(std::forward<Arg>(arg));
  }
  // clang-format on

  // forward get_absolute_coo
  template <class Arg>
  auto inline get_absolute_coo(Arg &&arg) const
      -> decltype(std::declval<ELP>().get_absolute_coo(std::declval<Arg &&>()))
  {
    return elp.get_absolute_coo(std::forward<Arg>(arg));
  }

  // forward get_dim
  CPP14CONSTEXPR auto inline get_dim() const -> decltype(std::declval<ELP>().get_dim())
  {
    return elp.get_dim();
  }

  // forward write_lp_to_stream
  inline std::ostream &write_lp_to_stream(std::ostream &os, bool const include_norm2 = true,
                                          bool const include_approx = true) const
  {
    elp.write_lp_to_stream(os, include_norm2, include_approx);
    if (include_approx)
    {
      os << " SimHash: " << sim_hashes;
    }
    return os;
  }

  // clang-format off
  // forward write_lp_rep_to_stream
  template <class T = ELP>
  inline std::ostream &write_lp_rep_to_stream(std::ostream &os) const
  {
    T_IS_ELP;
    static_assert(Has_ExposesInternalRep<T>::value, "");
    return elp.write_lp_rep_to_stream(os);
  }
  // clang-format on

  //
  // TODO: read_from_stream
  //

  void fill_with_zero() { elp.fill_with_zero(); }
  void make_negative() { elp.make_negative(); }
  bool is_zero() const { return elp.is_zero(); }

  AddBitApproximationToLP make_copy() const &
  {
    return static_cast<AddBitApproximationToLP>(elp.make_copy());
  }

  // TODO: More overloads?
  void sanitize() { elp.sanitize(); }

  // clang-format off
  template <class Arg>
  void sanitize(Arg &&arg) { elp.sanitize(std::forward<Arg>(arg)); }
  // clang-format on

  // forward get_norm2()
  inline auto get_norm2() const -> decltype(std::declval<ELP>().get_norm2())
  {
    return elp.get_norm2();
  }

  template <unsigned int level>
  inline auto get_norm2_at_level() const
      -> decltype(std::declval<ELP>().template get_norm2_at_level<level>())
  {
    return elp.template get_norm2_at_level<level>();
  }

  inline auto get_norm2_full() const -> decltype(std::declval<ELP>().get_norm2_full())
  {
    return elp.get_norm2_full();
  }

  // forward scalar products:
  template <class Arg>
  inline auto do_compute_sc_product(Arg &&arg) const
      -> decltype(std::declval<ELP>().do_compute_sc_product(std::declval<Arg>()))
  {
    return elp.do_compute_sc_product(std::forward<Arg>(arg));
  }

  // clang-format off
  template <unsigned int level, class Arg>
  inline auto do_compute_sc_product_at_level(Arg &&arg) const
      -> decltype( std::declval<ELP>().template do_compute_sc_product_at_level<level>(std::declval<Arg>()) )
  {
    return elp.template do_compute_sc_product_at_level<level>(std::forward<Arg>(arg));
  }
  // clang-format on

  template <class Arg>
  inline auto do_compute_sc_product_full(Arg &&arg) const
      -> decltype(std::declval<ELP>().do_compute_sc_product_full(std::declval<Arg>()))
  {
    return elp.do_compute_sc_product_full(std::forward<Arg>(arg));
  }
};

/**
  Static initializer for AddBitApproximation<ELP, CooSelection>:
  We just forward to the static initializer of ELP by wrapping around it.
  Note that static initialization of CooSelection (or rather of GlobalBitApproxData<CooSelection>)
  is not our job: According to SimHash.h, this is done by the main sieve.
*/

// static initializer
template <class ELP, class CooSelection>
class StaticInitializer<AddBitApproximationToLP<ELP, CooSelection>> final
    : public DefaultStaticInitializer<AddBitApproximationToLP<ELP, CooSelection>>
{
  StaticInitializer<ELP> const init_elp;

public:
  template <class InitArg>
  explicit StaticInitializer(InitArg &&init_arg) : init_elp(std::forward<InitArg>(init_arg))
  {
    static_assert(IsArgForStaticInitializer<mystd::decay_t<InitArg>>::value, "");
  }
};

}  // end namespace GaussSieve

#undef T_IS_ELP

#endif  // POINT_WITH_BITAPPROX_H
