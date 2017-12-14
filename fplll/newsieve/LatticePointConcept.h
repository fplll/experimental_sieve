#ifndef LATTICE_POINT_CONCEPT_H
#define LATTICE_POINT_CONCEPT_H

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"
#include "SieveUtility.h"
#include <gmpxx.h>

// clang-format status: done

namespace GaussSieve
{

/**
  GeneralLatticePoint defines an interface for the various classes of lattice points we use.
  Notably, we define a class template GeneralLatticePoint<Actual class> from which all
  lattice point classes inherit.
  To avoid overhead, we use CRTP rather than dynamic polymorphism, i.e.
  GeneralLatticePoint is templated by its child classes.

  Note that we have a lot of "default" implementations, which are defined depending on various
  lattice traits. Most of the actual code of these implementations is found in
  LatticePointGeneric.h

  GeneralLatticePoint also specifies which methods an implementation is allowed to override.
*/

// This class template stores the trait typedefs that the individual lattice point classes have.
// There MUST be a specialization for each lattice point class.
// The general LatticePointsTraits template must never be instantiated.
// Note that we can not put these traits into the lattice points classes directly, because that
// would cause circular references due to CRTP.

// Cf. ExacLatticePoint.h for an example of how a concrete specialization is supposed to look like.

template <class LatticePoint> struct LatticePointTraits
{
public:
  // clang-format off
  using Invalid = std::true_type;  // do not set this at all in an specialization.
  LatticePointTraits(...) = delete;
  // static_assert(false) is invalid due to subtleties of C++, even if it may work on some compilers

  using Trait_ScalarProductStorageType      = void;
  using Trait_ScalarProductStorageType_Full = void;
  using Trait_ExposesCoos                   = std::false_type;
  using Trait_CoordinateType                = void;
  using Trait_Coos_RW                       = std::false_type;
  using Trait_AbsoluteCooType               = void;
  using Trait_RepCooType                    = void;
  using Trait_ExposesInternalRep            = std::false_type;
  using Trait_InternalRepLinear             = std::false_type;
  using Trait_InternalRep_RW                = std::false_type;
  using Trait_InternalRepByCoos             = std::false_type;
  using Trait_InternalRepIsAbsolute         = std::false_type;
  using Trait_AbsoluteCoos                  = std::false_type;
  using Trait_CheapNorm2                    = std::false_type;
  using Trait_CheapNegate                   = std::false_type;
  using Trait_BitApprox                     = std::false_type;

  using Trait_ApproxLevel = std::integral_constant<unsigned int, 0>;
  using Trait_Leveled                       = std::false_type;
  // clang-format on
};

/**
  Lattice Point traits:

  Note: In the actual LatticePointTraits class, the traits are actually prefixed with Trait_
  (i.e. LatticePointTraits<LP> contains "using Trait_ExposesCoos = std::true_type/false_type" etc.)

  To retrieve a trait, use Has_TraitName for binary traits or Get_TraitName for non-binary traits.
  Do NOT access the Trait_* traits directly.

  TODO: Implications between traits may not be up-to-date (both in documentation and in code)

  ScalarProductStorageType: A type that can hold the result of a scalar product computation.
                            THIS TRAIT IS MANDATORY.
                            Note that the result from a scalar product computation might actually
                            differ. (due to delayed evaluation)

  ScalarProductStorageType_Full:  A type that can hold the result of a scalar product computation.
                                  In case the Point has approximations, will contain an approximate
                                  scalar product as well.
                                  Default: Same as ScalarProductStorageType

  ExposesCoos : Has an operator[].
                implied by existence of CoordinateType, InternalRepByCoos. Mandates Coordinate Type
  CoordinateType :  return type of operator[] is a (possibly const) reference to this, if available.
                    implies ExposesCoos, required if ExposesCoos is defined
  Coos_RW      : operator[] can return non-const reference.

  NOTE: ExposesCoos / CoordinateType is a purely syntactic trait.
        It does not tell anything about the semantics of [].

  AbsoluteCooType : return type of get_absolute_coo
                    defaults to CoordinateType

  ExposesInternalRep :  Indicates that get_internal_rep(i), get_internal_rep_size() exist.
                        get_internal_rep(i) may be read for 0<= i < get_internal_rep_size().
                        We are guaranteed that these entries determine the point.
                        implied by InternalRepVector_R, InternalRepVector_RW, InternalRepByCoos,
                        InternalRepIsAbsolute

  RepCooType      : return type of get_internal_rep (if available)
                    defaults to CoordinateType

  InternalRepLinear:    Indicates that the internal representation is linear.
                        implies Has_ExposesInternalRep
  InternalRep_RW:       Indicates that the internal representation may be written to.
                        After calling sanitize(), the class will be in a valid state.
                        implies Has_ExposesInternalRep

  InternalRepByCoos:  Indicates that the internal representation is accessible by operator[]
                      implies ExposesInternalRep, ExposesCoos

  InternalRepIsAbsolute:  Indicates that the internal representation gives absolute coordinates.
                          implies ExposesInternalRep

  AbsoluteCoos    : Indicates that [] gives absolute coordinates.

  CheapNorm2 : Set to true_type to indicate that get_norm2() is cheap.
               (typically, it's precomputed and stored with the point)

  BitApprox : Set to indicate that we have a bit-approximation.
              This gives the following promises to the user:
              There are public typedefs:
                  using SimHashBlock (equal to std::bitset<lenght> or compatible)
                  using SimHashes    (a container of SimHashBlocks)
              There are public member functions:

              SimHashBlock [const &] access_bitapproximation(unsigned int level);
              SimHashes take_bitapproximations() &&; (or possibly a lvalue-version)
              void update_bitapprox();

              access_bitapproximation(i) gives const-access to the i'th sim_hash.
              take_bitapproximations() moves the bitapproximations out of the point
              update_bitapprox() recomputes the bitapproximations.
              NOTE: For efficiency reasons, bitapproximations are NOT recomputed when modifying the
              point by += etc. It is currently the caller's job to trigger recomputation.
              (The reason is that recomputation is too slow)
              NOTE: This may be subject to change.

  CheapNegate: Set to true_type to indicate that negation needs no sanitize().

  // Currently unused:

  Leveled : Set to true_type to indicate that the point is a leveled object.
            Cf. Lazy.h for a more details on that.
            For lattice points, leveled objects correspond to lattice points with several layers
            of approximations.
            Note that, to simplify things, we always have leveled functions available, but by
            default, they only work for level 0.

  ApproxLevel: Set to std::integral_constant<unsigned int,...> to indicate the approximation level.
               Defaults to std::integral_constant<unsigned int,0>

               NOTE:  Several classes have a public static constexpr unsigned int member ApproxLevel
                      (without the Trait_ prefix) to indicate the corresponding approximation level.
                      Such classes include the lattice point class itself, but also classes storing
                      exact & approximated scalar product, etc.
                      Use ApproxLevelOf<Some_Class>::value to obtain Some_Class::ApproxLevel
                      (with a default of 0 if Some_Class::ApproxLevel does not exist)

  NOTE: Leveled and ApproxLevel do NOT relate to Bitapproximation.
*/

// forward declaration
template <class Implementation> class GeneralLatticePoint;

/**
  Trait getters

  Usage Has_CheapNegate<PlainLatticePoint>::value
  (as opposed to LatticePointTraits<PlainLatticePoint>::CheapNegate::value )

  The reason is that it is easier to read and that
  Has_CheapNegate<...> defaults to false, whereas
  LatticePointTraits<...>::... defaults to a compile-time error.
  Furthermore, Has_ExposesCoos<...> is aware of neccessary implications between traits.
  I.e. if a trait is implied by another trait, it is automagically set.
*/

// In order to define e.g. Has_CoosRW, we use 3 steps (some traits only require 2):
// First, we need to give the "map" LatP -> LatticePointTraits<LatP>::Trait_CoosRW a name.
// For this define template<class T> using Predicate_CoosRW = LatticePointTraits<T>::Trait_CoosRW;
// Second we use (my)std::is_detected to turn this into a "map" T -> true/false.
// More precisely, we define
// template<class T> using T_CoosRW = mystd::is_detected<Predicate_CoosRW,T>;
// (is_detected<Op,Args,...> equals true_type iff Op<Args,...> is valid, false_type otherwise)
// Third, we post-process T_CoosRW into the actual Get_CoosRW by taking dependencies into account
// Notably, if ExposesCoos && InternalRepByCoos && InternalRep_RW are all true, CoosRW is true as
// well.
// The intermediate steps 1 and 2 are contained in a helper namespace.

namespace TraitHelpers
{
// These are "predicates" expressed as templates Predicate<T> depending on a type T,
// where Predicate<T> is a valid expression iff the predicate holds true. This is used in SFINAE.
// With is_detected<Predicate,T>, such a predicate can be turned into std::true_type/std::false_type
template <class T> using IsTrueType = mystd::enable_if_t<std::is_same<std::true_type, T>::value>;
template <class T> using LatticePointPredicate = IsTrueType<typename T::LatticePointTag>;

// Declares alias for LatticePointTraits<T>::Trait_TraitName. Required to use is_detected.
#define GAUSS_SIEVE_OBTAIN_TRAIT_EXPRESSION(TraitName)                                             \
  template <class T> using Obtain_##TraitName = typename LatticePointTraits<T>::Trait_##TraitName

// Defines Predicate_TraitName as a check whether Trait_TraitName is set to true_type
#define GAUSS_SIEVE_BINARY_TRAIT_PREDICATE(TraitName)                                              \
  template <class T>                                                                               \
  using Predicate_##TraitName = IsTrueType<typename LatticePointTraits<T>::Trait_##TraitName>

// uses (my)std::is_detected to turn the latter into a true/false T_TraitName
#define GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(TraitName)                                          \
  GAUSS_SIEVE_BINARY_TRAIT_PREDICATE(TraitName);                                                   \
  template <class T> using T_##TraitName = mystd::is_detected<Predicate_##TraitName, T>

GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(ExposesCoos);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(Coos_RW);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(ExposesInternalRep);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(InternalRepLinear);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(InternalRep_RW);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(InternalRepByCoos);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(InternalRepIsAbsolute);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(AbsoluteCoos);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(CheapNorm2);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(CheapNegate);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(Leveled);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(AccessNorm2);
GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET(BitApprox);

// "Invalid" is set to true_type in the default instantiation of LatticePointTraits. This is used to
// detect whether we use a specialization.
template <class T> using Invalid_SieveTrait = IsTrueType<typename LatticePointTraits<T>::Invalid>;

// For Get_Trait<LatticePoint>
GAUSS_SIEVE_OBTAIN_TRAIT_EXPRESSION(ScalarProductStorageType);
GAUSS_SIEVE_OBTAIN_TRAIT_EXPRESSION(ScalarProductStorageType_Full);
GAUSS_SIEVE_OBTAIN_TRAIT_EXPRESSION(ApproxLevel);
GAUSS_SIEVE_OBTAIN_TRAIT_EXPRESSION(CoordinateType);
GAUSS_SIEVE_OBTAIN_TRAIT_EXPRESSION(AbsoluteCooType);
GAUSS_SIEVE_OBTAIN_TRAIT_EXPRESSION(RepCooType);

#undef GAUSS_SIEVE_BINARY_TRAIT_PREDICATE_GET
#undef GAUSS_SIEVE_BINARY_TRAIT_PREDICATE
#undef GAUSS_SIEVE_OBTAIN_TRAIT_PREDICATE
}
// clang-format off
template<class T> using IsALatticePoint                  = mystd::is_detected<TraitHelpers::LatticePointPredicate, T>;
template<class T> using DeclaresScalarProductStorageType = mystd::is_detected<TraitHelpers::Obtain_ScalarProductStorageType, T>;
template<class T> using DoesDeclareCoordinateType        = mystd::is_detected<TraitHelpers::Obtain_CoordinateType, T>;
template<class T> using HasNoLPTraits                    = mystd::is_detected<TraitHelpers::Invalid_SieveTrait, T>;
// clang-format on

// Make actual getter for TraitName with a Default:
#define GAUSS_SIEVE_MAKE_TRAIT_GETTER(TraitName, Default)                                          \
  template <class T>                                                                               \
  using Get_##TraitName = mystd::detected_or_t<Default, TraitHelpers::Obtain_##TraitName, T>
GAUSS_SIEVE_MAKE_TRAIT_GETTER(CoordinateType, void);
GAUSS_SIEVE_MAKE_TRAIT_GETTER(ScalarProductStorageType, void);
GAUSS_SIEVE_MAKE_TRAIT_GETTER(ScalarProductStorageType_Full, Get_ScalarProductStorageType<T>);
GAUSS_SIEVE_MAKE_TRAIT_GETTER(AbsoluteCooType, Get_CoordinateType<T>);
GAUSS_SIEVE_MAKE_TRAIT_GETTER(RepCooType, Get_CoordinateType<T>);

// "," in macro argument would cause trouble, so we just write it out:
template <class T>
using Get_ApproxLevel = mystd::detected_or_t<std::integral_constant<unsigned int, 0>,
                                             TraitHelpers::Obtain_ApproxLevel, T>;

#undef GAUSS_SIEVE_MAKE_TRAIT_GETTER

// finally, we define the actual Has_TraitName's. These just take T_TraitNames and
// apply some boolean formulae to enforce implied traits.

// TODO: Make sure these implications are correct and document them.
// TODO: Ensure that if we use implications to set Has_TraitName to true, then Trait_TraitName
//       is not manually set to false_type

// Alas, there is no way to format this nicely.
// We manually make the linebreaks at least internally consistent.
// clang-format off
template<class T> using Has_ExposesCoos =
    mystd::disjunction< TraitHelpers::T_ExposesCoos<T>,
                        TraitHelpers::T_InternalRepByCoos<T>,
                        mystd::negation< std::is_void<Get_CoordinateType<T>> >,
                        TraitHelpers::T_Coos_RW<T>,
                        TraitHelpers::T_AbsoluteCoos<T> >;

template<class T> using Has_Coos_RW =
    mystd::disjunction< TraitHelpers::T_Coos_RW<T>,
                        mystd::conjunction< TraitHelpers::T_InternalRepByCoos<T>,
                                            TraitHelpers::T_InternalRep_RW<T> > >;

template<class T> using Has_ExposesInternalRep =
    mystd::disjunction< TraitHelpers::T_ExposesInternalRep<T>,
                        TraitHelpers::T_InternalRepLinear<T>,
                        TraitHelpers::T_InternalRep_RW<T>,
                        TraitHelpers::T_InternalRepByCoos<T>,
                        TraitHelpers::T_InternalRepIsAbsolute<T> >;

template<class T> using Has_InternalRepLinear = TraitHelpers::T_InternalRepLinear<T>;

template<class T> using Has_InternalRep_RW =
    mystd::disjunction< TraitHelpers::T_InternalRep_RW<T>,
                        mystd::conjunction< TraitHelpers::T_InternalRepByCoos<T>,
                                            TraitHelpers::T_Coos_RW<T> > >;

template<class T> using Has_InternalRepByCoos = TraitHelpers::T_InternalRepByCoos<T>;

template<class T> using Has_InternalRepIsAbsolute =
    mystd::disjunction< TraitHelpers::T_InternalRepIsAbsolute<T>,
                        mystd::conjunction< TraitHelpers::T_InternalRepByCoos<T>,
                                            TraitHelpers::T_AbsoluteCoos<T> > >;

template<class T> using Has_AbsoluteCoos =
    mystd::disjunction< TraitHelpers::T_AbsoluteCoos<T>,
                        mystd::conjunction< TraitHelpers::T_InternalRepByCoos<T>,
                                            TraitHelpers::T_InternalRepIsAbsolute<T> > >;

template<class T> using Has_CheapNorm2  = TraitHelpers::T_CheapNorm2<T>;
template<class T> using Has_CheapNegate = TraitHelpers::T_CheapNegate<T>;
template<class T> using Has_Leveled     = TraitHelpers::T_Leveled<T>;
template<class T> using Has_BitApprox   = TraitHelpers::T_BitApprox<T>;

// Deprecated. Note is is and AND of things above, not an OR.
template<class LatP> using IsRepLinear_RW =
    mystd::conjunction< Has_InternalRepLinear<LatP>, Has_InternalRep_RW<LatP> >;
// clang-format on

// For approximations.
// TODO: Move this somewhere else. It's not really a lattice point trait.

// clang-format off
template<class T> using Obtain_ApproxLevel = std::integral_constant<unsigned int, T::ApproxLevel>;
template<class T> using ApproxLevelOf      = mystd::detected_or_t<std::integral_constant<unsigned int,0>,
								  Obtain_ApproxLevel, T>;
// clang-format on

// have template<class Impl = Latp> in the declarations below. This is to ensure
// only the default is ever used. See below for an explanation.
// clang-format off
#define IMPL_IS_LATP                                                                               \
  static_assert(std::is_same<Impl, LatP>::value, 						   \
	       "Using template member function with wrong type.")
// clang-format on

// unsure whether to use reinterpret_casts here...

// clang-format off
#define CREALTHIS static_cast<LatP const*>(this)
#define REALTHIS  static_cast<LatP*>(this)
// clang-format on

// clang-format off
template<class LatP>
class GeneralLatticePoint
// clang-format on
{
  static_assert(!HasNoLPTraits<LatP>::value, "Trait class not specialized.");
  static_assert(DeclaresScalarProductStorageType<LatP>::value,
                "Lattice Point class does not typedef its scalar product type");
  static_assert(!Has_ExposesCoos<LatP>::value || DoesDeclareCoordinateType<LatP>::value,
                "Lattice Point exposes coos, but does not tell its type.");

  friend LatP;  // makes children able to access private (in addition to protected) members.
                // since the constructor is private, this enforces correct usage.
                // (Note that it may prevent multi-level inheritance)
public:
  using ScalarProductStorageType            = Get_ScalarProductStorageType<LatP>;
  using ScalarProductStorageType_Full       = Get_ScalarProductStorageType_Full<LatP>;
  using CooType                             = Get_CoordinateType<LatP>;  // may be void
  using CoordinateType                      = Get_CoordinateType<LatP>;  // may be void
  using AbsCooType                          = Get_AbsoluteCooType<LatP>;
  using AbsoluteCooType                     = Get_AbsoluteCooType<LatP>;
  using RepCooType                          = Get_RepCooType<LatP>;
  static constexpr unsigned int ApproxLevel = Get_ApproxLevel<LatP>::value;
  static_assert((Has_Leveled<LatP>::value == true) || (ApproxLevel == 0),
                "Has ApproxLevel>0, but does not declare Trait_Leveled.");

private:
  // Empty base class, only callable from its friends (i.e. from LatP)
  explicit constexpr GeneralLatticePoint() = default;

public:
  // This is just to match the implementation of a typical instantiation.
  // Note the the deleted copy constructors and copy assignments prevents default copying
  // derived classes (since the (empty) base class must be copied as well).
  // clang-format off
  GeneralLatticePoint(GeneralLatticePoint const &)            = delete;
  GeneralLatticePoint(GeneralLatticePoint &&)                 = default;
  GeneralLatticePoint &operator=(GeneralLatticePoint const &) = delete;
  GeneralLatticePoint &operator=(GeneralLatticePoint &&)      = default;
protected:
  ~GeneralLatticePoint()                                      = default;
  // clang-format on
public:
  // This one should be overloaded by every derived class. Used for diagnostic.
  static std::string class_name() { return "General Lattice Point."; }

  // This may exists in an overload (if Has_ExposesCoo is true).
  template <class Arg> inline void operator[](Arg &&arg) = delete;

  // comparison with < or > are by length (by default).
  // Note that == or != are intendend to mean actual equality comparisons,
  // so P1 <= P2 <= P1 does not imply P1 == P2.

  // We require < and > to be strict weak orderings.
  // MAKE SURE THIS IS A STRICT WEAK ORDERING EVEN IF LatP HAS APPROXIMATIONS.

  // Be aware that P1.get_norm2() < P2.get_norm2() *may* give a different result than
  // P1 < P2 due to approximation errors.
  // (currently, these types of approximations are not used)

  // You should not overload <, >, <=, >=, !=
  // You may overload +=, -=, *=, unary-, ==

  // clang-format off
  template <class LatP2> inline bool operator< (LatP2 const &rhs) const;
  template <class LatP2> inline bool operator> (LatP2 const &rhs) const;
  template <class LatP2> inline bool operator<=(LatP2 const &rhs) const;
  template <class LatP2> inline bool operator>=(LatP2 const &rhs) const;
  // clang-format on
  // arithmetic:

  // Remark: The class Impl=LatP paramter serves the single purpose that the argument
  // to TEMPL_RESTRICT_DECL is not always false for any possible set of template arguments.
  // (Only the template parameters of the member count, not of the class, so LatP is a fixed type
  // is this context and not a template argument)
  // Otherwise, the SFINAE techniques behind it won't work.
  // More specifically, when the compiler is instantiation the class GeneralLatticePoint<LatP> for a
  // fixed LatP, it will already know LatP. If, based on that information, the compiler is able to
  // short-circuit the argument to TEMPL_RESTRICT_* to false / false_type, you may get an error.
  // This is prevented by the Impl=LatP, because LatP is only the default and it does not know
  // whether we stick to that.
  // Of course, we never use these templates with
  // Impl!=LatP and we usually even static_assert(Impl==LatP) inside the implementation.

  template <class LatP2, class Impl = LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline LatP &operator+=(LatP2 const &x2);  // default implementation asserts linearity

  template <class LatP2, class Impl = LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline LatP &operator-=(LatP2 const &x2);  // default implementation asserts linearity

  // add_multiply performs the same as +=(x2*multiplier), but is possibly more efficient.
  // sub_multiply performs the same as -=(x2*mutliplier), but is possibly more efficient.
  // The default implementation has the same restrictions as += resp. -=
  // TODO : may consider variant where the scalar product of *this and x2 is known.
  template <class LatP2, class Integer, class Impl = LatP,
            TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>, std::is_integral<Integer>)>
  inline void add_multiply(LatP2 const &x2, Integer const multiplier);

  template <class LatP2, class Integer, class Impl = LatP,
            TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>, std::is_integral<Integer>)>
  inline void sub_multiply(LatP2 const &x2, Integer multiplier);

  template <class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatP2>>)>
  inline bool operator!=(LatP2 &&x2) const
  {
    return !(CREALTHIS->operator==(std::forward<LatP2>(x2)));
  };

  template <class LatP2, class Impl = LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline bool operator==(LatP2 const &x2) const;  // default implementation asserts internal rep

  template <class Integer, class Impl = LatP, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline LatP &operator*=(Integer const multiplier);  // default implementation asserts linearity
  inline LatP &operator*=(mpz_class const &multiplier) = delete;  // not implemented yet

  inline LatP operator-() &&;  // unary-
  // operator+ is entirely defined in terms of +=
  // operator- is entirely defined in terms of -=
  // operator* is entirely defined in terms of *=
  // Definitions are out-of-class and not supposed to be overloaded.

  // get_internal_rep_size() returns the number of coordinates we get when using get_internal_rep
  // get_dim returns the (ambient) dimension the vector is supposed to represent.
  // By default, vec_size is the same as dim.

  // get_dim must be overloaded.
  // get_internal_rep_size may be overloaded.

  // Note: Overload may have different return type.
  template <class Impl = LatP>  // default implementation asserts Has_InternalRep
  inline auto get_internal_rep_size() const -> decltype(std::declval<Impl>().get_dim());

  template <class Arg, class Impl = LatP>
  inline RepCooType const &get_internal_rep(Arg &&arg) const;

  // clang-format off
  template <class Arg, class Impl = LatP>
  inline RepCooType &get_internal_rep(Arg &&arg);

  template <class Arg, class Impl = LatP>
  inline AbsoluteCooType get_absolute_coo(Arg &&arg) const;
  // clang-format on

  // Note that the overload shall NOT have void return type. It may be static / constexpr.
  void get_dim() const = delete;

  /**
    Used for output to stream. Note that operator<< calls this (or an overloaded version)
    Supposed to be overloaded for the specific class.
  */

  inline std::ostream &write_lp_to_stream(std::ostream &os, bool const include_norm2 = true,
                                          bool const include_approx = true) const;

  template <class Impl = LatP, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<Impl>)>
  inline std::ostream &write_lp_rep_to_stream(std::ostream &os) const;

  std::istream &read_from_stream(std::istream &is) = delete;

  /**
    Fills a lattice point with zeros.
    Note that freshly constructed lattice points may contain uninitialized values unless this
    function is used.

    The latter depends subtly on the constructors used (empty vs. default constructor
    and whether we use POD types) May be overloaded by Derived class.
  */
  template <class Impl = LatP> inline void fill_with_zero();

  /**
    Changes vector from v to -v. May be overloaded.
  */
  template <class Impl = LatP> inline void make_negative();

  /**
    Tests whether a lattice point is all-zero.
    May be overloaded.
  */
  template <class Impl = LatP> inline bool is_zero() const;

  /**
    Makes an (explicit) copy of the current point.

    Default: Copying components as a vector and sanitize additional data.
    (The latter probably means needless recomputing of data that the original already had)
    With the intended semantics of [] and entries, this results in a deep copy.

    May be overloaded.
  */

  // Calling this is probably an error.
  template <class Impl = LatP> inline LatP make_copy() const &;
  [[deprecated]] constexpr inline LatP make_copy() const && { return *CREALTHIS; }

  /**
    brings the lattice point into a defined state. Defaults to "do nothing".

    Note that this function must *not* be const, since it actually changes observable behaviour:
    For efficiency reasons, sanitize is not called upon write access to LatP[]; rather,
    sanitize is called from outside the class.
    As such, the user can leave LatP in an invalid state and sanitize remedies that.

    The second version takes norm2 as an argument (to avoid recomputing it).
  */

  void sanitize(){};
  inline void sanitize(ScalarProductStorageType const &norm2) { REALTHIS->sanitize(); }

  // must be overloaded if meaningful
  // (i.e. if ScalarProductStorageType_Full != ScalarProductStorageType)
  template <class Impl = LatP,
            TEMPL_RESTRICT_DECL(!(std::is_same<Get_ScalarProductStorageType<Impl>,
                                               Get_ScalarProductStorageType_Full<Impl>>::value))>
  inline void sanitize(ScalarProductStorageType_Full const &norm2_full) = delete;

  /**
    obtains the squared length.
    By default, we use the scalar product to compute it.
    Usually norm2 will be precomputed, so this function should be overridden by LatP.
    Overload may have different return type, but convertible to ScalarProductStorageType.
    If Trait_CheapNorm2 is set, this function *must* be overloaded.
  */

  inline ScalarProductStorageType get_norm2() const;

  /**
    The functions below are overloaded in ApproximatedPoint.h for points with (non-bit)
    Approximations.
    The _at_level variant gets the value at the given approximation level.
    The _full variant gets an object holding all levels.
    Note that for the plain and _exact variants, the return types might differ in the overloads,
    but need to be convertible to ScalarProductStorageType.
  */

  template <unsigned int level, class Impl = LatP>
  inline ScalarProductStorageType get_norm2_at_level() const
  {
    IMPL_IS_LATP;
    static_assert(level == 0, "Default only has level 0");
    return CREALTHIS->get_norm2();
  }

  template <class Impl = LatP> inline ScalarProductStorageType_Full get_norm2_full() const
  {
    static_assert(Has_Leveled<Impl>::value == false, "Need to overload");
    return CREALTHIS->get_norm2();
  }

  // don't call directly. We use compute_sc_product(x1,x2) for a more symmetric syntax.
  // However, out-of-class definitions get messy with overloading.

  // TODO: Make protected and befriend compute_* functions
  // TODO: Better std::forwarding (note const-ness)

  template <class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline ScalarProductStorageType do_compute_sc_product(LatP2 const &x2) const;

  template <unsigned int level, class Impl = LatP, class LatP2>
  inline ScalarProductStorageType do_compute_sc_product_at_level(LatP2 const &x2) const
  {
    IMPL_IS_LATP;
    static_assert(level == 0, "Default only has level 0");
    return CREALTHIS->do_compute_sc_product(x2);
  }
  template <class Impl = LatP, class LatP2>
  inline ScalarProductStorageType_Full do_compute_sc_product_full(LatP2 const &x2) const
  {
    IMPL_IS_LATP;
    static_assert(Has_Leveled<Impl>::value == false, "Need to overload");
    return CREALTHIS->do_compute_sc_product(x2);
  }

public:
  void access_bitapproximation(unsigned int level) = delete;  // need to overload
  void take_bitapproximations() &&                 = delete;  // need to overload
  void update_bitapprox()                          = delete;
};

/**
 Non-member functions
 */

// compute scalar products. These are the functions that users of lattice points should use.
// Users should only ever have to use the plain compute_sc_product and the _bitapprox versions.

// the other variants are only relevant when having approximations:
// Notably: compute_sc_product may return a Wrapper than encapsulates a lazy evalution.
// The _full variant computes an exact and an approximate scalar product and returns a type capable
// of holding both.
// The _at_level<level> variant computes the approximation at the given level (0 for exact)

template <class LP1, class LP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LP1>>,
                                                     IsALatticePoint<mystd::decay_t<LP2>>)>
inline auto compute_sc_product(LP1 &&lp1, LP2 &&lp2)
    -> decltype(std::declval<LP1>().do_compute_sc_product(std::declval<LP2>()))
{
  return std::forward<LP1>(lp1).do_compute_sc_product(std::forward<LP2>(lp2));
}

template <unsigned int level, class LP1, class LP2,
          TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LP1>>,
                               IsALatticePoint<mystd::decay_t<LP2>>)>
inline auto compute_sc_product_at_level(LP1 &&lp1, LP2 &&lp2)
    -> decltype(std::declval<LP1>().do_compute_sc_product_at_level<level>(std::declval<LP2>()))
{
  return std::forward<LP1>(lp1).template do_compute_sc_product_at_level<level>(
      std::forward<LP2>(lp2));
}

template <class LP1, class LP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LP1>>,
                                                     IsALatticePoint<mystd::decay_t<LP2>>)>
inline auto compute_sc_product_full(LP1 &&lp1, LP2 &&lp2)
    -> decltype(std::declval<LP1>().do_compute_sc_product_full(std::declval<LP2>()))
{
  return std::forward<LP1>(lp1).do_compute_sc_product_full(std::forward<LP2>(lp2));
}

// this function can be used to initialize an LP with container types that allow []-access.
// Note that there is an explicit static_cast to LP's entry types.
// In particular, this can convert mpz_t to mpz_class...
//
// Note: Return type does not take part in template argument deduction.
// Usage: make_from_any_vector<TargetType>(source_container, dim).

template <class LP, class SomeContainer, class DimType,
          TEMPL_RESTRICT_DECL2(IsALatticePoint<LP>, IsRepLinear_RW<LP>)>
LP make_from_any_vector(SomeContainer const &container, DimType dim)
{
  static_assert(DoesDeclareCoordinateType<LP>::value, "Not declaring coordinate types");
  using ET = Get_CoordinateType<LP>;
  DEBUG_TRACEGENERIC("generically converting vector to LP for" << LP::class_name())
  LP result(dim);
  //  auto dim = result.get_dim();
  for (uint_fast16_t i = 0; i < dim; ++i)
  {
    result[i] = static_cast<ET>(container[i]);
  }
  result.sanitize();
  return result;
}

// Same as above, but un-Z_NR's the container.

template <class LP, class SomeZNRContainer, class DimType,
          TEMPL_RESTRICT_DECL2(IsALatticePoint<LP>, IsRepLinear_RW<LP>)>
LP make_from_znr_vector(SomeZNRContainer const &container, DimType dim)
{
  static_assert(DoesDeclareCoordinateType<LP>::value, "Not declaring coordinate types");
  using ET = Get_CoordinateType<LP>;
  DEBUG_TRACEGENERIC("generically converting vector to LP and un-ZNRing for" << LP::class_name())
  LP result(dim);
  for (uint_fast16_t i = 0; i < dim; ++i)
  {
    result[i] = ConvertMaybeMPZ<ET>::convert_to_inttype(container[i].get_data());
  }
  result.sanitize();
  return result;
}

}  // end namespace GaussSieve

#include "LatticePointGeneric.h"

// TODO: Redeclare in LatticePointGeneric.h
// TODO: Rename these

// cleaning up internal macros.

#undef IMPL_IS_LATP
#undef CREALTHIS
#undef REALTHIS

#endif  // include guard
