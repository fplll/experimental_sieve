#ifndef LATTICE_POINT_CONCEPT_H
#define LATTICE_POINT_CONCEPT_H

// TODO: Remove the concept idea below, we replace it by CRTP (i.e. inheritance without virtual
// functions, using templates)

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include <cmath>
#include <cstdint>
#include <gmpxx.h>
#include <array>
#include "GlobalStaticData.h"

#include <bitset> //for approximation
#include <boost/dynamic_bitset.hpp> //for approximation

//#include "RelevantCoords.h"


// clang-format off

namespace GaussSieve{

/**
  GeneralLatticePoint defines an interface for the various classes of lattice points we use.
  Notably, we define a class template GeneralLatticePoint<Actual class> from which all
  lattice point classes inherit.
  To avoid overhead, we use CRTP rather than dynamic polymorphism, i.e.
  GeneralLatticePoint is templated by its child classes.

  Note that we have a lot of "default" implementations, which are defined depending on various
  lattice traits. Mosto of the actual code of these implementations is found in
  LatticePointGeneric.h

  GeneralLatticePoint also specifies which methods an implementation is allowed to override.
*/


// This class template stores the trait typedefs that the individual lattice point classes have
// There has to be a specialization for each lattice point class.
// The general template must never be instantiated.
// Note that we can not put these traits into the lattice points classes directly, because that
// would cause circular references due to CRTP.

// Cf. ExacLatticePoint.h for an example of how a concrete specialization is supposed to look like.

template<class LatticePoint> struct LatticePointTraits
{
  public:

  using Invalid=std::true_type; // do not set this at all in an specialization.
  LatticePointTraits(...) = delete;
  // static_assert(false) is invalid due to subtleties of C++, even if it may work on some compilers

  using Trait_ScalarProductStorageType = void;
  using Trait_ScalarProductStorageType_Full = void;
  using Trait_ExposesCoos = std::false_type;
  using Trait_CoordinateType = void;
  using Trait_Coos_RW = std::false_type;
  using Trait_AbsoluteCooType = void;
  using Trait_RepCooType = void;
  using Trait_ExposesInternalRep = std::false_type;
  using Trait_InternalRepLinear = std::false_type;
  using Trait_InternalRep_RW = std::false_type;
  using Trait_InternalRepByCoos = std::false_type;
  using Trait_InternalRepIsAbsolute = std::false_type;
  using Trait_AbsoluteCoos = std::false_type;
  using Trait_CheapNorm2 = std::false_type;
  using Trait_CheapNegate = std::false_type;
//  using Trait_Approximations = std::false_type;
  using Trait_Leveled = std::false_type;
  using Trait_BitApprox = std::false_type;
  using Trait_ApproxLevel = std::integral_constant<unsigned int, 0>;
};

/**
  Available traits
  (everything that needs to be set to true_type defaults to false, unless specified otherwise)
  These traits are used to selectively enable some meaningful default operations on lattice points.
  Note: In the actual LatticePointTraits class, the traits are actually prefixed with Trait_
  (i.e. LatticePointTraits<LP> contains using Trait_ExposesCoos etc.)
  To retrieve a trait, use Has_TraitName for binary traits.
  or Get_TraitName for non-binary traits

  ScalarProductStorageType: A type that can hold the result of a scalar product computation. Mandatory.
                           Note that the result from a scalar product computation might actually differ.
                           (due to delayed evaluation)

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

  RepCooType      : return type of get_internal_rep (if available)
                    defaults to CoordinateType

  ExposesInternalRep :  Indicates that get_internal_rep(i), get_internal_rep_size() exists and may be
                        read for 0<= i < get_internal_rep_size().
                        We are guaranteed that these entries determine the point.
                        implied by InternalRepVector_R, InternalRepVector_RW, InternalRepByCoos, InternalRepIsAbsolute

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
              There is a public member functions:
              get_bitapprox_norm2() : Returns the number of bits that we are using
              - defaults to get_dim()
              There is a member function:
              do_compute_sc_product_bitapprox(LatP const &): computes the bit-wise scalar product of
              bit-approximations of *this with LatP.
              Use compute_sc_product_bitapprox(LatP const & x1, LatP const &x2) to actually compute
              the bit-approx scalar product.


// Does not work ATM, might be needed later...
//  AccessNorm2:  Set to true_type to indicate that we we have an access_norm2() function to
//                return a const-reference to a precomputed norm2.

  CheapNegate: Set to true_type to indicate that negation needs no sanitize().

// TODO: maybe deprecate this trait.
//  Approximations: Set to true_type to indicate that the point has approximations.

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

    NOTE: Approximations and ApproxLevel do not relate to Bitapproximation.
*/

// forward declaration
template<class Implementation> class GeneralLatticePoint;

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

// These check for the existence of the trait in the Trait class:
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(LatticePointTag, std::true_type, IsALatticePoint);

CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, Trait_ScalarProductStorageType, DeclaresScalarProductStorageType);
CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, Trait_CoordinateType, DoesDeclareCoordinateType);

CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Invalid, std::true_type, HasNoLPTraits);

// should enclose these as static members of struct, declared as friend to Traits classes.
//namespace TraitHelpers // check whether trait exists and equals std::true_type.
//{
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_ExposesCoos, std::true_type, T_ExposesCoos);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_Coos_RW, std::true_type, T_Coos_RW);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_ExposesInternalRep, std::true_type, T_ExposesInternalRep);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_InternalRepLinear, std::true_type, T_InternalRepLinear);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_InternalRep_RW, std::true_type, T_InternalRep_RW);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_InternalRepByCoos, std::true_type, T_InternalRepByCoos);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_InternalRepIsAbsolute, std::true_type, T_InternalRepIsAbsolute);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_AbsoluteCoos, std::true_type, T_AbsoluteCoos);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_CheapNorm2, std::true_type, T_CheapNorm2);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_CheapNegate, std::true_type, T_CheapNegate);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_Leveled, std::true_type, T_Leveled);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_AccessNorm2, std::true_type, T_AccessNorm2);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Trait_BitApprox, std::true_type, T_BitApprox);

template<class T,class = int>
struct T_ApproxLevelOf
{
static constexpr unsigned int value = 0;
};

// the , inside the decltype is the comma operator for void (hence not overloaded)
// this template specialization is only valid if T::ApproxLevel exists and will be preferred over
// the general one above.
template<class T>
struct T_ApproxLevelOf<T, decltype( static_cast<void>(T::ApproxLevel), static_cast<int>(0))>
{
static constexpr unsigned int value = T::ApproxLevel;
};

template<class T>
using ApproxLevelOf = std::integral_constant<unsigned int,T_ApproxLevelOf<T>::value>;

//}
// Retrieving Traits with defaults:

// because the macro below does not like commas in arguments...
using Zero_Constant = std::integral_constant<unsigned int,0>;
MAKE_TRAIT_GETTER(LatticePointTraits, Trait_ApproxLevel, Zero_Constant, Get_ApproxLevel);
MAKE_TRAIT_GETTER(LatticePointTraits, Trait_CoordinateType, void, Get_CoordinateType);
MAKE_TRAIT_GETTER(LatticePointTraits, Trait_ScalarProductStorageType, void, Get_ScalarProductStorageType);
// ClassToCheck is the argument of the constructed Traits getter inside the macro def. This makes it
// default to GetScalarproductStorageType
MAKE_TRAIT_GETTER(LatticePointTraits, Trait_ScalarProductStorageType_Full,
  Get_ScalarProductStorageType<ClassToCheck>, Get_ScalarProductStorageType_Full);
MAKE_TRAIT_GETTER(LatticePointTraits, Trait_AbsoluteCooType,
  Get_CoordinateType<ClassToCheck>, Get_AbsoluteCooType);
MAKE_TRAIT_GETTER(LatticePointTraits, Trait_RepCooType,
  Get_CoordinateType<ClassToCheck>, Get_RepCooType);

//{
//using namespace TraitHelpers;
// These are what the rest of the code should be actually using:
template<class LatP> using Has_ExposesCoos = mystd::bool_constant<
  T_ExposesCoos<LatP>::value || T_InternalRepByCoos<LatP>::value ||
  (!std::is_void<Get_CoordinateType<LatP>>::value) || T_Coos_RW<LatP>::value || T_AbsoluteCoos<LatP>::value >;

template<class LatP> using Has_Coos_RW = mystd::bool_constant<
  T_Coos_RW<LatP>::value || (T_InternalRepByCoos<LatP>::value && T_InternalRep_RW<LatP>::value)>;

template<class LatP> using Has_ExposesInternalRep = mystd::bool_constant<
  T_ExposesInternalRep<LatP>::value || T_InternalRepLinear<LatP>::value || T_InternalRep_RW<LatP>::value ||
  T_InternalRepByCoos<LatP>::value || T_InternalRepIsAbsolute<LatP>::value >;

template<class LatP> using Has_InternalRepLinear = mystd::bool_constant<
  T_InternalRepLinear<LatP>::value>;

template<class LatP> using Has_InternalRep_RW = mystd::bool_constant<
  T_InternalRep_RW<LatP>::value || (T_InternalRepByCoos<LatP>::value && T_Coos_RW<LatP>::value)>;

template<class LatP> using Has_InternalRepByCoos = mystd::bool_constant<
  T_InternalRepByCoos<LatP>::value>;

template<class LatP> using Has_InternalRepIsAbsolute = mystd::bool_constant<
  T_InternalRepIsAbsolute<LatP>::value || ( T_InternalRepByCoos<LatP>::value && T_AbsoluteCoos<LatP>::value)>;

template<class LatP> using Has_AbsoluteCoos = mystd::bool_constant<
  T_AbsoluteCoos<LatP>::value || (T_InternalRepByCoos<LatP>::value && T_InternalRepIsAbsolute<LatP>::value)>;

template<class LatP> using Has_CheapNorm2 = mystd::bool_constant<
  T_CheapNorm2<LatP>::value>;

template<class LatP> using Has_CheapNegate = mystd::bool_constant<
  T_CheapNegate<LatP>::value>;

template<class LatP> using Has_Leveled = mystd::bool_constant<
  T_Leveled<LatP>::value>;

template<class LatP> using Has_AccessNorm2 = mystd::bool_constant<
  T_AccessNorm2<LatP>::value>;

template<class LatP> using Has_BitApprox = mystd::bool_constant<
  T_BitApprox<LatP>::value>;

// Deprecated
template<class LatP> using IsRepLinear_RW = mystd::bool_constant<
Has_InternalRepLinear<LatP>::value && Has_InternalRep_RW<LatP>::value>;

//}

#define IMPL_IS_LATP \
static_assert(std::is_same<Impl,LatP>::value,"Using template member function with wrong type.")

// unsure whether to use reinterpret_casts here...

#define CREALTHIS static_cast<LatP const*>(this)
#define REALTHIS  static_cast<LatP*>(this)

template<class LatP>
class GeneralLatticePoint
{
    static_assert(!HasNoLPTraits<LatP>::value, "Trait class not specialized.");
    static_assert(DeclaresScalarProductStorageType<LatP>::value,
                  "Lattice Point class does not typedef its scalar product type");
    static_assert(!Has_ExposesCoos<LatP>::value || DoesDeclareCoordinateType<LatP>::value,
                  "Lattice Point exposes coos, but does not tell its type.");

    friend LatP; // makes children able to access private (in addition to protected) members.
                 // since the constructor is private, this enforces correct usage.
                 // (Note that it may prevent multi-level inheritance)
    public:
    using ScalarProductStorageType      = Get_ScalarProductStorageType<LatP>;
    using ScalarProductStorageType_Full = Get_ScalarProductStorageType_Full<LatP>;
    // deprecated
    using CooType         = Get_CoordinateType<LatP>; //may be void
    using CoordinateType  = Get_CoordinateType<LatP>; //may be void
//    deprecated
    using AbsCooType      = Get_AbsoluteCooType<LatP>;
    using AbsoluteCooType = Get_AbsoluteCooType<LatP>;
    using RepCooType      = Get_RepCooType<LatP>;
    static constexpr unsigned int ApproxLevel = Get_ApproxLevel<LatP>::value;
    static_assert((Has_Leveled<LatP>::value == true) || (ApproxLevel==0 ),"Has ApproxLevel>0, but does not declare Trait_Leveled.");

    private:
    // Empty base class, only callable from its friends (i.e. from LatP)
    explicit constexpr GeneralLatticePoint()=default;
    public:

    // This is just to match the implementation of a typical instantiation.
    // Note the the deleted copy constructors and copy assignments prevents default copying
    // derived classes (since the (empty) base class must be copied as well).
    GeneralLatticePoint(GeneralLatticePoint const &other)=delete;
    GeneralLatticePoint(GeneralLatticePoint &&other)=default;
    GeneralLatticePoint& operator=(GeneralLatticePoint const & other) = delete;
    GeneralLatticePoint& operator=(GeneralLatticePoint && other) = default;
    protected:
    ~GeneralLatticePoint()=default;
    public:

// This one should be overloaded by every derived class. Used for diagnostic.
    static std::string class_name() {return "General Lattice Point.";};

    template<class Arg>
    void operator[](Arg &&arg) = delete; // This may exists in an overload (if Has_ExposesCoo is true).

// comparison with < or > are by length (by default).
// Note that == or != are intendend to mean actual equality comparisons,
// so P1 <= P2 <= P1 does not imply P1==P2.

// We require < and > to be strict weak orderings.
// MAKE SURE THIS IS A STRICT WEAK ORDERING EVEN IF LatP HAS APPROXIMATIONS.

// Be aware that P1.get_norm2() < P2.get_norm2() *may* give a different result than
// P1 < P2 due to approximation errors.

// You should not overload <, >, <=, >=, !=
// You may overload +=, -=, *=, unary-, ==

    inline bool operator< ( LatP const &rhs ) const;
    inline bool operator> ( LatP const &rhs ) const;
    inline bool operator<=( LatP const &rhs ) const;
    inline bool operator>=( LatP const &rhs ) const;
// arithmetic:

// Remark: The class Impl=LatP paramter serves the single purpose that the argument
// to TEMPL_RESTRICT_DECL is not always false for any possible set of template arguments.
// (Note that LatP is a fixed type is this context and not a template argument)
// Otherwise, the SFINAE magic behind it won't work.
// More specifically, when the compiler is instantiation the class GeneralLatticePoint<LatP> for a
// fixed LatP, it will already know LatP. If, based on that information, the compiler is able to
// short-circuit the argument to TEMPL_RESTRICT_* to false / false_type, you may get an error.
// This is prevented by the Impl=LatP, because LatP is only the default and it does not know whether
// we stick to that.
// Of course, we never use these templates with
// Impl!=LatP and we usually even static_assert(Impl==LatP) inside the implementation.

    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
    inline LatP& operator+=(LatP2 const &x2); // default implementation asserts linearity

    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
    inline LatP& operator-=(LatP2 const &x2); // default implementation asserts linearity

    template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatP2>>)>
    inline bool operator!=(LatP2 && x2) const {return !(CREALTHIS->operator==(std::forward<LatP2>(x2)));};

    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
    inline bool operator==(LatP2 const &x2) const; // default implementation asserts internal rep

    template<class Integer, class Impl=LatP, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
    inline LatP& operator*=(Integer const multiplier); // default implementation asserts linearity
    inline LatP& operator*=(mpz_class const &multiplier) = delete; // not implemented yet

    inline LatP operator-() &&; //unary-
    // operator+ is entirely defined in terms of +=
    // operator- is entirely defined in terms of -=
    // operator* is entirely defined in terms of *=
    // Definitions are out-of-class and not supposed to be overloaded.

    // get_internal_rep_size() returns the number of coordinates we get when using get_internal_rep
    // get_dim returns the (ambient) dimension the vector is supposed to represent.
    // By default, vec_size is the same as dim.

    // get_dim must be overloaded.
    // get_internal_rep_size may be overloaded.
    template<class Impl=LatP> // default implementation asserts Has_InternalRep
    inline auto get_internal_rep_size() const -> decltype( std::declval<Impl>().get_dim() ); // Note: Overload may have different return type.

    template<class Arg, class Impl=LatP>
    inline RepCooType const & get_internal_rep(Arg &&arg) const;

    template<class Arg, class Impl=LatP>
    inline RepCooType & get_internal_rep(Arg &&arg);

    template<class Arg, class Impl=LatP>
    inline AbsoluteCooType get_absolute_coo(Arg &&arg) const;

    void get_dim() const = delete; // Note that the overload shall NOT have void return type. It may be static / constexpr.
                                   // It's just not possible to specific it here w/o C++14 auto.

/**
  Used for output to stream. Note that operator<< calls this (or an overloaded version)
  Supposed to be overloaded for the specific class.
*/

    inline std::ostream& write_lp_to_stream(std::ostream &os, bool const include_norm2=true, bool const include_approx=true) const;

    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(Has_ExposesInternalRep<Impl>)>
    inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const;

    std::istream& read_from_stream(std::istream &is) = delete;

/**
  Fills a lattice point with zeros.
  Note that freshly constructed lattice points may contain uninitialized values unless this function is used.

  The latter depends subtly on the constructors used (empty vs. default constructor...)
  May be overloaded by Derived class.
*/
    template<class Impl=LatP> inline void fill_with_zero();

/**
  Changes vector from v to -v. May be overloaded.
*/
    template<class Impl=LatP> inline void make_negative();

/**
  Tests whether a lattice point is all-zero.
  May be overloaded.
*/
    template<class Impl=LatP> inline bool is_zero() const;

/**
  Makes an (explicit) copy of the current point.

  Default: Copying components as a vector and sanitize additional data.
  (The latter probably means needless recomputing of data that the original already had)
  With the intended semantics of [] and entries, this results in a deep copy.

  May be overloaded.
*/

    template<class Impl=LatP> inline LatP make_copy() const &;
    [[deprecated]] constexpr inline LatP make_copy() && {return *CREALTHIS;} // Calling this is probably an error.


/**
     brings the lattice point into a defined state. Defaults to "do nothing".

     Note that this function must *not* be const, since it actually changes observable behaviour:
     For efficiency reasons, sanitize is not called upon write access to LatP[]; rather,
     sanitize is called from outside the class.
     As such, the user can leave LatP in an invalid state and sanitize remedies that.

     The second version takes norm2 as an argument (to avoid recomputing it).
*/

// TODO: Arg type overloads for Full StorageType?
    void sanitize() {};
    inline void sanitize(ScalarProductStorageType const &norm2) { REALTHIS->sanitize(); };

    template<class Impl=LatP,TEMPL_RESTRICT_DECL(!(std::is_same<Get_ScalarProductStorageType<Impl>,Get_ScalarProductStorageType_Full<Impl>>::value))>
    inline void sanitize(ScalarProductStorageType_Full const &norm2_full) = delete; // must be overloaded if meaningful.


/**
     obtains the squared length.
     By default, we use the scalar product to compute it.
     Usually norm2 will be precomputed, so this function should be overridden by LatP.
     Overload may have different return type, but convertible to ScalarProductStorageType.
     If Trait_CheapNorm2 is set, this function *must* be overloaded.
*/

    inline ScalarProductStorageType get_norm2() const;

/** The functions below are overloaded in ApproximatedPoint.h for points with (non-bit) Approximations.
    The _at_level variant gets the value at the given approximation level.
    The _full variant gets an object holding all levels.
    Note that for the plain and _exact variants, the return types might differ in the overloads,
    but need to be convertible to ScalarProductStorageType.
*/

    template<unsigned int level, class Impl=LatP>
    inline ScalarProductStorageType get_norm2_at_level() const
    {
      IMPL_IS_LATP;
      static_assert(level==0,"Default only has level 0");
      return CREALTHIS->get_norm2();
    }

    template<class Impl=LatP>
    inline ScalarProductStorageType_Full get_norm2_full() const
    {
      static_assert(Has_Leveled<Impl>::value == false, "Need to overload");
      return CREALTHIS->get_norm2();
    }

    // don't call directly. We use compute_sc_product(x1,x2) for a more symmetric syntax.
    // However, out-of-class definitions get messy with overloading.

    // TODO: Make protected and befriend compute_* functions

    inline ScalarProductStorageType do_compute_sc_product(LatP const &x2) const;

    template<unsigned int level,class Impl=LatP>
    inline ScalarProductStorageType do_compute_sc_product_at_level(LatP const &x2) const
    {
      IMPL_IS_LATP;
      static_assert(level==0,"Default only has level 0");
      return CREALTHIS->do_compute_sc_product(x2);
    }
    template<class Impl=LatP>
    inline ScalarProductStorageType_Full do_compute_sc_product_full(LatP const &x2) const
    {
      IMPL_IS_LATP;
      static_assert(Has_Leveled<Impl>::value==false, "Need to overload");
      return CREALTHIS->do_compute_sc_product(x2);
    }
    public:

    template<class Impl=LatP>
    inline auto get_bitapprox_norm2() const -> decltype( std::declval<Impl>().get_dim() ); // Note: Overload may have different return type.

    // This function is only ever called from within a block protected by
    // if(templated constexpr resuting in false){ }. Lacking if constexpr, we bail out at runtime.
    // Use CPP17CONSTEXPRIF macro for the if.
    #if __if_constexpr
    template<class Impl=LatP, class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatP2>>)>
    inline int do_compute_sc_product_bitapprox(LatP2 const &) const = delete;
    #else
    template<class Impl=LatP, class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatP2>>)>
    inline int do_compute_sc_product_bitapprox(LatP2 const &) const { assert(false); }
    #endif

    #if __if_constexpr
    template<class Impl=LatP, class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatP2>>)>
    inline int do_compute_sc_product_bitapprox_2nd_order(LatP2 const &) const = delete;
    #else
    template<class Impl=LatP, class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatP2>>)>
    inline int do_compute_sc_product_bitapprox_2nd_order(LatP2 const &) const { assert(false); }
    #endif


    //FOR SIM-HASH
    #if __if_constexpr
      template<class Impl=LatP, class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<typename std::decay<LatP2>::type>)>
      inline int do_compute_sc_product_bitapprox_fixed(LatP2 const &) const = delete;
    #else
      template<class Impl=LatP, class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<typename std::decay<LatP2>::type>)>
      inline int do_compute_sc_product_bitapprox_fixed(LatP2 const &) const { assert(false); }
    #endif

    //FOR SIM-HASH 2nd order
    #if __if_constexpr
      template<class Impl=LatP, class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<typename std::decay<LatP2>::type>)>
      inline int do_compute_sc_product_bitapprox_fixed2(LatP2 const &) const = delete;
    #else
      template<class Impl=LatP, class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<typename std::decay<LatP2>::type>)>
      inline int do_compute_sc_product_bitapprox_fixed2(LatP2 const &) const { assert(false); }
    #endif
 };

 /**
  Non-member functions
  */

// compute scalar products. These are the functions that users of lattice points should use.
// Users should only ever have to use the plain compute_sc_product and the _bitapprox versions.

// the other variants are only relevant when having approximations:
// Notably: compute_sc_product may return a Wrapper than encapsulates a lazy evalution.
// The _full variant computes an exact and an approximate scalar product and returns a type capable of holding both.
// The _at_level<level> variant computes the approximation at the given level (0 for exact)

template<class LP1,class LP2, TEMPL_RESTRICT_DECL2(
  IsALatticePoint<mystd::decay_t<LP1>>,IsALatticePoint<mystd::decay_t<LP2>>)>
inline auto compute_sc_product(LP1 &&lp1, LP2 &&lp2)
-> decltype( std::declval<LP1>().do_compute_sc_product(std::declval<LP2>() ) )
{ return std::forward<LP1>(lp1).do_compute_sc_product(std::forward<LP2>(lp2)); }

template<unsigned int level, class LP1, class LP2, TEMPL_RESTRICT_DECL2(
  IsALatticePoint<mystd::decay_t<LP1>>,IsALatticePoint<mystd::decay_t<LP2>>)>
inline auto compute_sc_product_at_level(LP1 &&lp1, LP2 &&lp2)
-> decltype( std::declval<LP1>().do_compute_sc_product_at_level<level>(std::declval<LP2>() ))
{ return std::forward<LP1>(lp1).template do_compute_sc_product_at_level<level>(std::forward<LP2>(lp2)); }

template<class LP1, class LP2, TEMPL_RESTRICT_DECL2(
  IsALatticePoint<mystd::decay_t<LP1>>,IsALatticePoint<mystd::decay_t<LP2>>)>
inline auto compute_sc_product_full(LP1 &&lp1, LP2 &&lp2)
-> decltype( std::declval<LP1>().do_compute_sc_product_full(std::declval<LP2>() ) )
{ return std::forward<LP1>(lp1).do_compute_sc_product_full(std::forward<LP2>(lp2)); }

template<class LP1, class LP2, TEMPL_RESTRICT_DECL2(
IsALatticePoint<mystd::decay_t<LP1>>,IsALatticePoint<mystd::decay_t<LP2>>)>
inline auto compute_sc_product_bitapprox(LP1 &&lp1, LP2 &&lp2)
-> decltype( std::declval<LP1>().do_compute_sc_product_bitapprox(std::declval<LP2>() )  )
{ return std::forward<LP1>(lp1).do_compute_sc_product_bitapprox(std::forward<LP2>(lp2)); }

template<class LP1, class LP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LP1>,IsALatticePoint<LP2>)>
inline auto compute_sc_product_bitapprox_2nd_order(LP1 const &lp1, LP2 const &lp2)
-> decltype( std::declval<LP1>().do_compute_sc_product_bitapprox_2nd_order(std::declval<LP2>() )  )
{ return lp1.do_compute_sc_product_bitapprox_2nd_order(lp2); }


//FOR SIM-HASH
template<class LP1, class LP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LP1>,IsALatticePoint<LP2>)>
inline auto compute_sc_product_bitapprox_fixed(LP1 const &lp1, LP2 const &lp2)
-> decltype( std::declval<LP1>().do_compute_sc_product_bitapprox_fixed(std::declval<LP2>() )  )
{ return lp1.do_compute_sc_product_bitapprox_fixed(lp2); }

//FOR SIM-HASH 2nd order
template<class LP1, class LP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LP1>,IsALatticePoint<LP2>)>
inline auto compute_sc_product_bitapprox_fixed2(LP1 const &lp1, LP2 const &lp2)
-> decltype( std::declval<LP1>().do_compute_sc_product_bitapprox_fixed2(std::declval<LP2>() )  )
{ return lp1.do_compute_sc_product_bitapprox_fixed2(lp2); }





// this function can be used to initialize an LP with container types that allow []-access.
// Note that there is an explicit static_cast to LP's entry types.
// In particular, this can convert mpz_t to mpz_class...
//
// Note: Return type does not take part in template argument deduction.
// Usage: make_from_any_vector<TargetType>(source_container, dim).

template<class LP, class SomeContainer, class DimType, TEMPL_RESTRICT_DECL2(
         IsALatticePoint<LP>, IsRepLinear_RW<LP>)>
LP make_from_any_vector(SomeContainer const &container, DimType dim)
{
  static_assert(DoesDeclareCoordinateType<LP>::value, "Not declaring coordinate types");
  using ET = Get_CoordinateType<LP>;
  DEBUG_TRACEGENERIC("generically converting vector to LP for" << LP::class_name() )
  LP result(dim);
//  auto dim = result.get_dim();
  for(uint_fast16_t i =0; i<dim; ++i)
  {
    result[i] = static_cast<ET>( container[i] );
  }
  result.sanitize();
  return result;
}

// Same as above, but un-Z_NR's the container.

template<class LP, class SomeZNRContainer, class DimType, TEMPL_RESTRICT_DECL2(
         IsALatticePoint<LP>, IsRepLinear_RW<LP> )>
LP make_from_znr_vector(SomeZNRContainer const &container, DimType dim)
{
  static_assert(DoesDeclareCoordinateType<LP>::value, "Not declaring coordinate types");
  using ET = Get_CoordinateType<LP>;
  DEBUG_TRACEGENERIC("generically converting vector to LP and un-ZNRing for" << LP::class_name() )
  LP result(dim);
  for(uint_fast16_t i =0; i<dim; ++i)
  {
    result[i] = ConvertMaybeMPZ<ET>::convert_to_inttype( container[i].get_data() );
  }
  result.sanitize();
  return result;
}

/*

template<class LP>
std::istream & operator>> (std::istream & is, typename std::enable_if<IsALatticePoint<LP>::value, LP>::type &lp)
{
    lp.read_from_stream(is, IgnoreAnyArg{});
    return is;
}

template<class LP>
std::ostream & operator<< (std::ostream & os, typename std::enable_if<IsALatticePoint<LP>::value,LP>::type &lp )
{
    lp.write_to_stream(os,IgnoreAnyArg{});
    return os;
}

*/

} // end namespace

#include "LatticePointGeneric.h"


// cleaning up internal macros.

#undef IMPL_IS_LATP

#undef CREALTHIS
#undef REALTHIS

#endif

//clang-format on
