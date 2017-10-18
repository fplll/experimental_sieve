#ifndef LATTICE_POINT_CONCEPT_H
#define LATTICE_POINT_CONCEPT_H

// TODO: Remove the concept idea below, we replace it by CRTP (i.e. inheritance without virtual
// functions, using templates)

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include <cmath>
#include <cstdint>
#include <gmpxx.h>
#include "GlobalStaticData.h"

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
  using Invalid=std::true_type;
  LatticePointTraits(...) = delete;
  // void Invalid
  // static_assert(false) is invalid due to subtleties of C++, even if it may work on some compilers
};

/**
  Available traits
  (everything that needs to be set to true_type defaults to false, unless specified otherwise)
  These traits are used to selectively enable some meaningful default operations on lattice points.

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

  NOTE: ExposesCoos / CoordinateType is a purely syntactic trait.
        It does not tell anything about the semantics of [].

  AbsoluteCooType : return type of get_absolute_coo
                    defaults to CoordinateType

  RepCooType      : return type of get_internal_rep
                    defaults to CoordinateType

  ExposesInternalRep :  Indicates that get_internal_rep(i), get_internal_rep_size() exists and may be
                        read for 0<= i < get_internal_rep_size().
                        We are guaranteed that these entries determine the point.
                        implied by InternalRepVector_R, InternalRepVector_RW, InternalRepByCoos, InternalRepIsAbsolute

  InternalRepVector:    Indicates that the internal representation is linear.
                        implies HasInternalRep
  InternalRep_RW:       Indicates that the internal representation may be written to.
                        After calling sanitize(), the class will be in a valid state.
                        implies HasInternalRep

  InternalRepByCoos:  Indicates that the internal representation is accessible by operator[]
                      implies HasInternalRep, ExposesCoos

  InternalRepIsAbsolute:  Indicates that the internal representation gives absolute coordinates.
                          implies HasInternalRep

  CheapNorm2 : Set to true_type to indicate that get_norm2() is cheap.
               (typically, it's precomputed and stored with the point)

  CheapNegate: Set to true_type to indicate that negation needs no sanitize().

  HasApproximations: Set to true_type to indicate that the point has approximations.

DEPRECATED:

  CoordinateAccess : Set to true_type to indicate that the class exposes an operator[].
                     We may read LatticePoint[i] for 0 <= i < get_internal_rep_size().
                     These entries behave like coordinates.

  CoordinateVector : Set to true_type to indicate that we have the same guarantees as
                     CoordinateAccess (which still should be set) and that additionally, we may
                     also write to the coordinates of LatticePoint[i].
                     After calling sanitize(), the class will be in a valid state.

  AbsoluteCoos : Set to true_type to indicate that the class exposes an operator[] and that
                 the entries are absolute coordinates (as opposed to e.g. coordiates wrt. some basis
                 which is part of the class-wide data).
                 We assume that CoordinateAccess is set as well.


*/

/**
  Example:
  template<>
  struct / class LatticePointTraits<MyCoolLatticePointClass>
  {
    public: // automatic for structs
    using ScalarProductStorageType = long (say)
    using CoordinateVector = std::true_type
    using CoordinateType = int
    ...
  }
*/


template<class Implementation> class GeneralLatticePoint;



/**
  Trait getters
*/
// Usage IsNegateCheap<PlainLatticePoint>::value
// (as opposed to LatticePointTraits<PlainLatticePoint>::CheapNegate::value )

// The reason is that it is easier to read and that
// IsNegateCheap<...> defaults to false, whereas
// LatticePointTraits<...>::... defaults to a compile-time error.

// Trait getter names prefixed with T_ are internal and only serve to construct other getters
// via the template<class LatP> using ... typedefs below.


// These check for the existence of the trait in the Trait class:
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(LatticePointTag, std::true_type, IsALatticePoint);

CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, ScalarProductStorageType, DeclaresScalarProductStorageType);
CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, CoordinateType, DoesDeclareCoordinateType);

CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Invalid, std::true_type, HasNoLPTraits);

CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, ExposesCoos, std::true_type, T_ExposesCoos);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, ExposesInternalRep, std::true_type, T_InternalRep);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, InternalRepVector, std::true_type, T_RepVector_R);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, InternalRep_RW, std::true_type, T_Rep_RW);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, InternalRepByCoos, std::true_type, T_InternalRepByCoos);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, InternalRepIsAbsolute, std::true_type, T_AbsoluteCoos);

CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CoordinateVector, std::true_type, IsCooVector);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CoordinateAccess, std::true_type, HasCoos);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, AbsoluteCoos, std::true_type, CoosAreAbsolute);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CheapNorm2, std::true_type, T_Norm2Cheap);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CheapNegate, std::true_type, T_NegateCheap);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, HasApproximations, std::true_type, T_Approximation);

// Retrieving Traits with defaults:
MAKE_TRAIT_GETTER(LatticePointTraits, CoordinateType, void, GetCooType);
MAKE_TRAIT_GETTER(LatticePointTraits, ScalarProductStorageType, void, GetScalarProductStorageType);
// ClassToCheck is the argument of the constructed Traits getter inside the macro def. This makes it
// default to GetScalarproductStorageType
MAKE_TRAIT_GETTER(LatticePointTraits, ScalarProductStorageType_Full,
  typename GetScalarProductStorageType<ClassToCheck>::type, GetScalarProductStorageType_Full);
MAKE_TRAIT_GETTER(LatticePointTraits, AbsoluteCooType,
  typename GetCooType<ClassToCheck>::type, GetAbsoluteCooType);
MAKE_TRAIT_GETTER(LatticePointTraits, RepCooType,
  typename GetCooType<ClassToCheck>::type, GetRepCooType);


// These are what the rest of the code should be actually using:
template<class LatP> using DoesExposeCoos = std::integral_constant<bool,
  T_ExposesCoos<LatP>::value || T_InternalRepByCoos<LatP>::value ||
  (!std::is_same<typename GetCooType<LatP>::type, void>::value) >;

template<class LatP> using HasInternalRep = std::integral_constant<bool,
  T_InternalRep<LatP>::value || T_RepVector_R<LatP>::value || T_Rep_RW<LatP>::value ||
  T_InternalRepByCoos<LatP>::value || T_InternalRepByCoos<LatP>::value >;
template<class LatP> using IsRepLinear = std::integral_constant<bool,
  T_RepVector_R<LatP>::value>;
template<class LatP> using IsRepRW = std::integral_constant<bool,
  T_Rep_RW<LatP>::value>;
template<class LatP> using HasRepByCoos = std::integral_constant<bool,
  T_InternalRepByCoos<LatP>::value>;
template<class LatP> using HasAbsoluteRep = std::integral_constant<bool,
  T_AbsoluteCoos<LatP>::value>;
template<class LatP> using IsNorm2Cheap = std::integral_constant<bool,
  T_Norm2Cheap<LatP>::value>;
template<class LatP> using IsNegateCheap = std::integral_constant<bool,
  T_NegateCheap<LatP>::value>;
template<class LatP> using IsRepLinear_RW = std::integral_constant<bool,
  IsRepLinear<LatP>::value && IsRepRW<LatP>::value>;
template<class LatP> using HasApproximations = std::integral_constant<bool,
  T_Approximation<LatP>::value>;


#define IMPL_IS_LATP \
static_assert(std::is_same<Impl,LatP>::value,"Using template member function with wrong type")

// unsure whether to use reinterpret_casts here...

#define CREALTHIS static_cast<LatP const*>(this)
#define REALTHIS  static_cast<LatP*>(this)

template<class LatP>
class GeneralLatticePoint
{
    static_assert(!HasNoLPTraits<LatP>::value, "Trait class not specialized.");
    static_assert(DeclaresScalarProductStorageType<LatP>::value,
                  "Lattice Point class does not typedef its scalar product type");
    static_assert(!DoesExposeCoos<LatP>::value || DoesDeclareCoordinateType<LatP>::value,
                  "Lattice Point exposes coos, but does not tell its type.");

    friend LatP; // makes children able to access private (in addition to protected) members.
                 // since the constructor is private, this enforces correct usage.
                 // (Note that it may prevent multi-level inheritance)
    public:
    using ScalarProductStorageType = typename GetScalarProductStorageType<LatP>::type;
    using ScalarProductStorageType_Full = typename GetScalarProductStorageType_Full<LatP>::type;
    using CooType = typename GetCooType<LatP>::type; //may be void
    using AbsCooType = typename GetAbsoluteCooType<LatP>::type;
    using RepCooType = typename GetRepCooType<LatP>::type;

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
// Otherwise, the SFINAE magic behind it won't work. Of course, we never use these templates with
// Impl!=LatP and we may even static_assert(Impl==LatP) inside the implementation.

    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(
      IsALatticePoint<LatP2>, IsRepLinear_RW<Impl>, IsRepLinear<LatP2>)>
    inline LatP& operator+=(LatP2 const &x2);

    // Fun fact: TEMPL_RESTRICT_DECL2 misinterprets the second comma as a macro argument separator
    // so it has 3 arguments, the third one being garbage like " IsRepLinear<LatP2> > > "
    // However, due to way TEMPL_RESTRICT_DECL2 is defined, this actually still works!

    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(
      IsALatticePoint<LatP2>, MyNAND<IsRepLinear_RW<Impl> , IsRepLinear<LatP2> > )>
    LatP& operator+=(LatP2 const &x2) = delete;

    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(
      IsALatticePoint<LatP2>, IsRepLinear_RW<Impl>, IsRepLinear<LatP2>)>
    inline LatP& operator-=(LatP2 const &x2);
    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(
      IsALatticePoint<LatP2>, MyNAND<IsRepLinear_RW<Impl>, IsRepLinear<LatP2> >)>
    LatP& operator-=(LatP2 const &x2) = delete;

    template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<typename std::decay<LatP2>::type>)>
    inline bool operator!=(LatP2 && x2) const {return !(CREALTHIS->operator==(std::forward<LatP2>(x2)));};

    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(
      IsALatticePoint<LatP2>, HasInternalRep<Impl>, HasInternalRep<LatP2>)>
    inline bool operator==(LatP2 const &x2) const;
    template<class LatP2, class Impl=LatP, TEMPL_RESTRICT_DECL2(
      IsALatticePoint<LatP2>, MyNAND<HasInternalRep<Impl>, HasInternalRep<LatP2> >)>
    bool operator==(LatP2 const &x2) = delete;

    template<class Integer, class Impl=LatP, TEMPL_RESTRICT_DECL2(
      IsRepLinear_RW<Impl>, std::is_integral<Integer>)>
    inline LatP& operator*=(Integer const multiplier);
    template<class Integer, class Impl=LatP, TEMPL_RESTRICT_DECL2(
      std::is_integral<Integer>, MyNegation<IsRepLinear_RW<Impl> >)>
    LatP& operator*=(Integer const multiplier) = delete;
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
    template<class Impl=LatP,TEMPL_RESTRICT_DECL2(HasInternalRep<Impl>)>
    inline auto get_internal_rep_size() const -> decltype( std::declval<Impl>().get_dim() );
    template<class Arg, class Impl=LatP, TEMPL_RESTRICT_DECL2(HasRepByCoos<Impl>)>
    inline CooType const & get_internal_rep(Arg &&arg) const
    {
      return CREALTHIS->operator[](std::forward<Arg>(arg));
    }
    template<class Arg, class Impl=LatP, TEMPL_RESTRICT_DECL2(HasRepByCoos<Impl>, IsRepLinear_RW<Impl>)>
    inline CooType & get_internal_rep(Arg &&arg)
    {
      return REALTHIS->operator[](std::forward<Arg>(arg));
    }

    // decltype(auto) from C++14 would be a godsend...
    template<class Arg, class Impl=LatP, TEMPL_RESTRICT_DECL2(HasAbsoluteRep<Impl>)>
    inline auto get_absolute_coo(Arg &&arg) const
    -> typename std::remove_reference< decltype (std::declval<Impl>().get_internal_rep(std::forward<Arg>(arg)) )>::type
    {
      return CREALTHIS->get_internal_rep(std::forward<Arg>(arg));
    }
    template<class Arg, class Impl=LatP, TEMPL_RESTRICT_DECL2(MyNegation<HasAbsoluteRep<Impl> >)>
    int get_absolute_coo(Arg &&arg) const = delete; // needs to be overloaded

    void get_dim() const = delete; // Note that the overload shall NOT have void return type.
                                   // It's just not possible to specific it here w/o C++14 auto.

/**
  Used for output to stream. Note that operator<< calls this (or an overloaded version)
  Supposed to be overloaded for the specific class.
*/

//    MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE // This may be too strict.
    inline std::ostream& write_lp_to_stream(std::ostream &os, bool const include_norm2=true) const;

    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(HasInternalRep<Impl>)>
    inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const;

    std::istream& read_from_stream(std::istream &is) = delete;

/**
  Fills a lattice point with zeros.
  Note that freshly constructed lattice points may contain uninitialized values unless this function is used.

  The latter depends subtly on the constructors used (empty vs. default constructor...)
  May be overloaded by Derived class.
*/
    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(IsRepLinear_RW<Impl>)>
    inline void fill_with_zero();
    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(MyNegation<IsRepLinear_RW<Impl>>)>
    inline void fill_with_zero() = delete;
/**
  Changes vector from v to -v. May be overloaded.
*/
    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(IsRepLinear_RW<Impl>)>
    inline void make_negative();
    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(MyNegation<IsRepLinear_RW<Impl>>)>
    inline void make_negative() = delete;

/**
  Tests whether a lattice point is all-zero.
  May be overloaded.
*/
    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(IsRepLinear<Impl>)>
    inline bool is_zero() const;
    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(MyNegation<IsRepLinear<Impl>>)>
    inline bool is_zero() const = delete;

/**
  Makes an (explicit) copy of the current point.

  Default: Copying components as a vector and sanitize additional data.
  (The latter probably means needless recomputing of data that the original already had)
  With the intended semantics of [] and entries, this results in a deep copy.

  May be overloaded.
*/

    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(IsRepLinear_RW<Impl>)>
    inline LatP make_copy() const;
    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(MyNegation<IsRepLinear_RW<Impl> >)>
    LatP make_copy() const = delete;

/**
     brings the lattice point into a defined state. Defaults to "do nothing".

     Note that this function must *not* be const, since it actually changes observable behaviour:
     For efficiency reasons, sanitize is not called upon write access to LatP[]; rather,
     sanitize is called from outside the class.
     As such, the user can leave LatP in an invalid state and sanitize remedies that.

     The second version takes norm2 as an argument (to avoid recomputing it).
*/

// TODO: Arg type overloads
    void sanitize() {};
    void sanitize(ScalarProductStorageType const &norm2) { sanitize(); };


/**
     obtains the squared length.
     By default, we use the scalar product to compute it.
     Usually norm2 will be precomputed, so this function should be
     overridden by LatP.
*/

// TODO: This looks wrong regarding types.

    inline ScalarProductStorageType get_norm2() const;

/**   This function returns the exact norm, ignoring any approximations. */
    inline ScalarProductStorageType get_norm2_exact() const {return CREALTHIS->get_norm2(); }
    inline ScalarProductStorageType_Full get_norm2_full() const { return CREALTHIS->get_norm2(); }

    // don't call directly. We use compute_sc_product(x1,x2) for a more symmetric syntax.
    // However, out-of-class definitions get messy with overloading.

// TODO: friends

    inline ScalarProductStorageType do_compute_sc_product(LatP const &x2) const;

    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(MyNegation<HasApproximations<Impl> >)>
    inline ScalarProductStorageType do_compute_sc_product_exact(LatP const &x2) const
    {
      return CREALTHIS->do_compute_sc_product(x2);
    }
    template<class Impl=LatP, TEMPL_RESTRICT_DECL2(MyNegation<HasApproximations<Impl> >)>
    inline ScalarProductStorageType_Full do_compute_sc_product_full(LatP const &x2) const
    {
      return CREALTHIS->do_compute_sc_product(x2);
    }
 };

 /**
  Non-member functions
  */

template<class LP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LP>)>
inline typename LP::ScalarProductStorageType compute_sc_product(LP const &lp1, LP const &lp2)
{ return lp1.do_compute_sc_product(lp2); }

template<class LP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LP>)>
inline typename LP::ScalarProductStorageType compute_sc_product_exact(LP const &lp1, LP const &lp2)
{ return lp1.do_compute_sc_product_exact(lp2); }

template<class LP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LP>)>
inline typename LP::ScalarProductStorageType_Full compute_sc_product_full(LP const &lp1, LP const &lp2)
{ return lp1.do_compute_sc_product_full(lp2); }


#define FOR_LATTICE_POINT_LP \
template<class LP, typename std::enable_if<IsALatticePoint<LP>::value, int>::type=0>

#define FOR_LATTICE_POINTS_LP1_LP2 \
template<class LP1, class LP2, typename std::enable_if< \
         IsALatticePoint<LP1>::value && IsALatticePoint<LP2>::value,int>::type=0>


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
  using ET = typename GetCooType<LP>::type;
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
  using ET = typename GetCooType<LP>::type;
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
template<class Implementation>
bool GeneralLatticePoint<Implementation>::operator> (Implementation const & other) const
{
    return Implementation::get_norm2(*this) > other.get_norm2();
}
 */

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
#undef FOR_LATTICE_POINT_LP
#undef FOR_LATTICE_POINTS_LP1_LP2
#undef CREALTHIS
#undef REALTHIS

#endif

//clang-format on
