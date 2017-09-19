#ifndef LATTICE_POINT_CONCEPT_H
#define LATTICE_POINT_CONCEPT_H

// TODO: Remove the concept idea below, we replace it by CRTP (i.e. inheritance without virtual
// functions, using templates)

#include "DebugAll.h"
#include "SieveUtility.h"
#include "assert.h"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
#include <utility>
#include <gmpxx.h>

// clang-format off

namespace GaussSieve{

/**
  TODO: Write explanation, we use CRTP now.
*/

// This class template stores the typedefs that the individual lattice point classes have
// There has to be a specialization for each lattice point class.
// The general template must never be instantiated.
// Note that we can not put these traits into the lattice points classes directly, because that
// would cause circular references due to CRTP.

// Cf. ExacLatticePoint.h for an example of how a concrete specialization is supposed to look like.

template<class LatticePoint> struct LatticePointTraits
{
  public:
  using Invalid=std::true_type;
  // void Invalid
  // static_assert(false) is invalid due to subtleties of C++, even if it may work on some compilers
};

/**
  Available traits
  (everything that needs to be set to true_type defaults to false, unless specified otherwise)
  These traits are used to selectively enable some meaningful default operations on lattice points.

  AuxDataType : type of class-wide data that is required to interpret a given point
                (e.g. Dimension, custom memory allocator)
                Default: IgnoreAnyArg
                Needs to be initializable from an int (??? DimensionType)
  ScalarProductStorageType: A type that can hold the result of a scalar product computation. Mandatory.
                           Note that the result from a scalar product computation might actually differ.
                           (due to delayed evaluation)
  ScalarProductStorageType_Full:  A type that can hold the result of a scalar product computation.
                                  In case the Point has approximations, will contain an approximate
                                  scalar product as well.
                                  Default: Same as ScalarProductStorageType

  CoordinateType : return type of operator[] if available.

  CoordinateAccess : Set to true_type to indicate that the class exposes an operator[].
                     We may read LatticePoint[i] for 0 <= i < get_vec_size().
                     These entries behave like coordinates.

  CoordinateVector : Set to true_type to indicate that we have the same guarantees as
                     CoordinateAccess (which still should be set) and that additionally, we may
                     also write to the coordinates of LatticePoint[i].
                     After calling sanitize(), the class will be in a valid state.

  AbsoluteCoos : Set to true_type to indicate that the class exposes an operator[] and that
                 the entries are absolute coordinates (as opposed to e.g. coordiates wrt. some basis
                 which is part of the class-wide data).
                 We assume that CoordinateAccess is set as well.

  CheapNorm2 : Set to true_type to indicate that get_norm2() is cheap.
               (typically, it's precomputed and stored with the point)

  CheapNegate: Set to true_type to indicate that negation needs no sanitize().

  HasApproximations: Set to true_type to indicate that the point has approximations.
*/

/**
  Example:
  template<>
  struct / class LatticePointTraits<MyCoolLatticePointClass>
  {
    public: // automatic for structs
    using AuxDataType = Whatever
    using ScalarProductStorageType = long (say)
    using CoordinateVector = std::true_type
    using CoordingateType = int
    ...
  }
*/


template<class Implementation> class GeneralLatticePoint;
template<class T> class StaticInitializer;


/**
  Trait getters
*/
// Usage IsNegateCheap<PlainLatticePoint>::value
// (as opposed to LatticePointTraits<PlainLatticePoint>::CheapNegate::value )

// The reason is that it is easier to read and that
// IsNegateCheap<...> defaults to false, whereas
// LatticePointTraits<...>::... defaults to a compile-time error.


//CREATE_MEMBER_TYPEDEF_CHECK_CLASS(ScalarProductStorageType, DeclaresScalarProductReturnType);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(LatticePointTag, std::true_type, IsALatticePoint);
CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, ScalarProductStorageType, HasScalarProductReturnType);
MAKE_TRAIT_GETTER(LatticePointTraits, AuxDataType, IgnoreAnyArg, GetAuxDataType);
MAKE_TRAIT_GETTER(LatticePointTraits, ScalarProductStorageType, void, GetScalarProductStorageType);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Invalid, std::true_type, HasNoLPTraits);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CoordinateVector, std::true_type, IsCooVector);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CoordinateAccess, std::true_type, HasCoos);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, AbsoluteCoos, std::true_type, CoosAreAbsolute);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CheapNorm2, std::true_type, IsNorm2Cheap);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CheapNegate, std::true_type, IsNegateCheap);

CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, CoordinateType, DoesDeclareCoordinateType);

MAKE_TRAIT_GETTER(LatticePointTraits, CoordinateType, void, GetCooType);

// ClassToCheck is the argument of the constructed Traits getter inside the macro def.
MAKE_TRAIT_GETTER(LatticePointTraits, ScalarProductStorageType_Full,
  typename GetScalarProductStorageType<ClassToCheck>::type, GetScalarProductStorageType_Full);


#define MEMBER_ONLY_EXISTS_IF_COO_READ \
template<class Impl=LatP, typename std::enable_if<HasCoos<Impl>::value,int>::type = 0>

#define MEMBER_ONLY_EXISTS_IF_COO_READ_IMPL \
template<class Impl, typename std::enable_if<HasCoos<Impl>::value,int>::type>


#define MEMBER_ONLY_EXISTS_IF_COO_WRITE \
template<class Impl=LatP, typename std::enable_if<IsCooVector<Impl>::value,int>::type = 0>

#define MEMBER_ONLY_EXISTS_IF_COO_WRITE_IMPL \
template<class Impl, typename std::enable_if<IsCooVector<Impl>::value,int>::type>


#define MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE \
template<class Impl=LatP, typename std::enable_if<CoosAreAbsolute<Impl>::value,int>::type = 0>

#define MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE_IMPL \
template<class Impl, typename std::enable_if<CoosAreAbsolute<Impl>::value,int>::type>


#define IMPL_IS_LATP \
static_assert(std::is_same<Impl,LatP>::value,"Using template member function with wrong type")

// unsure whether to use reinterpret_casts here...

#define CREALTHIS static_cast<LatP const*>(this)
#define REALTHIS  static_cast<LatP*>(this)

template<class LatP>
class GeneralLatticePoint
{
    static_assert(!HasNoLPTraits<LatP>::value, "Trait class not specialized.");
    static_assert(HasScalarProductReturnType<LatP>::value,
                  "Lattice Point class does not typedef its scalar product type");
    public:
    friend LatP; // makes children able to access private (in addition to protected) members.
                 // since the constructor is private, this enforces correct usage.
                 // (Note that it may prevent multi-level inheritance)

    using AuxDataType = typename GetAuxDataType<LatP>::type;
    using ScalarProductStorageType = typename GetScalarProductStorageType<LatP>::type;
    using ScalarProductStorageType_Full = typename GetScalarProductStorageType_Full<LatP>::type;

    private:
    explicit constexpr GeneralLatticePoint()=default; //only callable from its friends
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

// You should not overload <, >, <=, >=, !=
// You may overload +=, -=, *=, unary-, ==

    inline bool operator< ( LatP const &rhs ) const;
    inline bool operator> ( LatP const &rhs ) const;
    inline bool operator<=( LatP const &rhs ) const;
    inline bool operator>=( LatP const &rhs ) const;
// arithmetic:
    template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value && IsCooVector<LatP>::value && HasCoos<LatP2>::value)>
    inline LatP& operator+=(LatP2 const &x2);

    template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value && IsCooVector<LatP>::value && HasCoos<LatP2>::value)>
    inline LatP& operator-=(LatP2 const &x2);

    template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<typename std::decay<LatP2>::type >::value)>
    inline bool operator!=(LatP2 && x2) const {return !(CREALTHIS->operator==(std::forward<LatP2>(x2)));};

    template<class LatP2, TEMPL_RESTRICT_DECL(IsALatticePoint<LatP2>::value && HasCoos<LatP>::value && HasCoos<LatP2>::value)>
    inline bool operator==(LatP2 const &x2) const;

    template<class Integer,
      TEMPL_RESTRICT_DECL(IsCooVector<LatP>::value && std::is_integral<Integer>::value)>
    inline LatP& operator*=(Integer const multiplier);
    inline LatP& operator*=(mpz_class const &multiplier) = delete; // not implemented yet

    inline LatP operator-() &&; //unary-
    // operator+ is entirely defined in terms of +=
    // operator- is entirely defined in terms of -=
    // operator* is entirely defined in terms of *=
    // Definitions are out-of-class and not supposed to be overloaded.


    // get_vec_size() returns the number of coordinates we get when using latp[i]
    // get_dim returns the (ambient) dimension the vector is supposed to represent.
    // By default, vec_size is the same as dim.

    // get_dim must be overloaded.
    // get_vec_size may be overloaded.
    MEMBER_ONLY_EXISTS_IF_COO_READ
    inline auto get_vec_size() const -> decltype( std::declval<Impl>().get_dim() );
    void get_dim() const = delete; // Note that the overload shall NOT have void return type.
                                   // It's just not possible to specific it here w/o C++14 auto.

/**
  Used for output to stream. Note that operator<< calls this (or an overloaded version)
  Supposed to be overloaded for the specific class.
*/

    MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE // This may be too strict.
    inline std::ostream& write_to_stream(std::ostream &os, bool const include_norm2=true) const;

    std::istream& read_from_stream(std::istream &is) = delete;

/**
  Fills a lattice point with zeros.
  Note that freshly constructed lattice points may contain uninitialized values unless this function is used.

  The latter depends subtly on the constructors used (empty vs. default constructor...)
  May be overloaded by Derived class.
*/

    MEMBER_ONLY_EXISTS_IF_COO_WRITE inline void fill_with_zero();

/**
  Changes vector from v to -v. May be overloaded.
*/
    MEMBER_ONLY_EXISTS_IF_COO_WRITE inline void make_negative();

/**
  Tests whether a lattice point is all-zero.
  May be overloaded.
*/
    MEMBER_ONLY_EXISTS_IF_COO_READ inline bool is_zero() const;

/**
  Makes an (explicit) copy of the current point.

  Default: Copying components as a vector and sanitize additional data.
  (The latter probably means needless recomputing of data that the original already had)
  With the intended semantics of [] and entries, this results in a deep copy.

  May be overloaded.
*/

    MEMBER_ONLY_EXISTS_IF_COO_WRITE inline LatP make_copy() const;

/**
     brings the lattice point into a defined state. Defaults to "do nothing".

     Note that this function must *not* be const, since it actually changes observable behaviour:
     For efficiency reasons, sanitize is not called upon write access to LatP[]; rather,
     sanitize is called from outside the class.
     As such, the user can leave LatP in an invalid state and sanitize remedies that.

     The second version takes norm2 as an argument (to avoid recomputing it).
*/
    void sanitize() {};
    void sanitize(ScalarProductStorageType const &norm2) { sanitize(); };


/**
     obtains the squared length.
     By default, we use the scalar product to compute it.
     Usually norm2 will be precomputed, so this function should be
     overridden by LatP.
*/

    inline ScalarProductStorageType get_norm2() const;

/**   This function returns the exact norm, ignoring any approximations. */
    inline ScalarProductStorageType get_norm2_exact() const {return CREALTHIS->get_norm2(); }
    inline ScalarProductStorageType_Full get_norm2_full() const { return CREALTHIS->get_norm2(); }


    // don't call directly. We use compute_sc_product(x1,x2) for a more symmetric syntax.
    // However, out-of-class definition get messy with overloading.
    MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE
    inline ScalarProductStorageType do_compute_sc_product(LatP const &x2) const;

    inline ScalarProductStorageType do_compute_sc_product_exact(LatP const &x2) const
    {
      return CREALTHIS->do_compute_sc_product(x2);
    }
    inline ScalarProductStorageType_Full do_compute_sc_product_full(LatP const &x2) const
    {
      return CREALTHIS->do_compute_sc_product(x2);
    }
 };

 /**
  Non-member functions
  */

template<class LP, TEMPL_RESTRICT_DECL(IsALatticePoint<LP>::value)>
inline typename LP::ScalarProductStorageType compute_sc_product(LP const &lp1, LP const &lp2);

template<class LP, TEMPL_RESTRICT_DECL(IsALatticePoint<LP>::value)>
inline typename LP::ScalarProductStorageType compute_sc_product_exact(LP const &lp1, LP const &lp2);

template<class LP, TEMPL_RESTRICT_DECL(IsALatticePoint<LP>::value)>
inline typename LP::ScalarProductStorageType_Full compute_sc_product_full(LP const &lp1, LP const &lp2);


// Initializer for static data.
// This is the default initializer, which does nothing.

template<class T>
class DefaultStaticInitializer
{
  public:
#ifndef DEBUG_SIEVE_LP_INIT
  explicit constexpr DefaultStaticInitializer() = default;
#else
  static unsigned int init_count; // counts the number of objects of this type that exist, essentially.
  static bool is_initialized(){ return init_count > 0; }; // Does an object exist?
  explicit DefaultStaticInitializer(){++init_count;};
  ~DefaultStaticInitializer()
  {
    assert(init_count>0);
    --init_count;
  }
#endif
};

#ifdef DEBUG_SIEVE_LP_INIT
  template<class T> unsigned int DefaultStaticInitializer<T>::init_count = 0;
#endif


// Non-member functions:

#define FOR_LATTICE_POINT_LP \
template<class LP, typename std::enable_if<IsALatticePoint<LP>::value, int>::type=0>

#define FOR_LATTICE_POINTS_LP1_LP2 \
template<class LP1, class LP2, typename std::enable_if< \
         IsALatticePoint<LP1>::value && IsALatticePoint<LP2>::value,int>::type=0>

// TODO: Assert Traits!

// default comparison is coordinate-wise.
// Note: We allow comparing different types of lattice points.
// This may or not work.



// dispatch to sub/subval function


FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 const &x1, LP2 const &x2){ return sub(x1,x2); }

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 && x1, LP2 const &x2)
{
  LP1 tmp = std::move(x1);
  tmp-=x2;
  return tmp;
}

FOR_LATTICE_POINT_LP
LP operator-(LP const &x1, LP &&x2)
{
  LP tmp = std::move(x2);
  tmp.make_negative();
  tmp+=x1;
  return tmp; //addval(tmp,x1);
}

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 && x1, LP2 && x2)
{
  LP1 tmp = std::move(x1);
  tmp-=std::move(x2);
  return tmp;
}

// dispatch to possible member function

FOR_LATTICE_POINT_LP
std::istream& operator>>(std::istream &is, LP &LatP)
{
  return LatP.read_from_stream(is);
}

FOR_LATTICE_POINT_LP
std::ostream& operator<<(std::ostream &os, LP const &LatP)
{
  return LatP.write_to_stream(os);
}




// implementations of generic functions for lattice points:

/*
template<class LP, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type = 0>
LP add(LP const &x1, LP const &x2)
{
  DEBUG_TRACEGENERIC( "generically adding" << x1.class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = x1.get_vec_size();
  auto const dim2 = x2.get_vec_size();
  assert( dim1 == dim2 );
  auto const real_dim1 = x1.get_dim();
  auto const real_dim2 = x2.get_dim();
  assert(real_dim1 == real_dim2);
  #endif
  auto const dim = x1.get_vec_size();
  auto const real_dim = x1.get_dim();
  LP NewLP(real_dim);
  for(uint_fast16_t i = 0; i < dim; ++i )
  {
    NewLP[i] = x1[i] +x2[i];
  }
  NewLP.sanitize();
  return NewLP;
}
*/


template<class LP, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type = 0>
LP sub(LP const &x1, LP const &x2)
{
  DEBUG_TRACEGENERIC("generically subtracting" << x1.class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = x1.get_vec_size();
  auto const dim2 = x2.get_vec_size();
  assert(dim1==dim2);
  auto const real_dim1 = x1.get_dim();
  auto const real_dim2 = x2.get_dim();
  assert(real_dim1 == real_dim2);
  #endif
  auto const dim = x1.get_vec_size();
  auto const real_dim = x1.get_dim();
  LP NewLP(real_dim);
  for(uint_fast16_t i=0; i < dim ; ++i )
  {
    NewLP[i] = x1[i] - x2[i];
  }
  NewLP.sanitize();
  return NewLP;
}


// this function can be used to initialize an LP with container types that allow []-access.
// Note that there is an explicit static_cast to LP's entry types.
// In particular, this can convert mpz_t to mpz_class...
//
// Note: Return type does not take part in template argument deduction.
// Usage: make_from_any_vector<TargetType>(source_container, dim).

template<class LP, class SomeContainer, class DimType, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type=0>
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

template<class LP, class SomeZNRContainer, class DimType, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type=0>
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
    static_assert(std::is_same< typename LP::AuxDataType, IgnoreAnyArg>::value == true, "This Lattice Point class requires auxiliary data for input");
    lp.read_from_stream(is, IgnoreAnyArg{});
    return is;
}

template<class LP>
std::ostream & operator<< (std::ostream & os, typename std::enable_if<IsALatticePoint<LP>::value,LP>::type &lp )
{
    static_assert(std::is_same< typename LP::AuxDataType, IgnoreAnyArg>::value == true, "This Lattice Point class requires auxiliary data for output");
    lp.write_to_stream(os,IgnoreAnyArg{});
    return os;
}

*/

} // end namespace

#include "LatticePointGeneric.h"


// cleaning up internal macros.
#undef MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE
#undef MEMBER_ONLY_EXISTS_IF_COO_READ
#undef MEMBER_ONLY_EXISTS_IF_COO_WRITE

#undef MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE_IMPL
#undef MEMBER_ONLY_EXISTS_IF_COO_READ_IMPL
#undef MEMBER_ONLY_EXISTS_IF_COO_WRITE_IMPL

#undef IMPL_IS_LATP
#undef FOR_LATTICE_POINT_LP
#undef FOR_LATTICE_POINTS_LP1_LP2
#undef CREALTHIS
#undef REALTHIS

#endif

//clang-format on
