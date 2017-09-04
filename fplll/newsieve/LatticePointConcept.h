
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
// would case circular references due to CRTP.
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
  ScalarProductReturnType: Return type of scalar products. Mandatory!

  CoordinateType : return type of operator[] if available.

  CoordinateAccess : Set to true_type to indicate that the class exposes an operator[].
                     We may read LatticePoint[i] for 0<= i < get_vec_size().
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

*/

/**
  Example:
  template<>
  struct / class LatticePointTraits<MyCoolLatticePointClass>
  {
    public: // automatic for structs
    using AuxDataType = Whatever
    using ScalarProductReturnType = long (say)
    using CoordinateVector = std::true_type
    using CoordingateType = int
    ...
  }
*/


//template<class Implementation> class ImplementationTraits;
template<class Implementation> class GeneralLatticePoint;

//CREATE_MEMBER_TYPEDEF_CHECK_CLASS(ScalarProductReturnType, DeclaresScalarProductReturnType);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(LatticePointTag, std::true_type, IsALatticePoint);
CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, ScalarProductReturnType, HasScalarProductReturnType);
MAKE_TRAIT_GETTER(LatticePointTraits, AuxDataType, IgnoreAnyArg, GetAuxDataType);
MAKE_TRAIT_GETTER(LatticePointTraits, ScalarProductReturnType, void, GetScPType);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, Invalid, std::true_type, HasNoLPTraits);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CoordinateVector, std::true_type, IsCooVector);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CoordinateAccess, std::true_type, HasCoos);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, AbsoluteCoos, std::true_type, CoosAreAbsolute);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CheapNorm2, std::true_type, IsNorm2Cheap);
CREATE_TRAIT_EQUALS_CHECK(LatticePointTraits, CheapNegate, std::true_type, IsNegateCheap);

CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, CoordinateType, DoesDeclareCoordinateType);

MAKE_TRAIT_GETTER(LatticePointTraits, CoordinateType, void, GetCooType);

// Usage IsNegateCheap<PlainLatticePoint>::value
// (as opposed to LatticePointTraits<PlainLatticePoint>::CheapNegate::value )
// also, IsNegateCheap<...> defaults to false
// whereas LatticePointTraits<...>::... defaults to a compile-time error.


#define MEMBER_ONLY_EXISTS_IF_COO_READ \
template<class Impl=LatP, typename std::enable_if<HasCoos<Impl>::value,int>::type = 0>

#define MEMBER_ONLY_EXISTS_IF_COO_WRITE \
template<class Impl=LatP, typename std::enable_if<IsCooVector<Impl>::value,int>::type = 0>

#define MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE \
template<class Impl=LatP, typename std::enable_if<CoosAreAbsolute<Impl>::value,int>::type = 0>

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
    using AuxDataType = typename GetAuxDataType<LatP>::type;
    using ScalarProductReturnType = typename GetScPType<LatP>::type;
    private:
    explicit constexpr GeneralLatticePoint()=default; //only callable from its friends
    public:


    // This is just to match the implementation of a typical instantiation.
    // Note the the deleted copy constructors and copy assignments prevents default copying
    // derived classes (since the base class must be copied as well).
    GeneralLatticePoint(GeneralLatticePoint const &other)=delete;
    GeneralLatticePoint(GeneralLatticePoint &&other)=default;
    GeneralLatticePoint& operator=(GeneralLatticePoint const & other) = delete;
    GeneralLatticePoint& operator=(GeneralLatticePoint && other) = default;
    protected:
    ~GeneralLatticePoint()=default;
    public:

// This one should be overloaded by every derived class.
// It is used to initalized the static data members.

// class_init may keep a counter of how may times it was called that is decremented by class_uninit.
// So for every call of class_init, there needs to be a call of class_uninit.
// If the counter is >0, calling class_init with different aux_data may fail.
// The return value of class_init indicates success. The return value of class_uninit indicates
// whether the counter is now at 0.
// Note: If there is no static data, the implementation may just always return true and not use
// a static counter at all.

    static bool class_init(AuxDataType const &aux_data) = delete;
    static bool class_uninit() = delete;

// This one should be overloaded by every derived class. Used for diagnostic.

    static std::string class_name() {return "General Lattice Point.";};

// comparison with < or > are by length (by default).
// Note that == or != are intendend to mean actual equality comparisons,
// so P1 <= P2 <= P1 does not imply P1==P2.

// These operations may be overloaded by the implementation.
// We require < and > to be strict weak orderings.

    bool operator< ( LatP const &x1 ) const
    {
      DEBUG_TRACEGENERIC("Generically comparing < for" << LatP::class_name() )
      return CREALTHIS->get_norm2() < x1.get_norm2();
    }

    bool operator> ( LatP const &x1 ) const
    {
      DEBUG_TRACEGENERIC("Generically comparing > for" << LatP::class_name() )
      return CREALTHIS->get_norm2() > x1.get_norm2();
    }

    bool operator<= ( LatP const &x1 ) const
    {
      DEBUG_TRACEGENERIC("Generically comparing <= for" << LatP::class_name() )
      return CREALTHIS->get_norm2() <= x1.get_norm2();
    }

    bool operator>= ( LatP const &x1 ) const
    {
      DEBUG_TRACEGENERIC("Generically comparing >= for" << LatP::class_name() )
      return CREALTHIS->get_norm2() >= x1.get_norm2();
    }

    // get_vec_size() returns the number of coordinates we get when using latp[i]
    // get_dim returns the (ambient) dimension the vector is supposed to represent.
    // By default, vec_size is the same as dim.
    // get_dim must be overloaded.


    MEMBER_ONLY_EXISTS_IF_COO_READ
    auto get_vec_size() const -> decltype( std::declval<Impl>().get_dim() )
    {
      DEBUG_TRACEGENERIC("Generically getting vec_size for" << LatP::class_name() )
      return CREALTHIS->get_dim();
    }

    void get_dim() const = delete; // Note that we have void return type. I really want C++14 auto.

/**
  Used for output to stream. Note that operator<< calls this (or an overloaded version)
  Supposed to be overloaded for the specific class.
*/

    MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE // This may be too strict.
    std::ostream& write_to_stream(std::ostream &os, bool include_norm2=true) const
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using generic writer for " << LatP::class_name() )
      auto const dim = CREALTHIS->get_vec_size();
      os << "[ "; // makes spaces symmetric
      for (uint_fast16_t i =0; i<dim; ++i)
        {
            os << CREALTHIS->operator[](i) << " ";
        }
        os <<"]";
      if(include_norm2)
      {
        os <<", norm2= " << CREALTHIS->get_norm2();
      }
      // No endl here (this is the caller's job).
      return os;
    }

    std::istream& read_from_stream(std::istream &is) = delete;

/**
  Fills a lattice point with zeros.
  Note that freshly constructed lattice points may contain uninitialized values.
  The latter depends subtly on the constructors used (empty vs. default constructor...)
  and may differ depending on whether nfixed==-1.

  Default version exists if Derived[i] is defined and writeable.
  Note that we additionally assume that we may assign 0 as in Derived[i] = 0;

  May be overloaded by Derived class.
*/

    MEMBER_ONLY_EXISTS_IF_COO_WRITE
    void fill_with_zero()
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using generic fill with zero for " << LatP::class_name() )
      auto const dim = CREALTHIS->get_vec_size();
      for (uint_fast16_t i=0;i<dim;++i)
      {
        REALTHIS->operator[](i) = 0;
      }
      REALTHIS->sanitize(0);
    }

    MEMBER_ONLY_EXISTS_IF_COO_WRITE
    void make_negative()
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using generic negation function for " << LatP::class_name() )
      auto const dim = CREALTHIS->get_vec_size();
      for (uint_fast16_t i=0;i<dim;++i)
      {
        REALTHIS->operator[](i) = - REALTHIS->operator[](i);
      }
      if( IsNegateCheap<Impl>::value) // constexpr if
      {
        return;
      }

      if( !IsNegateCheap<Impl>::value ) // constexpr if, really...
      {
        if(IsNorm2Cheap<Impl>::value)
        {
          REALTHIS->sanitize(CREALTHIS->get_norm2() );
        }
        else
        {
          REALTHIS->sanitize();
        }

      }

      return;
    }

/**
  Tests whether a lattice point is all-zero.
  Defaults to checking component-wise by doing != 0 - comparison.
  May be overloaded.
*/
    MEMBER_ONLY_EXISTS_IF_COO_READ
    bool is_zero() const
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using (possibly inefficient) test for zero for " << LatP::class_name() )
      // constexpr if, really...
      if (IsNorm2Cheap<LatP>::value)
      {
        return (CREALTHIS->get_norm2() == 0);
      }
      else
      {
        auto const dim = CREALTHIS->get_vec_size();
        for (uint_fast16_t i=0;i<dim;++i)
        {
          if(CREALTHIS->operator[](i) != 0)
          {
            return false;
          }
        }
      return true;
      }
    }

/**
  Makes an (explicit) copy of the current point.

  Default: Copying components as a vector and sanitize additional data.
  (The latter probably means needless recomputing of data that the original already had)
  With the intended semantics of [] and entries, this results in a deep copy.

  May be overloaded.
*/

    MEMBER_ONLY_EXISTS_IF_COO_WRITE
    LatP make_copy() const
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using generic copy for " << LatP::class_name() )
      auto const real_dim=CREALTHIS->get_dim(); // means ambient dimension.
      auto const dim = CREALTHIS->get_vec_size(); // number of coordinates stored. May be rank.
      LatP NewLP(real_dim);
      for (uint_fast16_t i=0; i<dim; ++i)
      {
        NewLP[i] = CREALTHIS->operator[](i);
      }
      if (IsNorm2Cheap<LatP>::value)
      {
        NewLP.sanitize(CREALTHIS->get_norm2() );
      }
      else
      {
        NewLP.sanitize();
      }
      return NewLP;
    }

    // brings the lattice point into a defined state. Defaults to "do nothing".

    // Note that this function must *not* be const, since it actually changes observable behaviour:
    // For efficiency reasons, sanitize is not called upon member access to LatP[]; rather,
    // sanitize is called from outside the class.
    // As such, the user can leave LatP in an invalid state and sanitize remedies that.

    // The second version takes norm2 as an argument (to avoid recomputing it).

    void sanitize() {};
    void sanitize(ScalarProductReturnType const &norm2) { sanitize(); };


/**
     obtains the squared length.
     By default, we use the scalar product to compute it.
     Usually norm2 will be precomputed, so this function should be
     overridden by LatP.
     Note that get_norm2() is const. If LatP decides to lazyly compute and store norm2 on demand,
     make norm2 a mutable (and be really wary of thread-safety...)
*/
    ScalarProductReturnType get_norm2() const
    {
      DEBUG_TRACEGENERIC("Generically computing norm2 for " << LatP::class_name() )

      // This function should not be called if IsNorm2Cheap is set.
      static_assert(IsNorm2Cheap<LatP>::value == false, "");
      return compute_sc_product(*(CREALTHIS),*(CREALTHIS) );
    }

/**
  multiplies by a scalar (which is of integral type or mpz_class)
*/
    template<class Integer, class Impl=LatP,
      typename std::enable_if<IsCooVector<Impl>::value && std::is_integral<Integer>::value,
      int>::type = 0>
    void scalar_multiply(Integer const multiplier)
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Generically scalar-multiplying for " << LatP::class_name() )
      auto const dim = CREALTHIS->get_vec_size();
      for(uint_fast16_t i=0;i<dim;++i)
      {
        REALTHIS->operator[](i) *= multiplier;
      }
      if(IsNorm2Cheap<LatP>::value) // constexpr if
      {
        REALTHIS->sanitize(CREALTHIS->get_norm2() * multiplier * multiplier );
      }
      else
      {
        REALTHIS->sanitize();
      }
    }


};

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

FOR_LATTICE_POINTS_LP1_LP2
bool operator==(LP1 const &x1, LP2 const &x2)
{
  DEBUG_TRACEGENERIC("Generically comparing norm2 for" << LP1::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = x1.get_vec_size();
  auto const dim2 = x2.get_vec_size();
  assert(dim1 == dim2);
  #endif // DEBUG_SIEVE_LP_MATCHDIM
  auto const dim = x1.get_vec_size();
  for(uint_fast16_t i=0;i<dim;++i)
  {
    if (x1[i] != x2[i])
    {
      return false;
    }
  }
  return true;
}

FOR_LATTICE_POINTS_LP1_LP2
bool operator!=(LP1 const &x1, LP2 const &x2)
{
  return !(x1==x2);
}


// dispatch to add / addval function.
// This is done as non-member functions to simplify cases where we add different types of points.

// addval(x1,x2) : x1+=x2;
// add(x1,x2) : x1 + x2;

// x1 += x2
FOR_LATTICE_POINTS_LP1_LP2
LP1& operator+=(LP1 &x1, LP2 const &x2) { return addval(x1, x2); }

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator+(LP1 const &x1, LP2 const &x2) { return add(x1,x2); }

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator+(LP1 && x1, LP2 const &x2)
{
LP1 tmp = std::move(x1);
return addval(tmp,x2);
}

// We don't want to return LP2 here...
// If this causes trouble, the caller should change the order of arguments.
FOR_LATTICE_POINT_LP
LP operator+(LP const &x1, LP && x2)
{
LP tmp = std::move(x2);
return addval(tmp,x1);
}

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator+(LP1 &&x1, LP2 && x2)
{
LP1 tmp = std::move(x1);
return addval(tmp,std::move(x2));
}

// Note: We could define general add in terms of addval as well.
// (i.e. define + as copy and += also for the non-move case)
// Doing it as above might be slightly faster due to memory-access and copying additional data.

// unary minus

FOR_LATTICE_POINT_LP
LP operator-(LP &&x1)
{
  LP tmp = std::move(x1);
  x1.make_negative();
  return x1;
}

template<class LP, class Integer,
  typename std::enable_if<IsALatticePoint<LP>::value &&
  (std::is_integral<Integer>::value || std::is_same<Integer,mpz_class>::value),
  int>::type=0>
LP operator*(LP const &x1, Integer const multiplier)
{
//  assert(false);
  LP tmp = x1.make_copy();
  tmp.scalar_multiply(multiplier);
  return tmp;
}

template<class LP, class Integer,
  typename std::enable_if<IsALatticePoint<LP>::value &&
  (std::is_integral<Integer>::value || std::is_same<Integer,mpz_class>::value),
  int>::type=0>
LP operator*(LP &&x1, Integer const multiplier)
{
  LP tmp = std::move(x1);
  tmp.scalar_multiply(multiplier);
  return tmp;
}

// dispatch to sub/subval function

FOR_LATTICE_POINTS_LP1_LP2
LP1& operator-=(LP1 &x1, LP2 const &x2){ return subval(x1,x2); }

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 const &x1, LP2 const &x2){ return sub(x1,x2); }

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 && x1, LP2 const &x2)
{
  LP1 tmp = std::move(x1);
  subval(tmp,x2);
  return tmp;
}

FOR_LATTICE_POINT_LP
LP operator-(LP const &x1, LP &&x2)
{
  LP tmp = std::move(x2);
  tmp.make_negative();
  addval(tmp,x1);
  return tmp; //addval(tmp,x1);
}

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 && x1, LP2 && x2)
{
  LP1 tmp = std::move(x1);
  subval(tmp,std::move(x2));
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

// x1+=x2 for x1, x2 of the same type
template<class LP, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type = 0>
LP& addval(LP &x1, LP const &x2)
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
//  auto const real_dim = x1.get_dim();
//  LP NewLP(real_dim);
  for(uint_fast16_t i = 0; i < dim; ++i )
  {
    x1[i] += x2[i];
  }
  x1.sanitize();
  return x1;
}

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

template<class LP, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type = 0>
LP& subval(LP &x1, LP const &x2)
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

  for(uint_fast16_t i=0; i < dim ; ++i )
  {
    x1[i] -= x2[i];
  }
  x1.sanitize();
  return x1;
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
    result[i] = static_cast<ET>( container[i].get_data() );
  }
  result.sanitize();
  return result;
}

template<class LP, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type = 0>
typename GetCooType<LP>::type compute_sc_product(LP const &lp1, LP const &lp2)
{
  DEBUG_TRACEGENERIC("Generically computing scalar product for" << LP::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = lp1.get_dim();
  auto const dim2 = lp2.get_dim();
  assert(dim1 == dim2 );
  #endif // DEBUG_SIEVE_LP_MATCHDIM
  using ET = typename GetCooType<LP>::type;
  auto const dim = lp1.get_dim();
  ET result = 0; // assumes that ET can be initialized from 0...
  for(uint_fast16_t i=0; i<dim; ++i)
  {
    result += lp1[i]*lp2[i];
  }
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



// cleaning up internal macros.
#undef MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE
#undef MEMBER_ONLY_EXISTS_IF_COO_READ
#undef MEMBER_ONLY_EXISTS_IF_COO_WRITE
#undef IMPL_IS_LATP
#undef FOR_LATTICE_POINT_LP
#undef FOR_LATTICE_POINTS_LP1_LP2
#undef CREALTHIS
#undef REALTHIS

#endif

//clang-format on