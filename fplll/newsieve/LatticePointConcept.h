
#ifndef LATTICE_POINT_CONCEPT_H
#define LATTICE_POINT_CONCEPT_H

//TODO: Remove the concept idea below, we replace it by CRTP (i.e. inheritance without virtual functions, using templates)

#include "SieveUtility.h"
#include <iostream>
#include <string>
#include "DebugAll.h"
#include "assert.h"

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

  AuxDataType : type of class-wide data that is required to interpret a given point
                (e.g. Dimension, custom memory allocator)
                Default: IgnoreAnyArg
  ScalarProductReturnType: Return type of scalar products. Mandatory!
  CoordinateAccess : Set to true_type to indicate that the class exposes an operator[].
  CoordinateVector : Set to true_type to indicate that the class exposes an operator[] and further
                     that the class is a "vector with extra data". This means we have
                     get_dim() coordinate entries and write operations to these behave linearly.
                     After calling sanitize(), the class will be in a valid state.
  AbsoluteCoos : Set to true_type to indicate that the class exposes an operator[] and that
                 the entries are absolute coordinates.
  CoordinateType : return type of operator[] if available.
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
CREATE_TRAIT_CHECK_CLASS(LatticePointTraits, CoordinateType, DoesDeclareCoordinateType);
MAKE_TRAIT_GETTER(LatticePointTraits, CoordinateType, void, GetCooType);

#ifdef MEMBER_ONLY_EXISTS_IF_COOS_ARE_ABSOLUTE
#error Macros clash
#endif // MEMBER_ONLY_EXISTS_IF_COOS_ARE_ABSOLUTE
#ifdef MEMBER_ONLY_EXISTS_IF_IS_COO_VECTOR
#error Macros clash
#endif
#ifdef IMPL_IS_LATP
#error Macros clash
#endif

#define MEMBER_ONLY_EXISTS_IF_COOS_ARE_ABSOLUTE \
template<class Impl=LatP, typename std::enable_if<CoosAreAbsolute<Impl>::value,int>::type = 0>

#define MEMBER_ONLY_EXISTS_IF_IS_COO_VECTOR \
template<class Impl=LatP, typename std::enable_if<IsCooVector<Impl>::value,int>::type = 0>

#define IMPL_IS_LATP \
static_assert(std::is_same<Impl,LatP>::value,"Using template member function with wrong type")

template<class LatP>
class GeneralLatticePoint
{
    static_assert(!HasNoLPTraits<LatP>::value, "Trait class not specialized.");
//    Does not work: We cannot access typedefs of our children directly (hence the traits class)
//    static_assert(IsALatticePoint<LatP>::value,"Could not recognize lattice point class.");
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
    //This is just to match the implementation of a typical instantiation.
    GeneralLatticePoint(GeneralLatticePoint const &other)=delete;
    GeneralLatticePoint(GeneralLatticePoint &&other)=default;
    GeneralLatticePoint& operator=(GeneralLatticePoint const & other) = delete;
    GeneralLatticePoint& operator=(GeneralLatticePoint && other) = default;
    ~GeneralLatticePoint()=default;
    static void class_init(AuxDataType const &aux_data) = delete;
    static constexpr std::string class_name() {return "General Lattice Point";};



//    bool operator< ( LatP const &x1 ) const
//    { return static_cast<LatP const*>(this)->get_norm2() < x1.get_norm2(); }


    MEMBER_ONLY_EXISTS_IF_COOS_ARE_ABSOLUTE
    std::ostream& write_to_stream(std::ostream &os) const
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using generic writer for " << LatP::class_name() )
      auto const dim = static_cast<LatP const*>(this)->get_dim();
      os << "[";
      for (unsigned int i =0; i<dim; ++i)
        {
            os << static_cast<LatP const*>(this)->operator[](i) << " ";
        }
        os <<"]" << std::endl;
      return os;
    }

/**
  Fills a lattice point with zeros.
  Note that freshly constructed lattice points may contain uninitialized values.
  The latter depends subtly on the constructors used (empty vs. default constructor...)
  and may differ depending on whether nfixed==-1.
*/

    MEMBER_ONLY_EXISTS_IF_IS_COO_VECTOR
    void fill_with_zero()
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using generic fill with zero for " << LatP::class_name() )
      auto const dim = static_cast<LatP*>(this)->get_dim();
      for (unsigned int i=0;i<dim;++i)
      {
        static_cast<LatP*>(this)->operator[](i) = 0;
      }
      static_cast<LatP*>(this)->sanitize();
    }

/**
  Tests whether a lattice point is all-zero
*/
    MEMBER_ONLY_EXISTS_IF_IS_COO_VECTOR
    bool is_zero()
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using (possibly inefficient) test for zero for " << LatP::class_name() )
      auto const dim = static_cast<LatP*>(this)->get_dim();
      for (unsigned int i=0;i<dim;++i)
      {
        if(static_cast<LatP*>(this)->operator[](i) != 0)
        {
          return false;
        }
      }
      return true;
    }

/**
  Makes an (explicit) copy of the current point;
*/

// This assumes that NewLP[i] actually allows assignment. Be aware that Z_NR has no move-assignment.
    MEMBER_ONLY_EXISTS_IF_IS_COO_VECTOR
    LatP make_copy() const
    {
      IMPL_IS_LATP;
      DEBUG_TRACEGENERIC("Using generic copy for " << LatP::class_name() )
      auto const dim=static_cast<LatP const*>(this)->get_dim();
      LatP NewLP(dim);
      for (unsigned int i=0; i<dim; ++i)
      {
        NewLP[i] = static_cast<LatP*>(this)->operator[](i);
      }
      NewLP.sanitize();
      return NewLP;
    }


    //Note:
    //Once the derived class Implementation defines a member function with the same name f as a function below,
    //then overload resolution to a call LP.f(blah) will *first* look at whatever is declared by Implementation.
    //If any viable function is found, GeneralLatticePoints' functions are no longre looked at.
    //This is exactly what we want, s.t. derived functions can use types that are convertible from AuxDataType instead.
    //
    //A using GeneralLatticePoint - directive in the derived class changes that behavior, so you may not want to do that.

//    LatP make_copy(AuxDataType const & aux_data={}) = delete;
//    ScalarProductReturnType get_norm2(AuxDataType const & aux_data={}) const = delete;
//    unsigned int get_dim(AuxDataType const &aux_data={}) const = delete;
//    std::istream & read_from_stream(std::istream &is = std::cin, AuxDataType const &aux_data={})=delete;
//    std::ostream & write_to_stream(std::ostream &os = std::cout, AuxDataType const &aux_data={})=delete; //=delete;
//    void increment_by(GeneralLatticePoint const &how_much, AuxDataType const &aux_data={})=delete;
//    void decrement_by(GeneralLatticePoint const &how_much, AuxDataType const &aux_data={})=delete;
//    bool operator< (Implementation const & other) const;// = delete;
//    bool operator> (Implementation const & other) const = delete;
    void sanitize() {};
    unsigned int get_dim() const = delete; //note: using unsigned int to return here!

//    Implementation& operator+(Implementation const &second) = delete;
//    Implementation& operator+(Implementation &&second) = delete;
//    Implementation& operator-(Implementation const &second) = delete;
//    Implementation& operator-(Implementation &&second) = delete;
//    Implementation& operator-() = delete;
};

//template<class Implementation>
//ScalarProductReturnType
//GeneralLatticePoint<Implementation>::get_norm2(AuxDataType const &aux_data) const
//{
//    ScalarProductReturnType res;
//    for(unsigned int i=0; i < Implementation::get_dim();++i)
//    {
//
//    }
//}

// dispatch to add function.

#ifdef FOR_LATTICE_POINT_LP
#error Macros clash
#endif // FOR_LATTICE_POINT_LP
#ifdef FOR_LATTICE_POINTS_LP1_LP2
#error Macros clash
#endif

#define FOR_LATTICE_POINT_LP \
template<class LP, typename std::enable_if<IsALatticePoint<LP>::value, int>::type=0>

#define FOR_LATTICE_POINTS_LP1_LP2 \
template<class LP1, class LP2, typename std::enable_if< \
         IsALatticePoint<LP1>::value && IsALatticePoint<LP2>::value,int>::type=0>


FOR_LATTICE_POINTS_LP1_LP2
LP1 operator+(LP1 const &x1, LP2 const &x2) { return add(x1,x2); }

FOR_LATTICE_POINT_LP
LP operator+(LP && x1, LP const &x2) { return add(std::move(x1),x2); }

FOR_LATTICE_POINT_LP
LP operator+(LP const &x1, LP && x2) { return add(std::move(x2),x1); }

// dispatch to sub function

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 const &x1, LP2 const &x2){ return sub(x1,x2); }

FOR_LATTICE_POINT_LP
LP operator-(LP && x1, LP const &x2) { return sub(std::move(x1),x2); }

FOR_LATTICE_POINT_LP
LP operator-(LP const &x1, LP &&x2) { return sub(x1,std::move(x2)); }

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

template<class LP, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type = 0>
LP add(LP const &x1, LP const &x2)
{
  DEBUG_TRACEGENERIC( "generically adding" << x1.class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = x1.get_dim();
  auto const dim2 = x2.get_dim();
  assert( dim1 == dim2 );
  #endif
  auto const dim = x1.get_dim();
  LP NewLP(dim);
  for(unsigned int i = 0; i < dim; ++i )
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
  auto const dim1 = x1.get_dim();
  auto const dim2 = x2.get_dim();
  assert(dim1==dim2);
  #endif
  auto const dim = x1.get_dim();
  LP NewLP(dim);
  for(unsigned int i=0; i < dim ; ++i )
  {
    NewLP[i] = x1[i] - x2[i];
  }
  NewLP.sanitize();
  return NewLP;
}

// this function can be used to initialize an LP with container types that allow []-access.
// Note that there is an explicit static_cast to LP's entry types.
// In particular, this can convert mpz_t to mpz_class...

template<class LP, class SomeContainer, int fixeddim, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::value=0>
LP make_from_any_vector(SomeContainer const &container, Dimension<fixeddim> dim)
{
  static_assert(DoesDeclareCoordinateType<LP>::value, "Not declaring coordinate types");
  using ET = typename GetCooType<LP>::type;
  DEBUG_TRACEGENERIC("generically converting vector to LP for" << LP::class_name() )
  LP NewLP(dim);
  for(unsigned int i =0; i<dim; ++i)
  {
    NewLP[i] = static_cast<ET>( container[i] );
  }
  NewLP.sanitize();
  return NewLP;
}

// Same as above, but un-Z_NR's the container.

template<class LP, class SomeZNRContainer, int fixeddim, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::value=0>
LP make_from_znr_vector(SomeZNRContainer const &container, Dimension<fixeddim> dim)
{
  static_assert(DoesDeclareCoordinateType<LP>::value, "Not declaring coordinate types");
  using ET = typename GetCooType<LP>::type;
  DEBUG_TRACEGENERIC("generically converting vector to LP and un-ZNRing for" << LP::class_name() )
  LP NewLP(dim);
  for(unsigned int i =0; i<dim; ++i)
  {
    NewLP[i] = static_cast<ET>( container[i].get_data() );
  }
  NewLP.sanitize();
  return NewLP;
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
  for(unsigned int i=0; i<dim; ++i)
  {
    result += lp1[i]*lp2[i];
  }
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


}

// cleaning up macros, since there is no namespace for them.

#undef MEMBER_ONLY_EXISTS_IF_COOS_ARE_ABSOLUTE
#undef MEMBER_ONLY_EXISTS_IF_IS_COO_VECTOR
#undef IMPL_IS_LATP
#undef FOR_LATTICE_POINT_LP
#undef FOR_LATTICE_POINTS_LP1_LP2

#endif

//clang-format on
