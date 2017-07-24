// clang-format off

#ifndef LATTICE_POINT_CONCEPT_H
#define LATTICE_POINT_CONCEPT_H

//TODO: Remove the concept idea below, we replace it by CRTP (i.e. inheritance without virtual functions, using templates)

#include "SieveUtility.h"
#include <iostream>
#include <string>
#include "DebugAll.h"
#include "assert.h"


/**

TODO: Rewrite explanation, we use CRTP instead.

Lattice points represent points in the lattice.

We treat LatticePoints as a concept, in the sense that a lattice point is any class that implements the following interface:
(or a subset thereof)
(This is WIP, not sure what should be mandatory)
//Note that functions might also have a different interface with arguments that are convertible from AuxDataType


class LatticePoint{
public: //typedefs
using LatticePoint_Tag = true_type //required, because we static_assert that this is there or use it for SFINAE
using AuxDataType = ... //optional public typedef that specifies some class, defaults to IgnoreAnyArg
using ScalarProductReturnType //type returned by scalar products and norms

public: //member functions:
LatticePoint make_copy(AuxDataType& const aux_data); //actually makes a copy
LatticePoint(LatticePoint &&); //move constructor
LatticePoint(LatticePoint const &)=delete //copy constructor deleted. Always make copies explicitly
LatticePoint& operator=(LatticePoint const &old)=delete;
LatticePoint& operator=(LatticePoint && old);
ScalarProductReturnType get_norm2(AuxDataType const &aux_data);
void negate(AuxDataType); //negate function.
explicit LatticePoint(Dimension<nfixed> dim, AuxDataType const& aux_data); //creates an unitialized point of dimension dim.
void make_zero(AuxDataType);


}

//non-member functions:
LatticePoint add(LatticePoint const & A, LatticePoint const & B, A::AuxDataType);
LatticePoint subtract(LatticePoint const &A, LatticePoint const & B, A::AuxDataType);
A::ScalarProductReturnType compute_sc_product(LatticePoint const &A, LatticePoint const B, AuxDataType const &aux_data;
bool compare_abs_sc_product(LatticePoint const &A, LatticePoint const & B, ScalarProductReturnType target, AuxDataType const &aux_data); //compares whether the absolute value of the scalar product is at least target. This function might err.
bool compare_sc_product(LatticePoint const &A, LatticePoint const & B, ScalarProductReturnType target, AuxDataType const &aux_data); // Checks whether <A,B> > t

*/

namespace GaussSieve{

//This class template stores the typedefs that the individual lattice point classes have
//There has to be a specialization for each lattice point class
//class ImplementationTraitsBase
//{
//    public:
//    using AuxDataType = IgnoreAnyArg;
//};

// This class template stores the typedefs that the individual lattice point classes have
// There has to be a specialization for each lattice point class.
// The general template must never be instantiated.
template<class LatticePoint> struct LatticePointTraits
{
  public:
  using Invalid=std::true_type;
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

template<class LatP>
class GeneralLatticePoint
{
    static_assert(!HasNoLPTraits<LatP>::value);
//    static_assert(IsALatticePoint<LatP>::value,"Could not recognize lattice point class.");
    static_assert(HasScalarProductReturnType<LatP>::value, "Lattice Point class does not typedef its scalar product type");
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

template<class LP, typename std::enable_if<
         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
         int>::type = 0


         >
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
    NewLP[i].add(x1[i],x2[i]);
  }
  NewLP.sanitize();
  return NewLP;
}



template<class LP1, class LP2,
         typename std::enable_if<
          IsALatticePoint<LP1>::value && IsALatticePoint<LP2>::value
         ,int>::type=0>
LP1 operator+(LP1 const &x1, LP2 const &x2){ return add(x1,x2); }

/*
template<class LP,
         typename std::enable_if<
          IsALatticePoint<LP>::value
         ,int>::type=0>
LP operator+(LP && x1, LP const &x2) { return add(std::move(x1),x2); }

template<class LP,
         typename std::enable_if<
         IsALatticePoint<LP>::value
         ,int>::type=0>
LP operator+(LP const &x1, LP && x2) { return add(std::move(x2),x1); }
*/

//template<class LP,
//         typename std::enable_if<
//         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
//         int>::type = 0
//        >
//LP operator+(LP &&x1, LP const &x2)
//{
//  DEBUG_TRACEGENERIC( "generically adding; consume first" << x1.class_name() )
//  #ifdef DEBUG_SIEVE_LP_MATCHDIM
//  auto const dim1 = x1.get_dim();
//  auto const dim2 = x2.get_dim();
//  assert( dim1 == dim2 );
//  #endif
//  auto const dim = x1.get_dim();
//  LP NewLP(std::move(x1));
//  for(unsigned int i = 0; i < dim; ++i )
//  {
//    NewLP[i].add(NewLP[i],x2[i]); //there is no increment function in Z_NR...
//  }
//  NewLP.sanitize();
//  return NewLP;
//}
//
//template<class LP,
//         typename std::enable_if<
//         IsALatticePoint<LP>::value && IsCooVector<LP>::value,
//         int>::type = 0
//        >
//LP operator+(LP const x1, LP &&x2)
//{
//  DEBUG_TRACEGENERIC( "generically adding; consume second" << x1.class_name() )
//  #ifdef DEBUG_SIEVE_LP_MATCHDIM
//  auto const dim1 = x1.get_dim();
//  auto const dim2 = x2.get_dim();
//  assert( dim1 == dim2 );
//  #endif
//  auto const dim = x2.get_dim();
//  LP NewLP(std::move(x2));
//  for(unsigned int i = 0; i < dim; ++i )
//  {
//    NewLP[i].add(NewLP[i],x1[i]); //there is no increment function in Z_NR...
//  }
//  NewLP.sanitize();
//  return NewLP;
//}
//



//This code fails

//template<class Implementation>
//bool GeneralLatticePoint<Implementation>::operator< (Implementation const & other) const
//{
//  return (static_cast<Implementation const&>(*this)).get_norm2() < other.get_norm2();
//}

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
#endif

//clang-format on
