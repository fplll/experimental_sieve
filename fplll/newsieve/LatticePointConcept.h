#ifndef LATTICE_POINT_CONCEPT_H
#define LATTICE_POINT_CONCEPT_H

//TODO: Remove the concept idea below, we replace it by CRTP (i.e. inheritance without virtual functions, using templates)

#include "Utility.h"
#include <iostream>
#include <string>


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



//This class template stores the typedefs that the individual lattice point classes have
//There has to be a specialization for each lattice point class
class ImplementationTraitsBase
{
    public:
    using AuxDataType = IgnoreAnyArg;
};

template<class Implementation> class ImplementationTraits;
template<class Implementation> class GeneralLatticePoint;

CREATE_MEMBER_TYPEDEF_CHECK_CLASS(ScalarProductReturnType, DeclaresScalarProductReturnType);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(LatticePointTag, std::true_type, IsALatticePoint);

template<class Implementation>
class GeneralLatticePoint
{
    public:
    //using LatticePointTag = std::true_type;
    using AuxDataType = typename ImplementationTraits<Implementation>::AuxDataType;
    static_assert(DeclaresScalarProductReturnType<ImplementationTraits<Implementation>>::value, "Lattice Point class does not typedef its scalar product type");
    using ScalarProductReturnType = typename ImplementationTraits<Implementation>::ScalarProductReturnType;
    explicit constexpr GeneralLatticePoint()=default; //only ever called from its children
    GeneralLatticePoint(GeneralLatticePoint const &other)=delete; //This is just to match the implementation of a typical instantiation. Only the default constructor is ever used.
    GeneralLatticePoint(GeneralLatticePoint &&other)=default;
    GeneralLatticePoint& operator=(GeneralLatticePoint const & other) = delete;
    GeneralLatticePoint& operator=(GeneralLatticePoint && other) = default;
    ~GeneralLatticePoint()=default;

    //The calling syntax for the constructor of the derived object(i.e. Implementation) is supposed to be
    //Implementation(DIM dim, AUX aux={}),
    //where DIM is convertible from int
    //and AUX is convertible from Implementation::AuxDataType

    //Note:
    //Once the derived class Implementation defines a member function with the same name f as a function below,
    //then overload resolution to a call LP.f(blah) will *first* look at whatever is declared by Implementation.
    //If any viable function is found, GeneralLatticePoints' functions are no longre looked at.
    //This is exactly what we want, s.t. derived functions can use types that are convertible from AuxDataType instead.
    //
    //A using GeneralLatticePoint - directive in the derived class changes that behavior, so you may not want to do that.

    Implementation make_copy(AuxDataType const & aux_data={}) = delete;
    ScalarProductReturnType get_norm2(AuxDataType const & aux_data={}) const = delete;
    unsigned int get_dim(AuxDataType const &aux_data={}) const = delete;
    std::istream & read_from_stream(std::istream &is = std::cin, AuxDataType const &aux_data={})=delete;
    std::ostream & write_to_stream(std::ostream &os = std::cout, AuxDataType const &aux_data={}); //=delete;
    static std::string class_name() {return "General Lattice Point";};
    void increment_by(GeneralLatticePoint const &how_much, AuxDataType const &aux_data={})=delete;
    void decrement_by(GeneralLatticePoint const &how_much, AuxDataType const &aux_data={})=delete;
    bool operator< (Implementation const & other) const;
    bool operator> (Implementation const & other) const;
    void sanitize();
    Implementation& operator+(Implementation const &second) = delete;
    Implementation& operator+(Implementation &&second) = delete;
    Implementation& operator-(Implementation const &second) = delete;
    Implementation& operator-(Implementation &&second) = delete;
    Implementation& operator-() = delete;
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


//This code fails
/*
template<class Implementation>
bool GeneralLatticePoint<Implementation>::operator< (Implementation const & other) const
{
    return Implementation:: this->get_norm2() < other.get_norm2();
}

template<class Implementation>
bool GeneralLatticePoint<Implementation>::operator> (Implementation const & other) const
{
    return Implementation::get_norm2(*this) > other.get_norm2();
}
 */



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



#endif
