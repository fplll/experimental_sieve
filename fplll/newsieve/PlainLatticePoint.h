#ifndef PLAIN_LATTICE_POINT_H
#define PLAIN_LATTICE_POINT_H

#include "LatticePointConcept.h"
#include "Utility.h"

//most simple Lattice Point that just wraps around a vector / array of ET's.

template<class ET, int nfixed> class PlainLatticePoint;

template<class ET,int nfixed>
class ImplementationTraits< PlainLatticePoint<ET,nfixed> > : public ImplementationTraitsBase
{
    public:
    using AuxDataType = Dimension<nfixed>;
    using ScalarProductReturnType = ET;
};

template<class ET, int nfixed> class PlainLatticePoint;

template<class ET, int nfixed>
class PlainLatticePoint : public GeneralLatticePoint< PlainLatticePoint<ET,nfixed> >
{

    public:
    using LatticePointTag = std::true_type;
    using AuxDataType = typename ImplementationTraits<PlainLatticePoint>::AuxDataType;
    using ScalarProductReturnType = typename ImplementationTraits<PlainLatticePoint>::ScalarProductReturnType;
    static void class_init(AuxDataType const & aux_data)
    {
//        dim = aux_data;
    };

};



#endif
