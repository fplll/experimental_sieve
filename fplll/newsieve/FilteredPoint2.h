// clang-format off

//
//  FilteredPoint.h
//
//
//  Created by Elena on 06/03/17.
//
//

#ifndef _FilteredPoint2_h
#define _FilteredPoint2_h

#include "DebugAll.h"
#include "SieveUtility.h"
#include "ExactLatticePoint.h"

namespace GaussSieve{

template <class ET, int nfixed> class ExactLatticePoint;

template <class ET, int nfixed, class EntryType> class FilteredPointPointer;

// Template parameters are:
//  ET: entry type
//  nfixed: indicates whether the dimension is fixed or not
//  SC: scalar-product type

template <class ET, int nfixed, class SC>
class FilteredPointPointer
{
public:
    
    using StoredPoint = ExactLatticePoint<ET,nfixed>;

    FilteredPointPointer()=delete;
    FilteredPointPointer(const FilteredPointPointer &Point) = delete; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    FilteredPointPointer(FilteredPointPointer &&Point) = default ;
    FilteredPointPointer(StoredPoint const* x, SC sc)
    {
        this->point = x;
        this->sc_prod = sc;
    }


    FilteredPointPointer& operator=(FilteredPointPointer const &that) =delete;
    FilteredPointPointer& operator=(FilteredPointPointer && that) =default;


    ~FilteredPointPointer() {}

    inline StoredPoint const* get_point() const {return this->point;}
    inline SC get_sc_prod() const {return sc_prod;}
    //inline bool get_sign() const {return minus;}


private:
    //members
    StoredPoint const* point;

    SC sc_prod;

    // true if sc_prod is correct for point
    // false if for -point
    // not used now
    //bool minus;



};

}

#endif