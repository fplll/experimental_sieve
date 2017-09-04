// clang-format off

//
//  FilteredPoint.h
//
//
//  Created by Elena on 06/03/17.
//
//

#ifndef _FilteredPoint_h
#define _FilteredPoint_h

#include "DebugAll.h"
#include "SieveUtility.h"
#include "PlainLatticePoint.h"


template <class ET, int nfixed> class PlainLatticePoint;

template <class ET, int nfixed, class SC> class FilteredPoint;

// Template parameters are:
//  ET: entry type
//  nfixed: indicates whether the dimension is fixed or not
//  SC: scalar-product type

template <class ET, int nfixed, class SC>
class FilteredPoint
{
public:
    
    using StoredPoint = PlainLatticePoint<ET, nfixed>;

    FilteredPoint()=delete;
    FilteredPoint(const FilteredPoint &Point) = delete; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    FilteredPoint(FilteredPoint &&Point) = default ;
    FilteredPoint(StoredPoint x, SC sc, bool sign)
    {
        this->point = x;
        this->sc_prod = sc;
        this->minus = sign;
    }


    FilteredPoint& operator=(FilteredPoint const &that) =delete;
    FilteredPoint& operator=(FilteredPoint && that) =default;


    ~FilteredPoint() {}


    inline StoredPoint get_point() const {return this->point;}
    inline SC get_sc_prod() const {return sc_prod;}
    inline bool get_sign() const {return minus;}


private:
    //members
    StoredPoint point;

    // always positive
    SC sc_prod;

    //true if sc_prod is correct for point
    // false if for -point
    bool minus;



};

#endif
