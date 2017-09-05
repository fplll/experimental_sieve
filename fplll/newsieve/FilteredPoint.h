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
#include "ExactLatticePoint.h"

namespace GaussSieve{

template <class ET, int nfixed> class ExactLatticePoint;

template <class ET, int nfixed, class EntryType> class FilteredPoint;

// Template parameters are:
//  ET: entry type
//  nfixed: indicates whether the dimension is fixed or not
//  SC: scalar-product type

template <class ET, int nfixed, class SC>
class FilteredPoint
{
public:
    
    using StoredPoint = ExactLatticePoint<ET,nfixed>;

    FilteredPoint()=delete;
    FilteredPoint(const FilteredPoint &Point) = delete; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    FilteredPoint(FilteredPoint &&Point) = default ;
    FilteredPoint(StoredPoint x, SC sc)
    {
        //Store a pointer to point to avoid copying
        this->point = x.make_copy();
        this->sc_prod = sc;
    }


    FilteredPoint& operator=(FilteredPoint const &that) =delete;
    FilteredPoint& operator=(FilteredPoint && that) =default;


    ~FilteredPoint() {}


    //TODO: DO WE NEED TO COPY?
    inline StoredPoint get_point() const {return this->point.make_copy();}
    inline SC get_sc_prod() const {return sc_prod;}
    //inline bool get_sign() const {return minus;}


private:
    //members
    StoredPoint point;

    SC sc_prod;

    // true if sc_prod is correct for point
    // false if for -point
    // not used now
    //bool minus;



};

}

#endif
