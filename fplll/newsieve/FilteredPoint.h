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

#include "sieve_common.h"
#include "LatticePoint2.h"


template <class ET, int nfixed, class SC> class FilteredPoint;

template <class ET, int nfixed, class SC>
class FilteredPoint
{
    public:

    FilteredPoint()=default;
    FilteredPoint(const FilteredPoint &Point) = default; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    FilteredPoint(FilteredPoint &&Point) = default ;
    FilteredPoint(ApproxLatticePoint<ET, nfixed> x, SC sc, bool sign)
    {
        this->point = x;
        this->sc_prod = sc;
        this->minus = sign;
    }


    FilteredPoint& operator=(FilteredPoint const &that) =default;
    FilteredPoint& operator=(FilteredPoint && that) =default;


    ~FilteredPoint() {}


    inline PlainLatticePoint<ET, nfixed>  getApproxVector() const {return this->point;}
    inline SC get_sc_prod() const {return sc_prod;}
    inline bool get_sign() const {return minus;}


private:
    //members
    PlainLatticePoint<ET, nfixed> point;

    // always positive
    SC sc_prod;

    //true is sc_prod is correct for point
    // false if for -point
    bool minus;



};

#endif
