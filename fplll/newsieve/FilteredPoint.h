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

// using namespace LatticeApproximations;

/*  Never put "using namespace" declarations at file scope into header files.
    The issue is that any file that #includes this, has the using namespace declaration in effect...
    This is equivalent to not using namespaces at all. -- Gotti
*/

template <class ET, class SC> class FilteredPoint;

//template <class ET, bool insideMTList=false, int n_fixed=-1>
template <class ET, class SC>
class FilteredPoint
{
    public:

    FilteredPoint()=default;
    FilteredPoint(const FilteredPoint &Point) = default; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    FilteredPoint(FilteredPoint &&Point) = default ;
    FilteredPoint(ApproxLatticePoint<ET> x, SC sc, bool sign)
    {
        this->point = x;
        this->sc_prod = sc;
        this->minus = sign;
    }

    
    /*
    FilteredPoint(ApproxLatticePoint<ET> x, LatticeApproximations::ApproxTypeNorm2 sc)
    {
        this->point = x;
        this->sc_prod = sc;
    }
    

    FilteredPoint(ApproxLatticePoint<ET> x, float sc)
    {
        this->point = x;
        this->sc_prod = sc;
    }
    */


    //FilteredPoint(ApproxLatticePoint x, ApproxLatticePoint p)


    FilteredPoint& operator=(FilteredPoint const &that) =default;
    FilteredPoint& operator=(FilteredPoint && that) =default;


    ~FilteredPoint() {}


    inline ApproxLatticePoint<ET>  getApproxVector() const {return this->point;}
    inline SC get_sc_prod() const {return sc_prod;}
    inline bool get_sign() const {return minus;}


private:
    //members
    ApproxLatticePoint<ET> point;
    
    // always positive
    SC sc_prod;
    
    //true is sc_prod is correct for point
    // false if for -point
    bool minus;



};

#endif
