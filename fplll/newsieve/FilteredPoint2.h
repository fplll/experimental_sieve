// clang-format off

//
//  FilteredPoint.h
//
//
//  Created by Elena on 06/03/17.
//
//

#ifndef FILTERED_POINT2_H
#define FILTERED_POINT2_H

#include "DefaultIncludes.h"
#include "DebugAll.h"
#include "SieveUtility.h"
#include "ExactLatticePoint.h"
//#include "HashedLatticePoint.h"
#include "GaussListBitapprox.h"
#include "BitApproximationNew.h"
#include "Typedefs.h"

namespace GaussSieve{

//template <class ET, int nfixed> class ExactLatticePoint;
//template <class ET, int nfixed> class HashedLatticePoint;

template <class SieveTraits, class ET> class FilteredPoint;

// Template parameters are:
//  ET: entry type
//  nfixed: indicates whether the dimension is fixed or not
//  SC: scalar-product type

template <class SieveTraits, class ET>
class FilteredPoint
{
  public:

//#ifndef USE_LSH
    //using SimHashes   = typename SimHashGlobalDataType::SimHashes;
    //using StoredPoint = typename SieveTraits::GaussList_StoredPoint;
    using StoredData    = STNode< SieveTraits>;
//#else
 //   using StoredPoint = HashedLatticePoint<ET,nfixed>;
//#endif

    FilteredPoint()=delete;
    FilteredPoint(const FilteredPoint &Point) = delete; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    FilteredPoint(FilteredPoint &&Point) = default ;


    FilteredPoint(StoredData const * pt, ET sc_pr)
    {
      this->point   = pt;
      this->sc_prod = sc_pr;
    }

    FilteredPoint& operator=(FilteredPoint const &that) =delete;
    FilteredPoint& operator=(FilteredPoint && that) =default;


    ~FilteredPoint() {}

    inline StoredData const& get_point() const {return *point;}
    inline ET  get_sc_prod() const {return sc_prod;}
    //inline StoredData& get_ptr_toexact() const {return ptr_to_exact;}
    //inline bool get_sign() const {return minus;}


private:
    //members
    //SimHashes const* approx_point; //may be copy?
    //StoredPoint* ptr_to_exact;

    StoredData const * point;
    ET sc_prod;





};

}

#endif
