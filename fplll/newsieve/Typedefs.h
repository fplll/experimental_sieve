#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include "LatticePointsNew.h"
#include "MyLatticePointClass.cpp"

//various typedef declarations that control the types used by our classes.
//This typedef defines the return type that the Sampler should have.

namespace GaussSieve{
// unfortunately, trigonometric functions to compute pi don't have constexpr variants on all compilers we want to support, so we just define pi directly
long double constexpr   pi_long     = 3.14159265358979323846264338327950288419716939937510L;
double constexpr        pi_double   = 3.14159265358979323846264338327950288419716939937510;
long double constexpr   pi          = 3.14159265358979323846264338327950288419716939937510L;
};


namespace GaussSieve{

    template<class ET,bool MT, int nfixed> using GaussSampler_ReturnType = MyLatticePoint<ET,nfixed>;

    template<class ET,bool MT, int nfixed> using GaussList_ReturnType = MyLatticePoint<ET, nfixed>;
    template<class ET,bool MT, int nfixed> using GaussList_StoredPoint = MyLatticePoint<ET,nfixed>;
    //template<class ET,bool MT, int nfixed> using GaussList_DereferencesTo = ApproximateLatticePoint<ET,nfixed>;

    template<class ET,bool MT, int nfixed> using GaussQueue_ReturnType = GaussSampler_ReturnType<ET,MT,nfixed>;
    template<class ET,bool MT, int nfixed> using GaussQueue_DataType   = GaussQueue_ReturnType<ET,MT,nfixed>;

    template<class ET,bool MT, int nfixed> using FastAccess_Point      = MyLatticePoint<ET, nfixed>; //for a small number of lattice points that we need to access very often.

};


#endif
