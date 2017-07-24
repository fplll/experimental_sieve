// This header defines (global) traits. We might want to encapsulte these in some trait classes that
// get forwarded, eventually.
// This header should not depend on anything.

// Consider renaming the file to avoid clashes.

#ifndef GAUSS_SIEVE_TYPEDEFS_H
#define GAUSS_SIEVE_TYPEDEFS_H

namespace GaussSieve
{
// unfortunately, trigonometric functions to compute pi don't have constexpr variants on all
// compilers we want to support, so we just define pi directly
long double constexpr pi_long = 3.14159265358979323846264338327950288419716939937510L;
double constexpr pi_double    = 3.14159265358979323846264338327950288419716939937510;
long double constexpr pi      = 3.14159265358979323846264338327950288419716939937510L;

// forward-declarations:
template <class ET, int nfixed> class MyLatticePoint;

// various typedef declarations that control the types used by our classes.

// lines are too long, clang-format destroys vertical alignment
// clang-format off

template <class ET, bool MT, int nfixed> using GaussSampler_ReturnType = MyLatticePoint<ET, nfixed>;

template <class ET, bool MT, int nfixed> using GaussList_ReturnType    = MyLatticePoint<ET, nfixed>;
template <class ET, bool MT, int nfixed> using GaussList_StoredPoint   = MyLatticePoint<ET, nfixed>;

template <class ET, bool MT, int nfixed> using GaussQueue_ReturnType   = GaussSampler_ReturnType<ET, MT, nfixed>;
template <class ET, bool MT, int nfixed> using GaussQueue_DataType     = GaussQueue_ReturnType<ET, MT, nfixed>;

// for a small number of lattice points that we need to access very often.
template <class ET, bool MT, int nfixed> using FastAccess_Point = MyLatticePoint<ET, nfixed>;

// clang-format on

};

#endif
