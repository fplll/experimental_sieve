//This is the header file for the Gauss Sieve (both single-threaded and multi-threaded version).
//Note:
//With C++ lacking a good "static_if" / "constexpr if" that is broadly supported by compilers (as of early 2017), we use some pre-processor shenanigans to avoid code duplication:

//TODO : explain

#define SIEVE_GAUSS_SINGLE_THREADED

#ifndef SIEVE_GAUSS_H
#define SIEVE_GAUSS_H

//SIEVE_GAUSS_SINGLE_THREADED : if defined by user, means that we request the single-threaded implementation
//SIEVE_GAUSS_MULTI_THREADED  : if defined by user, means that we request the multi -threaded implementation
//If neither is set, we default to both.

//SIEVE_GAUSS_DEFAULT_THREADED: if set to true, we default to multi-threaded. Set by ourself.

//GAUSS_SIEVE_IS_MULTI_THREADED : macro constant (true/false) that is internally used. Within SieveJoint.*, indicates whether we currently generate code for
//the single- or multi-threaded version.

//default values: instantiate both, default to multi-threaded
#if !defined (SIEVE_GAUSS_SINGLE_THREADED) && !defined (SIEVE_GAUSS_MULTI_THREADED)
#define SIEVE_GAUSS_SINGLE_THREADED
#define SIEVE_GAUSS_MULTI_THREADED
#define SIEVE_GAUSS_DEFAULT_THREADED true
#elif defined(SIEVE_GAUSS_MULTI_THREADED)
#define SIEVE_GAUSS_DEFAULT_THREADED true
#elif defined(SIEVE_GAUSS_SINGLE_THREADED)
#define SIEVE_GAUSS_DEFAULT_THREADED false
#endif

template<class ET, bool MultiThreaded, int nfixed=-1>
class Sieve;

template<class ET>
using SieveGauss = Sieve<ET,SIEVE_GAUSS_DEFAULT_THREADED>;

template<class ET>
using SieveST = Sieve<ET,false>;

template<class ET>
using SieveMT = Sieve<ET,true>;


#ifdef  GAUSS_SIEVE_IS_MULTI_THREADED
#undef  GAUSS_SIEVE_IS_MULTI_THREADED
#endif

#ifdef SIEVE_GAUSS_SINGLE_THREADED
#define GAUSS_SIEVE_IS_MULTI_THREADED false
#include "SieveJoint.h"
//#include "SieveJoint.cpp"
//#include "SieveST.cpp"
#undef GAUSS_SIEVE_IS_MULTI_THREADED
#endif // SIEVE_GAUSS_SINGLE_THREADED

#ifdef SIEVE_GAUSS_MULTI_THREADED
#define GAUSS_SIEVE_IS_MULTI_THREADED true
#include "SieveJoint.h"
//#include "SieveJoint.cpp"
//#include "SieveMT.cpp"
#undef GAUSS_SIEVE_IS_MULTI_THREADED
#endif // SIEVE_GAUSS_MULTI_THREADED


//global declarations go here.



#endif // SIEVE_GAUSS_H
