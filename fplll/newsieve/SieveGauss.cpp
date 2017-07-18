#define SIEVE_GAUSS_SINGLE_THREADED

#ifndef SIEVE_GAUSS_CPP
#define SIEVE_GAUSS_CPP

#include "SieveGauss.h" //This calls SieveJoint.h, possiby twice. We also set some preprocessor macros.

#ifdef  GAUSS_SIEVE_IS_MULTI_THREADED
#undef  GAUSS_SIEVE_IS_MULTI_THREADED
#endif


#ifdef SIEVE_GAUSS_SINGLE_THREADED
#define GAUSS_SIEVE_IS_MULTI_THREADED false
//#include "SieveJoint.h"
#include "SieveJoint.cpp"
#include "SieveST.cpp"
#undef GAUSS_SIEVE_IS_MULTI_THREADED
#endif // SIEVE_GAUSS_SINGLE_THREADED

#ifdef SIEVE_GAUSS_MULTI_THREADED
#define GAUSS_SIEVE_IS_MULTI_THREADED true
//#include "SieveJoint.h"
#include "SieveJoint.cpp"
#include "SieveMT.cpp"
#undef GAUSS_SIEVE_IS_MULTI_THREADED
#endif // SIEVE_GAUSS_MULTI_THREADED

#endif
