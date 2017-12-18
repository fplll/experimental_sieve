// Don't mess with include order:
// clang-format off

#ifndef SIEVE_GAUSS_CPP  // include guard for this file
#define SIEVE_GAUSS_CPP

#define SIEVE_GAUSS_SINGLE_THREADED

#include "fplll.h"

#include "SieveGauss.h" //This calls SieveJoint.h, possiby twice. We also set some preprocessor macros.
// This will then implicitly call the required other *.h files

#ifdef  GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
#undef  GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
#endif


#ifdef SIEVE_GAUSS_SINGLE_THREADED
#define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED false
//#include "SieveJoint.h"
#include "SieveJoint_impl.h"
#include "SieveST.cpp"
#undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
#endif // SIEVE_GAUSS_SINGLE_THREADED

#ifdef SIEVE_GAUSS_MULTI_THREADED
#define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED true
//#include "SieveJoint.h"
#include "SieveJoint_impl.h"
#include "SieveMT.cpp"
#undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
#endif // SIEVE_GAUSS_MULTI_THREADED

// We now include the implementation files after all header files...

#include "DefaultTermConds_impl.h"
#include "GaussQueue_impl.h"
#include "GPVSampler_impl.h"
#include "Sampler_impl.h"

#include "UniformSampler_impl.h"
#include "GPVSamplerExtended_impl.h"


#endif

//clang-format on
