/**
  Include this file to use the k-tuple GaussSieve for arbitrary template parameters.
*/

/**
  Note: If you edit this file, you will also have to edit Sieve.h and SieveGauss.cpp
*/

#ifndef SIEVE_GAUSS_H
#define SIEVE_GAUSS_H

/**
The macros
SIEVE_GAUSS_SINGLE_THREADED
and
SIEVE_GAUSS_SINGLE_THREADED
determine whether compile for single-threaded or multi-threaded mode or both.
If either macro is defined, we support the corresponding mode.
If neither is set, we fall back to a default (which is "both")
*/

// temporary, for as long as we do not support multithreaded variant.
#define SIEVE_GAUSS_SINGLE_THREADED
#ifdef SIEVE_GAUSS_MULTI_THREADED
#error not supported yet
#endif


// If neither of the above macros is set, we set a default value (compile both).
// We also set a default value for the template argument:
// SIEVE_GAUSS_DEFAULT_THREADED is set to true or false and determines the default template parameter
// (which is multi-threaded iff we compiled it)

#if !defined (SIEVE_GAUSS_SINGLE_THREADED) && !defined (SIEVE_GAUSS_MULTI_THREADED)
  #define SIEVE_GAUSS_SINGLE_THREADED
  #define SIEVE_GAUSS_MULTI_THREADED
  #define SIEVE_GAUSS_DEFAULT_THREADED true
#elif defined(SIEVE_GAUSS_MULTI_THREADED)
  #define SIEVE_GAUSS_DEFAULT_THREADED true
#elif defined(SIEVE_GAUSS_SINGLE_THREADED)
  #define SIEVE_GAUSS_DEFAULT_THREADED false
#endif

// We forward declare some templates and set default parameters.

#include "DefaultIncludes.h"
#include "Typedefs.h"

// we now include the "main (internal) header file SieveJoint.h for the Gauss Sieve.
// This file may be included twice with different values
// Note that this header file requires GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED to be set to either true or false.
// SieveJoint.h and SieveJoint_impl.h are designed in a way that allows them to be included twice
// with a different value of this macro.

#ifdef  GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
  #error Must never happen
#endif


#ifdef SIEVE_GAUSS_SINGLE_THREADED
  #define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED false
  #include "SieveJoint.h"
  #undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
#endif // SIEVE_GAUSS_SINGLE_THREADED

#ifdef SIEVE_GAUSS_MULTI_THREADED
  #define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED true
  #include "SieveJoint.h"
  #undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
#endif // SIEVE_GAUSS_MULTI_THREADED

#include "Sampler.h"
#include "TerminationConditions.h"
#include "DefaultTermConds.h"

namespace GaussSieve
{
template<class SieveTraits, bool MT> class Sieve;

template<class CoefficientType>
using TupleSieve = Sieve<DefaultSieveTraits<CoefficientType, SIEVE_GAUSS_DEFAULT_THREADED, -1>, SIEVE_GAUSS_DEFAULT_THREADED>;

#ifdef SIEVE_GAUSS_SINGLE_THREADED
template<class CoefficientType>
using SieveST = Sieve<DefaultSieveTraits<CoefficientType, false, -1> , false>;
#endif
#ifdef SIEVE_GAUSS_MULTI_THREADED
template<class SieveTraits>
using SieveMT = Sieve<DefaultSieveTraits,CoefficientType, true, -1>, true>;
#endif

}  // end namespace GaussSieve


// #if 0 decativated in Sieve.h, activated in the others
// (We still indent)

 /**
    Single-threaded-only implementation files
 */
  #ifdef SIEVE_GAUSS_SINGLE_THREADED
    #define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED false
    #include "SieveJoint_impl.h"
    #include "SieveST_impl.h"
    #include "SieveST2_impl.h"
    #include "SieveST3_impl.h"
    #undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
  #endif // SIEVE_GAUSS_SINGLE_THREADED

  /**
    Multi-threaded-only implementation files
  */

  #ifdef SIEVE_GAUSS_MULTI_THREADED
    #define GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED true
    #include "SieveJoint_impl.h"
    #include "SieveMT.cpp"
    #undef GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED
  #endif // SIEVE_GAUSS_MULTI_THREADED

  /**
    Implementation files for either case:
  */

  #include "DefaultTermConds_impl.h"
  #include "GaussQueue_impl.h"
  #include "GPVSampler_impl.h"
  #include "Sampler_impl.h"
  #include "UniformSampler_impl.h"
  #include "GPVSamplerExtended_impl.h"

// #endif // of if 0 block for the Sieve.h file

// Sieve.cpp has explicit instantiation, Sieve.h has extern template.
// We have nothing: We implicitly instatiatiate exactly what we need.

#if 0
namespace GaussSieve
{

#ifdef SIEVE_GAUSS_SINGLE_THREADED
  extern template class Sieve<DefaultSieveTraits<long,      false, -1>, false>;
// extern template class Sieve<DefaultSieveTraits<double,    false, -1>, false>;
  extern template class Sieve<DefaultSieveTraits<mpz_class, false, -1>, false>;


// testing only
  extern template class Sieve<DefaultSieveTraits<int32_t, false, -1, fplll::ZZ_mat<mpz_t>>,false>;



  extern template class Sampler<DefaultSieveTraits<long,      false, -1>, false>;
//  extern template class Sampler<DefaultSieveTraits<double,    false, -1>, false>;
  extern template class Sampler<DefaultSieveTraits<mpz_class, false, -1>, false>;

  extern template class TerminationCondition<DefaultSieveTraits<long,      false, -1>, false>;
//  extern template class TerminationCondition<DefaultSieveTraits<double,    false, -1>, false>;
  extern template class TerminationCondition<DefaultSieveTraits<mpz_class, false, -1>, false>;
#endif

#ifdef SIEVE_GAUSS_MULTI_THREADED
  extern template class Sieve<DefaultSieveTraits<long,      true, -1>, false>;
//  extern template class Sieve<DefaultSieveTraits<double,    true, -1>, false>;
  extern template class Sieve<DefaultSieveTraits<mpz_class, true, -1>, false>;

  extern template class Sampler<DefaultSieveTraits<long,      true, -1>, true>;
//  extern template class Sampler<DefaultSieveTraits<double,    true, -1>, true>;
  extern template class Sampler<DefaultSieveTraits<mpz_class, true, -1>, true>;

  extern template class TerminationCondition<DefaultSieveTraits<long,      true, -1>, true>;
//  extern template class TerminationCondition<DefaultSieveTraits<double,    true, -1>, true>;
  extern template class TerminationCondition<DefaultSieveTraits<mpz_class, true, -1>, true>;
#endif

}  // end namespace GaussSieve
#endif


#endif  // main include guard for Sieve.h
