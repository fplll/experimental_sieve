/**
  This file collects and specifies all debugging preprocessor macros that
  enable checks.
*/

#ifndef DEBUG_ALL_H
#define DEBUG_ALL_H

// clang-format kills the indenting for #defines
// clang-format off

// enables all silent checks
#define DEBUG_SIEVE_SILENT_ALL

#ifdef DEBUG_SIEVE_SILENT_ALL

// tests that LP classes were initialized prior to using.
  #define DEBUG_SIEVE_LP_INIT

// tests that operations on LP have matching dimensions
  #define DEBUG_SIEVE_LP_MATCHDIM

// tests that the RNG engine is called with valid parameters
  #define DEBUG_SIEVE_MTPRNG

// Tests that accesses to the mu and g-matrices only access entries that are potentially non-zero.
// While code may reasonably access such entries that are always zero, we give an error, because
// it's more likely to be a bug rather than intended.
  #define DEBUG_SIEVE_LOWERTRIANGULAR_MUG

#endif  // DEBUG_SIEVE_SILENT_ALL


// Outputs late-initializations of modules from GaussSieve
#define DEBUG_SIEVE_INITIALIZATIONS

#ifdef DEBUG_SIEVE_INITIALIZATIONS
  #define DEBUG_SIEVE_TRACEINITIATLIZATIONS(STR) std::cout << STR << std::endl << std::flush;
#else
  #define DEBUG_SIEVE_TRACEINITIATLIZATIONS(STR)
#endif


// verbosely trace calls of generic functions for LPs.
// Disabled by default:

// #define DEBUG_SIEVE_LP_TRACEGENERIC
#ifdef DEBUG_SIEVE_LP_TRACEGENERIC
  #define DEBUG_TRACEGENERIC(STR) std::cout << STR << std::endl << std::flush;
#else  // DEBUG_SIEVE_LP_TRACEGENERIC
  #define DEBUG_TRACEGENERIC(STR)
#endif  // DEBUG_SIEVE_LP_TRACEGENERIC

// clang-format on

#endif
