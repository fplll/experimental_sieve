/**
  This file collects and specifies all debugging preprocessor macros that
  enable checks.
*/

#ifndef DEBUG_ALL_H
#define DEBUG_ALL_H

// enables all silent checks
#define DEBUG_SIEVE_SILENT_ALL

#ifdef DEBUG_SIEVE_SILENT_ALL

// tests that LP classes were initialized prior to using.
  #define DEBUG_SIEVE_LP_INIT

// tests that operations on LP have matching dimensions
  #define DEBUG_SIEVE_LP_MATCHDIM

#endif // DEBUG_SIEVE_SILENT_ALL

// verbosely trace calls of generic functions for LPs.
// Disabled by default

// #define DEBUG_SIEVE_LP_TRACEGENERIC
#ifdef DEBUG_SIEVE_LP_TRACEGENERIC
  #define DEBUG_TRACEGENERIC( STR ) std::cout << STR << std::endl << std::flush;
#else // DEBUG_SIEVE_LP_TRACEGENERIC
  #define DEBUG_TRACEGENERIC( STR )
#endif // DEBUG_SIEVE_LP_TRACEGENERIC


#endif
