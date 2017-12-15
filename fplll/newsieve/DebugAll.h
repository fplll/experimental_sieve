/**
  This file collects and specifies all debugging preprocessor macros that
  enable checks / verbose debug output.
*/

#ifndef DEBUG_ALL_H
#define DEBUG_ALL_H

/**
  We recognize the following debug symbols:
  DEBUG_SIEVE_SILENT_ALL
      Turns on all other debug symbols that are silent (i.e. not not normally produce output) and do
      not change the behaviour. Turns on some asserts / static_asserts / exceptions
  DEBUG_SIEVE_LP_INIT
      Turns on checks that Lattice points were initialized prior to using
  DEBUG_OTHER_STATIC_INIT
      Turns on checks that non-lattice point objects were initialized with a static initializer
      prior to using.
      TODO: Not used yet (but it should)
  DEBUG_SIEVE_LP_MATCHDIM
      Turns on checks that (binary) operations on lattice points have matching dimensions
  DEBUG_SIEVE_MTPRNG_THREAD_RANGE
      Tests that the RNG engine is called with a valid thread-number parameter
  DEBUG_SIEVE_LOWERTRIANGULAR_MUG
      Tests that accesses to the mu- and g-matrices stored in the GSO only access entries that are
      potentially non-zero. (While code may reasonably access other entries, we give an error,
      because it's likely a bug and we do not guarantee that the GSO object works for those)

  DEBUG_SIEVE_ALL_APPROX_AT_SAME_LEVEL    [Lazy evaluation only, unused]
      Ensures that operations on leveled objects do not lose levels.
      ( E.g. assume the result on X+Y if a lazy evaluation object that first computes the result
        approximately (e.g. if X+Y is used in a comparison), then exactly.
        Then some operation (X+Y)*Z might not be able to produce another lattice evaluation object,
        because (X+Y)*Z is not well-defined for the approximations to X/Y/Z, but is well-defined
        for the exact X/Y/Z. Without this debug symbol, the approximations are silently dropped.
        With it set, we give an error. )
      Used to trace such missing implementations.

  DEBUG_SIEVE_STANDALONE_MODULES_ALL
      Turns on all other DEBUG_SIEVE_STANDALONE_* symbols.
      These enable modules which use a (back-) pointer to Sieve<...> to be used standalone without
      an associated sieve. Note that this changes semantics and disables some features.
      *Only use in unit testing!*
  DEBUG_SIEVE_STANDALONE_SAMPLER
      For a standalone sampler
  DEBUG_SIEVE_STANDALONE_QUEUE
      For a standalone queue

  DEBUG_SIEVE_LAZY_ALL
      Activates other (very verbose) debug symbols for lazy evaluations
      Lazy evaluation is currently unused and these are subject to change
      Only use if you are developing the lazy evaluations module
  DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS  [Lazy evaluations only, unused
      Logs constructions of lazy objects to stdout.
  DEBUG_SIEVE_LAZY_TRACE_EVALS
      Logs evaluations of lazy objects to stdout.

  DEBUG_SIEVE_SHOW_INITIALIZATIONS
      Logs (global) initializations to stdout. For debugging StaticInitializers, but also others.
      Has relatively moderate output size.
      -> Determines whether DEBUG_SIEVE_TRACEINITIATLIZATIONS(Output string) is a no-op

  DEBUG_SIEVE_LP_TRACEGENERIC
      Logs usage of generic implementations of lattice points.
      Extremely verbose output.
      -> Determines whether DEBUG_TRACEGENERIC(output string) is a no-op

TODO: Some of these need revision / usage clarification.
Their use is not very consistent.
*/

//#define DEBUG_SIEVE_LAZY_ALL

// clang-format kills the indentation for #defines, which does not help here.
// clang-format off


#ifdef DEBUG_SIEVE_SILENT_ALL
  #define DEBUG_SIEVE_LP_INIT
  #define DEBUG_SIEVE_LP_MATCHDIM
  #define DEBUG_SIEVE_MTPRNG_THREAD_RANGE
  #define DEBUG_SIEVE_LOWERTRIANGULAR_MUG
  #define DEBUG_OTHER_STATIC_INIT
//  #define DEBUG_SIEVE_ALL_APPROX_AT_SAME_LEVEL
#endif  // DEBUG_SIEVE_SILENT_ALL

#ifdef DEBUG_SIEVE_STANDALONE_MODULES_ALL
  #define DEBUG_SIEVE_STANDALONE_SAMPLER
  #define DEBUG_SIEVE_STANDALONE_QUEUE
#endif

#ifdef DEBUG_SIEVE_LAZY_ALL
  #define DEBUG_SIEVE_LAZY_TRACE_CONSTRUCTIONS
  #define DEBUG_SIEVE_LAZY_TRACE_EVALS
#endif

// Outputs late-initializations of modules from GaussSieve
//#define DEBUG_SIEVE_INITIALIZATIONS

#ifdef DEBUG_SIEVE_SHOW_INITIALIZATIONS
  #define DEBUG_SIEVE_TRACEINITIATLIZATIONS(STR) std::cout << STR << std::endl << std::flush;
#else
  #define DEBUG_SIEVE_TRACEINITIATLIZATIONS(STR)
#endif

#ifdef DEBUG_SIEVE_LP_TRACEGENERIC
  #define (STR) std::cout << STR << std::endl << std::flush;
#else  // DEBUG_SIEVE_LP_TRACEGENERIC
  #define DEBUG_TRACEGENERIC(STR)
#endif  // DEBUG_SIEVE_LP_TRACEGENERIC

// clang-format on

#endif
