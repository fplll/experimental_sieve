#ifndef TERM_COND_NEW_H
#define TERM_COND_NEW_H

/**
  This file describes the interface we use for termination conditions that tell the sieve algorithm
  when to stop.
  Since we want to support arbitrary user-provided termination conditions, we implement this via
  dynamic (runtime) polymorphism.
  the main sieve stores a pointer to a TerminationCondition<SieveTraits, MT> (or, rather a class
  derived from it).

  The termination condition itself is initialized with a pointer to the associated sieve, which
  allows to query the status of the sieve.

  Similar considerations as for the sampler apply. Notably, the actual termination condition object
  (if user-provided) may be constructed prior to the construction of the sieve and outlive it. For
  that reason, we perfom most initialization not inside the constructor but in a separate
  init(Sieve *ptr_to_sieve) routine that is called when the sieve starts.
*/

/**

  base class for Termination Condition.
  Any valid type of (possibly user-defined) termination condition needs to be a (probably template)
  class derived from this

  template<class SieveTraits, bool MT>
  class MyTerminationCondition : public TerminationCondition<ET,SieveTraits>
  {
  ... additional stuff
  }
  The derived class *needs* to overwrite check() and write a custom destructor (which may be empty).

  termination_condition_type must not be overwritten for user-defined classes, unless one makes
  according changes to the dumping / reading routines of the Sieve class.

  semantics:
  init(sieveptr) is called if the GaussSieve is started, which calls custom_init().
  custom_init() is supposed to be overloaded by the derived class to perform its custom
  initialization.
  Since we may eventually support dumping / reading, we do not guarantee that init() is only called
  once.

  check() returns 1 to indicate that sieve is considered finished, 0 otherwise. It is called
  regularly from the sieve. This needs to be thread-safe if derived from
  TerminationCondition<SieveTraits,true>.

  check_vec(TermCond_QueryType const &lattice_point) has the same semantics as check, but is called
  only whenever a new lattice point is point. The data type used for the lattice point is
  SieveTraits::TermCond_QueryType.

  The return type is an int rather than a bool to support extending to other return values
  that may indicate something to the sieve.

  Since a termination condition might depend on, say, the number of iterations / list size, we call
  check() regularly. To possibly limit that, you can overload dependency_type().
  This must returns a class-specific constant to indicate to the sieve when it should be called
  (e.g. if it returns TermCondDependencyType::new_shortest_vector, we might not call check() or
  check_vec() unless we find a new shortest vector).
  *** The sieve is free to ignore this information, so you need to implement check() ***

  NOTE: For custom termination conditions making use of statistics, assume that check() is only
  called rarely and that statistics collected might be lagging behind (due to e.g. potential
  synchronization effects in multithreading)
  In particular, only ever compare >= or <= and never ==.
*/

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"

namespace GaussSieve{

// forward declarations

template <class SieveTraits, bool MT> class TerminationCondition;
template <class SieveTraits, bool MT> class NeverTerminationCondition;
template <class SieveTraits, bool MT> class LengthTerminationCondition;
template <class SieveTraits, bool MT> class MinkowskiTerminationCondition;
template <class SieveTraits, bool MT> std::ostream & operator<<(std::ostream &os,TerminationCondition<SieveTraits,MT>* const term_cond); //printing
template <class SieveTraits, bool MT> std::istream & operator>>(std::istream &is,TerminationCondition<SieveTraits,MT>* const term_cond); //reading

// types of Termination Condition recognized. Note that we may have user-provided Termination
// conditions.
enum class TerminationConditionType
{
  user_defined        = 0,
  never_terminate     = 1,
  length_condition    = 2,
  minkowski_condition = 3
};

enum class TermCondDependencyType
{
  any,
  new_vector_only,
  new_shortest_vector
};

template <class SieveTraits, bool MT> class Sieve;

template <class SieveTraits, bool MT> class TerminationCondition
{

private:  // shorthand typedefs to avoid "typename ..."
  using TermCond_QueryType = typename SieveTraits::TermCond_QueryType;

public:

  // TODO: Remove pointer
  friend std::ostream & operator<< <SieveTraits,MT>(std::ostream &os,TerminationCondition<SieveTraits,MT> *const term_cond);
  friend std::istream & operator>> <SieveTraits,MT>(std::istream &is,TerminationCondition<SieveTraits,MT> *const term_cond);

  // associates the termination condition with the calling sieve
  void init(Sieve<SieveTraits,MT> * const ptr_to_caller)
  {
    sieveptr = ptr_to_caller;
    this->custom_init();  // dispatch to virtual custom_init, which is probably overloaded.
  }

  virtual void custom_init() { }  // defaults to "do nothing"

  // check() is called regularly from the sieve to query whether it should terminate.
  virtual int check() = 0;

  // variant that is called whenever a new lattice point is found.
  virtual int check_vec(TermCond_QueryType const &lattice_point) { return check(); }
  virtual ~TerminationCondition()=0;  // needs to be virtual

  // used to possibly optimize the frequency with which check() / check_vec() is called.
  virtual TermCondDependencyType dependency_type() const { return TermCondDependencyType::any; }

  // run-time type information. This may be needed / helpful to implement duming / reading back
  // dumps to let the sieve know the type of termination condition.
  virtual TerminationConditionType termination_condition_type() const
  {
    return TerminationConditionType::user_defined;
  }

private:
    virtual std::ostream & dump_to_stream(std::ostream &os)   {return os; }  //implementation of << operator.
    virtual std::istream & read_from_stream(std::istream &is) {return is; }  //implementation of >> operator.

    protected:
    Sieve<SieveTraits, MT> *sieveptr;  // holds an (observing) pointer to the caller sieve or
                                       // nullptr if unassociated
};

template <class SieveTraits,bool MT> TerminationCondition<SieveTraits,MT>::~TerminationCondition() {} //actually needed, even though destructor is pure virtual as the base class destructor is eventually called implicitly.

template<class SieveTraits,bool MT> std::ostream & operator<<(std::ostream &os,TerminationCondition<SieveTraits,MT>* const term_cond){return term_cond->dump_to_stream(os);};
template<class SieveTraits,bool MT> std::istream & operator>>(std::istream &is,TerminationCondition<SieveTraits,MT>* const term_cond){return term_cond->read_from_stream(is);};

}  // end namespace GaussSieve

#endif  // include guards
// clang-format on
