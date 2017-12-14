#ifndef TERM_COND_NEW_H
#define TERM_COND_NEW_H

// TODO : Change to this version later.
// Considerations :
// pure virtual function version breaks dumping / reading (without some RTT)
// templated version locks into static typing.
// However, type of term-cond might depend on user input.

// solution : provide limited form of RTT. Require user to supply a (possibly unitialized) pointer to correct derived class if reading from dump.


// base class for Termination Condition.
// Any valid type of (possibly user-defined) termination condition needs to be a template derived from this class, i.e.
//
// template<class ET, bool MT>
// class MyTerminationCondition : public TerminationCondition<ET,MT>
//{
// ... additional stuff
//}
// The derived class *needs* to overwrite check, check_vec and write a custom destructor.
// The derived class may want to overwrite dump_to_stream, read_from_stream, is_simple, init.
// termination_condition_type must not be overwritten for user-defined classes, unless one makes according changes to the dumping / reading routines of the Sieve class.
//
// semantics:
// init(sieve) is called if the GaussSieve is started. Note that init() may be called multiple times if the sieve is suspended and parameters change. You may assume that the main parameters do not change without init being called again.
// check(sieve) returns 1 if sieve is considered finished, 0 otherwise. Needs to be thread-safe if MT==true.
// check_vec(sieve,norm2) has the same semantics as check, but is called whenever a new lattice point of norm^2 norm2 is found. Needs to be thread-safe if MT==true.
//
// is_simple needs to be either {return true;} or {return false;} (i.e. essentially a class-dependant constant). If it returns true, we assume that check() is a function of length only. In this case, we only call check_vec() whenever a new shortest (so far) vector is found.
// Note: Even if is_simple returns true, check() may still be called rarely.
//
// Note : Output of check functions is int rather than bool to enable future extensions to output return values indicating "suspend" or "dump"
//
// The stream operators are for dumping / reading. To enable polymorphism, these are wrappers to dump_to_stream resp. read_from_stream.
// We assume that if we dump a Termination condition T1 to a filestream and read it into T2, then T2 will be in the same status as T1.
// You may need to overwrite dump_to_stream / read_from_stream to enable that.
// Reading in needs to work on already initialised data (i.e. has to overwrite previous values).
// After reading in, we are guaranteed that init() is called before we do any check() or check_vec().
// For custom termination conditions making use of statistics, assume that check() is only called rarely and that statistics collected might be lagging behind (due to inter-thread synchronisation buffers).
// In particular, only ever compare >= or <= and never ==.

#include "DefaultIncludes.h"

//#include "SieveGauss.h"
#include "SieveUtility.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"

namespace GaussSieve{

template <class SieveTraits, bool MT> class TerminationCondition;
template <class SieveTraits, bool MT> class NeverTerminationCondition;
template <class SieveTraits, bool MT> class LengthTerminationCondition;
template <class SieveTraits, bool MT> class MinkowskiTerminationCondition;
template <class SieveTraits, bool MT> std::ostream & operator<<(std::ostream &os,TerminationCondition<SieveTraits,MT>* const term_cond); //printing
template <class SieveTraits, bool MT> std::istream & operator>>(std::istream &is,TerminationCondition<SieveTraits,MT>* const term_cond); //reading (also used by constructor from istream)

// types of Termination Condition recognized. Note that we may have user-provided Termination
// conditions.
enum class TerminationConditionType
{
  user_defined        = 0,
  never_terminate     = 1,
  length_condition    = 2,
  minkowski_condition = 3
};

template <class SieveTraits, bool MT> class Sieve;

template <class SieveTraits, bool MT> class TerminationCondition
{
    public:
    friend std::ostream & operator<< <SieveTraits,MT>(std::ostream &os,TerminationCondition<SieveTraits,MT>* const term_cond);
    friend std::istream & operator>> <SieveTraits,MT>(std::istream &is,TerminationCondition<SieveTraits,MT>* const term_cond);
    virtual void init(Sieve<SieveTraits,MT> * const sieve) {};     //TODO: Fix const-correctness. Problem is with cbegin() from main_list, really...
    virtual int check(Sieve<SieveTraits,MT> * const sieve) = 0;
    virtual int check_vec(Sieve<SieveTraits,MT> * const sieve, typename SieveTraits::LengthType const & length2) = 0;
    virtual ~TerminationCondition()=0; //needs to be virtual
    virtual bool is_simple() const {return false;};
    virtual TerminationConditionType  termination_condition_type() const {return TerminationConditionType::user_defined;};    //run-time type information.
                                                                    //This is used to determine how to interpret a dump file.
                                                                    //0 means user-defined (which is the default).
                                                                    //Other values mean that the GaussSieve dumping routine is aware of the type, simplifying the syntax for dumping / reading.
    //TODO : Write some explanation how to do that.

    private:
    virtual std::ostream & dump_to_stream(std::ostream &os)  {return os;};    //implementation of << operator.
    virtual std::istream & read_from_stream(std::istream &is){return is;};    //implementation of >> operator.

};

template <class SieveTraits,bool MT> TerminationCondition<SieveTraits,MT>::~TerminationCondition() {} //actually needed, even though destructor is pure virtual as the base class destructor is eventually called implicitly.

template<class SieveTraits,bool MT> std::ostream & operator<<(std::ostream &os,TerminationCondition<SieveTraits,MT>* const term_cond){return term_cond->dump_to_stream(os);};
template<class SieveTraits,bool MT> std::istream & operator>>(std::istream &is,TerminationCondition<SieveTraits,MT>* const term_cond){return term_cond->read_from_stream(is);};

}  // end namespace GaussSieve

#endif  // include guards
// clang-format on
