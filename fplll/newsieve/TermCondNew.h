#ifndef TERM_COND_NEW_H
#define TERM_COND_NEW_H

//TODO : Change to this version later.
//Considerations :
//pure virtual function version breaks dumping / reading (without some RTT)
//templated version locks into static typing. However, type of term-cond might depend on user input.

//solution : provide limited form of RTT. Require user to supply a (possibly unitialized) pointer to correct derived class if reading from dump.


//base class for Termination Condition.
//Any valid type of (possibly user-defined) termination condition needs to be a template derived from this class, i.e.
//
//template<class ET, bool MT>
//class MyTerminationCondition : public TerminationCondition<ET,MT>
//{
// ... additional stuff
//}
//The derived class *needs* to overwrite check, check_vec and write a custom destructor.
//The derived class may want to overwrite dump_to_stream, read_from_stream, is_simple, init.
//termination_condition_type must not be overwritten for user-defined classes, unless one makes according changes to the dumping / reading routines of the Sieve class.
//
//semantics:
//init(sieve) is called if the GaussSieve is started. Note that init() may be called multiple times if the sieve is suspended and parameters change. You may assume that the main parameters do not change without init being called again.
//check(sieve) returns 1 if sieve is considered finished, 0 otherwise. Needs to be reentrant if MT==true.
//check_vec(sieve,norm2) has the same semantics as check, but is called whenever a new lattice point of norm^2 norm2 is found. Needs to be reentrant if MT==true.
//
//is_simple needs to be either {return true;} or {return false;} (i.e. essentially a class-dependant constant). If it returns true, we assume that check() is a function of length only. In this case, we only call check_vec() whenever a new shortest (so far) vector is found.
//Note: Even if is_simple returns true, check() may still be called rarely.
//
//Note : Output of check functions is int rather than bool to enable future extensions to output return values indicating "suspend" or "dump"
//
//The stream operators are for dumping / reading. To enable polymorphism, these are wrappers to dump_to_stream resp. read_from_stream.
//We assume that if we dump a Termination condition T1 to a filestream and read it into T2, then T2 will be in the same status as T1.
//You may need to overwrite dump_to_stream / read_from_stream to enable that.
//Reading in needs to work on already initialised data (i.e. has to overwrite previous values).
//After reading in, we are guaranteed that init() is called before we do any check() or check_vec().
//For custom termination conditions making use of statistics, assume that check() is only called rarely and that statistics collected might be lagging behind (due to inter-thread synchronisation buffers).
//In particular, only ever compare >= or <= and never ==.

template<class ET, bool MT, int nfixed=-1> class TerminationCondition;
template<class ET, bool MT, int nfixed=-1> class NeverTerminationCondition;
template<class ET, bool MT, int nfixed=-1> class LengthTerminationCondition;
template<class ET, bool MT, int nfixed=-1> class MinkowskiTerminationCondition;
template<class ET, bool MT, int nfixed=-1> ostream & operator<<(ostream &os,TerminationCondition<ET,MT,nfixed>* const term_cond); //printing
template<class ET, bool MT, int nfixed=-1> istream & operator>>(istream &is,TerminationCondition<ET,MT,nfixed>* const term_cond); //reading (also used by constructor from istream)

enum class TerminationConditionType
{
    user_defined = 0,
    never_terminate = 1,
    length_condition = 2,
    minkowski_condition =3
};


#include "SieveGauss.h"
#include "Utility.cpp"


template<class ET,bool MT, int nfixed> class TerminationCondition
{
    public:
    friend ostream & operator<< <ET,MT,nfixed>(ostream &os,TerminationCondition<ET,MT,nfixed>* const term_cond);
    friend istream & operator>> <ET,MT,nfixed>(istream &is,TerminationCondition<ET,MT,nfixed>* const term_cond);
    virtual void init(Sieve<ET,MT,nfixed> * const sieve) {};     //TODO: Fix const-correctness. Problem is with cbegin() from main_list, really...
    virtual int check(Sieve<ET,MT,nfixed> * const sieve) = 0;
    virtual int check_vec(Sieve<ET,MT,nfixed> * const sieve, ET const & length2) = 0;
    virtual ~TerminationCondition()=0; //needs to be virtual
    virtual bool is_simple() const {return false;};
    virtual TerminationConditionType  termination_condition_type() const {return TerminationConditionType::user_defined;};    //run-time type information.
                                                                    //This is used to determine how to interpret a dump file.
                                                                    //0 means user-defined (which is the default).
                                                                    //Other values mean that the GaussSieve dumping routine is aware of the type, simplifying the syntax for dumping / reading.
    //TODO : Write some explanation how to do that.

    private:
    virtual ostream & dump_to_stream(ostream &os)  {return os;};    //implementation of << operator.
    virtual istream & read_from_stream(istream &is){return is;};    //implementation of >> operator.

};

template <class ET,bool MT,int nfixed> TerminationCondition<ET,MT,nfixed>::~TerminationCondition() {} //actually needed, even though destructor is pure virtual as the base class destructor is eventually called implicitly.

template<class ET,bool MT,int nfixed> ostream & operator<<(ostream &os,TerminationCondition<ET,MT,nfixed>* const term_cond){return term_cond->dump_to_stream(os);};
template<class ET,bool MT,int nfixed> istream & operator>>(istream &is,TerminationCondition<ET,MT,nfixed>* const term_cond){return term_cond->read_from_stream(is);};



//default Termination Conditions here:
template<class ET,bool MT, int nfixed> class NeverTerminationCondition : public TerminationCondition<ET,MT,nfixed> //never terminate
{
    public:
    virtual int check(Sieve<ET,MT,nfixed> * const sieve) override                              {return 0;};
    virtual int check_vec(Sieve<ET,MT,nfixed> * const sieve, ET const & length2) override      {return 0;};
    virtual bool is_simple() const override                                             {return true;};
    virtual ~NeverTerminationCondition() {};
    virtual TerminationConditionType  termination_condition_type() const override       {return TerminationConditionType::never_terminate;};    //run-time type information.
};

template<class ET,bool MT,int nfixed> class LengthTerminationCondition : public TerminationCondition<ET,MT,nfixed> //Length Termination Condition
{
    public:
    LengthTerminationCondition(ET const & init_target_length) : target_length(init_target_length) {};
    virtual int check(Sieve<ET,MT,nfixed> * const sieve) override                               {return (sieve -> get_best_length2()<=target_length)?1:0;};
    virtual int check_vec(Sieve<ET,MT,nfixed> * const sieve, ET const & length2) override       {return (length2<=target_length)?1:0;};
    virtual bool is_simple() const override                                                     {return true;};
    virtual TerminationConditionType  termination_condition_type() const override               {return TerminationConditionType::length_condition;};    //run-time type information.
    virtual ~LengthTerminationCondition() {};
    class bad_dumpread_LengthTermCond:public std::runtime_error
    {
        public: bad_dumpread_LengthTermCond():runtime_error("Dump read failed for LengthTerminationCondition") {}
    }; //exception indicating that read from dump failed.
    virtual ostream & dump_to_stream(ostream &os) override                              {os << "Target Length=" << target_length << endl; return os;}
    virtual istream & read_from_stream(istream &is) override
    {
        if(!GaussSieve::string_consume(is,"Target Length=")) throw bad_dumpread_LengthTermCond();
        is >> target_length;
        return is;
    };
    private:
    ET target_length;
};

template<class ET,bool MT, int nfixed> class MinkowskiTerminationCondition : public TerminationCondition<ET,MT,nfixed> //Length Termination Condition
{
    public:
    MinkowskiTerminationCondition() : target_length() {}; //unitialised. We are guaranteed that init() is run before use.
    virtual int check(Sieve<ET,MT,nfixed>  * const sieve) override                          {return (sieve -> get_best_length2()<=target_length)?1:0;};
    virtual int check_vec(Sieve<ET,MT,nfixed>  * const sieve, ET const & length2) override  {return (length2<=target_length)?1:0;};
    virtual bool is_simple() const override                                                 {return true;};
    virtual TerminationConditionType  termination_condition_type() const override           {return TerminationConditionType::minkowski_condition;};    //run-time type information.
    virtual ~MinkowskiTerminationCondition() {};
    virtual void init(Sieve<ET,MT,nfixed> * const sieve) override                           {target_length = GaussSieve::compute_mink_bound(sieve->get_original_basis());};
    private:
    ET target_length;
};


#endif


//Note : Igonore idea below, because we want run-time typing.
/* Concept "TerminationCondition":

A termination condition is a template<class ET> class TerminationCondition with at least the following public members: (see DummyTerminationCondition for an example)

- public member typedef simple, set to either std::true_type or std::false_type (Note : May need to #include <type_traits> )
- public default constructor
- public member function templates
template<bool MT> void init(GaussSieve<ET,MT> const & sieve);
template<bool MT> int check(GaussSieve<ET,MT> const & sieve);
template<bool MT> int check_vec(GaussSieve<ET,MT> const & sieve, ET const & norm2);

- non-member function templates
template<class ET> ostream & operator<<(ostream &os,TerminationCondition<ET> const &term_cond);
template<class ET> istream & operator>>(istream &is,TerminationCondition<ET> &term_cond);

It needs to have the following semantics:

init(...) is called if the GaussSieve is started. You may assume that the parameters do no longer change. Note that init() may be called multiple times if the sieve is suspended and parameters change.
check(...) returns 1 if sieve is considered finished, 0 otherwise. Needs to be reentrant if MT==true.
check_vec(...) has the same semantics as check, but is called whenever a new lattice point of norm^2 norm2 is found. Needs to be reentrant if MT==true.

If simple==true_type, we assume that check() is a function of length only. In this case, we only call check_vec() if a new shortest (so far) vector is found.
check() may be called rarely.

Note : Output of check functions is int rather than bool to enable future extensions to output return values indicating "suspend" or "dump"

The stream operators are for dumping / reading.
We assume that if we dump a Termination condition T1 to a filestream and read it into T2, then T2 will be in the same status as T1.
TODO: reading is unused so far.
*/

/*
    example dummy Termination condition below
*/

//template<class ET> ostream & operator<<(ostream &os,DummyTerminationCondition<ET> const &term_cond); //printing
//template<class ET> istream & operator>>(istream &is,DummyTerminationCondition<ET> &term_cond); //reading (also used by constructor from istream)
//template<class ET>
//class DummyTerminationCondition //This class serves only to exemplify the interface, lacking support for "Concepts" in any compiler execpt GCC (experimental) atm.
//{
//    public:
//    template<class ET> ostream & operator<<(ostream &os,TerminationConditions<ET> const &term_cond); //printing
//    template<class ET> istream & operator>>(istream &is,TerminationConditions<ET> &term_cond); //reading (also used by constructor from istream)
//    using simple = std::true_type;
//
//    DummyTerminationCondition() = default;
//
//    template<bool MT>
//    void init(GaussSieve<ET,MT> const & sieve) = delete;
//
//    template<bool MT>
//    int check(GaussSieve<ET,MT> const & sieve) = delete;
//
//    template<bool MT>
//    int check_vec(GaussSieve<ET,MT> const & sieve, ET const & length) = delete;
//}
