#ifndef DEFAULT_TERM_CONDS_H
#define DEFAULT_TERM_CONDS_H

#include "TermCondNew.h" //consider merging with that file
#include <iostream>

namespace GaussSieve{
inline    fplll::Z_NR<mpz_t> compute_mink_bound(fplll::ZZ_mat<mpz_t> const & basis);                                        //computes a meaningful Minkowski bound for the length of the shortest vector
}

namespace GaussSieve{

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
    virtual inline int check(Sieve<ET,MT,nfixed> * const sieve) override;
    virtual int check_vec(Sieve<ET,MT,nfixed> * const sieve, ET const & length2) override       {return (length2<=target_length)?1:0;};
    virtual bool is_simple() const override                                                     {return true;};
    virtual TerminationConditionType  termination_condition_type() const override               {return TerminationConditionType::length_condition;};    //run-time type information.
    virtual ~LengthTerminationCondition() {};
    class bad_dumpread_LengthTermCond:public std::runtime_error
    {
        public: bad_dumpread_LengthTermCond():runtime_error("Dump read failed for LengthTerminationCondition") {}
    }; //exception indicating that read from dump failed.
    virtual std::ostream & dump_to_stream(std::ostream &os) override                              {os << "Target Length=" << target_length << std::endl; return os;}
    virtual std::istream & read_from_stream(std::istream &is) override
    {
        if(string_consume(is,"Target Length=")) throw bad_dumpread_LengthTermCond();
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
    virtual inline int check(Sieve<ET,MT,nfixed>  * const sieve) override; //{return (sieve -> get_best_length2()<=target_length)?1:0;};
    virtual int check_vec(Sieve<ET,MT,nfixed>  * const sieve, ET const & length2) override  {return (length2<=target_length)?1:0;};
    virtual bool is_simple() const override                                                 {return true;};
    virtual TerminationConditionType  termination_condition_type() const override           {return TerminationConditionType::minkowski_condition;};    //run-time type information.
    virtual ~MinkowskiTerminationCondition() {};
    inline virtual void init(Sieve<ET,MT,nfixed> * const sieve) override; //                           {target_length = GaussSieve::compute_mink_bound(sieve->get_original_basis());};
    private:
    ET target_length;
};

}

#endif
