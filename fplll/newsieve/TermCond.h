#ifndef TERM_COND_H
#define TERM_COND_H

#error "no longer used"

/* Concept "TerminationCondition":

A termination condition is a class template
template<class ET>
class TerminationCondition{
template<MT>
init(GaussSieve<ET,MT>
}

that implements (at least) the following methods:





*/


template<class ET>
class TerminationConditions;

class bad_dumpread_TermCond:public std::runtime_error
{
   public: bad_dumpread_TermCond():runtime_error("Dump read failed for Termination Condition") {}
}; //exception indicating that read from dump failed.

template<class ET> ostream & operator<<(ostream &os,TerminationConditions<ET> const &term_cond); //printing
template<class ET> istream & operator>>(istream &is,TerminationConditions<ET> &term_cond); //reading (also used by constructor from istream)
template<class ET>
class TerminationConditions
{
    friend ostream & operator<< <ET>(ostream &os,TerminationConditions const &term_cond); //printing
    friend istream & operator>> <ET>(istream &is,TerminationConditions &term_cond);
public:
    TerminationConditions() : default_condition(true) {}; //if this is set, we are to ignore all other values.
    explicit TerminationConditions(istream &is):default_condition(true)
    {
        is>> (*this);
    };
    TerminationConditions(TerminationConditions const &old)=default;
    TerminationConditions(TerminationConditions && old)=default;
    TerminationConditions & operator=(TerminationConditions const &old)=default;
    TerminationConditions & operator=(TerminationConditions && old)=default;
    ~TerminationConditions()= default;
    enum class TerminationType //different basic types of termination check
    {
        CheckCollisions=1,
        CheckLength=2,
        CheckListSize=3
    };
    bool do_we_check_collisions() const                     {return do_we_check_collisions_;};
    bool do_we_check_length() const                         {return do_we_check_length_;};
    bool do_we_check_list_size() const                      {return do_we_check_list_size_;};
    bool do_we_use_default_condition() const                {return default_condition;};
    unsigned long int get_allowed_collisions() const        {return allowed_collisions_;};
    ET get_target_length() const                            {return target_length_;};
    unsigned long int get_allowed_list_size() const         {return allowed_list_size_;};

    void set_allowed_collisions(unsigned long const &colls)
    {   allowed_collisions_=colls;
        do_we_check_collisions_=true;
        default_condition = false;
        return;
    };

    void set_allowed_list_size(unsigned long const &maxsize)
    {
        allowed_list_size_=maxsize;
        do_we_check_list_size_=true;
        default_condition = false;
        return;
    };

    void set_target_length(ET const &new_target_length)
    {
        target_length_=new_target_length;
        do_we_check_length_=true;
        default_condition = false;
        return;
    };

private:
    bool do_we_check_collisions_;
    bool do_we_check_length_;
    bool do_we_check_list_size_;
    bool default_condition;
    unsigned long int allowed_collisions_;
    unsigned long int allowed_list_size_;
    ET target_length_;
}; //end of termination condition class

template<class ET> ostream & operator<<(ostream &os,TerminationConditions<ET> const &term_cond) //printing
{

    os << "Default_Conditions=" << term_cond.default_condition << endl;
    if(!term_cond.default_condition)
    {
        os << "Check Collisions=" << term_cond.do_we_check_collisions_ << endl;
        if(term_cond.do_we_check_collisions_)
        {
            os << "Number=" << term_cond.allowed_collisions_ << endl;
        }
        os << "Check List Size=" << term_cond.do_we_check_list_size_ << endl;
        if(term_cond.do_we_check_list_size_)
        {
            os << "Number=" << term_cond.allowed_list_size_ << endl;
        }
        os << "Check Target Length=" << term_cond.do_we_check_length_ << endl;
        if(term_cond.do_we_check_length_)
        {
            os << "Target Length=" << term_cond.target_length_ << endl;
        }
    }
    return os;
}
template<class ET> istream & operator>>(istream &is,TerminationConditions<ET> &term_cond)
{
    bool do_we_check_collisions_;
    bool do_we_check_length_;
    bool do_we_check_list_size_;
    bool default_condition;
    unsigned long int allowed_collisions_;
    unsigned long int allowed_list_size_;
    ET target_length_;
//We should probably throw an exception rather than return is.
    if (!GaussSieve::string_consume(is,"Default_Conditions=")) throw bad_dumpread_TermCond();
    is >> term_cond.default_condition;
    if(!term_cond.default_condition)
    {
        if(!GaussSieve::string_consume(is,"Check Collisions=")) throw bad_dumpread_TermCond();
        is>> term_cond.do_we_check_collisions_;
        if(term_cond.do_we_check_collisions_)
        {
            if(!GaussSieve::string_consume(is,"Number=")) throw bad_dumpread_TermCond();
            is >> term_cond.allowed_collisions;
        }
        if(!GaussSieve::string_consume(is,"Check List Size=")) throw bad_dumpread_TermCond();
        is>> term_cond.do_we_check_list_size_;
        if(term_cond.do_we_check_list_size_)
        {
            if(!GaussSieve::string_consume(is,"Number=")) throw bad_dumpread_TermCond();
            is>> term_cond.allows_list_size_;
        }
        if(!GaussSieve::string_consume(is,"Check Target Length=")) throw bad_dumpread_TermCond();
        is>>term_cond.do_we_check_length_;
        if(term_cond.do_we_check_length_)
        {

        }

    }
    return is;

} //reading (also used by constructor from istream)


#endif
