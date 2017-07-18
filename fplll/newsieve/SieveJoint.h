/*
------DO NOT INCLUDE THIS FILE MANUALLY.--------
------     USE SieveGauss.h INSTEAD     --------
*/

//SIEVE_JOINT_H_ST and SIEVE_JOINT_H_MT are separate include guards.
//They are set if we already included this file (for the ST / MT case) respectively.
//There is also an include guard SIEVE_JOINT_H, which is set _after_ we included this file at least once already.
//Use this to condition on the second pass.

#ifndef GAUSS_SIEVE_IS_MULTI_THREADED
    #error wrong usage of SieveJoint.h
#endif

#undef DO_INCLUDE_SIEVE_JOINT_H
#if GAUSS_SIEVE_IS_MULTI_THREADED == false

#ifndef SIEVE_JOINT_H_ST
#define SIEVE_JOINT_H_ST
#define DO_INCLUDE_SIEVE_JOINT_H
#endif

#elif GAUSS_SIEVE_IS_MULTI_THREADED == true
#ifndef SIEVE_JOINT_H_MT
#define SIEVE_JOINT_H_MT
#define DO_INCLUDE_SIEVE_JOINT_H
#endif
#endif

#ifdef DO_INCLUDE_SIEVE_JOINT_H

//
//end of (most) preprocessor stuff
//


#ifndef SIEVE_JOINT_H
/*

THE FOLLOWING PARTS ARE ONLY PARSED ONCE BY THE COMPILER,
EVEN IF WE USE BOTH MULTI-THREADED AND SINGLE-THREADED VARIANTS.

DECLARATIONS MAY GO HERE OR TO SieveGauss.h.

HELPER CLASSES AND FUNCTIONS WHICH ARE NOT TEMPLATES WITH TEMPLATE PARAMETER
GAUSS_SIEVE_IS_MULTI_THREADED
NEED TO GO HERE OR TO SieveGauss.h:
*/

/*
  Forward Declarations
  These go before includes to allow cyclic dependencies
*/

template<class ET, bool MultiThreaded, int nfixed> class Sieve;


/* TODO: Move these into the files where they are used... */

//template<class ET,bool MultiThreaded> class CompareQueue;

/*INCLUDES */

#include <iostream>
#include <type_traits>
#include <sys/stat.h>
#include <fstream>
#include <string>
#include <exception>
//#include "TermCond.h"
#include "GaussQueue.h"
//#include "FilteredPoint.h"
#include "TermCondNew.h"
//#include "LatticePoint.h"
#include "LatticePointsNew.h"
//#include "PointList.h"
#include "PointListNew.h"
#include "Utility.h"
#include "Typedefs.h"

#endif //end of ONLY-ONCE part

/*
EVERYTHING BELOW HERE IS POTENTIALLY INCLUDED TWICE.
TEMPLATES WITH TEMPLATE ARGUMENT GAUSS_SIEVE_IS_MULTI_THREADED
GO HERE.
*/

//The following may be included once or twice (with different values for GAUSS_SIEVE_IS_MULTI_THREADED)

//TODO: Move to where it is actually used.

//template<class ET>
//class CompareQueue<ET, GAUSS_SIEVE_IS_MULTI_THREADED>{
//    public:
//     bool operator() (const FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> & el1, const FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> & el2) const
//     {
//        return el1.get_sc_prod() > el2.get_sc_prod();  // inner products are in decreasing order
//    }
//};

#ifndef GAUSS_SIEVE_IS_MULTI_THREADED
#error Something very bad just happened
#endif


template<class ET, int nfixed> class Sieve<ET, GAUSS_SIEVE_IS_MULTI_THREADED,nfixed >
{
public:
    /*DATA TYPES*/
    //using LPType           = ExactLatticePoint<ET,nfixed>;
    using FastAccess_Point = GaussSieve::FastAccess_Point<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>;
    using MainQueueType    = GaussQueue<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>; //FIXME
    using MainListType     = GaussListNew<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>;
    using LatticeBasisType = ZZ_mat<typename ET::underlying_data_type>;
    //using SamplerType      = KleinSampler<typename ET::underlying_data_type, FP_NR<double>> *; //TODO : Should be a class with overloaded operator() or with a sample() - member.;
    //using FilteredListType = std::vector<FilteredPoint<ET, float>>; //queue is also fine for our purposes; scalar products are not of type ET, two-templates; float for now; may be changed.


//    using FilteredListType = std::list<FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2>>;

    //using BlockDivisionType  = std::array< ApproxLatticePoint<ET,GAUSS_SIEVE_IS_MULTI_THREADED> , 100> ;
//    using LengthDivisionType = std::array< LatticeApproximations::ApproxTypeNorm2, 100 >;

    //stores the last element of each block for filtered_list
//    using FilterDivisionType = std::array<typename FilteredListType::iterator, 100 >;

    //number of elements per block in filtered_list // COULD BE LONG?
//    using FilterNumOfElems = std::array<int, 100 >;



//    using AppendixType =  std::priority_queue<FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2>, std::vector<FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> >,  CompareQueue<ET,GAUSS_SIEVE_IS_MULTI_THREADED> >;



    //map where a key is a pair <length_of_list_element, inner-product>
//    using FilteredListType2 = std::map<pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2>, FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2>, CompareFilteredPoint>;
//    using FilteredListTypeP = std::map<pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2>, FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> *, CompareFilteredPoint>;
    //using FilteredListTest  =  std::map<char,int,classcomp>;
    using TermCondType     = TerminationCondition<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed> *;

public:
    /*FRIENDS */
    friend GaussQueue<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>;
    /*CONSTRUCTORS / DESTRUCTORS */
    Sieve() = delete;
    Sieve(Sieve const &old ) = delete;
    Sieve(Sieve &&old) = delete;
    Sieve & operator=(Sieve const & old)=delete;
    Sieve & operator=(Sieve &&old) = delete; //neither movable nor copyable. (not movable due to mutexes)

    #if GAUSS_SIEVE_IS_MULTI_THREADED == true
    explicit Sieve(LatticeBasisType B, unsigned int k=2, unsigned int num_threads=0, TermCondType const termcond = nullptr, unsigned int verbosity_=2, int seed_sampler = 0);
    #else
    explicit Sieve(LatticeBasisType B, unsigned int k=2, TermCondType const termcond = nullptr, unsigned int verbosity_=2, int seed_sampler = 0);
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
    //explicit Sieve(std::string const &infilename); //read from dump file.
    ~Sieve();
    static bool constexpr class_multithreaded =  GAUSS_SIEVE_IS_MULTI_THREADED;
    //class_multithreaded is for introspection, is_multithreaded is what the caller wants (may differ if we dump and re-read with different params)

    //void run_sieve(int k); //runs k-sieve

    //LPType get_SVP() = delete;  //obtains Shortest vector and it's length. If sieve has not yet run, start it. Not yet implemented.

    void run();                 //runs the sieve specified by the parameters. Dispatches to the corresponding k-sieve

    void run_2_sieve(); //actually runs the Gauss Sieve with k=2
    //void run_3_sieve(); //actually runs the Gauss Sieve with k=3
    //void run_k_sieve(); //runs Gauss Sieve with arbitrary k

    #if GAUSS_SIEVE_IS_MULTI_THREADED == true
//    void sieve_2_thread(int const thread_id);   //function for worker threads
//    void sieve_3_thread(int const thread_id);
//    void sieve_k_thread(int const thread_id);
    #else
    void sieve_2_iteration (FastAccess_Point &p); //one run through the main_list (of 2-sieve)
    //void sieve_3_iteration (LatticePoint<ET> &p); //one run through the main_list (of 3-sieve)
    //void sieve_3_iteration_new (LatticePoint<ET> &p); //new run through the main_list (of 3-sieve) usign map for filtered_list
    //void sieve_3_iteration_test (LatticePoint<ET> &p);
    //void sieve_k_iteration (LatticePoint<ET> &p);
    #endif


    void print_status(int verb = -1, std::ostream &out = cout) {dump_status_to_stream(out,verb);};      //prints status to out. verb overrides the verbosity unless set to -1.
    void dump_status_to_file(std::string const &outfilename, bool overwrite = false);                   //dumps to file (verbosity overridden to 3)
    void dump_status_to_stream(ostream &of, int verb=-1);       //dumps to stream. Can be read back if verb>= 3. Otherwise, verbosity determines what is output.

//getter / setter functions

    int get_verbosity() const                                   {return verbosity;};                //non-thread-safe
    void set_verbosity(int new_verbosity)                       {verbosity=new_verbosity;return;};  //non-thread-safe
    unsigned int get_lattice_rank() const                       {return lattice_rank;};             //non-thread-safe
    Dimension<nfixed> get_ambient_dimension() const             {return ambient_dimension;};        //non-thread-safe
    unsigned int get_k() const                                  {return sieve_k;};                  //non-thread-safe
    void set_k(unsigned int new_k)                              {sieve_k=new_k;return;};            //non-thread-safe
    bool is_multithreaded_wanted() const                        {return multi_threaded_wanted;};    //Note: No setter
    #if GAUSS_SIEVE_IS_MULTI_THREADED == true
    void set_num_threads(unsigned int t);                                                            //non-thread safe, only call while suspended. In SieveMT.cpp
    unsigned int get_num_threads() const                        {return num_threads_wanted;};
    #else
    static unsigned int constexpr get_num_threads()             {return 1;};
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED

    bool update_shortest_vector_found(FastAccess_Point const & newvector);
    LatticeBasisType const & get_original_basis()               {return original_basis;};
//  LPType get_shortest_vector_found() const                    {return shortest_vector_found;}; //TODO: Thread-safety
    FastAccess_Point const & get_shortest_vector_found();
//  ET get_best_length2() const                                 {return get_shortest_vector_found().norm2;};
    ET get_best_length2();                                       //{return (main_list.cbegin())->get_details().norm2;}; //TODO: Change to above
    bool check_whether_sieve_is_running() const                 {return (sieve_status==SieveStatus::sieve_status_running);};
    unsigned long int get_number_of_collisions() const          {return number_of_collisions;};
    unsigned long int get_number_of_points_sampled() const      {return number_of_points_sampled;};
    unsigned long long get_number_of_points_constructed() const {return number_of_points_constructed;};
    unsigned long int get_current_list_size() const             {return current_list_size;};
    unsigned long int get_current_queue_size()                  {return main_queue.size();}; //TODO : fix const-correctness
    unsigned long long get_number_of_scprods_level1() const     {return number_of_total_scprods_level1;};
    //TODO:DO the same for all levels
    void set_termination_condition(TermCondType termcond)       {term_cond = termcond;}; //TODO: If we default - initialize (and own in this case), may need to delete previous value.
private:

//Use termination Condition to check whether we are done, based on statistics so far.
    bool check_if_done();
    //ET ComputeMinkowski2Bound(); //computes Minkowski bound for the square(!) length. May need to round upwards.

//Note: The member fields of Sieve denote the (global) "internal" status of the sieve during a run or execution.
//It should be possible to dump the status to harddisk and resume from dump using that information.
//It should also be possible to suspend the run of the sieve, change (certain) parameters (like k!) and resume.

//main data that is changing.

    MainListType main_list;
    //MainListType3 main_list_test;
    MainQueueType main_queue;
//    FilteredListType filtered_list;
//    FilteredListType2 filtered_list2;
//    FilteredListTypeP filtered_listp;

//information about lattice and algorithm we are using

    LatticeBasisType original_basis;
    unsigned int lattice_rank;
    Dimension<nfixed> ambient_dimension; //consider merging these into a latticespec struct.
    bool multi_threaded_wanted;
    #if GAUSS_SIEVE_IS_MULTI_THREADED == true
    unsigned int num_threads_wanted;        //number of threads that we spawn
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
    unsigned int sieve_k; //parameter k of the sieve currently running.
    //SamplerType sampler; //TODO: Thread-safety. Move control to queue.
    int verbosity;       //ranged from 0 to 3 (0 : silent, 1 : errors only, 2 : more output, 3 : debug

//public:  //made public to avoid complicated (due to template hack) friend - declaration.
//            //TODO : Change term-cond to a user-provided function.
//    TerminationConditions<ET> term_cond;
private:
    TermCondType term_cond;
    enum class SieveStatus
    {
        sieve_status_error  =  -1,      //indicates an error (add error codes as neccessary)
        sieve_status_init   =  1,       //we have initialized data (and may yet initialize some more, but sieve has not started
        sieve_status_running=  2,       //sieve is currently running
        sieve_status_suspended=3,       //sieve is currently suspended. Useful for dumping / cleanup of internal data structures.
        sieve_status_finished=100       //sieve has finished
    } sieve_status; //thread safety?
    FastAccess_Point *shortest_vector_found; //including its length //TODO: Thread-safety

//statistics
#if GAUSS_SIEVE_IS_MULTI_THREADED == false
    unsigned long int number_of_collisions;
    unsigned long int number_of_points_sampled;
    unsigned long long int number_of_points_constructed; //sampling  + succesful pairs
    unsigned long int current_list_size;
    unsigned long long int number_of_total_scprods;
    unsigned long long int number_of_total_scprods_level1; //for k=2, case
    unsigned long long int number_of_total_scprods_level2; //for k=2,3 cases
    unsigned long long int number_of_total_scprods_level3; //for k=2,3,4 cases
    unsigned long long int number_of_exact_scprods;
    unsigned long long int number_of_mispredictions; //could not reduce in spite of approximation saying so.
#else //note: we might collect statistics per-thread and merge occasionally. This means these statistics might be inaccurate.
    atomic_ulong number_of_collisions;
    atomic_ulong number_of_points_sampled;
    atomic_ullong number_of_points_constructed;
    atomic_ulong current_list_size;
    atomic_ullong number_of_scprods;
    atomic_ullong number_of_total_scprods_level1;
    atomic_ullong number_of_exact_scprods;
    atomic_ullong number_of_mispredictions;
#endif // GAUSS_SIEVE_IS_MULTI_THREADED

#if GAUSS_SIEVE_IS_MULTI_THREADED==true
    GarbageBin<typename MainListType::DataType> * garbage_bins; //dynamically allocated array of garbage bins.
    std::mutex dump_mutex;
    std::mutex shortest_vector_mutex;
#endif // GAUSS_SIEVE_IS_MULTI_THREADED

//TODO: total time spent?
};

template class Sieve<Z_NR<long>,GAUSS_SIEVE_IS_MULTI_THREADED,-1>;
template class Sieve<Z_NR<double>,GAUSS_SIEVE_IS_MULTI_THREADED,-1>;
template class Sieve<Z_NR<mpz_t>,GAUSS_SIEVE_IS_MULTI_THREADED,-1>;

/*DUMPING / READING ROUTINES */

//Note: Actually, we want an unformatted binary dump. Unfortunately, the underlying FPLLL types support only
// << and >> operations with formated input / output. So we will do with formatted input / output for now.
//The main issue here is that mpz_t only provides formatted stream - I/O or unformatted I/O via old-style C FILE* interfaces, neither of which is what we really want.

// We assume that << writing data to a filestream and >> reading it back will give back the same data.

// This might cause problems if e.g.:
// separation chars (whitespace) are used in a data field
// locales are different
// some subobject changes format flags
// output loses data (e.g. rounding of floats)

//helper functions:

/*
Reads length(str) chars from stream is, expecting them to equal str. If what is read differs we output false. If verbose, we also display an error.
*/


#define SIEVE_JOINT_H
#endif // DO_INCLUDE_SIEVE_JOINT_H
