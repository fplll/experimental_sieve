/*
------DO NOT INCLUDE THIS FILE MANUALLY.--------
------     USE SieveGauss.h INSTEAD     --------
*/

// SIEVE_JOINT_H_ST and SIEVE_JOINT_H_MT are separate include guards.
// They are set if we already included this file (for the ST / MT case) respectively.
// There is also an include guard SIEVE_JOINT_H, which is set _after_ we included this file at least
// once already.
// Use this to condition on the second pass.

// clang-format off
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
// clang-format on

#ifdef DO_INCLUDE_SIEVE_JOINT_H  // serves as the actual include guard

//
// end of (most) preprocessor stuff
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

/*INCLUDES */

// clang-format does for some strange reason not reorder includes here...
#include "DefaultIncludes.h"
#include "GlobalStaticData.h"
#include "GaussQueue.h"
#include "DefaultTermConds_impl.h"  // TODO: Change
#include "GaussListBitapprox.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "LatticeBases.h"
#include "Statistics.h"
#include "SimHash.h"
#include "FilteredPoint.h"

namespace GaussSieve
{
template <class SieveTraits, bool MT> class Sieve;
}

#endif  // end of ONLY-ONCE part

/*
EVERYTHING BELOW HERE IS POTENTIALLY INCLUDED TWICE.
TEMPLATES WITH TEMPLATE ARGUMENT GAUSS_SIEVE_IS_MULTI_THREADED
GO HERE.
*/

// The following may be included once or twice (with different values for
// GAUSS_SIEVE_IS_MULTI_THREADED)

// TODO: Move to where it is actually used.

// template<class ET>
// class CompareQueue<ET, GAUSS_SIEVE_IS_MULTI_THREADED>{
//    public:
//     bool operator() (const FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> & el1, const
//     FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> & el2) const
//     {
//        return el1.get_sc_prod() > el2.get_sc_prod();  // inner products are in decreasing order
//    }
//};

#ifndef GAUSS_SIEVE_IS_MULTI_THREADED
#error Something very bad just happened
#endif

namespace GaussSieve
{

template <class SieveTraits> class Sieve<SieveTraits, GAUSS_SIEVE_IS_MULTI_THREADED>
{
public:
  /*DATA TYPES*/

  using FastAccess_Point = typename SieveTraits::FastAccess_Point;
  using MainQueueType    = GaussQueue<SieveTraits, GAUSS_SIEVE_IS_MULTI_THREADED>;
  using MainListType     = GaussListWithBitApprox<SieveTraits, GAUSS_SIEVE_IS_MULTI_THREADED>;
  using LatticeBasisType = SieveLatticeBasis<SieveTraits, GAUSS_SIEVE_IS_MULTI_THREADED>;
  using InputBasisType   = typename SieveTraits::InputBasisType;
  using DimensionType    = typename SieveTraits::DimensionType;
  using LengthType       = typename SieveTraits::LengthType;
  using CoordinateSelectionUsed = typename SieveTraits::CoordinateSelectionUsed;
  using SimHashGlobalData       = GlobalBitApproxData<CoordinateSelectionUsed>;

  // TODO: Remove this typedef
  using Filtered_Point = FilteredPoint2<SieveTraits, GAUSS_SIEVE_IS_MULTI_THREADED>;

  using FilteredListType            = std::vector<Filtered_Point>;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;
  using SieveStatistics = GaussSieveStatistics<SieveTraits, GAUSS_SIEVE_IS_MULTI_THREADED>;

  template <class, bool> friend class GaussSieveStatistics;

  // Termination condition may be some class *derived from* this.
  using TermCondType     = TerminationCondition<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>;

public:
  /*FRIENDS */
  friend GaussQueue<SieveTraits, GAUSS_SIEVE_IS_MULTI_THREADED>;

  /*CONSTRUCTORS / DESTRUCTORS */

  // neither movable nor copyable. (not movable due to mutexes)
  Sieve()                         = delete;
  Sieve(Sieve const &)            = delete;
  Sieve(Sieve &&)                 = delete;
  Sieve &operator=(Sieve const &) = delete;
  Sieve &operator=(Sieve &&)      = delete;

#if GAUSS_SIEVE_IS_MULTI_THREADED == true
  explicit Sieve InputBasisType const &B, unsigned int const k=2,
                  unsigned int const num_threads=0, TermCondType * const termcond = nullptr,
                  unsigned int const verbosity_=2, int seed_sampler = 0);
#else
  explicit Sieve( InputBasisType const & B, unsigned int const k=2,
                    TermCondType * const termcond = nullptr, unsigned int const verbosity_=2,
                    int seed_sampler = 0);
#endif // GAUSS_SIEVE_IS_MULTI_THREADED

  //explicit Sieve(std::string const &infilename); // read from dump file (NOT IMPLEMENTED)
  ~Sieve();
  static bool constexpr class_multithreaded =  GAUSS_SIEVE_IS_MULTI_THREADED;
  //class_multithreaded is for introspection, is_multithreaded is what the caller wants (may differ if we dump and re-read with different params)


  void run();         //runs the sieve specified by the parameters. Dispatches to the corresponding k-sieve

  void run_2_sieve(); //calls sieve_2_iteration for 2-reduction until the termination conditions are satisfied
  void run_3_sieve(); //calls sieve_3_iteration for 3-reduction until the termination conditions are satisfied
  //void run_k_sieve(); //runs Gauss Sieve with arbitrary k

#if GAUSS_SIEVE_IS_MULTI_THREADED == true
//  void sieve_2_thread(int const thread_id);   //function for worker threads
//  void sieve_3_thread(int const thread_id);
//  void sieve_k_thread(int const thread_id);
#else
  void sieve_2_iteration (FastAccess_Point &p); //one run through the main_list (of 2-sieve)
  template<class LHS, class RHS>
  bool check2red(LHS &&p1, RHS &&p2, int &scalar);

  void sieve_3_iteration (FastAccess_Point &p); //one run through the main_list (of 3-sieve)
  //void sieve_k_iteration (LatticePoint<ET> &p);
#endif
  
  /* PRINTING / DUMPING ROUTINES */
  
  //prints status to out. verb overrides the verbosity unless set to -1.
  void print_status(int verb = -1, std::ostream &out = std::cout) {dump_status_to_stream(out,verb);};
  void dump_status_to_file(std::string const &outfilename, bool overwrite = false);                   //dumps to file (verbosity overridden to 3)
  void dump_status_to_stream(std::ostream &of, int verb=-1);       //dumps to stream. Can be read back if verb>= 3. Otherwise, verbosity determines what is output.


  /* GETTERS / SETTERS */

  int get_verbosity() const                                   {return verbosity;};                //non-thread-safe
  void set_verbosity(int const new_verbosity)                 {verbosity=new_verbosity;return;};  //non-thread-safe
  unsigned int get_lattice_rank() const                       {return lattice_rank;};             //non-thread-safe
  DimensionType get_ambient_dimension() const                 {return ambient_dimension;};        //non-thread-safe

#ifdef PROGRESSIVE
  uint_fast16_t get_progressive_rank() const                  {return progressive_rank;};
  void increase_progressive_rank();
  void set_target_list_size(double const target)              {target_list_size = target; return;};
  double get_target_list_size() const                         {return target_list_size;};
  bool check_if_enough_short_vectors();
#endif

  unsigned int get_k() const                                  {return sieve_k;};                  //non-thread-safe
  void set_k(unsigned int const new_k)                        {sieve_k=new_k;return;};            //non-thread-safe
  bool is_multithreaded_wanted() const                        {return multi_threaded_wanted;};    //Note: No setter

#if GAUSS_SIEVE_IS_MULTI_THREADED == true
  void set_num_threads(unsigned int t);                                                            //non-thread safe, only call while suspended. In SieveMT.cpp
  unsigned int get_num_threads() const                        {return num_threads_wanted;};
#else
  static unsigned int constexpr get_num_threads()             {return 1;};
#endif // GAUSS_SIEVE_IS_MULTI_THREADED

  bool update_shortest_vector_found(FastAccess_Point const & newvector);
  LatticeBasisType const & get_basis()                        {return lattice_basis;};
  InputBasisType const & get_input_basis()                    {return original_basis;};
  FastAccess_Point const & get_shortest_vector_found();
  LengthType get_best_length2()                               {return shortest_vector_found->get_norm2();};
  bool check_whether_sieve_is_running() const                 {return (sieve_status==SieveStatus::sieve_status_running);};
  unsigned long get_final_list_size()                         {return main_list.size();};

  // STAT_MARK

  inline void set_termination_condition(TermCondType * const termcond)
  {
    if(term_cond_owned && (term_cond!=nullptr)) delete term_cond;
    term_cond_owned = false;
    term_cond = termcond;
  }
private:

//Use termination Condition to check whether we are done, based on statistics so far.
  bool check_if_done();

//Note: The member fields of Sieve denote the (global) "internal" status of the sieve during a run or execution.
//It should be possible to dump the status to harddisk and resume from dump using that information.
//It should also be possible to suspend the run of the sieve, change (certain) parameters (like k!) and resume.

  DimensionType ambient_dimension; //consider merging these into a latticespec struct.
  
  GlobalStaticDataInitializer global_static_data;
  StaticInitializer<SimHashGlobalData> static_init_sim_hash_global_data;
  StaticInitializer<FastAccess_Point> static_init_fast_access_point;


  InputBasisType original_basis;
  LatticeBasisType lattice_basis; //coverted variant of original_basis

  //main data that is changing.
  MainListType main_list;
  MainQueueType main_queue;
  
  unsigned int lattice_rank;

#ifdef PROGRESSIVE
  uint_fast16_t progressive_rank;
#endif

  bool multi_threaded_wanted;
#if GAUSS_SIEVE_IS_MULTI_THREADED == true
  unsigned int num_threads_wanted;        //number of threads that we spawn
#endif // GAUSS_SIEVE_IS_MULTI_THREADED


#ifdef PROGRESSIVE
  double target_list_size;
#endif
  unsigned int sieve_k; //parameter k of the sieve currently running.
  int verbosity;       //ranged from 0 to 3 (0 : silent, 1 : errors only, 2 : more output, 3 : debug


private:
  bool term_cond_owned;  // Do we own the following pointer. Required for proper clean-up.
  TermCondType *term_cond;  // consider using smart pointers
  enum class SieveStatus
  {
    sieve_status_error     =  -1,      //indicates an error (add error codes as neccessary)
    sieve_status_init      =  1,       //we have initialized data (and may yet initialize some more, but sieve has not started
    sieve_status_running   =  2,       //sieve is currently running
    sieve_status_suspended =3,       //sieve is currently suspended. Useful for dumping / cleanup of internal data structures.
    sieve_status_finished  =100       //sieve has finished
  } sieve_status; // thread safety?

  // Note: This is a pointer for now due to initialization order issues.
  // TODO: Create RAII class (this would actually solve several problems...)
  // Note: The latter was actually done...

  FastAccess_Point *shortest_vector_found; // including its length //TODO: Thread-safety
public:  // switched to public to get list-sizes
  SieveStatistics statistics;

#if GAUSS_SIEVE_IS_MULTI_THREADED==true
  GarbageBin<typename MainListType::DataType> * garbage_bins; //dynamically allocated array of garbage bins.
  std::mutex dump_mutex;
  std::mutex shortest_vector_mutex;
#endif // GAUSS_SIEVE_IS_MULTI_THREADED

}; // end of class definition

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
}

#define SIEVE_JOINT_H
#endif // DO_INCLUDE_SIEVE_JOINT_H

// clang-format on
