// clang-format off

//SIEVE_JOINT_CPP_ST and
//SIEVE_JOINT_CPP_MT are separate include guards for the single and multi-threaded pass over this file.
//SIEVE_JOINT_CPP is set *after* we have passed over this file once. Use it to condition on the second run.

#undef DO_INCLUDE_SIEVE_JOINT_CPP
#ifndef GAUSS_SIEVE_IS_MULTI_THREADED
    #error wrong usage of SieveJoint.cpp --1
#endif

#if GAUSS_SIEVE_IS_MULTI_THREADED == false
    #if !defined(SIEVE_GAUSS_SINGLE_THREADED)
        #error wrong usage of SieveJoint.cpp -- 2
    #endif

    #ifndef SIEVE_JOINT_CPP_ST
        #define SIEVE_JOINT_CPP_ST
        #define DO_INCLUDE_SIEVE_JOINT_CPP
    #endif
#elif GAUSS_SIEVE_IS_MULTI_THREADED == true
    #if !defined(SIEVE_GAUSS_MULTI_THREADED)
        #error wrong usage of SieveJoint.cpp -- 3
    #endif
    #ifndef SIEVE_JOINT_CPP_MT
        #define SIEVE_JOINT_CPP_MT
        #define DO_INCLUDE_SIEVE_JOINT_CPP
    #endif
#endif

#ifdef DO_INCLUDE_SIEVE_JOINT_CPP

//actual code starts here.
//Be aware that code may be read twice.

#ifndef SIEVE_JOINT_CPP //code in this block only read once.

//#include "Sampler_impl.h"
//#include "GaussQueue_impl.h"
//#include "PointListNew.cpp"

//End of things included only once.
#endif // SIEVE_JOINT_CPP

//#define SIEVE_GAUSS_H //This means that SieveGauss.h is already included. SieveGauss really is a wrapper around SieveJoint.h, and here we have to call SieveJoint.h direcly.
//#include "SieveJoint.h" //will include it with GAUSS_SIEVE_IS_MULTI_THREADED set correctly. This file is actually allowed to include SieveJoint.h directly.

#define SIEVE_FILE_ID "kTuple-Sieve dump file"
//version string for dump file
#define SIEVE_VER_STR "Version TEST1"

namespace GaussSieve{

template<class SieveTraits>
void Sieve<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>::dump_status_to_file(std::string const &outfilename, bool overwrite)
{
    assert(!check_whether_sieve_is_running());
    if(verbosity>=2)
    {
        std::cout << "Dumping to file " << outfilename << " ..." << std::endl;
    }
    if(!overwrite)
    {
        //checks if file exists
        struct stat buffer;
        if (::stat(outfilename.c_str(), &buffer) == 0)
        {
            if (verbosity>=1)
            {
                std::cerr << "Trying to dump to existing file without overwrite flag set. Aborting dump." << std::endl;
            }
            return;
        }
    }
    std::ofstream of(outfilename,std::ofstream::out | std::ofstream::trunc);
    dump_status_to_stream(of, 3); //verbosity set to 3 to dump everything.
    if(verbosity>=2)
    {
        std::cout << "Dump finished." << std::endl;
    }
}

template<class SieveTraits>
void Sieve<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>::dump_status_to_stream(std::ostream &of, int verb)
{
    using std::endl;
    int howverb = ( (verb==-1) ? verbosity : verb);
    if(howverb>=2) of << SIEVE_FILE_ID << endl;
    if(howverb>=2) of << SIEVE_VER_STR << endl;
    if(howverb>=2) of << "--Params--" << endl;
    if(howverb>=2) of << "Multithreaded version=" << class_multithreaded << endl;
    if(howverb>=1) of << "Multithreaded is wanted" << multi_threaded_wanted << endl;
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    if(howverb>=1) of << "Number of Threads=" << num_threads_wanted << endl;
    #endif
    if(howverb>=2) of << "k=" << sieve_k << endl;
    if(howverb>=2) of << "verbosity=" << verbosity << endl;
    if(howverb>=1) of << "sieve_status=" << static_cast<int>(sieve_status) << endl;
    if(howverb>=2) of << "Lattice rank=" << lattice_rank << endl;
    if(howverb>=2) of << "Ambient dimension=" << ambient_dimension << endl;
    if(howverb>=2) of << "Termination Conditions:" << endl;
    if(howverb>=2) of << term_cond;
    if(howverb>=2) of << "Sampler:"<< endl;
    //if(howverb>=2) of << "Sampler Initialized" << static_cast<bool>(sampler!=nullptr) << endl;
    //if(sampler!=nullptr)
    //{
    //    if(howverb>=3) cerr << "Note : Dumping of internal data of sampler not yet supported" << endl;
    //    //dump internals of sampler?
    //}
    if(howverb>=3) of << "Original Basis:" << endl;
    if(howverb>=3) of << original_basis;
    if(howverb>=2) of << "--End of Params--" << endl << endl;

    // STAT_MARK

    if(howverb>=1) of << "--Statistics--" << endl;
    statistics.dump_status_to_stream(of,howverb);

    if(howverb>=1) {
        of << "Best vector found so far=";
        shortest_vector_found->write_lp_to_stream(of,true);
        of << endl;
    }


    //of << "Best vector found so far=" << shortest_vector_found << endl; //TODO : Display length seperately

    //TODO: Check output
    // if(howverb>=1) {of << "sv is: "; main_list.cbegin().get_exact_point().printLatticePoint();} //TODO: Remove (shortest_vector_found above should rather do this).
	// to check the ratio between the shortest and the longest vectors in the list
//    if(howverb>=1) {
//		auto it = main_list.cend();
//		// operator--() is needed for tests
//		--(it);
//		of << "longest vector is: "; it.get_exact_point().printLatticePoint();
//    }
    if(howverb>=1) of << "--End of Statistics--" << endl << endl;
    if(howverb>=3)
    {
        of << "--Main List--" << endl;

        //TODO: ACTUALLY OUTPUT THE LIST

//        for(auto it = main_list.cbegin(); it!=main_list.cend(); ++it)
//        {
//            of << (*it);
//        }
        of << "--End of Main List--" << endl << endl;
        of << "--Main Queue--" << endl;
//  for(auto it = main_queue.begin();it!=main_queue.end();++it)
        std::cerr << "Dumping of main queue not supported yet.";
        of << "--End of Main Queue--";
        {

        }
    }
}


// constructor
template<class SieveTraits>
#if GAUSS_SIEVE_IS_MULTI_THREADED==true
Sieve<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>::Sieve(
    InputBasisType const & B, unsigned int k, unsigned int num_threads,
    TermCondType * const termcond, unsigned int verbosity_, int seed_sampler):
#else
Sieve<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>::Sieve(
    InputBasisType const & B, unsigned int k,
    TermCondType * const termcond, unsigned int verbosity_, int seed_sampler):
#endif
/**
  Initializer list:
*/
    ambient_dimension(B.get_cols()), //Note : this means that rows of B form the basis.
    global_static_data(ambient_dimension),
    static_init_sim_hash_global_data(global_static_data, seed_sampler),
    static_init_fast_access_point(global_static_data),
    original_basis(B),
    lattice_basis(B,global_static_data),
    main_list(global_static_data),
    main_queue(this,global_static_data),
    lattice_rank(B.get_rows()),
    multi_threaded_wanted(GAUSS_SIEVE_IS_MULTI_THREADED),
    #if GAUSS_SIEVE_IS_MULTI_THREADED == true
    num_threads_wanted(num_threads),
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
    #ifdef USE_LSH
    hash_tables(),
    number_of_hash_tables(0),
    #endif
    sieve_k(k),
    verbosity(verbosity_),
    term_cond_owned(termcond==nullptr),
    term_cond(termcond),
    sieve_status(SieveStatus::sieve_status_init),
    shortest_vector_found(nullptr), // NOTE: Static data in class not initialized! //TODO:True???
    statistics(this)

#if GAUSS_SIEVE_IS_MULTI_THREADED==true
    ,garbage_bins(nullptr)
#endif // GAUSS_SIEVE_IS_MULTI_THREADED
{
// BIG TODO: Move these into their respective classes!

    if (SieveTraits::get_nfixed!=-1)
    {
        assert(B.get_cols() == SieveTraits::get_nfixed );
    }
#if GAUSS_SIEVE_IS_MULTI_THREADED==true
    if(num_threads_wanted==0) //0 means we take a meaningful default, which is given by thread::hardware_concurrency
      num_threads_wanted = std::max(std::thread::hardware_concurrency(),static_cast<unsigned int>(1)); //Note: hardware_concurrency might return 0 for "unknown".
#endif
    //unsigned int n = lattice_rank;
    //auto it = main_list.before_begin();
    //assert(main_list.empty()); We don't have a function to check that yet...
    if (verbosity>=2) {std::cout <<"Initializing list with original basis..." << std::endl;}



    auto it = main_list.cbegin();
    for (unsigned int i=0; i<lattice_rank; ++i)
    {
        it = main_list.insert_before(it, static_cast<typename SieveTraits::GaussList_StoredPoint> (
                                     lattice_basis.get_basis_vector(i).make_copy() ) );
        ++it;
    }

//    statistics.increment_current_list_size_by(lattice_rank); // TODO: Let list manage that itself.

#ifdef USE_LSH
    if (verbosity>=2) {std::cout <<"Initializing LSH..." << std::endl;}
    hash_tables.initialize_hash_tables(ambient_dimension);
    for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it)
    {
            hash_tables.add_to_all_hash_tables(*it);
    }
    //TODO
    //number_of_hash_tables = hash_tables.get_num_of_tables();
  //if(verbosity>=2)    {hash_tables.print_all_tables();};
#endif
//    #if GAUSS_SIEVE_IS_MULTI_THREADED == false

#ifdef PROGRESSIVE
  assert(lattice_rank > 0);
  // Note : The +1 is to ensure correctness if lattice_rank == 1 -- Gotti
  progressive_rank = ((lattice_rank+1) / 2); //TODO: to adjust
  std::cout << "set progressive_rank to " <<progressive_rank << std::endl;
#endif
  if(verbosity>=2)    {std::cout << "Sorting ...";}
  main_list.sort();
    //for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it) {std::cout << (*it).get_norm2() << std::endl;}; //check for sort()

  if(verbosity>=2)    {std::cout << "is finished." << std::endl;}


    //FIXME: Initialize shortest vector

  shortest_vector_found = new FastAccess_Point (main_list.cbegin()->make_copy());
  std::cout << "shortest_vector_found is initialized " << std::endl << std::flush;

//    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
    //TODO : enable sorting for multithreaded case.
#if GAUSS_SIEVE_IS_MULTI_THREADED==true
  garbage_bins = new GarbageBin<typename MainListType::DataType>[num_threads_wanted]; //maybe init later.
#endif
  assert(main_queue.sampler!=nullptr);
  main_queue.sampler->init(this, lattice_basis);

  std::cout << "sampler is initialized " << std::endl << std::flush;

// This should be moved to Statistics.h -- Gotti

#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED

  unsigned int size_of_stat_arrays = GaussSieve::SimHash::num_of_levels * GaussSieve::SimHash::sim_hash_len+1;

  for (unsigned int lvl=0; lvl<SimHash::num_of_levels; ++lvl)
  {
    statistics.red_stat[lvl].resize(size_of_stat_arrays);
    statistics.no_red_stat[lvl].resize(size_of_stat_arrays);

    statistics.red_stat_innloop[lvl].resize(size_of_stat_arrays);
    statistics.no_red_stat_innloop[lvl].resize(size_of_stat_arrays);
  }
#endif

/*
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX
  statistics.no_red_stat.resize(this->ambient_dimension+1);
  statistisc.red_stat.resize(this->ambient_dimension+1);
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_2ND_ORDER
  statistisc.no_red_stat2.resize(2*this->ambient_dimension+1);
  statistisc.red_stat2.resize(2*this->ambient_dimension+1);
#endif
#endif
 */
};

template<class SieveTraits>
Sieve<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>::~Sieve()
{
#if GAUSS_SIEVE_IS_MULTI_THREADED==true
  delete[] garbage_bins;
#endif
  delete shortest_vector_found;
};

template<class SieveTraits>
bool Sieve<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>::check_if_done()
{
  return (term_cond->check(this) != 0)?true:false;
};

#ifdef PROGRESSIVE
//for progressive sieving
//TODO: currently the function checks only for the size of the list
//      theoretically, it should check whether we have enough *short* vectors
template<class SieveTraits>
bool Sieve<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>::check_if_enough_short_vectors()
{
  // check if the current list is long enough and contains enough short vectors

  //TODO:
  double log_bound = 1.23* double(this->get_progressive_rank())*(this->get_target_list_size());
  unsigned long int bound = pow(2, log_bound);
  return (this->statistics.get_current_list_size() > bound);
}
#endif // PROGRESSIVE

#ifdef PROGRESSIVE
template<class SieveTraits>
void Sieve<SieveTraits,GAUSS_SIEVE_IS_MULTI_THREADED>::increase_progressive_rank()
{
  assert(this->progressive_rank < this->get_lattice_rank());
  ++(this->progressive_rank);
  if (this->get_progressive_rank() == this->get_lattice_rank())
  {
    std::cout << "From now on we are full-rank" << std::endl;
  }
  main_queue.sampler->set_progressive_rank(progressive_rank);
}
#endif // PROGRESSIVE

} // end namespace

#define SIEVE_JOINT_CPP
#endif

//clang-format on
