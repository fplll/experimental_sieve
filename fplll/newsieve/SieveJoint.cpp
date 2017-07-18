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

#include "Sampler.cpp"
#include "GaussQueue.cpp"
#include "PointListNew.cpp"
#include "TermCondNew.cpp"

//End of things included only once.
#endif // SIEVE_JOINT_CPP

//#define SIEVE_GAUSS_H //This means that SieveGauss.h is already included. SieveGauss really is a wrapper around SieveJoint.h, and here we have to call SieveJoint.h direcly.
//#include "SieveJoint.h" //will include it with GAUSS_SIEVE_IS_MULTI_THREADED set correctly. This file is actually allowed to include SieveJoint.h directly.

#define SIEVE_FILE_ID "kTuple-Sieve dump file"
//version string for dump file
#define SIEVE_VER_STR "Version TEST1"

template<class ET,int nfixed>
void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::dump_status_to_file(std::string const &outfilename, bool overwrite)
{
    assert(!check_whether_sieve_is_running());
    if(verbosity>=2)
    {
        cout << "Dumping to file " << outfilename << " ..." << endl;
    }
    if(!overwrite)
    {
        //checks if file exists
        struct stat buffer;
        if (::stat(outfilename.c_str(), &buffer) == 0)
        {
            if (verbosity>=1)
            {
                cerr << "Trying to dump to existing file without overwrite flag set. Aborting dump." << endl;
            }
            return;
        }
    }
    std::ofstream of(outfilename,std::ofstream::out | std::ofstream::trunc);
    dump_status_to_stream(of, 3); //verbosity set to 3 to dump everything.
    if(verbosity>=2)
    {
        cout << "Dump finished." << endl;
    }
}

template<class ET,int nfixed>
void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::dump_status_to_stream(ostream &of, int verb)
{
    int howverb = verb==-1 ? verbosity : verb;
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
    if(howverb>=1) of << "--Statistics--" << endl;
    if(howverb>=1) of << "Number of collisions=" << number_of_collisions << endl;
    if(howverb>=1) of << "Number of points Sampled=" << number_of_points_sampled << endl;
    if(howverb>=1) of << "Number of points Constructed=" << number_of_points_constructed << endl;
    //if(howverb>=1) of << "Number of approx. scalar products=" << number_of_scprods << endl;
    if(howverb>=1) of << "Number of exact scalar products=" << number_of_exact_scprods << endl;
    if(howverb>=1) of << "Number of mispredictions=" << number_of_mispredictions << endl;
    if(howverb>=1) of << "Current List Size=" << get_current_list_size() << endl;
    if(howverb>=1) of << "Current Queue Size="<< get_current_queue_size()<< endl;
    if(howverb>=1) {
        of << "Best vector found so far=";
        shortest_vector_found->write_to_stream(of);
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
        cerr << "Dumping of main queue not supported yet.";
        of << "--End of Main Queue--";
        {

        }
    }
}


template<class ET, int nfixed> //ET : underlying entries of the vectors. Should be a Z_NR<foo> - type.
#if GAUSS_SIEVE_IS_MULTI_THREADED==true
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::Sieve(LatticeBasisType B, unsigned int k, unsigned int num_threads, TermCondType const termcond, unsigned int verbosity_, int seed_sampler):
#else
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::Sieve(LatticeBasisType B, unsigned int k, TermCondType const termcond, unsigned int verbosity_, int seed_sampler):
#endif
    main_list(),
    main_queue(this),
//    filtered_list(),
    original_basis(B),
    lattice_rank(B.get_rows()),
    ambient_dimension(B.get_cols()), //Note : this means that rows of B form the basis.
    multi_threaded_wanted(GAUSS_SIEVE_IS_MULTI_THREADED),
    #if GAUSS_SIEVE_IS_MULTI_THREADED == true
    num_threads_wanted(num_threads),
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
    sieve_k(k),
    verbosity(verbosity_),
    term_cond(termcond),
    sieve_status(SieveStatus::sieve_status_init),
    shortest_vector_found(nullptr), //TODO : initialize to meaningful value, e.g. any vector of B.
    number_of_collisions(0),
    number_of_points_sampled(0),
    number_of_points_constructed(0),
    current_list_size(0),
    number_of_total_scprods_level1(0),
    number_of_total_scprods_level2(0),
    number_of_total_scprods_level3(0),
    number_of_exact_scprods(0),
    number_of_mispredictions(0)
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    ,garbage_bins(nullptr)
    #endif // GAUSS_SIEVE_IS_MULTI_THREADED

{
    if (nfixed!=-1)
    {
        assert(B.get_cols() == nfixed );
    }
    FastAccess_Point::class_init(ambient_dimension);
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    if(num_threads_wanted==0) //0 means we take a meaningful default, which is given by thread::hardware_concurrency
    num_threads_wanted = std::max(std::thread::hardware_concurrency(),static_cast<unsigned int>(1)); //Note: hardware_concurrency might return 0 for "unknown".
    #endif
    //unsigned int n = lattice_rank;
    //auto it = main_list.before_begin();
    //assert(main_list.empty()); We don't have a function to check that yet...
    if (verbosity>=2) {cout <<"Initializing list with original basis..." << endl;}
    auto it = main_list.cbegin();
    for (unsigned int i=0; i<lattice_rank; ++i)
    {
        FastAccess_Point tmppoint(original_basis[i]);
        it = main_list.insert_before(it, std::move(tmppoint));
        ++it;
//        ExactLatticePoint<ET,nfixed> * new_basis_vector = new ExactLatticePoint<ET,nfixed> ( conv_matrixrow_to_lattice_point<ET,nfixed> (original_basis[i]));
//        main_list.insert_before(it,  static_cast<CompressedPoint<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed> >(new_basis_vector) );
    }
    current_list_size+=lattice_rank;
//    #if GAUSS_SIEVE_IS_MULTI_THREADED == false
    if(verbosity>=2)    {cout << "Sorting ...";}
        main_list.sort();
    for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it) {cout << (*it).get_norm2() << endl;}; //check for sort()
        
    if(verbosity>=2)    {cout << "is finished." << endl;}


    //FIXME: Initialize shortest vector

    shortest_vector_found = new FastAccess_Point (main_list.cbegin()->make_copy());
    cout << "shortest_vector_found is initialized " << endl << flush;
    
//    #endif // GAUSS_SIEVE_IS_MULTI_THREADED
    //TODO : enable sorting for multithreaded case.
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    garbage_bins = new GarbageBin<typename MainListType::DataType>[num_threads_wanted]; //maybe init later.
    #endif
    assert(main_queue.sampler!=nullptr);
    main_queue.sampler->init(this);
    
    cout << "sampler is initialized " << endl << flush;

};

template<class ET,int nfixed>
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::~Sieve()
{
    #if GAUSS_SIEVE_IS_MULTI_THREADED==true
    delete[] garbage_bins;
    delete shortest_vector_found;
    #endif
};

template<class ET, int nfixed>
bool Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED,nfixed>::check_if_done()
{
    return (term_cond->check(this) != 0)?true:false;
};


#define SIEVE_JOINT_CPP
#endif
