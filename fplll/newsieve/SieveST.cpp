#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(SIEVE_GAUSS_SINGLE_THREADED) || GAUSS_SIEVE_IS_MULTI_THREADED == true
#error SieveST.cpp included with wrong macro settings
#endif

namespace GaussSieve{

//We may always assumed that SieveJoint.cpp is prepended before this file


//class next_block;
//using namespace LatticeApproximations; // to be able to use ApproxTypeNorm2 to store inner-produces scaled by length

//class main_list;

template<class SieveTraits>
bool Sieve<SieveTraits,false>::update_shortest_vector_found(FastAccess_Point const & newvector)
{
    if(newvector < (*shortest_vector_found) )
    {
        delete shortest_vector_found;
        shortest_vector_found = new FastAccess_Point (newvector.make_copy());
        return true;
    }
    return false;
}

template<class SieveTraits>
typename Sieve<SieveTraits,false>::FastAccess_Point const & Sieve<SieveTraits,false>::get_shortest_vector_found()
{
    return *shortest_vector_found;
}

template<class SieveTraits>
typename Sieve<SieveTraits,false>::EntryType Sieve<SieveTraits,false>::get_best_length2()
{
    return shortest_vector_found->get_norm2();
}

template<class SieveTraits> void Sieve<SieveTraits,false>::run()
{
    if (verbosity >=2) std::cout << "the shortest vector in the input basis has norm2 = " << get_shortest_vector_found().get_norm2() << std::endl;
    //int MaxIteration = 8000;
    if (term_cond==nullptr)
    {
        if (verbosity >=1) std::cerr << "No termination condition set. Aborting.";
        return;
    }
    if(sieve_k <=1)
    {
        if(verbosity >=1) std::cerr << "sieve_k <=1. Aborting.";
        return;
    }
    sieve_status=SieveStatus::sieve_status_running;
    term_cond->init(this); //initialisation of termination conditions.
    if (verbosity >=2) std::cout << "start " << sieve_k << "-Sieve" << std::endl;

    //dispatch

    switch (sieve_k)
    {
        case 2: run_2_sieve(); break;
        case 3: run_3_sieve(); break;
        //default:run_k_sieve(); break;
        default: assert(false);
    }
    sieve_status = SieveStatus::sieve_status_finished;
}

template<class SieveTraits> void Sieve<SieveTraits,false>::run_2_sieve()
{
//    typename SieveTraits::GaussQueue_ReturnType p;
  int i=0;

  std::cout << "start 2-sieve " << std::endl;

  //set target list-size for progressive sieving
#ifdef PROGRESSIVE
  set_target_list_size(list_size_k2);
#endif

  while (!check_if_done() )
  {

#ifdef PROGRESSIVE
    if ((progressive_rank < lattice_rank)  && check_if_enough_short_vectors() )
    {
      increase_progressive_rank();
    }
#endif
  //convert here???

//        GaussSieve::FastAccess_Point<ET, false, nfixed> p_converted (std::move(p));
    typename SieveTraits::FastAccess_Point p = main_queue.true_pop(); // may need conversion.
  //typename SieveTraits::FastAccess_Point p = static_cast<typename SieveTraits::GaussList_StoredPoint>(main_queue.true_pop());
//        Sieve<ET,false,nfixed>::sieve_2_iteration(p_converted);
//        std::cout << p << std::endl << std::flush;

#ifdef USE_LSH
    hash_sieve_2_iteration(p);
#else
    sieve_2_iteration(p);
#endif

    ++i;
    if (( i % 1000 == 0) && (verbosity >=2))
    {
    // STAT_MARK
      std::cout << "[" << i << "]"
      << " |L|=" << statistics.get_current_list_size()
      << " |Q|=" << main_queue.size()
      << " #samples = " << statistics.get_number_of_points_sampled()
      << " |sv|= " <<  get_best_length2() << std::endl << std::flush;
    }
  }
}


template<class SieveTraits> void Sieve<SieveTraits,false>::run_3_sieve()
{
  int i=0;

  std::cout << "start 3-sieve " << std::endl;

#ifdef PROGRESSIVE
  set_target_list_size(list_size_k3);
#endif

  while (!check_if_done() )
  {
#ifdef PROGRESSIVE
    if ((progressive_rank < lattice_rank)  && check_if_enough_short_vectors() )
    {
      increase_progressive_rank();
    }
#endif
    typename SieveTraits::FastAccess_Point p = main_queue.true_pop();
    
    sieve_3_iteration(p);
    
    ++i;
    if (( i % 1000 == 0) && (verbosity >=2))
    {
      std::cout << "[" << i << "]"
      << " |L|=" << statistics.get_current_list_size()
      << " |Q|=" << main_queue.size() // STAT_MARK
      << " #samples = " << statistics.get_number_of_points_sampled()
      << " |FL|= " << statistics.get_filtered_list_size()
      << " |sv|= " << get_best_length2() <<  std::endl;
    }
  }
}


} // end namespace

// TODO: #include placement

//#include "HyperplaneLSH.h"
#include "SieveST2.cpp"
#include "SieveST3.cpp"
//#include "SieveSTk.cpp"

#endif

