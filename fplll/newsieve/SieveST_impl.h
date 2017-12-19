// clang-format on

#ifndef SIEVE_GAUSS_SINGLE_THREADED_CPP
#define SIEVE_GAUSS_SINGLE_THREADED_CPP

#if !defined(SIEVE_GAUSS_SINGLE_THREADED) || GAUSS_SIEVE_COMPILE_FOR_MULTI_THREADED == true
#error SieveST_impl.h included with wrong macro settings
#endif

namespace GaussSieve
{

// We may always assumed that SieveJoint.cpp is prepended before this file
// This file collects utility functions that are called from sieving routines

/*
  Checks if newvector is shorter that the shortest_vector_found
  If yes, updates the shortest_vector_found
 */
template <class SieveTraits>
bool Sieve<SieveTraits, false>::update_shortest_vector_found(FastAccess_Point const &newvector)
{
  if (newvector < (*shortest_vector_found))
  {
    delete shortest_vector_found;
    shortest_vector_found = new FastAccess_Point(newvector.make_copy());
    return true;
  }
  return false;
}

/*
  Returns the shortest vector found so far
 */
template <class SieveTraits>
typename Sieve<SieveTraits, false>::FastAccess_Point const &
Sieve<SieveTraits, false>::get_shortest_vector_found()
{
  return *shortest_vector_found;
}

/*
 Returns the squared length of shortest vector found so far
 */
/*
template <class SieveTraits>
typename Sieve<SieveTraits, false>::LengthType Sieve<SieveTraits, false>::get_best_length2()
{
  return shortest_vector_found->get_norm2();
}
*/
/*
 Initializes termination conditionsx
 Depending on k={2,3} from user's input, calls the corresponding sieving iteration
*/
template <class SieveTraits> void Sieve<SieveTraits, false>::run()
{
  if (verbosity >= 2)
    std::cout << "the shortest vector in the input basis has norm2 = "
              << get_shortest_vector_found().get_norm2() << std::endl;
  if (term_cond == nullptr)
  {
    if (verbosity >= 1)
      std::cerr << "No termination condition set. Aborting.";
    return;
  }
  if (sieve_k <= 1)
  {
    if (verbosity >= 1)
      std::cerr << "sieve_k <=1. Aborting.";
    return;
  }
  sieve_status = SieveStatus::sieve_status_running;
  term_cond->init(this);  // initialisation of termination conditions.
  if (verbosity >= 2)
    std::cout << "start " << sieve_k << "-Sieve" << std::endl;

  // dispatch
  switch (sieve_k)
  {
  case 2:
//    run_2_sieve();
    break;
  case 3:
//    run_3_sieve();
    break;
  // default:run_k_sieve(); break;
  default:
    assert(false);
  }
  sieve_status = SieveStatus::sieve_status_finished;
}

/*
    runs 2-Sieve; after every print_step_2sieve sieve iterations, prints statistics
 */

// CHANGE TO VECTOR:
/*
template <class SieveTraits> void Sieve<SieveTraits, false>::run_2_sieve()
{
  int i = 0;

// set target list-size for progressive sieving
#ifdef PROGRESSIVE
  set_target_list_size(list_size_k2);
#endif

  while (!check_if_done())
  {
#ifdef PROGRESSIVE
    if ((progressive_rank < lattice_rank) && check_if_enough_short_vectors())
    {
      increase_progressive_rank();
    }
#endif

    // pop p from  main_queue
    typename SieveTraits::FastAccess_Point p = main_queue.true_pop();

    // checks if p participates in 2-reduction, inserts p into main_list
    sieve_2_iteration(p);
    ++i;
    if ((i % SieveTraits::print_step_2sieve == 0) && (verbosity >= 2))
    {
      std::cout << "[" << i << "]"
                << " |L|=" << statistics.get_current_list_size()
                << " |Q|=" << main_queue.size()  // STAT_MARK
                << " #samples = " << statistics.get_number_of_points_sampled()
                << " |sv|= " << get_best_length2() << std::endl;
    }
  }
}
*/

/*
 runs 3-Sieve; after every print_step_2sieve sieve iterations, prints statistics
 */


 // CHANGE TO VECTOR

 /*
template <class SieveTraits> void Sieve<SieveTraits, false>::run_3_sieve()
{
  int i = 0;

#ifdef PROGRESSIVE
  set_target_list_size(list_size_k3);
#endif

  while (!check_if_done())
  {
#ifdef PROGRESSIVE
    if ((progressive_rank < lattice_rank) && check_if_enough_short_vectors())
    {
      increase_progressive_rank();
    }
#endif
    // pop p from  main_queue
    typename SieveTraits::FastAccess_Point p = main_queue.true_pop();

    // checks if p participates in 2-reduction, inserts p into main_list
    sieve_3_iteration(p);

    ++i;
    if ((i % SieveTraits::print_step_3sieve == 0) && (verbosity >= 2))
    {
      std::cout << "[" << i << "]"
                << " |L|=" << statistics.get_current_list_size()
                << " |Q|=" << main_queue.size()  // STAT_MARK
                << " #samples = " << statistics.get_number_of_points_sampled()
                << " |sv|= " << get_best_length2() << std::endl;
    }
  }
}
*/

}  // end namespace

// TODO: #include placement
#endif
