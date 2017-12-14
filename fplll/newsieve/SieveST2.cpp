// clang-format on

/*
  MAIN ROUTINES FOR 2-GAUSS SIEVE
*/

namespace GaussSieve
{

/**
 The function takes two arguents p1, p2, each of which can be either a point (e.g.,FastAccess_Point)
 or an iterator (e.g., over the main_list). This is done becase the function is called twice with
 reversed type of RHS and LSH, thus allowing us to assume ||p1|| >= ||p2||

 The function Checks whether we can perform a 2-reduction between p1, p2. It modifies scalar setting
 it to the value s.t. || p1 - scalar * p2 || < ||p1||
 */

template <class SieveTraits>
template <class LHS, class RHS>
bool Sieve<SieveTraits, false>::check2red(LHS &&p1, RHS &&p2, int &scalar)
{
  // check if SimHash xor of p1 and p2 is promissing
  // do not perform exact check is the SimHash check returns false
  if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
          p1, p2, SieveTraits::threshold_lvls_2sieve_lb, SieveTraits::threshold_lvls_2sieve_ub))
  {
    return false;
  }

  statistics.increment_number_of_scprods_level1();

  using std::round;
  using std::abs;

  using LengthType = typename SieveTraits::LengthType;

  LengthType sc_prod =
      compute_sc_product(turn_maybe_iterator_to_point(p1), turn_maybe_iterator_to_point(p2));

  // if (2|<p1, p1>|-||p2||^2 > 0), can reduce
  LengthType abs_2scprod = abs(sc_prod * 2);

  if (abs_2scprod <= turn_maybe_iterator_to_point(p2).get_norm2())
  {
    return false;
  }

  double const mult =
      convert_to_double(sc_prod) / convert_to_double(turn_maybe_iterator_to_point(p2).get_norm2());
  // TODO: Check over- / underflows.
  scalar = round(mult);
  return true;
}


template <class SieveTraits>
void Sieve<SieveTraits, false>::sieve_2_iteration(typename SieveTraits::FastAccess_Point &p)
{
  using std::abs;
  assert(!p.is_zero());

  // it_comparison_flip stores the point where the list elements become larger than p
  auto it_comparison_flip = main_list.cend();

// go here if p is changed (alternatively, we can p into the queue and pop it in the next run)
start_over:
  
  double approx_norm2 = convert_to_double(p.get_norm2());
  // holds a double-approximation to approx-norm2. Note that the main list stores such an
  // approximation directly inside the nodes. Using these double-approximations should be
  // better in the multithreaded case. For the single-threaded case, it might lead to
  // better cache usage in higher dimensions. Please test whether t makes a difference.

  for (auto it = main_list.cbegin(); it != main_list.cend(); ++it)  // while p keeps changing
  {

    // if (p < *it)
    if (approx_norm2 < it.get_approx_norm2())
    {
      it_comparison_flip = it;
      break;
    }

    int scalar = 0;

    // checks if || p - scalar*it || < ||p||
    // scalar is modified by check2red
    if (check2red(p, it, scalar))
    {
      assert(scalar != 0);
      // p -= (*it) * scalar;
      p.sub_multiply(*it, scalar);
      if (p.is_zero())
      {
        statistics.increment_number_of_collisions();
        return;
      }
      p.update_bitapprox();
      goto start_over;
    }
  }

  // p no longer changes. it_comparison_flip is iterator to first (shortest) element in the list
  // that is longer than p.
  // If no such element exists, it_comparison_flip refers to after-the-end.

  for (auto it = it_comparison_flip; it != main_list.cend();)  // ++it inside body of loop
  {
    int scalar = 0;
    if (check2red(it, p, scalar))
    {

      // TODO: We can move the collision counting to the very beginning of the first loop.
      // and steal *it here for the new point.
      // This requires the sampler to never output 0.

      // true_pop_point erases it, sio the place to insert p is incremented to avoid segfaults
      if (it == it_comparison_flip)
      {
        ++it_comparison_flip;
      }

      auto v_new = main_list.true_pop_point(it);  // implicitly performs ++it.
      // v_new = (*it) - (p * scalar);
      v_new.sub_multiply(p, scalar);

      if (v_new.is_zero())  // this only happens if the list contains a non-trivial multiple of p.
      {
        statistics.increment_number_of_collisions();
      }
      else
      {
        main_queue.push(std::move(v_new));
      }
      continue;  // This increments the iterator in the sense that its point to the next element now
    }
    else  // no reduction.
    {
      ++it;
    }
  }

  if (update_shortest_vector_found(p))
  {
    if (verbosity >= 2)
    {
      std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
    }
  }

  // insert p into the list
  main_list.insert_before(it_comparison_flip, std::move(p));
}

}  // End namespace
