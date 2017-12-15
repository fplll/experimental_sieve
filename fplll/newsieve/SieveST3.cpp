/*
 MAIN ROUTINES FOR 3-GAUSS SIEVE
 */

namespace GaussSieve
{

/*
  main 3-sieve iteration
  for the input point p, we check if ||p +/- x1 +/- x2|| < max {||p||, ||x1||, ||x2||}
  for all x1, x2 from main_list. In case, max is ||p||, p is modified and the iteration
  is re-started with this new p (via goto). If max is x1 or x2, the max is modified,
  removed from main_list and pushed into main_queue.
  This routines also checks for collisions (0-vector after any of the modifications).

*/

template <class SieveTraits>
void Sieve<SieveTraits, false>::sieve_3_iteration(typename SieveTraits::FastAccess_Point &p)
{
  using std::abs;
  using std::round;

  // TODO: may be make  filtered_list a member of the Sieve class

  static FilteredListType filtered_list;
  filtered_list.reserve(SieveTraits::filtered_list_size_max);

start_over:
  // it_comparison_flip stores the point where the list elements become larger than p
  auto it_comparison_flip = main_list.cend();
  filtered_list.clear();

  // ||p|| >= ||x1|| >= ||x2||
  for (auto it_x1 = main_list.cbegin(); it_x1 != main_list.cend(); ++it_x1)
  {
    if (p < (*it_x1))  // TODO: use approx_norm2_p (may be)
    {
      it_comparison_flip = it_x1;
      break;  // we proceed to the case where p is no longer the largest point of the triple
    }

    // if the sim_hash - scalar product between p and it_x1 is bad, don't bother with this x1:
    if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
            p, it_x1, SieveTraits::threshold_lvls_3sieve_lb_out,
            SieveTraits::threshold_lvls_3sieve_ub_out))
    {
      continue;
    }

    LengthType const sc_prod_px1 = compute_sc_product(p, *it_x1);
    statistics.increment_number_of_scprods_level1();

    // check for 2-reduction. We already computed the (exact) scalar product anyway.
    LengthType const abs_sc_prod_px1 = abs(sc_prod_px1);

    // cond_x1 = 2*|<p,x_1>| - ||x_1||^2.
    // If >0, we can perform 2-reduction. It is useful to keep this term, since we will re-use it
    // in the exact check for 3-reduction
    LengthType const cond_x1 = 2 * abs_sc_prod_px1 - it_x1->get_norm2();

    if (cond_x1 > 0)  // We can perform 2-reduction, changing p:
    {
      double const mult = convert_to_double(sc_prod_px1) / convert_to_double(it_x1->get_norm2());
      int const scalar  = round(mult);
      p.sub_multiply(*it_x1, scalar);
      if (p.is_zero())  // might move this to after start_over.
      {
        statistics.increment_number_of_collisions();
        return;
      }
      // we start the current iteration all over:
      // This should be faster than main_queue.push(std::move(p)); return;
      statistics.increment_number_of_2reds();
      p.update_bitapprox();
      goto start_over;
    }

    // could not perform 2-reduction, but possibly 3-reduction
    // compute scaled inner-product: <p, x1> / ( ||p||^2 * ||x1||^2)
    // the result should always be positive
    double const sc_prod_px1_normalized =
        convert_to_double(sc_prod_px1) * convert_to_double(sc_prod_px1) /
        (convert_to_double(p.get_norm2()) * convert_to_double(it_x1->get_norm2()));

    // If the scalar product is too small, we cannot perform 3-reduction, so we take the next x1
    if (sc_prod_px1_normalized < SieveTraits::x1x2_target)
    {
      continue;  // for loop over it_x1;
    }

    // From here : x1 is a candidate for 3-reduction and will eventually be put into filtered_list.
    //             To avoid checking the triple (p, x1, x1), we only append to filtered_list after
    //             we iterate over candidates for x2.
    bool const sign_px1 = (sc_prod_px1 > 0);
    for (auto &filtp_x2 : filtered_list)  // Note that we know ||p|| >= ||*it_x1|| >= ||x2||
    {

      // (approximate) check if x1 and x2 are close enough to participate in 3-reduction
      // we can remove it and perform the exact check only, however, it seems to speed-up the search
      if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
              it_x1, filtp_x2.sim_hashes, SieveTraits::threshold_lvls_3sieve_lb_inn,
              SieveTraits::threshold_lvls_3sieve_ub_inn))
      {
        continue;
      }
      LengthType sc_prod_x1x2 = (filtp_x2.sign_flip == sign_px1)
                                    ? compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact))
                                    : -compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact));
      statistics.increment_number_of_scprods_level2();
      //      cond_x1 == 2|<p,x_1>| - ||x_1||^2.
      // filp_x2.cond == 2|<p,x_2>| - ||x_2||^2, which we computed and stored in a previous
      // iteration. The condition is equivalent to ||p+/-x_1 +/- x_2||^2 < ||p||^2,
      // where we have a minus sign @x1 iff sign_px1 == true
      // and                        @x2 iff filtp_x2.sign_flip == true
      if (2 * sc_prod_x1x2 < cond_x1 + filtp_x2.cond)  // perform 3-reduction:
      {
        // LengthType const debug_test = p.get_norm2();
        if (sign_px1)
        {
          p -= *it_x1;
        }
        else
        {
          p += *it_x1;
        }
        if (filtp_x2.sign_flip)
        {
          p -= *(filtp_x2.ptr_to_exact);
        }
        else
        {
          p += *(filtp_x2.ptr_to_exact);
        }
        // assert(p.get_norm2() <= debug_test);  // make sure we are making progress.
        if (p.is_zero())
        {
          statistics.increment_number_of_collisions();
          return;
        }
        p.update_bitapprox();
        statistics.increment_number_of_3reds();
        goto start_over;
      }
    }
    filtered_list.emplace_back(it_x1, sign_px1, cond_x1);
  }  // end of first part of for-loop where p is the largest of the triple

  /**********

  **********/

  // p no longer changes now. it_comparison_flip is iterator to first (shortest) element in the list
  // that is longer than p. If no such element exists, it_comparison_flip refers to after-the-end.
  for (auto it_x1 = it_comparison_flip; it_x1 != main_list.cend();)  // ++it inside loop body.
  {
    // if <p,x1> is bad, don't bother with x1
    if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
            p, it_x1, SieveTraits::threshold_lvls_3sieve_lb_out,
            SieveTraits::threshold_lvls_3sieve_ub_out))
    {
      ++it_x1;
      continue;
    }
    // x1 is somewhat promising according to the above check
    LengthType const sc_prod_px1 = compute_sc_product(p, *it_x1);
    statistics.increment_number_of_scprods_level1();
    bool const sign_px1 = (sc_prod_px1 > 0);

    // check for 2-reduction. We already computed the (exact) scalar product anyway.
    LengthType twice_abs_sc_prod_px1 = 2 * abs(sc_prod_px1);
    // cond_x1_store = 2*|<p,x_1>| - ||x_1||^2.
    // cond_x1_p     = 2*|<p,x_1>| - ||p||^2.
    LengthType const cond_x1_p = twice_abs_sc_prod_px1 - p.get_norm2();
    if (cond_x1_p > 0)  // In this case, we can perform 2-reduction, changing x1
    {
      double const mult = convert_to_double(sc_prod_px1) / convert_to_double(p.get_norm2());
      int const scalar  = round(mult);
      assert(scalar != 0);

      // true_pop_point erases it, sio the place to insert p is incremented to avoid segfaults
      if (it_x1 == it_comparison_flip)
      {
        ++it_comparison_flip;
      }

      auto v_new = main_list.true_pop_point(it_x1);  // also performs ++it_x1 !
      v_new.sub_multiply(p, scalar);
      statistics.increment_number_of_2reds();
      if (v_new.is_zero())
      {
        statistics.increment_number_of_collisions();
      }
      else
      {
        main_queue.push(std::move(v_new));
      }
      continue;  // with next it_x1
    }

    // no 2-reduction possible, consider 3-reductions
    // compute scaled inner-product: <p, x1> / ( ||p||^2 * ||x1||^2)
    // the result should always be positive
    double const sc_prod_px1_normalized =
        convert_to_double(sc_prod_px1) * convert_to_double(sc_prod_px1) /
        (convert_to_double(p.get_norm2()) * convert_to_double(it_x1->get_norm2()));

    // If the scalar product is too small, 3-reduction is unlikely to happen, so we take the next x1
    if (sc_prod_px1_normalized < SieveTraits::x1x2_target)
    {
      ++it_x1;
      continue;  // for loop over it_x1;
    }

    // From here : x1 is a candidate for 3-reduction and will eventually be put into filtered_list.
    //             To avoid checking the triple (p, x1, x1), we only append to filtered_list after
    //             we iterate over candidates for x2.
    for (auto const &filtp_x2 : filtered_list)
    {

      if (!check_simhash_scalar_product<typename SieveTraits::CoordinateSelectionUsed>(
              it_x1, filtp_x2.sim_hashes, SieveTraits::threshold_lvls_3sieve_lb_inn,
              SieveTraits::threshold_lvls_3sieve_ub_inn))
      {
        continue;
      }

      // Note that we know ||p|| < ||x1|| and ||x1|| >= ||x2||, so x1 is the maximum.
      LengthType sc_prod_x1x2 = (filtp_x2.sign_flip == sign_px1)
                                    ? compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact))
                                    : -compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact));
      statistics.increment_number_of_scprods_level2();
      // The correct condition here has cond_x1_p = 2|<p,x_1>| - ||p||^2.
      // This differs from the case above, because now x_1 is larger than p.
      if (2 * sc_prod_x1x2 < cond_x1_p + filtp_x2.cond)  // perform 3-reduction on x1
      {
        statistics.increment_number_of_3reds();
        // LengthType const debug_test = it_x1->get_norm2();
        if (it_x1 == it_comparison_flip)
        {
          ++it_comparison_flip;
        }
        auto v_new = main_list.true_pop_point(it_x1);  // also performs ++it_x1 !
        // Note: sign_px1 says whether we need to change x1 (i.e.
        //       we need to consider p+/- v_new. We instead flip the global sign
        //       and look at v_new +/- p
        if (sign_px1)
        {
          v_new -= p;
        }
        else
        {
          v_new += p;
        }
        // If sign_px1 == true, we need to invert the sign of x2 because of the global sign flip.
        if (filtp_x2.sign_flip != sign_px1)
        {
          v_new -= *(filtp_x2.ptr_to_exact);
        }
        else
        {
          v_new += *(filtp_x2.ptr_to_exact);
        }
        // assert(v_new.get_norm2() < debug_test);  // make sure we are making progress.
        if (v_new.is_zero())
        {
          statistics.increment_number_of_collisions();
        }
        else
        {
          main_queue.push(std::move(v_new));
        }
        // we break the inner loop over filtp_x2 and continue the it_x1-loop here
        goto end_of_x1_loop;
      }
      // No 3-reduction for this x2
    }
    // No 3-reduction for for this x1 for any x2
    twice_abs_sc_prod_px1 -= it_x1->get_norm2();
    filtered_list.emplace_back(it_x1, sign_px1, std::move(twice_abs_sc_prod_px1));
    ++it_x1;
  end_of_x1_loop:;
  }  // end of second part of for loop

  // put p into the main_list
  assert(!(p.is_zero()));
  if (update_shortest_vector_found(p))
  {
    if (verbosity >= 2)
    {
      std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
    }
  }

  // could be done between the two for-loops, but it would require making a copy of p
  main_list.insert_before(it_comparison_flip, std::move(p));
}

}  // namespace GaussSieve
