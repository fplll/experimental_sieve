// TODO: Remove define
#define filtered_list_size_max  1000

namespace GaussSieve{

/*
template<class SieveTraits>
bool Sieve<SieveTraits,false>::check3red_approx(typename SieveTraits::FastAccess_Point const &p1,
                typename SieveTraits::FastAccess_Point const &p2)

{

  uint_fast32_t approx_scprod = 0;

  //std::cout<< approx_scprod << std::endl;
  for (unsigned int lvl = 0; lvl<SimHash::num_of_levels; ++lvl)
  {
    statistics.increment_number_of_approx_scprods_level1();
    approx_scprod += static_cast<uint_fast32_t> (compute_sc_product_bitapprox_level(p1, p2, lvl));
    if (approx_scprod >= SimHash::sim_hash_len/2 + SimHash::threshold_lvls_3sieve[lvl] ||
        approx_scprod <= SimHash::sim_hash_len/2 - SimHash::threshold_lvls_3sieve[lvl] )
        {
          continue;
            //approx_scprod+=static_cast<uint_fast32_t>(compute_sc_product_bitapprox_level(p1, p2, lvl));
        }
    else
    {
      //std::cout<< approx_scprod << std::endl;
      return false;

    }
  }

  return true;
}
*/

template<class SieveTraits>
template<class IT>
bool Sieve<SieveTraits,false>::check_sc_prod_outer (typename SieveTraits::FastAccess_Point const &x1, IT &&x2,
                    typename SieveTraits::LengthType & sc_prod_x1x2)
{


  if (!check_simhash_scalar_product<typename SieveTraits::SimHashGlobalDataType>(
                                            x1, x2,
                                            SieveTraits::threshold_lvls_3sieve_lb_out,
                                            SieveTraits::threshold_lvls_3sieve_ub_out))
  {
    return false;
  }


  using std::abs;


  sc_prod_x1x2 = compute_sc_product(x1, turn_maybe_iterator_to_point(x2));
  statistics.increment_number_of_scprods_level1();


  double sc_prod_px1_norm = convert_to_double( sc_prod_x1x2)*convert_to_double(sc_prod_x1x2 ) /
                  ( convert_to_double (x1.get_norm2()) *
                  convert_to_double( turn_maybe_iterator_to_point(x2).get_norm2() )) ;

  if (abs(sc_prod_px1_norm) > SieveTraits::x1x2_target)
    return true;

  return false;
}

template<class SieveTraits>
template<class IT>
bool Sieve<SieveTraits,false>::check_sc_prod_inner (Filtered_Point const &x1, IT &&x2,
                    typename SieveTraits::LengthType & sc_prod_x1x2)
{

  //TODO: make the approx check
  // currently: type mismatch
  /*
  if (!check_simhash_scalar_product<typename SieveTraits::SimHashGlobalDataType>(
                                            *(x1.get_point()), x2,
                                            SieveTraits::threshold_lvls_3sieve_lb_inn,
                                            SieveTraits::threshold_lvls_3sieve_lb_inn))
  {
    return false;
  }
  */

  using std::abs;

  sc_prod_x1x2 = compute_sc_product( *(x1.get_point().ptr_to_exact), turn_maybe_iterator_to_point(x2));
  statistics.increment_number_of_scprods_level1();

  double x1_norm2 = convert_to_double( (x1.get_point().ptr_to_exact)->get_norm2() );
  double x2_norm2 = convert_to_double( turn_maybe_iterator_to_point(x2).get_norm2() );


  double sc_prod_px1_norm = convert_to_double( sc_prod_x1x2)*convert_to_double(sc_prod_x1x2 ) /
                  ( x1_norm2 * x2_norm2) ;

  if (abs(sc_prod_px1_norm) > SieveTraits::x2x3_target)
    return true;

  return false;
}

    // The function checks if
    //    || x1 \pm x2 \pm x3 || < || x1 ||
    // The first arguments is assumed to have the largest norm
    // The last two arguments are modified. They return the correct signs, i.e.
    // x_new = x1 + sgn1*x2 + x3*sng2
    // p_is_max is true if x1 ==p, in which case x3_X stores <x3,x2>
    // otherwise x3_X stores <x3, x1>
    //  EXACT CHECKS ONLY
template<class SieveTraits>
template<class ARG1, class ARG2>
bool Sieve<SieveTraits,false>::check_triple (ARG1 &&x1, ARG2 &&x2,
                 Filtered_Point const &x3,
                 typename SieveTraits::LengthType const &x1x2,
                 typename SieveTraits::LengthType const &x3_X,
                 int &sgn2, int &sgn3, bool p_is_max)


{
  //retrieve the signs of the 3 inner-products;
  //they should satisfy (<p, x1>*<p, x2>*<x1, x2>) < 0
  //check for all 4 combinations that do not lead to a reduction
  typename SieveTraits::LengthType x3x1;
  typename SieveTraits::LengthType x3x2;
  using std::abs;

  if (x1x2 * x3_X * x3.get_sc_prod() >0)
  {
    return false;
  }
  if (!p_is_max) {
    x3x1 = x3_X;
    x3x2 = x3.get_sc_prod();
  }
  else
  {
    x3x1 = x3.get_sc_prod();
    x3x2 = x3_X;
  }

  LengthType x1_norm2 = turn_maybe_iterator_to_point(x1).get_norm2();
  LengthType x2_norm2 = turn_maybe_iterator_to_point(x2).get_norm2();
  LengthType x3_norm2 = (x3.get_point().ptr_to_exact)->get_norm2();

  if (x3x1 <0 && x1x2<0 && x3x2<0 &&
      x2_norm2 + x3_norm2 < 2 * ( abs(x3x1 + x1x2 + x3x2) ) )
  {

    sgn2 = 1;
    sgn3 = 1;
        //std::cout << "sgns: case 1" << std::endl;
    return true;
  }

    /* bool f = ((x ^ y) < 0); // true iff x and y have opposite signs*/

    //if (x3x1 <0 && !((x1x2^x3x2) <0) &&
  if (x3x1 <0 && x1x2>0 && x3x2 >0 &&
      x2_norm2 + x3_norm2 < 2 * ( -x3x1 + x1x2 + x3x2 ) )
  {

    sgn2 = -1;
    sgn3 = 1;
        //std::cout << "sgns: case 2" << std::endl;
    return true;

  }

  if (x3x1 >0 && x1x2<0 && x3x2>0 &&
      x2_norm2 + x3_norm2 < 2 * ( x3x1 - x1x2 + x3x2 ) )
  {
    sgn2 = 1;
    sgn3 = -1;
    //std::cout << "sgns: case 3" << std::endl;
    return true;
  }

  if (x3x1 >0 && x1x2>0 && x3x2<0 &&
      x2_norm2 + x3_norm2 < 2 * (  x3x1 + x1x2 - x3x2 ) )
  {

    sgn2 = -1;
    sgn3 = -1;
        //std::cout << "sgns: case 4" << std::endl;
    return true;
  }

  return false;
}


template<class SieveTraits>
void Sieve<SieveTraits,false>::sieve_3_iteration (typename SieveTraits::FastAccess_Point &p)
{
//  std::cout << "Starting Iteration" << std::endl << std::flush;
  using std::abs;
  using std::round;

  // TODO: Reuse memory in different iterations. Consider making this static or a member
  // of the Sieve class

  thread_local static FilteredListType filtered_list;
  filtered_list.reserve(filtered_list_size_max);

start_over:
  auto it_comparison_flip = main_list.cend(); //to store the point where the list elements become larger than p.
  filtered_list.clear();

//  auto it = main_list.cbegin();
  double approx_norm2_p = convert_to_double(p.get_norm2());

  // this part of the outer loop is for the case where p is the largest of the triple
  for (auto it_x1 = main_list.cbegin(); it_x1 != main_list.cend(); ++it_x1)
  {
    if (p  < (*it_x1) )  // TODO: use approx_norm2_p
    {
      it_comparison_flip = it_x1;
      break; // we proceed to the case where p is no longer the largest point of the triple
    }

    // if the sim_hash - scalar product between p and it_x1 is bad, don't bother with this x1:
    if (!check_simhash_scalar_product<typename SieveTraits::SimHashGlobalDataType>(
                                                 p, it_x1,
                                                 SieveTraits::threshold_lvls_3sieve_lb_out,
                                                 SieveTraits::threshold_lvls_3sieve_ub_out))
    {
      continue;
    }
    // TODO: Check whether that computation is actually necessary.
    LengthType const sc_prod_px1 = compute_sc_product(p, *it_x1);
    bool const sign_px1 = (sc_prod_px1 > 0);

    // check for 2-reduction. We already computed the (exact) scalar product anyway.
    LengthType const abs_sc_prod_px1 = abs(sc_prod_px1);
    // cond_x1 = 2*|<p,x_1>| - ||x_1||^2.
    // If >0, we can perform 2-reduction. It is useful to keep this term, since we will re-use it.
    LengthType const cond_x1 = 2*abs_sc_prod_px1 - it_x1->get_norm2();
    if (cond_x1 > 0 )  // We can perform 2-reduction, changing p:
    {
      double const mult = convert_to_double(sc_prod_px1) / convert_to_double(it_x1->get_norm2());
      int const scalar = round(mult);
      p.sub_multiply(*it_x1, scalar);
      if (p.is_zero())  // might move this to after start_over.
      {
        statistics.increment_number_of_collisions();
        return;
      }
      // we changed p and have to start the current iteration all over:
      // Note that this is faster than main_queue.push(std::move(p)); return;
      p.update_bitapprox();
      goto start_over;
    }
    // could not perform 2-reduction, but possibly 3-reduction:
    double const sc_prod_px1_normalized = convert_to_double(sc_prod_px1) * convert_to_double(sc_prod_px1)
                  / ( convert_to_double(p.get_norm2()) * convert_to_double(it_x1->get_norm2()) );
    assert(sc_prod_px1_normalized >= 0); //  old code had an abs here.

    // If the scalar product is too small, we cannot perform 3-reduction, so we take the next x1
    if (sc_prod_px1_normalized < SieveTraits::x1x2_target)
    {
      continue;  // for loop over it_x1;
    }
    // From here : x1 is a candidate for 3-reduction and will eventually be put into filtered_list.
    //             To avoid checking the triple (p, x1, x1), we only append to filtered_list after
    //             we iterate over candidates for x2.
    for (auto & filtp_x2 : filtered_list) // Note that we know ||p|| >= ||*it_x1|| >= ||x2||
    {
      // Note: We do not check approximately, just like the old code below.
      LengthType sc_prod_x1x2 = (filtp_x2.sign_flip==sign_px1)
                                    ?  compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact))
                                    : -compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact));
      // Recall cond_x1 == 2|<p,x_1> - ||x_1||^2.
      // filp_x2.cond == 2|<p,x_2>| - ||x_2||^2, which we computed and stored in a previous
      // iteration. The condition is equivalent to ||p+/-x_1 +/- x_2||^2 < ||p||^2,
      // where we have a minus sign @x1 iff sign_px1 == true
      // and                        @x2 iff filtp_x2.sign_flip == true
      if (sc_prod_x1x2 < cond_x1 + filtp_x2.cond)  // perform 3-reduction:
      {
        LengthType const debug_test = p.get_norm2();
        if (sign_px1)           { p-=*it_x1; }
        else                    { p+=*it_x1; }
        if (filtp_x2.sign_flip) { p-=*(filtp_x2.ptr_to_exact); }
        else                    { p+=*(filtp_x2.ptr_to_exact); }
        assert(p.get_norm2() <= debug_test);  // make sure we are making progress.
        if (p.is_zero())  { statistics.increment_number_of_collisions(); return; }
        p.update_bitapprox();
        goto start_over;
      }
    }
    filtered_list.emplace_back( it_x1, sign_px1, cond_x1 );
  }  // end of first part of for-loop where p is the largest of the triple

/**********

**********/

//  std::cout << "Entering part 2" << std::endl << std::flush;

  // p no longer changes now. it_comparison_flip is iterator to first (shortest) element in the list
  // that is longer than p. If no such element exists, it_comparison_flip refers to after-the-end.

  for( auto it_x1 = it_comparison_flip; it_x1 != main_list.cend(); )  // ++it inside loop body.
  {
    // if <p,x1> is bad, don't bother with x1
    if (!check_simhash_scalar_product<typename SieveTraits::SimHashGlobalDataType>(
                                                 p, it_x1,
                                                 SieveTraits::threshold_lvls_3sieve_lb_out,
                                                 SieveTraits::threshold_lvls_3sieve_ub_out))
    {
      ++it_x1;
      continue;
    }
    // x1 is at least somewhat promising:
    LengthType const sc_prod_px1 = compute_sc_product(p, *it_x1);
    bool const sign_px1 = (sc_prod_px1 > 0);

    // check for 2-reduction. We already computed the (exact) scalar product anyway.
    LengthType twice_abs_sc_prod_px1 = 2*abs(sc_prod_px1);
    // cond_x1_store = 2*|<p,x_1>| - ||x_1||^2.
    // cond_x1_p     = 2*|<p,x_1>| - ||p||^2.
    LengthType const cond_x1_p     = twice_abs_sc_prod_px1 - p.get_norm2();
    if (cond_x1_p > 0)  // In this case, we can perform 2-reduction, changing x1
    {
      double const mult = convert_to_double(sc_prod_px1) / convert_to_double(p.get_norm2());
      int const scalar = round(mult);
      assert(scalar!=0);
      if (it_x1 == it_comparison_flip) { ++it_comparison_flip; }
      auto v_new = main_list.true_pop_point(it_x1);  // also performs ++it_x1 !
      v_new.sub_multiply(p, scalar);
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
    // no 2-reduction possible, consider 3-reductions:
    double const sc_prod_px1_normalized = convert_to_double(sc_prod_px1) * convert_to_double(sc_prod_px1)
                  / ( convert_to_double(p.get_norm2()) * convert_to_double(it_x1->get_norm2()) );
    assert(sc_prod_px1_normalized >= 0); //  old code had an abs here.
    // If the scalar product is too small, we cannot perform 3-reduction, so we take the next x1
    if (sc_prod_px1_normalized < SieveTraits::x1x2_target)
    {
      ++it_x1;
      continue;  // for loop over it_x1;
    }
    // From here : x1 is a candidate for 3-reduction and will eventually be put into filtered_list.
    //             To avoid checking the triple (p, x1, x1), we only append to filtered_list after
    //             we iterate over candidates for x2.
    for (auto const & filtp_x2 : filtered_list)
    {

      // Note that we know ||p|| < ||x1|| and ||x1|| >= ||x2||, so x1 is the maximum.
      LengthType sc_prod_x1x2 = (filtp_x2.sign_flip == sign_px1)
                                    ?  compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact))
                                    : -compute_sc_product(*it_x1, *(filtp_x2.ptr_to_exact));
      // Note that the correct condition here has cond_x1_p = 2|<p,x_1>| - ||p||^2.
      // This differs from the case above, because now x_1 is larger than p.
      if (sc_prod_x1x2 < cond_x1_p + filtp_x2.cond)  // perform 3-reduction on x1
      {
        LengthType const debug_test = it_x1->get_norm2();
        if (it_x1 == it_comparison_flip) { ++it_comparison_flip; }
        auto v_new = main_list.true_pop_point(it_x1);  // also performs ++it_x1 !
        // Note: sign_px1 says whether we need to change x1 (i.e.
        //       we need to consider p+/- v_new. We instead flip the global sign
        //       and look at v_new +/- p
        if (sign_px1)           { v_new-=p; }
        else                    { v_new+=p; }
        // If sign_px1 == true, we need to invert the sign of x2 because of the global sign flip.
        if (filtp_x2.sign_flip != sign_px1) { v_new-=*(filtp_x2.ptr_to_exact); }
        else                                { v_new+=*(filtp_x2.ptr_to_exact); }
        assert(v_new.get_norm2() < debug_test);  // make sure we are making progress.
        if (p.is_zero())
        {
          statistics.increment_number_of_collisions();
        }
        else
        {
          main_queue.push(std::move(v_new));
        }
        // we need to break the inner loop over filtp_x2 and continue the it_x1-loop here
        goto end_of_x1_loop;
      }
      // No 3-reduction for this x2
    }
    // No 3-reduction for for this x1 for any x2
    twice_abs_sc_prod_px1 -= it_x1 -> get_norm2();
    // twice_abs_sc_prod_px1 now holds 2|<p,x_1> - ||x_1||^2
    filtered_list.emplace_back( it_x1, sign_px1, std::move(twice_abs_sc_prod_px1) );
    ++it_x1;
    end_of_x1_loop: ;
  }  // end of second part of for loop

//  std::cout << "End of both loops" << std::endl << std::flush;

  // put p into the main_list, delayed until now to use std::move()

  assert(!(p.is_zero()));
  if(update_shortest_vector_found(p))
  {
    if(verbosity>=2)
    {
      std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
    }
  }
  main_list.insert_before(it_comparison_flip,std::move(p));



/*
  //if (p.is_zero() )
  //{
  //  return; //TODO: Ensure sampler does not output 0 (currently, it happens).
  //}

  auto it_comparison_flip = main_list.cend(); //to store the point where the list elements become larger than p.
 thread_local static FilteredListType filtered_list;

  //double px1_target  = .1111; // TO ADJUST

  int scalar = 0; //for 2-reduction

  filtered_list.clear();
  filtered_list.reserve(filtered_list_size_max);

  auto it = main_list.cbegin();
  double approx_norm2_p = convert_to_double(p.get_norm2());

  while (it != main_list.cend()) // TODO: Turn into for loop over it for clarity
  {
    if (p  < (*it) )  // TODO: use approx_norm2_p
    {
      it_comparison_flip = it;
      break;
    }

    //
    //check for 2-reduction
    //

    // TODO: This computes bitapprox_scalar product TWICE

    if ( check2red(p, it, scalar) )
    {
      assert(scalar!=0);
      p-= (*it) * scalar;
      approx_norm2_p = convert_to_double(p.get_norm2());
      p.update_bitapprox();

      if (p.is_zero() )
      {
        statistics.increment_number_of_collisions();
      }
      else
      {
        main_queue.push(std::move(p)); // TODO: Change to goto start (which is faster)
      }
      return;
    }

//    statistics.set_filtered_list_size(0);


    LengthType sc_prod_px1;
    if (check_sc_prod_outer(p, it, sc_prod_px1))
    {
      for (auto & filtered_list_point: filtered_list)
      {

        LengthType sc_prod_x1x2;
        if (check_sc_prod_inner(filtered_list_point, it, sc_prod_x1x2))
        {
          int sgn2 = 1;
          int sgn3 = 1;
          //check if || p \pm x1 \pm x2 || < || p ||
          // ! check_triple assumes that the first argument has the largest norm
          if ( check_triple(p, it, filtered_list_point, sc_prod_px1, sc_prod_x1x2, sgn2, sgn3, true) )
          {
            //LengthType pnorm_old = p.get_norm2();
            //TODO:
            p += (*it)*sgn2;
            p += *(filtered_list_point.get_point().ptr_to_exact)* sgn3;
            p.update_bitapprox();
            if (p.is_zero() )
            {
              statistics.increment_number_of_collisions();
            }
            else
            {
              main_queue.push(std::move(p));
            }
            return;
          } //if(check_triple())
        }
      } // end for loop for inner point
      //typename SieveTraits::FlilteredPointType new_filtered_point((*it).make_copy(), sc_prod_px1);
      Filtered_Point new_filtered_point( &(*(it.true_star())), sc_prod_px1);
      filtered_list.push_back(std::move(new_filtered_point));

    } //if (check_sc_prod(p, *it, sc_prod_px1))

    ++it;
  } //while-loop (for-loop over it, p largest point)

*/

/*
  main_list.insert_before(it_comparison_flip,p.make_copy());

  if(update_shortest_vector_found(p))
  {
    if(verbosity>=2)
    {
      std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
    }
  }
*/

/*
  //now p is not the largest
  //it_comparison_flip points to the next after p list-element
  it = it_comparison_flip;

  while (it !=main_list.cend())
  {

    //if true, do not put into the filtered_list
    bool x1_reduced = false;


    //
    //check for 2-reduction
    //

    if ( check2red(it, p, scalar) )
    {
      assert(scalar!=0);  // should not be 0 in any case
      typename SieveTraits::FastAccess_Point v_new = (*it) - (p*scalar);


      if (v_new.is_zero() )
      {
        statistics.increment_number_of_collisions();
      }

      main_queue.push(std::move(v_new));
      it = main_list.erase(it);
      continue;
    }

    //
    // 3-reduction
    //
    // Now x1 is the largest

    LengthType sc_prod_px1;
    if (check_sc_prod_outer(p, it, sc_prod_px1))
    {
      for (auto & filtered_list_point: filtered_list)
      {
        LengthType sc_prod_x1x2;
        if (check_sc_prod_inner(filtered_list_point, it, sc_prod_x1x2))
        {
          int  sgn2 = 1;
          int  sgn3 = 1;
          // ! check_triple assumes that the first argument has the largest norm
          if ( check_triple( it, p, filtered_list_point, sc_prod_px1, sc_prod_x1x2, sgn2, sgn3, false) )
          {
            typename SieveTraits::FastAccess_Point v_new =(*it) + p*sgn2 + *(filtered_list_point.get_point().ptr_to_exact) * sgn3;
            if (v_new.is_zero() )
            {
              statistics.increment_number_of_collisions();
            }
            main_queue.push(std::move(v_new));
            it = main_list.erase(it);
            x1_reduced = true;
            break; //for-loop over the filtered_list
          }
        }
      } //for-loop

      if (!x1_reduced)
      {
        //typename SieveTraits::FlilteredPointType new_filtered_point((*it).make_copy(), sc_prod_px1);
        Filtered_Point new_filtered_point( &(*(it.true_star())), sc_prod_px1);
        filtered_list.push_back(std::move(new_filtered_point));
      }


    } //if (check_sc_prod_outer(p, it, sc_prod_px1))



    if (!x1_reduced)
    {
      ++it;
    }
//        if (filtered_list.size()>0) {
//            std::cout << "filtered.size() = " << filtered_list.size() << std::endl;
//        }

  } // 'lower' while-loop

//  statistics.set_filtered_list_size(filtered_list.size());
  filtered_list.clear();
*/

}


} //namespace GaussSieve


/*Unused*/


//The function checks if ||x1 + scalar* x2|| < ||x1||
// the value <x1,x2> is provided as input
// scalar is modified
template<class SieveTraits, class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
bool check_2red_with_scprod (typename SieveTraits::FastAccess_Point const &x1,
                             typename SieveTraits::FastAccess_Point const &x2,
                             typename SieveTraits::LengthType const &x1x2, Integer & scalar)
{
  using std::abs;
  using std::round;

  typename SieveTraits::LengthType const abs_2scprod = abs(x1x2 * 2);

  // check if |2 * <x1, x2>| <= |x2|^2. If yes, no reduction
  if (abs_2scprod <= x2.get_norm2())
  {
    return false;
  }

  //compute the multiple mult s.t. res = x1 \pm mult* x2;
  double const mult = convert_to_double( x1x2 ) / convert_to_double( x2.get_norm2() );
  // TODO: Check over- / underflows.
  scalar =  round (mult);
  return true;
}

/*
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
      unsigned int lvl = 0;
      uint_fast32_t approx_scprod = static_cast<uint_fast32_t>(compute_sc_product_bitapprox_level(p, *it, lvl));
      statistics.no_red_stat_innloop[lvl][static_cast<uint_fast32_t>(approx_scprod)]++;
#endif
 */

// OLD CODE
       //Collecting stats
          /*
  #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
          double sc_prod_x1x2_norm = convert_to_double( sc_prod_x1x2 )* convert_to_double( sc_prod_x1x2 )  /
          ( convert_to_double ( filtered_list_point.get_point().get_norm2()) * convert_to_double( (*it).get_norm2() )) ;
          if (std::abs(sc_prod_x1x2_norm) > px1_target)
          {

            uint_fast32_t approx_scprod = static_cast<uint_fast32_t> (compute_sc_product_bitapprox_level(filtered_list_point.get_point(), *it, 0));
            statistics.red_stat_innloop[0][approx_scprod]++;
            for (unsigned int lvl = 1; lvl<SimHash::num_of_levels; ++lvl)
            {
              if (approx_scprod >= SimHash::sim_hash_len/2 + SimHash::threshold_lvls_3sieve[lvl-1] ||
                  approx_scprod <= SimHash::sim_hash_len/2 - SimHash::threshold_lvls_3sieve[lvl-1] )
              {
                approx_scprod+=static_cast<uint_fast32_t>(compute_sc_product_bitapprox_level(filtered_list_point.get_point(), *it, lvl));
                statistics.red_stat_innloop[lvl][approx_scprod]++;
              }
            }

          }
          else
          {
            uint_fast32_t approx_scprod = static_cast<uint_fast32_t> (compute_sc_product_bitapprox_level(filtered_list_point.get_point(), *it, 0));
            statistics.no_red_stat_innloop[0][approx_scprod]++;
            for (unsigned int lvl = 1; lvl<SimHash::num_of_levels; ++lvl)
            {
              if (approx_scprod >= SimHash::sim_hash_len/2 + SimHash::threshold_lvls_3sieve[lvl-1] ||
                  approx_scprod <= SimHash::sim_hash_len/2 - SimHash::threshold_lvls_3sieve[lvl-1] )
              {
                approx_scprod+=static_cast<uint_fast32_t>(compute_sc_product_bitapprox_level(filtered_list_point.get_point(), *it, lvl));
                statistics.no_red_stat_innloop[lvl][approx_scprod]++;
              }
            }

          }

  #endif
           */
