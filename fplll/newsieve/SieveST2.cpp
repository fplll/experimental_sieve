// clang-format off






namespace GaussSieve{

/*
 Assume ||p1|| > ||p2||
 */

 /**
  Checks whether we can perform a 2-reduction. Modifies scalar.
 */
template<class SieveTraits, class Integer, typename std::enable_if<
  std::is_integral<Integer>::value
  ,int>::type =0 >
bool check2red (typename SieveTraits::FastAccess_Point const &p1,
                typename SieveTraits::FastAccess_Point const &p2,
                Integer & scalar)
{
  //assert(!p2.is_zero()); <-TODO: ERRS NOW
  using EntryType = typename SieveTraits::EntryType;
  using std::abs;
  using std::round;
    /*
    ET sc_prod, abs_2scprod;
    sc_prod= compute_sc_product(p1, p2);
    abs_2scprod.mul_ui(sc_prod,2);
    abs_2scprod.abs(abs_2scprod);
    */

    EntryType const sc_prod = compute_sc_product(p1,p2);
    EntryType const abs_2scprod = abs(sc_prod * 2);

    // check if |2 * <p1, p2>| <= |p2|^2. If yes, no reduction
    if (abs_2scprod <= p2.get_norm2())
        return false;
    //compute the multiple mult s.t. res = p1 \pm mult* p2;
    //mult = round ( <p1, p2> / ||p2||^2 )

    double const mult = convert_to_double( sc_prod ) / convert_to_double( p2.get_norm2() );
    // TODO: Check over- / underflows.
    scalar =  round (mult);
    return true;

}
  /**
   Checks whether we can perform a 2-reduction. Modifies scalar.
   Contrary to the above function, it does not assume that p1 is max, but deduces it from p_is_max
   */
  
template<class SieveTraits, class Integer, typename std::enable_if<
    std::is_integral<Integer>::value
    ,int>::type =0 >
    bool check2red (typename SieveTraits::FastAccess_Point const &p1,
                    typename SieveTraits::FastAccess_Point const &p2,
                    Integer & scalar, bool& p_is_max)
{
  #ifdef USE_APPROXPOINT
    using EntryType = typename GaussSieve::EMVScalar;
    EntryType const sc_prod = compute_sc_product_approx(p1.access_approx(), p2.access_approx());
  #else
    using EntryType = typename SieveTraits::EntryType;
    EntryType const sc_prod = compute_sc_product(p1,p2);
  #endif
  
    using std::abs;
    using std::round;

  
    EntryType const abs_2scprod = abs(sc_prod * 2);


    if (p1.get_norm2() > p2.get_norm2() && abs_2scprod > p2.get_norm2() )
    {
        p_is_max = true;
        double const mult = convert_to_double( sc_prod ) / convert_to_double( p2.get_norm2() );
        scalar =  round (mult);
        return true;
    }


    if (p1.get_norm2() < p2.get_norm2() && abs_2scprod > p1.get_norm2() )
    {
        p_is_max = false;
        double const mult = convert_to_double( sc_prod ) / convert_to_double( p1.get_norm2() );
        scalar =  round (mult);
        return true;
    }


    return false;
}

// Note: scalar changed from EntryType to int.
template<class LatticePoint, class Integer, typename std::enable_if<
  (IsALatticePoint<LatticePoint>::value) && (std::is_integral<Integer>::value)
  ,int>::type =0 >
LatticePoint perform2red (LatticePoint const &p1, LatticePoint const &p2, Integer const scalar)
{
//  typename SieveTraits::FastAccess_Point res;
//  res = scalar_mult<ET,nfixed>(p2, scalar);
//  res = sub(p1, res);
  return p1 - (p2*scalar);
}


template<class SieveTraits> void Sieve<SieveTraits,false>::sieve_2_iteration (typename SieveTraits::FastAccess_Point &p)
{
  if (p.is_zero() ) return; //TODO: Ensure sampler does not output 0 (currently, it happens).
  bool loop = true;

  auto it_comparison_flip=main_list.cend(); //used to store the point where the list elements become larger than p.
//  auto it = main_list.cbegin();

  while (loop) //while p keeps changing
  {
    loop = false;
    for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it)
    {
      if (p  < (*it) )
      {
        it_comparison_flip = it;
        break;
      }

      ++number_of_scprods_level1;
      int scalar;
      if ( check2red<SieveTraits>(p, *it, scalar) )
      {
        assert(scalar!=0);
        p-= (*it) * scalar; //The efficiency can be improved here, but it does not matter, probably.
        loop = true;
        break;
      }

    }

  }

//p no longer changes. it_comparison_flip is iterator to first (shortest) element in the list that is longer than p.
//If no such element exists, it_comparison_flip refers to after-the-end.
    if (p.is_zero() )
    {
      //std::cout << "collision is found " << std::endl;
        number_of_collisions++;
        return;
    }

// Changed, conversion done in insert_before.
/*
    GaussSieve::FastAccess_Point<ET,false,nfixed> p_copy = p.make_copy();

    //convert FastAccess_Point to GaussList_StoredPoint
    GaussSieve::GaussList_StoredPoint<ET, false, nfixed> p_converted (std::move(p_copy));

    main_list.insert_before(it_comparison_flip, std::move(p_converted));
*/

    main_list.insert_before(it_comparison_flip, p.make_copy() );
    ++current_list_size; // TODO: Manage by list and / or guard by DEBUGS.

//    it = it_comparison_flip;

//    while( it!=main_list.cend() )
    for(auto it = it_comparison_flip; it!=main_list.cend(); ) //++it done in body of loop
    {

      ++number_of_scprods_level1;

      int scalar;
      if ( check2red<SieveTraits>(*it, p, scalar) )
      {
//        GaussSieve::FastAccess_Point<ET,false,nfixed> v_new;
//        v_new = perform2red(*it, p, scalar );

// TODO: We can move the collision counting to the very beginning of the first loop.
// and steal *it here for the new point.
// This requires the sampler to never output 0.

        typename SieveTraits::FastAccess_Point v_new = (*it) - (p*scalar);

        //cout << "new v of norm = " << v_new.get_norm2() << endl << flush;

        if (v_new.is_zero() ) // this only happens if the list contains a non-trivial multiple of p.
        {
          //std::cout << "collision on v_new is found " << std::endl;
          number_of_collisions++;
//          ++it;
//          continue;
        }
        main_queue.push(std::move(v_new));

        // This increments the iterator in the sense that its point to the next element now,
        // effectively doubling as a ++it;
        it = main_list.erase(it);
        --current_list_size;

      }
      else // no reduction.
      {
        ++it;
      }
    }


    if(update_shortest_vector_found(p))
     {
         if(verbosity>=2)
         {
             std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
         }
     }
}

#ifdef USE_LSH

    template<class SieveTraits>
    void Sieve<SieveTraits,false>::hash_sieve_2_iteration (typename SieveTraits::FastAccess_Point &p)
    {
      auto constexpr number_of_hash_tables = SieveTraits::number_of_hash_tables;
        start_of_function:
        if (p.is_zero() )
        {
          ++number_of_collisions;
         return; //TODO: Ensure sampler does not output 0 (currently, it happens).
        }
        //std::cout << "p = " << p.get_norm2() << std::endl;


        for (int t=0 ; t < number_of_hash_tables; ++t)
        {

          HashTableType* hash_table = hash_tables.get_ith_hash_table(t);

          int hash_value = hash_table->hash(p);
          //TODO: just stream candidates

          for (auto it = hash_table->candidates(hash_value).cbegin(); it!=hash_table->candidates(hash_value).cend();)
          {

                int scalar;
                bool p_is_max;

                if ( check2red<SieveTraits>(p, it->get_point(), scalar, p_is_max) )
                {
                    assert(scalar!=0);
                    /*
                    if (scalar == 0)
                    {
                        std::cout << p.get_norm2() << " " << it->get_point().get_norm2() << std::endl<<std::flush;
                        assert(false);
                    }
                    */
                    if (p_is_max)
                    {
                        //std::cout << "reduce p" << std::endl << std::flush;
                        p-= it->get_point() * scalar;
                        goto start_of_function;
                    }
                    else
                    {
                        //std::cout << "reduce v" << std::endl << std::flush;
                        typename SieveTraits::FastAccess_Point v_new = it->get_point() - (p*scalar);

                        /*
                        if (v_new.get_norm2() > (it->get_point()).get_norm2())
                        {
                            std::cout << "bug in v_new reduction : " << std::endl << std::flush;
                            std::cout << scalar << " " << p.get_norm2() << " " << it->get_point().get_norm2() << std::endl <<std::flush;
                            assert(false);
                        }
                        */

                        //std::cout<<"v_new = " << &v_new << std::endl;
                        main_queue.push(std::move(v_new));
                        //std::cout<<"v_new is in the queue" << &v_new << std::endl;
                        //std::cout << "about to erase " << &(it->get_point()) << std::endl;


                        it = hash_tables.remove_from_all_hash_tables(&(it->get_point()), t);
                        //std::cout << "finished erase " << std::endl;



                    }
                }
                else
                {
                    ++it;
                }
            }//for-loop over a bucket


        }//iteration over hash-tables

        /*
        if (p.is_zero() )
        {
            number_of_collisions++;
            return;
        }
        */
      hash_tables.add_to_all_hash_tables(p);

        if(update_shortest_vector_found(p))
        {
            if(verbosity>=2)
            {
                std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
            }
        }
    }

#endif // USE_LSH

/*
        OLD IMPLEMENTATION
*/


/*
 if (p.access_exact_point_r().norm2==0) return; //cerr << "Trying to reduce 0"; //TODO: Ensure sampler does not output 0 (currently, it happens).
 ApproximateLatticePoint<ET,nfixed> p_approx (p.access_approximation_r(), n); //makes a copy.
 ET p_exact_norm = p.get_exact_norm2();
 ExactLatticePoint<ET,nfixed> p_exact = p.get_exact_point();

 bool loop = true;

 ET scalar; //reduction multiple output by check2red_new

 typename MainListType::Iterator it_comparison_flip=main_list.cend(); //used to store the point where the list elements become larger than p.

 while (loop) //while p keeps changing
 {
  loop = false;
  for (auto it = main_list.cbegin(); it!=main_list.cend(); ++it)
  {
    if (p_exact_norm < it.get_true_norm2())
    {
      it_comparison_flip = it;
      break;
    }
    ++number_of_scprods;
    bool const predict =GaussSieve::compare_abs_approx_scalar_product(p_approx,*it,it->get_norm2_mantissa(), it->get_norm2_exponent() -1,n); //the -1 acts a a constant factor 1/2, related to our target scalar product.
    //        bool predict = LatticeApproximations::Compare_Sc_Prod(p_approx,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );
    if(!predict) continue;

    ++number_of_exact_scprods;
    if ( GaussSieve::check2red_exact<ET,nfixed>(p_exact, it.dereference_exactly_r(), scalar) )
    {
      p = GaussSieve::perform2red_exact_to_compressed<ET,false,nfixed>(p_exact,it.dereference_exactly_r(),scalar);
      p_approx.replace_by(p.access_approximation_r(), n);
      p_exact  = p.get_exact_point();
      p_exact_norm = p_exact.norm2;
      //p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
      //update the approximation of f
      //if (p.norm2!=0) p_approx = static_cast< ApproxLatticePoint<ET,false> >(p);
      loop = true;
      break;
    }
    else
    {
      ++number_of_mispredictions;
    }
  }
}

 //p no longer changes. it_comparison_flip is iterator to first (shortest) element in the list that is longer than p.
 //If no such element exists, it_comparison_flip refers to after-the-end.

 if (p_exact_norm == 0) //essentially means that p was already inside the list.
	{
 //cout << "p has norm 2" << endl;
 number_of_collisions++;
 return;
	}

 //insert p into main_list;
 main_list.insert_before(it_comparison_flip,p.deep_copy_compressed_point()); //FIXME: This should cause an error without .deep_copy, but doesn't. Something is wrong!
 ++current_list_size;
 if(update_shortest_vector_found(p_exact))
 {
 if(verbosity>=2)
 {
 cout << "New shortest vector found. Norm2 = " << get_best_length2() << endl;
 }
 }

 for(auto it = it_comparison_flip; it!=main_list.cend(); ) //go through rest of the list.
 //We know that every element current_list_point=*it is at least as long as p, so we reduce x using p.
 {
 ++number_of_scprods;
 bool const predict = GaussSieve::compare_abs_approx_scalar_product(p_approx,*it,p_approx.get_norm2_mantissa(),p_approx.get_norm2_exponent()-1,n);
 if(!predict){++it;continue;} //if prediction is bad, don't even bother to reduce.
 //We believe we can reduce *it
 ExactLatticePoint<ET,nfixed> current_list_point = it.dereference_exactly_r();
 ++number_of_exact_scprods;
 if (GaussSieve::check2red_exact(current_list_point, p_exact, scalar)) //We can reduce *it.
 {
 //create new list point
 CompressedPoint<ET,false,nfixed> reduced_point = GaussSieve::perform2red_exact_to_compressed<ET,false,nfixed>(current_list_point, p_exact, scalar);
 //if (!predict) cerr << "Misprediction 2" << endl;
 //cout << "v was found" <<  endl;

 if (reduced_point.get_exact_norm2() == 0) //Note : this cannot happen unless the list contains a non-trivial multiple of p (because the collision would have triggered before).
 {
 number_of_collisions++;
 ++it;
 continue; //was a break. Changed to ++it;continue; -- Gotti
 }

 main_queue.push(std::move(reduced_point));
 it = main_list.erase(it); //this also moves it forward
 --current_list_size;
 }
 else
 {
 ++number_of_mispredictions;
 //	prev = it1;
 ++it;
 }
 }
 */
// print for debugging
//for (it1 = main_list.begin(); it1!=main_list.end(); ++it1) {
//	(*it1).printLatticePoint();
//}

//clang-format on

} // End namespace
