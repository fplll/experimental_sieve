#define filtered_list_size_max  1000

namespace GaussSieve{
  
  
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
  

template<class SieveTraits>
bool Sieve<SieveTraits,false>::check_sc_prod (typename SieveTraits::FastAccess_Point const &x1,
                    typename SieveTraits::FastAccess_Point const &x2,
                    typename SieveTraits::EntryType & sc_prod_x1x2)
{
  //#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
  //  if(!check3red_approx(x1, x2)) return false;
  //#endif
  
  using std::abs;
  
  // ! targets are squared
  //double px1_target = 0.1024;
  double px1_target  = .1111; // TO ADJUST
  //double px1_target = 0.123;
  
  sc_prod_x1x2 = compute_sc_product(x1, x2);
  statistics.increment_number_of_scprods_level1();
  
  double sc_prod_px1_norm = convert_to_double( sc_prod_x1x2)*convert_to_double(sc_prod_x1x2 ) /
                  ( convert_to_double ( x1.get_norm2()) * convert_to_double( x2.get_norm2() )) ;
                  
  if (abs(sc_prod_px1_norm) >px1_target)
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
template<class SieveTraits>
bool check_triple (typename SieveTraits::FastAccess_Point  const &x1,
                 typename SieveTraits::FastAccess_Point const &x2,
                 typename SieveTraits::FlilteredPointType const &x3,
                 typename SieveTraits::EntryType const &x1x2,
                 typename SieveTraits::EntryType const &x3_X,
                 int &sgn2, int &sgn3, bool p_is_max)


{
  //retrieve the signs of the 3 inner-products;
  //they should satisfy (<p, x1>*<p, x2>*<x1, x2>) < 0
  //check for all 4 combinations that do not lead to a reduction
  typename SieveTraits::EntryType x3x1;
  typename SieveTraits::EntryType x3x2;
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


  if (x3x1 <0 && x1x2<0 && x3x2<0 &&
      x2.get_norm2() + (x3.get_point()).get_norm2() <
      2 * ( abs(x3x1 + x1x2 + x3x2) ) )
  {

    sgn2 = 1;
    sgn3 = 1;
        //std::cout << "sgns: case 1" << std::endl;
    return true;
  }

    /* bool f = ((x ^ y) < 0); // true iff x and y have opposite signs*/

    //if (x3x1 <0 && !((x1x2^x3x2) <0) &&
  if (x3x1 <0 && x1x2>0 && x3x2 >0 &&
      x2.get_norm2() + (x3.get_point()).get_norm2() <
      2 * ( -x3x1 + x1x2 + x3x2 ) )
  {

    sgn2 = -1;
    sgn3 = 1;
        //std::cout << "sgns: case 2" << std::endl;
    return true;

  }

  if (x3x1 >0 && x1x2<0 && x3x2>0 &&
      x2.get_norm2() + (x3.get_point()).get_norm2() < 2 * ( x3x1 - x1x2 + x3x2 ) )
  {
    sgn2 = 1;
    sgn3 = -1;
    //std::cout << "sgns: case 3" << std::endl;
    return true;
  }

  if (x3x1 >0 && x1x2>0 && x3x2<0 &&
      (x2.get_norm2() + (x3.get_point()).get_norm2() < 2 * (  x3x1 + x1x2 - x3x2 )) )
  {

    sgn2 = -1;
    sgn3 = -1;
        //std::cout << "sgns: case 4" << std::endl;
    return true;
  }

  return false;
}

    /*
     runs 3-reduction
     TODO: input should be GaussQueue_ReturnType (same for sieve_2_iteration)
     */

template<class SieveTraits> void Sieve<SieveTraits,false>::sieve_3_iteration (typename SieveTraits::FastAccess_Point &p)
{
  
  if (p.is_zero() )
  {
    return; //TODO: Ensure sampler does not output 0 (currently, it happens).
  }
  
  //double px1_target  = .1111; // TO ADJUST
  int scalar=0; //for 2-reduction

  auto it_comparison_flip=main_list.cend(); //to store the point where the list elements become larger than p.

  FilteredListType filtered_list;

  auto it = main_list.cbegin();

  while (it!=main_list.cend())
  {
    if (p  < (*it) )
    {
      it_comparison_flip = it;
      break;
    }

    //EntryType sc_prod_px1 = compute_sc_product(p,*it);

    //
    //check for 2-reduction
    //
    if ( check2red(p, *it, scalar) )
    {
      assert(scalar!=0); //should not be 0 in any case
      p-= (*it) * scalar;

      if (p.is_zero() )
      {
        statistics.increment_number_of_collisions();
      }
      else
      {
        main_queue.push(std::move(p));
      }
      return;
    }

    statistics.set_filtered_list_size(0);

    filtered_list.reserve(filtered_list_size_max);

    EntryType sc_prod_px1;
    if (check_sc_prod(p, *it, sc_prod_px1))
    {
        
      for (auto & filtered_list_point: filtered_list)
      {

        //compute exact sc_prod right away
        EntryType sc_prod_x1x2 = compute_sc_product(*it, filtered_list_point.get_point());
        statistics.increment_number_of_scprods_level2();
        
        int sgn2, sgn3;

        //check if || p \pm x1 \pm x2 || < || p ||
        // ! check_triple assumes that the first argument has the largest norm
        if ( check_triple<SieveTraits> ( p, *it, filtered_list_point, sc_prod_px1, sc_prod_x1x2, sgn2, sgn3, true) )
        {
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
          unsigned int lvl = 0;
          uint_fast32_t approx_scprod = static_cast<uint_fast32_t>(compute_sc_product_bitapprox_level(p, *it, lvl));
          statistics.red_stat_innloop[lvl][static_cast<uint_fast32_t>(approx_scprod)]++;
#endif
          //TODO:  RETRIEVE ||p|| from the sc_prods
          //EntryType pnorm_old = p.get_norm2();
          p += (*it)*sgn2 + (filtered_list_point).get_point() * sgn3;
          //p +=  (*it)*sgn2 + *(filtered_list_point).get_point() * sgn3;

          //FOR DEBUGGING
          /*
          if (p.get_norm2() > pnorm_old)
          {
            std::cout << "bug in computing p " << std::endl;
            assert(false);
          }
          */
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
      //typename SieveTraits::FlilteredPointType new_filtered_point((*it).make_copy(), sc_prod_px1);
      typename SieveTraits::FlilteredPointType new_filtered_point(&(*it), sc_prod_px1);
      filtered_list.push_back(std::move(new_filtered_point));
      
    } //if (check_sc_prod(p, *it, sc_prod_px1))
    else
    {
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
      unsigned int lvl = 0;
      uint_fast32_t approx_scprod = static_cast<uint_fast32_t>(compute_sc_product_bitapprox_level(p, *it, lvl));
      statistics.no_red_stat_innloop[lvl][static_cast<uint_fast32_t>(approx_scprod)]++;
#endif
    }
    
    ++it;
  } //while-loop

  main_list.insert_before(it_comparison_flip,p.make_copy());
  statistics.increment_current_list_size();
  //std::cout << "list_size = " <<current_list_size << std::endl;
  if(update_shortest_vector_found(p))
  {
    if(verbosity>=2)
    {
      std::cout << "New shortest vector found. Norm2 = " << get_best_length2() << std::endl;
    }
  }


  //now p is not the largest
  //it_comparison_flip points to the next after p list-element
  it = it_comparison_flip;
  //for(auto it = it_comparison_flip; it!=main_list.cend(); )
  while (it !=main_list.cend())
  {

    //if true, do not put into the filtered_list
    //if true due to 2-reduction, re-compute the sc_product
    bool x1_reduced = false;

    //EntryType sc_prod_px1 = compute_sc_product(p,*it);
    //statistics.increment_number_of_scprods_level2();

        
    //
    //check for 2-reduction
    //

    if ( check2red(*it, p, scalar) )
    {
      assert(scalar!=0); //should not be 0 in any case
      typename SieveTraits::FastAccess_Point v_new = (*it) - (p*scalar);


      if (v_new.is_zero() )
      {
        statistics.increment_number_of_collisions();
      }

      main_queue.push(std::move(v_new));

      it = main_list.erase(it);
      statistics.decrement_current_list_size();

      //it was increased by the erase; we may already reach the end, so the code below will segfalut without the if-cond below
      if (it == main_list.cend())
      {
        return;
      }
    }
        
    //
    // 3-rediction
    //
    // Now x1 is the largest      

    EntryType sc_prod_px1;
    if (check_sc_prod(p, *it, sc_prod_px1))
    {
      
      for (auto & filtered_list_point: filtered_list)
      {
        
        EntryType sc_prod_x1x2 = compute_sc_product(*it, (filtered_list_point).get_point());
        statistics.increment_number_of_scprods_level2();
        
        int  sgn2, sgn3;

        // ! check_triple assumes that the first argument has the largest norm
        if ( check_triple<SieveTraits> ( *it, p, filtered_list_point, sc_prod_px1, sc_prod_x1x2, sgn2, sgn3, false) )
        {
          
          typename SieveTraits::FastAccess_Point v_new =(*it) + p*sgn2 + filtered_list_point.get_point() * sgn3;
                    //typename SieveTraits::FastAccess_Point v_new =(*it) + p*sgn2 + *(filtered_list_point).get_point() * sgn3;

          if (v_new.is_zero() )
          {
            statistics.increment_number_of_collisions();
          }
          //FOR DEBUG
          /*
          if (v_new.get_norm2() > (*it).get_norm2())
          {
              std::cout << "bug in computing v_new" << std::endl;
              assert(false);
          }
          */
          main_queue.push(std::move(v_new));
          it = main_list.erase(it);

          statistics.decrement_current_list_size();

          x1_reduced = true;

          break; //for-loop over the filtered_list
        }
      } //for-loop

      if (!x1_reduced)
      {
        //typename SieveTraits::FlilteredPointType new_filtered_point((*it).make_copy(), sc_prod_px1);
        typename SieveTraits::FlilteredPointType new_filtered_point(&(*it), sc_prod_px1);
        filtered_list.push_back(std::move(new_filtered_point));
      }
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
      unsigned int lvl = 0;
      uint_fast32_t approx_scprod = static_cast<uint_fast32_t>(compute_sc_product_bitapprox_level(p, *it, lvl));
      statistics.red_stat_innloop[lvl][static_cast<uint_fast32_t>(approx_scprod)]++;
#endif
  
    } //if (check_sc_prod<SieveTraits>(p, *it, sc_prod_px1))
    else
    {
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
      unsigned int lvl = 0;
      uint_fast32_t approx_scprod = static_cast<uint_fast32_t>(compute_sc_product_bitapprox_level(p, *it, lvl));
      statistics.no_red_stat_innloop[lvl][static_cast<uint_fast32_t>(approx_scprod)]++;
#endif
    }
    if (!x1_reduced)
    {
      ++it;
    }

        /*
        if (filtered_list.size()>0) {
            std::cout << "filtered.size() = " << filtered_list.size() << std::endl;
        }
         */
  } // 'lower' while-loop

  statistics.set_filtered_list_size(filtered_list.size());
  filtered_list.clear();
}


} //namespace GaussSieve


/*Unused*/


//The function checks if ||x1 + scalar* x2|| < ||x1||
// the value <x1,x2> is provided as input
// scalar is modified
template<class SieveTraits, class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
bool check_2red_with_scprod (typename SieveTraits::FastAccess_Point const &x1,
                             typename SieveTraits::FastAccess_Point const &x2,
                             typename SieveTraits::EntryType const &x1x2, Integer & scalar)
{
  using std::abs;
  using std::round;
  
  typename SieveTraits::EntryType const abs_2scprod = abs(x1x2 * 2);
  
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

/* OLD CODE */

       /*Collecting stats*/
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