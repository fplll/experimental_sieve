

namespace GaussSieve{

// The function checks if
//    || x1 \pm x2 \pm x3 || < || x1 ||
// The first arguments is assumed to have the largest norm
// The last two arguments are modified. They return the correct signs, i.e.
// x1 := x1 + sgn1*x2 + x3*sng2
// p_is_max is true if x1 ==p, in which case x3_X stores <x3,x2>
// otherwise x3_X stores <x3, x1>


    int filtered_list_size_max = 1000;

template<class SieveTraits>
bool check_3red (typename SieveTraits::FastAccess_Point  const &x1,
                 typename SieveTraits::FastAccess_Point const &x2,
                 typename SieveTraits::FlilteredPointType const &x3,
                 typename SieveTraits::EntryType const &x1x2,
                 typename SieveTraits::EntryType const &x3_X,
                 bool p_is_max, int &sgn2, int &sgn3)
{
        //retrieve the signs of the 3 inner-products;
        //they should satisfy (<p, x1>*<p, x2>*<x1, x2>) < 0
        //check for all 4 combinations that do not lead to a reduction

    typename SieveTraits::EntryType x3x1;
    typename SieveTraits::EntryType x3x2;
    
    if (p_is_max) {
        x3x1 = x3_X;
        x3x2 = x3.get_sc_prod();
    }
    else
    {
        x3x1 = x3.get_sc_prod();
        x3x2 = x3_X;
    }
    
    if (x1x2 * x3x1 * x1x2 >0)
        return false;

    if (x3x1 <0 && x1x2<0 && x3x2<0 &&
        x2.get_norm2() + (*x3.get_point()).get_norm2() <
        2 * ( abs(x3x1) + abs(x1x2) + abs(x3x2) ) )
    {
        
        sgn2 = 1;
        sgn3 = 1;
        //std::cout << "sgns: case 1" << std::endl;
        return true;
    }
    
    if (x3x1 <0 && x1x2>0 && x3x2>0 &&
        x2.get_norm2() + (*x3.get_point()).get_norm2() <
        2 * ( x3x1 + x1x2 + x3x2 ) )
    {
        
        sgn2 = -1;
        sgn3 = 1;
        //std::cout << "sgns: case 2" << std::endl;
        return true;
            
    }
    
    if (x3x1 >0 && x1x2<0 && x3x2>0 &&
        x2.get_norm2() + (*x3.get_point()).get_norm2() <
        2 * ( x3x1 - x1x2 + x3x2 ) )
    {
        
        sgn2 = 1;
        sgn3 = -1;
        //std::cout << "sgns: case 3" << std::endl;
        return true;
    }
    
    if (x3x1 >0 && x1x2>0 && x3x2<0 &&
        x2.get_norm2() + (*x3.get_point()).get_norm2() <
        2 * (  x3x1 + x1x2 - x3x2 ) )
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
    if (p.is_zero() ) return; //TODO: Ensure sampler does not output 0 (currently, it happens).
    
    
    
    // ! targets are squared
    double px1  = .102; // TO ADJUST
    double x1x2 = .102;

    int scalar; //for 2-reduction

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
        
        //
        //check for 2-reduction
        //
        if ( GaussSieve::check2red<SieveTraits>(p, *it, scalar) )
        {
            assert(scalar!=0); //should not be 0 in any case
            p-= (*it) * scalar;

            if (p.is_zero() )
            {
                number_of_collisions++;
            }
            else
            {
                main_queue.push(std::move(p));
            }
            return;
        }
        
        filtered_list_size = 0;
        
        filtered_list.reserve(filtered_list_size_max);

        //
        //compare <p, x1> with px1
        //
        //TODO: Change the computations below to smth. faster/better/cleverer/etc.
        EntryType sc_prod_px1 = compute_sc_product(p,*it);
        double sc_prod_px1_norm = convert_to_double( sc_prod_px1)*convert_to_double(sc_prod_px1 ) /
                        ( convert_to_double ( p.get_norm2()) * convert_to_double( (*it).get_norm2() )) ;
        if (std::abs(sc_prod_px1_norm) > px1)
        {
            //This is a fast iteration accodring to the Internet
            //use 'auto &' to take a reference (copy-constructor is deleted)
            for (auto & filtered_list_it: filtered_list)
            {
                

                //EntryType sc_prod_x1x2 = compute_sc_product(*it, (filtered_list_it).get_point());
                EntryType sc_prod_x1x2 = compute_sc_product(*it, *(filtered_list_it).get_point());

                int sgn2, sgn3;
                
                //check if || p \pm x1 \pm x2 || < || p ||
                // ! check_3red assumes that the first argument has the largest norm
                if ( check_3red<SieveTraits> ( p, *it, filtered_list_it, sc_prod_px1, sc_prod_x1x2, true, sgn2, sgn3) )
                {
                    
                    //TODO:  RETRIEVE ||p|| from the sc_prods
                    
                   
                    //p += (*it)*sgn2 + (filtered_list_it).get_point() * sgn3;
                    p +=  (*it)*sgn2 + *(filtered_list_it).get_point() * sgn3;
                    
                    //FOR DEBUGGING
                    /* 
                    EntryType pnorm_old = p.get_norm2();
                    if (p.get_norm2() > pnorm_old)
                    {
                        assert(false);
                    }
                    */
                    
                    if (p.is_zero() )
                    {
                        number_of_collisions++;
                    }
                    else
                    {
                        main_queue.push(std::move(p));
                    }
                    return;
                }
            }

            //typename SieveTraits::FlilteredPointType new_filtered_point((*it).make_copy(), sc_prod_px1);
            typename SieveTraits::FlilteredPointType new_filtered_point(&(*it), sc_prod_px1);
            filtered_list.push_back(std::move(new_filtered_point));
            
        }

        ++it;
    } //while-loop

    main_list.insert_before(it_comparison_flip,p.make_copy());
    ++current_list_size;
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
        bool x1_reduced = false;
        
        
        //
        //check for 2-reduction
        //

        if ( GaussSieve::check2red<SieveTraits>(*it, p, scalar) )
        {
            assert(scalar!=0); //should not be 0 in any case
            typename SieveTraits::FastAccess_Point v_new = (*it) - (p*scalar);
            
          
            if (v_new.is_zero() )
            {
                number_of_collisions++;
            }
            
            main_queue.push(std::move(v_new));
            
            it = main_list.erase(it);
            --current_list_size;
            
            
            //it was increased by the erase; we may already reach the end, so the code below will segfalut
            if (it == main_list.cend())
                return;
            
        }
        
        
        //
        // 3-rediction
        //
        // Now x1 is the largest
        EntryType sc_prod_px1 = compute_sc_product(p,*it);
        double sc_prod_px1_norm = convert_to_double( sc_prod_px1 )* convert_to_double( sc_prod_px1 )  /
        ( convert_to_double ( p.get_norm2()) * convert_to_double( (*it).get_norm2() )) ;
        
        if (std::abs(sc_prod_px1_norm) > px1)
        {
            
            for (auto & filtered_list_it: filtered_list)
            {
                
                //EntryType sc_prod_x1x2 = compute_sc_product(*it, (filtered_list_it).get_point());
                EntryType sc_prod_x1x2 = compute_sc_product(*it, *(filtered_list_it).get_point());
                int  sgn2, sgn3;
                
                // ! check_3red assumes that the first argument has the largest norm
                if ( check_3red<SieveTraits> ( *it, p, filtered_list_it, sc_prod_px1, sc_prod_x1x2, false, sgn2, sgn3) )
                {
                    //std::cout <<  sgn2 << " " << sgn3 << " " << std::endl;
                    
                    
                    //typename SieveTraits::FastAccess_Point v_new =(*it) + p*sgn2 + (filtered_list_it).get_point() * sgn3;
                    typename SieveTraits::FastAccess_Point v_new =(*it) + p*sgn2 + *(filtered_list_it).get_point() * sgn3;
                    
                    if (v_new.is_zero() )
                    {
                        number_of_collisions++;
                    }
                    
                    //FOR DEBUGGINS
                    /*
                    if (v_new.get_norm2() > (*it).get_norm2())
                    {
                        std::cout << "bug in computing v_new" << std::endl;
                        assert(false);
                    }
                    */
                    main_queue.push(std::move(v_new));
                    it = main_list.erase(it);
                    
                    --current_list_size;
                    
                    x1_reduced = true;
                    
                    break; //for-loop over the filtered_list
                }
            } //if-cond
            
            
            if (!x1_reduced)
            {
                //typename SieveTraits::FlilteredPointType new_filtered_point((*it).make_copy(), sc_prod_px1);
                typename SieveTraits::FlilteredPointType new_filtered_point(&(*it), sc_prod_px1);
                filtered_list.push_back(std::move(new_filtered_point));
            }
            
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
    
    filtered_list_size =filtered_list.size();
    filtered_list.clear();

}


} //namespace
