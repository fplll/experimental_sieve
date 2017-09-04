

namespace GaussSieve{
    
// The function checks if
//    || p \pm x1 \pm x2 || < || p ||
// The first arguments is assumed to have the largest norm
// The last three arguments are modified. They return the correct signs
    
template<class SieveTraits >
bool check_3red (typename SieveTraits::FastAccess_Point  const &p,
                 typename SieveTraits::FastAccess_Point const &x1,
                 typename SieveTraits::Filtered_Point const &x2,
                 typename SieveTraits::EntryType const &px1,
                 typename SieveTraits::EntryType const &x1x2,
                 int &sgn1, int &sgn2, int &sgn3)
{
        //retrieve the signs of the 3 inner-products;
        //they should satisfy (<p, x1>*<p, x2>*<x1, x2>) < 0
        //check for all 4 combinations that do not lead to a reduction
        if ( (x2.get_sc_prod() >0 && px1>0 && x1x2>0) ||
             (x2.get_sc_prod() <0 && px1<0 && x1x2>0) ||
             (x2.get_sc_prod() <0 && px1>0 && x1x2<0) ||
             (x2.get_sc_prod() >0 && px1<0 && x1x2<0) )
            return false;
    
    using EntryType = typename SieveTraits::EntryType;
    
    //check if the there is a reduction assumming all <,> are negative -> all sgns are '+'
    //no reduction if ||x_1 ||^2 +||x_2 ||^2 >= 2* ( <p, x1> + <p, x2>+ <x1,x2>)
    if ( x1.get_norm2() + x2.get_point().get_norm2 >=
        2 * ( abs(x2.get_sc_prod()) + abs(px1) + abs(x1x2)  ) )
        return false;
    
    // all good sign-configurations: (p+x1+x2), (p+x1-x2), (p-x1+x2), (-p+x1+x2)
    if (x2.get_sc_prod() <0 && px1<0 && x1x2<0)
    {
        sgn1 = 1;
        sgn2 = 1;
        sgn3 = 1;
    }
    
    if (x2.get_sc_prod() <0 && px1>0 && x1x2>0)
    {
        sgn1 = 1;
        sgn2 = -1;
        sgn3 = 1;
    }
    
    if (x2.get_sc_prod() >0 && px1<0 && x1x2>0)
    {
        sgn1 = 1;
        sgn2 = 1;
        sgn3 = -1;
    }
    
    if (x2.get_sc_prod() >0 && px1>0 && x1x2<0)
    {
        sgn1 = -1;
        sgn2 = 1;
        sgn3 = 1;
    }

    return true;
}
    
    /*
     runs 3-reduction
     TODO: input should be GaussQueue_ReturnType (same for sieve_2_iteration)
     */

template<class SieveTraits> void Sieve<SieveTraits,false>::sieve_3_iteration (typename SieveTraits::FastAccess_Point &p)
{
    if (p.is_zero() ) return; //TODO: Ensure sampler does not output 0 (currently, it happens).
    
    double px1  = .31; // TO ADJUST
    double x1x2 = .31;
    
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
        if ( GaussSieve::check2red<SieveTraits>(*it, p, scalar) )
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
        
        
        //
        //compare <p, x1> with px1
        //
        //TODO: Change the computations below to smth. faster/better/cleverer/etc.
        EntryType sc_prod_px1 = compute_sc_product(p,*it);
        double sc_prod_px1_norm = convert_to_double( sc_prod_px1 ) / (convert_to_double ( p.get_norm2() * (*it.get_norm2())));
        
        if (std::abs(sc_prod_px1_norm) > px1)
        {
            //This is a fast iteration accodring to the Internet
            for (auto filtered_list_it: filtered_list)
            {
                //check if || p \pm x1 \pm x2 || < || p ||
        
                EntryType sc_prod_x1x2 = compute_sc_product(*it, (*filtered_list_it).get_point());
                
                int sgn1, sgn2;
                
                if ( check_3red ( p, *it, *filtered_list_it, *sc_prod_px1,  *sc_prod_x1x2, sgn1, sgn2) )
                {
                    
                }
            }
        }
        
        
    } //while-loop
}






} //namespace