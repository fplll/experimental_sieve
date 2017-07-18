/* DO NOT INCLUDE THIS FILE DIRECTLY
*/


template<class ET>
void Sieve<ET,false>::sieve_3_iteration (LatticePoint<ET> &p)
{
    if (p.norm2==0) return;
    ApproxLatticePoint<ET,false> pApprox (p);

    //cout << "------------" << endl;
    //cout << "Run iteration on p of approx norm " << pApprox.get_approx_norm2() << endl;

    int const n = get_ambient_dimension();
    typename MainListType::Iterator it_comparison_flip=main_list.cend();
    ET scalar;

    //if abs( <p, x1> / (|p||x1|) ) < px1bound, do not put it in the filtered_list
    float  px1bound = 0.31; // TO ADJUST


    float x1x2=.31;

    // length_factor determines the difference (multiplicatively) between the legnth of the i-th and (i+1)-st blocks
    float length_factor = 1.3; //TO ADJUST

    //number of blocks
    int NumOfBlocks = 0;

    auto it = main_list.cbegin();

    typename MainListType::Iterator first_element_of_the_block = main_list.cbegin();
    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_the_current_block =  first_element_of_the_block->get_approx_norm2();
    LatticeApproximations::ApproxTypeNorm2 max_length_of_the_current_block = floor(length_factor * assumed_norm_of_the_current_block + 1);



    //BlockDivisionType BlockPointers;
    //ApproxLatticePoint<ET,false> first_block_element =*it;
    //BlockPointers[NumOfBlocks] = &first_block_element;
    //BlockPointers[NumOfBlocks] = *it;

    LengthDivisionType BlockDivision;
    BlockDivision[NumOfBlocks] = it->get_approx_norm2();


    // the same number of appendices as the number of blocks
    std::array<AppendixType, 100> appendices;

    FilteredListType filtered_list;
    FilterDivisionType last_elements;
    FilterNumOfElems num_of_elements_in_filter;
    num_of_elements_in_filter.fill(0);

    last_elements[0] = filtered_list.begin();

    while (it!=main_list.cend())
    {

        //cout << "1: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //check if p is still of max. length
        if (p.norm2 < it.get_true_norm2())
        {
            //cout << "reached the position to insert p" << endl;
            it_comparison_flip = it;
            break;
        }


         //--------------------------------2-red-------------------------------

        //check 2-red
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );

        if(predict){


        ++number_of_exact_scprods;

         // preform 2-reduction
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            {
                p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);

                //put p back into the queue and break
                if (p.norm2!=0)
                    main_queue.push(p);
                //cout << "2-reduced p, break all the loops" << endl;
                return;
            }
            else
                ++number_of_mispredictions;
        }


        //--------------------------------3-red-------------------------------


        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
                BlockDivision[NumOfBlocks] = assumed_norm_of_the_current_block;

                //cout << "enter the next block" << endl;

                //cout << "1: filtered_list is:" << endl;
                //for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                //    cout << it1->get_sc_prod() << endl;

                //cout << "max_length_of_the_current_block = " << max_length_of_the_current_block << endl;


                //in case non of the list-elements from the new block will end up in filtered_list, the last_element of this block is the last_element of the prev. block
                last_elements[NumOfBlocks] = last_elements[NumOfBlocks-1];
                /*
                last_elements[NumOfBlocks] = filtered_list.end();
                if (num_of_elements_in_filter[NumOfBlocks-1]==0)
                {
                    if (NumOfBlocks==1)
                        last_elements[NumOfBlocks-1] = filtered_list.end();
                    else
                       last_elements[NumOfBlocks-1] = last_elements[NumOfBlocks-2];
                }
                 */

        }


        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);

        // scaling factor for <px1>: ||p|| * ||x_1||
        //float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        float scale = sqrt ( (float)pApprox.get_approx_norm2() * (float)(it->get_approx_norm2()) );

        //cout  << "true_inner_product_px1 / scale " << (float)true_inner_product_px1 / scale << endl;

        //if true_inner_product_px1 is large enough, put it in filtered_list
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {

            //cout << "true_inner_product_px1 = " << true_inner_product_px1 << endl;

            typename FilteredListType::iterator it_filter =  filtered_list.begin();
            typename FilteredListType::iterator it_comparison_flip_filter = filtered_list.end();

            int block_counter = 0;
            bool position_is_found = false;

            //cout << "filtered_list is:" << endl;
            //for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
            //        cout << it1->get_sc_prod() << endl;


            //loop over all elements x2 from the filtered list to 3-reduc (p, x1, x2)
            while ( it_filter != filtered_list.cend())
            {

                // OR take the length from BlockDivision[appendixCounter]
                LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = (it_filter->getApproxVector()).get_approx_norm2();



                //in res_upper, we store the upper bound on the inner-product of px2 for x2 from filtered_list
                float res_upper = 6000;

                //Compute_one_third_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, res_upper);

                //cout << "res_upper = " << res_upper << endl;

                //check if we even need to iterate over the current block of filtered_list
                if (it_filter->get_sc_prod() > res_upper)
                {

                    //check if the top element of the relevant appendix is larger than res_upper

                    /*
                    bool merge = false;
                    if (!appendices[block_counter].empty())
                    {
                        cout << "appendix.top() has norm = " << (appendices[block_counter].top()).get_sc_prod() << endl;
                        if ( (appendices[block_counter].top()).get_sc_prod() > res_upper)
                            merge = true;
                    }
                    */


                    ApproxTypeNorm2 true_inner_product_x1x2;

                    typename FilteredListType::iterator it_block =  it_filter;

                    int inblock_counter = 0;


                    while(it_block!=filtered_list.end() && it_block->get_sc_prod() > res_upper && inblock_counter < num_of_elements_in_filter[block_counter])
                    {

                            //cout << "it_block->get_sc_prod() = " << it_block->get_sc_prod() << endl;

                            //check if we found the position to insert *it
                            if (it_block->get_sc_prod() < abs(true_inner_product_px1) && !position_is_found && block_counter == NumOfBlocks)
                            {
                                //cout << "position is found " << endl;
                                it_comparison_flip_filter = it_block;
                                position_is_found = true;
                            }

                            //------------------check for reduction-------------------

                            predict = LatticeApproximations::Compare_Sc_Prod_3red(it_block->getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);

                            //cout << "REDUCTION CHECK" << endl;

                            //sgn1 and sgn2 store the signs s.t. || p + sgn1*x1 + sgn*x2 || < || p ||
                            // these values are computed in check3red_signs
                            int sgn1 = 0;
                            int sgn2 = 0;


                            //in fistered_list we store abs(<p, x1>); retrieve the true sign of x1 stored in filtered_list
                            int sgn_px1 = 1;
                            int sgn_px2 = 1 ;
                            int sgn_x1x2 = 1;


                            if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                            if (!it_block->get_sign()) sgn_px2 = -1;
                            if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;

                            //check if reduction
                            if (GaussSieve::check3red_signs(p, *(it.access_details()), (it_block->getApproxVector()).get_details(), sgn_px1, sgn_px2, sgn_x1x2, sgn1, sgn2))
                            {

                                //perfrom actual reduction
                                p = GaussSieve::perform3red(p, *(it.access_details()), (it_block->getApproxVector()).get_details(), sgn1, sgn2);

                                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);

                                if(p.norm2!=0)
                                    main_queue.push(p);


                                //cout << "reduced p: RE-START with p on norm2 = " << p.norm2 << endl;
                                //assert(false);


                                filtered_list.clear();
                                return;
                            }

                            //else
                                //cout << "check3red is false" << endl;


                            //--------------------- end of reduction's check ----------------------

                            //make merge with an element from the queue
                            if(!appendices[block_counter].empty() && it_block->get_sc_prod() < (appendices[block_counter].top()).get_sc_prod() )
                            {

                                FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> poped_from_queue =  appendices[block_counter].top();
                                appendices[block_counter].pop();



                                //Assume insert means inseart_before (TO CHECK)
                                filtered_list.insert(it_block,poped_from_queue);
                                num_of_elements_in_filter[block_counter]++;

                                //cout << "insert from appendix the element " << poped_from_queue.get_sc_prod() << endl;
                                //cout << "last_elements[block_counter]->get_sc_prod() " << last_elements[block_counter]->get_sc_prod() << endl;
                                //cout << "it_block->get_sc_prod() " << it_block->get_sc_prod() << endl;

                                //never happens
                                /*
                                if ( last_elements[block_counter]->get_sc_prod() > poped_from_queue.get_sc_prod() )
                                {
                                    last_elements[block_counter] = it_block; //TODO: TO CHECK.
                                    cout << "upd last_elements from queue to " <<  last_elements[block_counter]->get_sc_prod() <<  endl;
                                    assert(false);

                                }
                                */


                            }

                            //check if the new top element of the relevant appendix is larger than res_upper

                            /*
                            merge = false;
                            if (!appendices[appendixCounter].empty())
                            {
                                if ( (appendices[appendixCounter].top().getApproxVector()).get_approx_norm2() > res_upper )
                                merge = true;
                            }
                            */

                            inblock_counter++;
                            //cout << "inblock_counter " << inblock_counter << endl;
                            ++it_block;

                    } //end of iteration over the block


                    //insert from the appendix if the last element of the block is larger
                    /*
                    if ( !appendices[block_counter].empty() && (appendices[block_counter].top()).get_sc_prod() < last_elements[block_counter]->get_sc_prod())
                    {
                        FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> poped_from_queue =  appendices[block_counter].top();
                        appendices[block_counter].pop();


                        //cout << "about to append from appendix "  << endl;

                        //BUG HERE: it_block does not always point to the end of the block
                        filtered_list.insert(it_block,poped_from_queue);
                        num_of_elements_in_filter[block_counter]++;
                        last_elements[block_counter] = --it_block; //TODO: TO CHECK.
                        //cout << "last_elements[block_counter] is upd to " << last_elements[block_counter]->get_sc_prod() << endl;
                        ++it_block;
                    }
                     */

                } // end of if (it_filter->get_sc_prod() > res_upper)

                //cout << "last_elements[block_counter]->get_sc_prod() " << last_elements[block_counter]->get_sc_prod() << endl;
                it_filter = last_elements[block_counter];

                //if(it_filter!=filtered_list.end())
                //   cout << "change it_filter to " << last_elements[block_counter]->get_sc_prod() << endl;

                block_counter++;

                //cout << "block_counter = " << block_counter << endl;
                if(it_filter!=filtered_list.end())
                    it_filter++;

            } //end of the loop over filtered_list

            //cout << "NumOfBlocks = " << NumOfBlocks << endl;
            //cout << "num_of_elements_in_filter[NumOfBlocks] = " << num_of_elements_in_filter[NumOfBlocks] << endl;


            bool px1_sign;
            if (true_inner_product_px1>0)
                px1_sign = true;
            else
                px1_sign = false;

            FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);

            //cout << "position_is_found = " << position_is_found << endl;

            //is the insertion position was reached during the filter_list traversal, insert into filtered_list; otherwise, push into the queue
            if (position_is_found || num_of_elements_in_filter[NumOfBlocks] == 0 )
            {
                filtered_list.insert(it_comparison_flip_filter,new_filtered_point);
                num_of_elements_in_filter[NumOfBlocks]++;

                //cout << "num_of_elements_in_filter[NumOfBlocks] = " << num_of_elements_in_filter[NumOfBlocks] << endl;

                //if it was the first insertion, assign last_elements[NumOfBlocks]
                if (num_of_elements_in_filter[NumOfBlocks]==1)
                {
                    last_elements[NumOfBlocks] = --it_comparison_flip_filter;
                    //cout << "upd last_elements as num_of_elements_in_filter[NumOfBlocks] == 0 to " << last_elements[NumOfBlocks]->get_sc_prod() <<  endl;
                }

                else if ( last_elements[NumOfBlocks]->get_sc_prod() > abs(true_inner_product_px1) )
                {
                    last_elements[NumOfBlocks] = --it_comparison_flip_filter;
                    //cout << "upd last_elements to " <<  last_elements[NumOfBlocks]->get_sc_prod() <<  endl;
                }
                //else
                //    cout << "do not change last_elements " << endl;
            }
            /*
            else if (last_elements[NumOfBlocks]->get_sc_prod() > abs(true_inner_product_px1) )
            {
                filtered_list.insert(filtered_list.end(),new_filtered_point);
                last_elements[NumOfBlocks] = --filtered_list.end();
            }
            */
            else
                appendices[NumOfBlocks].push(new_filtered_point);



            //cout << "filtered_list.size = " << filtered_list.size() << endl;
            //cout << "appendices[NumOfBlocks].size = " << appendices[NumOfBlocks].size() << endl;

            //if (!appendices[NumOfBlocks].empty())
            //    cout << "appendix[NumOfBlocks].top() is " << (appendices[NumOfBlocks].top()).get_sc_prod() << endl;

            /*
            if (filtered_list.size() > 10 && NumOfBlocks > 1)
            {
                for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                    cout << it1->get_sc_prod() << endl;

                cout << "Number of elements in filtered_list per block: " << endl;
                int i=0;
                while (num_of_elements_in_filter[i]>0)
                {
                    cout << num_of_elements_in_filter[i] << " ";
                    i++;
                }

                assert(false);
            }
            */

        } //if (px1 is larger enough)

        ++it;
    }

    //INSERT p
    main_list.insert_before(it_comparison_flip,p);
    ++current_list_size;
    //cout << "list_size = " <<current_list_size << endl;
    if(update_shortest_vector_found(p))
    {
        if(verbosity>=2)
        {
            cout << "New shortest vector found. Norm2 = " << get_best_length2() << endl;
        }
    }

    /*
    cout << "p is inserted" << endl;
    cout << "sizes of appendices are : "<< endl;
    int i=0;
    while (i<=NumOfBlocks)
    {
        cout << appendices[i].size() << " ";
        i++;
    }
    */

     it =it_comparison_flip; //points to the next after p list-element

     //assumed_norm_of_the_current_block = it->get_approx_norm2();
     //max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
     //NumOfBlocks++;

     bool reduced_x1 = false; //if true, do not increase ++it and do not include into filtered_list

     //now we are reducing *it
     while (it!=main_list.cend())
     {
        reduced_x1 = false;
        //cout << "2: consider a list element of approx norm = " << it->get_approx_norm2() << endl;


        //--------------------------------2-red-------------------------------
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );

        // preform 2-reduction
        if(predict){


        ++number_of_exact_scprods;
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
        {
                LatticePoint<ET> reduced = GaussSieve::perform2red(p, *(it.access_details()), scalar);

                //put reduced back into the queue
                if (reduced.norm2!=0)
                    main_queue.push(reduced);


                //cout << "put x1 of size " << reduced.norm2 << " into the queue " << endl;

                it = main_list.erase(it); //also makes ++it
                --current_list_size;
                reduced_x1 = true;
        }
        else
                ++number_of_mispredictions;
        }


        //--------------------------------3-red-------------------------------

        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
                BlockDivision[NumOfBlocks] = assumed_norm_of_the_current_block;

                //cout << "enter the next block" << endl;

                //cout << "filtered_list is:" << endl;
                //for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                //    cout << it1->get_sc_prod() << endl;

                //cout << "max_length_of_the_current_block = " << max_length_of_the_current_block << endl;


                //in case non of the list-elements from the new block will end up in filtered_list, the last_element of this block is the last_element of the prev. block
                last_elements[NumOfBlocks] = last_elements[NumOfBlocks-1];
                /*
                last_elements[NumOfBlocks] = filtered_list.end();
                if (num_of_elements_in_filter[NumOfBlocks-1]==0)
                {
                    if (NumOfBlocks==1)
                        last_elements[NumOfBlocks-1] = filtered_list.end();
                    else
                       last_elements[NumOfBlocks-1] = last_elements[NumOfBlocks-2];
                }
                 */

        }

        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);

        // scaling factor for <px1>: ||p|| * ||x_1||
        //float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        float scale = sqrt ( (float)pApprox.get_approx_norm2() * (float)(it->get_approx_norm2()) );


        //if true_inner_product_px1 is large enough, put it in filtered_list
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {
            //cout  << "true_inner_product_px1 " << true_inner_product_px1 << endl;
            typename FilteredListType::iterator it_filter =  filtered_list.begin();
            typename FilteredListType::iterator it_comparison_flip_filter = filtered_list.end();

            int block_counter = 0;
            bool position_is_found = false;

            //------for debug------

            /*
            cout << "2: filtered_list is:" << endl;
            for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                    cout << it1->get_sc_prod() << endl;

            cout << "2: Number of elements in filtered_list per block: " << endl;
            int i=0;
            while (i<=NumOfBlocks)
            {
                cout << num_of_elements_in_filter[i] << " ";
                i++;
            }
            cout <<endl;
            */
            //----------------------


            while (!reduced_x1 && it_filter != filtered_list.cend())
            {

                LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = (it_filter->getApproxVector()).get_approx_norm2();



                //in res_upper, we store the upper bound on the inner-product of px2 for x2 from filtered_list
                float res_upper = 6000;

                //Compute_one_third_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, res_upper);

                //cout << "res_upper = " << res_upper << endl;
                //cout << "it_filter->get_sc_prod() = " << it_filter->get_sc_prod() << endl;

                //check if we even need to iterate over the current block of filtered_list
                if (it_filter->get_sc_prod() > res_upper)
                {


                    ApproxTypeNorm2 true_inner_product_x1x2;

                    typename FilteredListType::iterator it_block =  it_filter;

                    int inblock_counter = 0;

                    while(it_block!=filtered_list.end() && it_block->get_sc_prod() > res_upper && inblock_counter < num_of_elements_in_filter[block_counter])
                    {

                            //cout << "it_block->get_sc_prod() = " << it_block->get_sc_prod() << endl;

                            //check if we found the position to insert *it
                            if (it_block->get_sc_prod() < abs(true_inner_product_px1) && !position_is_found && block_counter == NumOfBlocks)
                            {
                                //cout << "position is found " << endl;
                                it_comparison_flip_filter = it_block;
                                position_is_found = true;
                            }


                            //make merge with an element from the queue
                            if(!appendices[block_counter].empty() && it_block->get_sc_prod() < (appendices[block_counter].top()).get_sc_prod() )
                            {

                                FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> poped_from_queue =  appendices[block_counter].top();
                                appendices[block_counter].pop();



                                //Assume insert means inseart_before (TO CHECK)
                                filtered_list.insert(it_block,poped_from_queue);
                                num_of_elements_in_filter[block_counter]++;

                                //cout << "num_of_elements_in_filter[block_counter] = " << num_of_elements_in_filter[block_counter] << endl;

                                //cout << "insert from appendix the element " << poped_from_queue.get_sc_prod() << endl;
                                //cout << "last_elements[block_counter]->get_sc_prod() " << last_elements[block_counter]->get_sc_prod() << endl;
                                //cout << "it_block->get_sc_prod() " << it_block->get_sc_prod() << endl;

                            }

                            //------------------check for reduction-------------------

                            predict = LatticeApproximations::Compare_Sc_Prod_3red(it_block->getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);

                            //cout << "REDUCTION CHECK" << endl;

                            //sgn1 and sgn2 store the signs s.t. || p + sgn1*x1 + sgn*x2 || < || p ||
                            // these values are computed in check3red_signs
                            int sgn1 = 0;
                            int sgn2 = 0;


                            //in fistered_list we store abs(<p, x1>); retrieve the true sign of x1 stored in filtered_list
                            int sgn_px1 = 1;
                            int sgn_px2 = 1 ;
                            int sgn_x1x2 = 1;


                            if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                            if (!it_block->get_sign()) sgn_px2 = -1;
                            if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;

                            //check if reduction
                            if (GaussSieve::check3red_signs(*(it.access_details()), p, (it_block->getApproxVector()).get_details(), sgn_px1, sgn_px2, sgn_x1x2, sgn1, sgn2))
                            {

                                //perfrom actual reduction
                                LatticePoint<ET> reduced = GaussSieve::perform3red(*(it.access_details()), p, (it_block->getApproxVector()).get_details(), sgn1, sgn2);

                                if(reduced.norm2!=0)
                                        main_queue.push(reduced);

                                it = main_list.erase(it); //also makes ++it
                                --current_list_size;

                                //cout << "reduced x1" << endl;
                                //assert(false);

                                //break both while-loops and set the flag not to insert *it into filtered_list
                                reduced_x1 = true;
                                break;
                            }

                            //else
                                //cout << "check3red is false" << endl;


                            //--------------------- end of reduction's check ----------------------


                            inblock_counter++;
                            ++it_block;

                    } //end of iteration over the block


                    //insert from the appendix if it_block now points to a smaller element
                    /*
                    if ( ( !appendices[block_counter].empty() && it_block!=filtered_list.end() && (appendices[block_counter].top()).get_sc_prod() < it_block->get_sc_prod() ) )
                    {
                        FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> poped_from_queue =  appendices[block_counter].top();
                        appendices[block_counter].pop();


                        cout << "about to append from appendix "  << endl;

                        //BUG HERE: it_block does not always point to the end of the block
                        filtered_list.insert(it_block,poped_from_queue);
                        num_of_elements_in_filter[block_counter]++;
                        last_elements[block_counter] = --it_block; //TODO: TO CHECK.
                        cout << "num_of_elements_in_filter[block_counter] " << num_of_elements_in_filter[block_counter] << endl;
                        cout << "last_elements[block_counter] is upd to " << last_elements[block_counter]->get_sc_prod() << endl;
                        ++it_block;

                        cout << "2: filtered_list is:" << endl;
                        for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                            cout << it1->get_sc_prod() << endl;
                    }
                     */

                } // end of if (it_filter->get_sc_prod() > res_upper)

                //cout << "last_elements[block_counter]->get_sc_prod() " << last_elements[block_counter]->get_sc_prod() << endl;
                it_filter = last_elements[block_counter];

                //if(it_filter!=filtered_list.end())
                //   cout << "change it_filter to " << last_elements[block_counter]->get_sc_prod() << endl;

                block_counter++;

                //cout << "block_counter = " << block_counter << endl;
                if(it_filter!=filtered_list.end())
                    it_filter++;

            } //end of the loop over filtered_list

            if (!reduced_x1)
            {
                bool px1_sign;
                if (true_inner_product_px1>0)
                    px1_sign = true;
                else
                    px1_sign = false;

                FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);

                //cout << "position_is_found = " << position_is_found << endl;

                //is the insertion position was reached during the filter_list traversal, insert into filtered_list; otherwise, push into the queue
                if (position_is_found || num_of_elements_in_filter[NumOfBlocks] == 0 )
                {
                    filtered_list.insert(it_comparison_flip_filter,new_filtered_point);
                    num_of_elements_in_filter[NumOfBlocks]++;

                    //cout << "num_of_elements_in_filter[NumOfBlocks] = " << num_of_elements_in_filter[NumOfBlocks] << endl;

                    //if it was the first insertion, assign last_elements[NumOfBlocks]
                    if (num_of_elements_in_filter[NumOfBlocks]==1)
                    {
                        last_elements[NumOfBlocks] = --it_comparison_flip_filter;
                        //cout << "upd last_elements as num_of_elements_in_filter[NumOfBlocks] == 0 to " << last_elements[NumOfBlocks]->get_sc_prod() <<  endl;
                    }

                    else if ( last_elements[NumOfBlocks]->get_sc_prod() > abs(true_inner_product_px1) )
                    {
                        last_elements[NumOfBlocks] = --it_comparison_flip_filter;
                        //cout << "upd last_elements to " <<  last_elements[NumOfBlocks]->get_sc_prod() <<  endl;
                    }
                    //else
                    //    cout << "do not change last_elements " << endl;
                }
                /*
                else if (last_elements[NumOfBlocks]->get_sc_prod() > abs(true_inner_product_px1) )
                {
                    filtered_list.insert(filtered_list.end(),new_filtered_point);
                    last_elements[NumOfBlocks] = --filtered_list.end();
                }
                */
                else
                    appendices[NumOfBlocks].push(new_filtered_point);
                }

        } //if (px1 is larger enough)

        if (!reduced_x1)
        ++it;

    }

    /*
    if (filtered_list.size()>60)
    {

         cout << "filtered_list is:" << endl;
         for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
         cout << it1->get_sc_prod() << endl;

         cout << "Number of elements in filtered_list per block: " << endl;
         int i=0;
         while (i<=NumOfBlocks)
         {
         cout << num_of_elements_in_filter[i] << " ";
         i++;
         }

        assert(false);

    }
    */


     filtered_list.clear();
}


//
// Without appendies; filtered_list is not sorted. iterate over all filtered_list
//

template<class ET>
void Sieve<ET,false>::sieve_3_iteration_test (LatticePoint<ET> &p)
{
    if (p.norm2==0) return;
    ApproxLatticePoint<ET,false> pApprox (p);

    //cout << "------------" << endl;
    //cout << "Run iteration on p of approx norm " << pApprox.get_approx_norm2() << endl;

    int const n = get_ambient_dimension();
    typename MainListType::Iterator it_comparison_flip=main_list.cend();
    ET scalar;

    //if abs( <p, x1> / (|p||x1|) ) < px1bound, do not put it in the filtered_list
    float  px1bound = 0.31; // TO ADJUST


    float x1x2=.31;

    // length_factor determines the difference (multiplicatively) between the legnth of the i-th and (i+1)-st blocks
    float length_factor = 1.3; //TO ADJUST

    //number of blocks
    int NumOfBlocks = 0;

    auto it = main_list.cbegin();

    typename MainListType::Iterator first_element_of_the_block = main_list.cbegin();
    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_the_current_block =  first_element_of_the_block->get_approx_norm2();
    LatticeApproximations::ApproxTypeNorm2 max_length_of_the_current_block = floor(length_factor * assumed_norm_of_the_current_block + 1);



    //BlockDivisionType BlockPointers;
    //ApproxLatticePoint<ET,false> first_block_element =*it;
    //BlockPointers[NumOfBlocks] = &first_block_element;
    //BlockPointers[NumOfBlocks] = *it;

    LengthDivisionType BlockDivision;
    BlockDivision[NumOfBlocks] = it->get_approx_norm2();


    // the same number of appendices as the number of blocks
    //std::array<AppendixType, 100> appendices;

    FilteredListType filtered_list;
    FilterDivisionType last_elements;
    FilterNumOfElems num_of_elements_in_filter;
    num_of_elements_in_filter.fill(0);

    //last_elements[0] = filtered_list.begin();

    while (it!=main_list.cend())
    {

        //cout << "1: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //check if p is still of max. length
        if (p.norm2 < it.get_true_norm2())
        {
            //cout << "reached the position to insert p" << endl;
            it_comparison_flip = it;
            break;
        }


         //--------------------------------2-red-------------------------------

        //check 2-red
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );

        if(predict){


        ++number_of_exact_scprods;

         // preform 2-reduction
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            {
                p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);

                //put p back into the queue and break
                if (p.norm2!=0)
                    main_queue.push(p);
                //cout << "2-reduced p, break all the loops" << endl;
                return;
            }
            else
                ++number_of_mispredictions;
        }


        //--------------------------------3-red-------------------------------


        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
                BlockDivision[NumOfBlocks] = assumed_norm_of_the_current_block;

                //cout << "enter the next block" << endl;

                //cout << "1: filtered_list is:" << endl;
                //for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                //    cout << it1->get_sc_prod() << endl;

                //cout << "max_length_of_the_current_block = " << max_length_of_the_current_block << endl;


                //in case non of the list-elements from the new block will end up in filtered_list, the last_element of this block is the last_element of the prev. block
                //last_elements[NumOfBlocks] = last_elements[NumOfBlocks-1];
                /*
                last_elements[NumOfBlocks] = filtered_list.end();
                if (num_of_elements_in_filter[NumOfBlocks-1]==0)
                {
                    if (NumOfBlocks==1)
                        last_elements[NumOfBlocks-1] = filtered_list.end();
                    else
                       last_elements[NumOfBlocks-1] = last_elements[NumOfBlocks-2];
                }
                 */

        }


        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);

        // scaling factor for <px1>: ||p|| * ||x_1||
        //float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        float scale = sqrt ( (float)pApprox.get_approx_norm2() * (float)(it->get_approx_norm2()) );

        //cout  << "true_inner_product_px1 / scale " << (float)true_inner_product_px1 / scale << endl;

        //if true_inner_product_px1 is large enough, put it in filtered_list
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {

            //cout << "true_inner_product_px1 = " << true_inner_product_px1 << endl;

            typename FilteredListType::iterator it_filter =  filtered_list.begin();
            //typename FilteredListType::iterator it_comparison_flip_filter = filtered_list.end();

            int block_counter = 0;

            //cout << "filtered_list is:" << endl;
            //for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
            //        cout << it1->get_sc_prod() << endl;


            //loop over all elements x2 from the filtered list to 3-reduc (p, x1, x2)
            while ( it_filter != filtered_list.cend())
            {

                // OR take the length from BlockDivision[appendixCounter]
                LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = (it_filter->getApproxVector()).get_approx_norm2();



                //in res_upper, we store the upper bound on the inner-product of px2 for x2 from filtered_list
                //float res_upper = 6000;

                //Compute_one_third_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, res_upper);

                //cout << "res_upper = " << res_upper << endl;

                //check if we even need to iterate over the current block of filtered_list
                //if (it_filter->get_sc_prod() > res_upper)
                //{


                    ApproxTypeNorm2 true_inner_product_x1x2;

                    typename FilteredListType::iterator it_block =  it_filter;

                    int inblock_counter = 0;


                    //cout << "num_of_elements_in_filter[block_counter] = " << num_of_elements_in_filter[block_counter] << endl;

                    //while(it_block!=filtered_list.end() && it_block->get_sc_prod() > res_upper && inblock_counter < num_of_elements_in_filter[block_counter])
                    //{

                            //cout << "it_block->get_sc_prod() = " << it_block->get_sc_prod() << endl;

                            //------------------check for reduction-------------------

                            predict = LatticeApproximations::Compare_Sc_Prod_3red(it_block->getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);

                            //cout << "REDUCTION CHECK" << endl;

                            //sgn1 and sgn2 store the signs s.t. || p + sgn1*x1 + sgn*x2 || < || p ||
                            // these values are computed in check3red_signs
                            int sgn1 = 0;
                            int sgn2 = 0;


                            //in fistered_list we store abs(<p, x1>); retrieve the true sign of x1 stored in filtered_list
                            int sgn_px1 = 1;
                            int sgn_px2 = 1 ;
                            int sgn_x1x2 = 1;


                            if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                            if (!it_block->get_sign()) sgn_px2 = -1;
                            if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;

                            //check if reduction
                            if (GaussSieve::check3red_signs(p, *(it.access_details()), (it_block->getApproxVector()).get_details(), sgn_px1, sgn_px2, sgn_x1x2, sgn1, sgn2))
                            {

                                //perfrom actual reduction
                                p = GaussSieve::perform3red(p, *(it.access_details()), (it_block->getApproxVector()).get_details(), sgn1, sgn2);

                                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);

                                if(p.norm2!=0)
                                    main_queue.push(p);


                                //cout << "reduced p: RE-START with p on norm2 = " << p.norm2 << endl;
                                //assert(false);


                                filtered_list.clear();
                                return;
                            }

                            //else
                            //    cout << "check3red is false" << endl;



                            //inblock_counter++;
                            //cout << "inblock_counter " << inblock_counter << endl;
                            //++it_block;

                    //} //end of iteration over the block


                //} // end of if (it_filter->get_sc_prod() > res_upper)

                //cout << "last_elements[block_counter]->get_sc_prod() " << last_elements[block_counter]->get_sc_prod() << endl;
                //it_filter = last_elements[block_counter];


                //if(it_filter!=filtered_list.end())
                //   cout << "change it_filter to " << last_elements[block_counter]->get_sc_prod() << endl;

                //block_counter++;

                //cout << "block_counter = " << block_counter << endl;
                //if(it_filter!=filtered_list.end())
                    it_filter++;

            } //end of the loop over filtered_list

            //cout << "NumOfBlocks = " << NumOfBlocks << endl;
            //cout << "num_of_elements_in_filter[NumOfBlocks] = " << num_of_elements_in_filter[NumOfBlocks] << endl;


            bool px1_sign;
            if (true_inner_product_px1>0)
                px1_sign = true;
            else
                px1_sign = false;

            FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
            filtered_list.insert(filtered_list.end(),new_filtered_point);
            num_of_elements_in_filter[NumOfBlocks]++;


        } //if (px1 is larger enough)

        ++it;
    }

    //INSERT p
    main_list.insert_before(it_comparison_flip,p);
    ++current_list_size;
    //cout << "list_size = " <<current_list_size << endl;
    if(update_shortest_vector_found(p))
    {
        if(verbosity>=2)
        {
            cout << "New shortest vector found. Norm2 = " << get_best_length2() << endl;
        }
    }

    /*
    cout << "p is inserted" << endl;
    cout << "sizes of appendices are : "<< endl;
    int i=0;
    while (i<=NumOfBlocks)
    {
        cout << appendices[i].size() << " ";
        i++;
    }
    */

     it =it_comparison_flip; //points to the next after p list-element

     //assumed_norm_of_the_current_block = it->get_approx_norm2();
     //max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
     //NumOfBlocks++;

     bool reduced_x1 = false; //if true, do not increase ++it and do not include into filtered_list

     //now we are reducing *it
     while (it!=main_list.cend())
     {
        reduced_x1 = false;
        //cout << "2: consider a list element of approx norm = " << it->get_approx_norm2() << endl;


        //--------------------------------2-red-------------------------------
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );

        // preform 2-reduction
        if(predict){


        ++number_of_exact_scprods;
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
        {
                LatticePoint<ET> reduced = GaussSieve::perform2red(p, *(it.access_details()), scalar);

                //put reduced back into the queue
                if (reduced.norm2!=0)
                    main_queue.push(reduced);


                //cout << "put x1 of size " << reduced.norm2 << " into the queue " << endl;

                it = main_list.erase(it); //also makes ++it
                --current_list_size;
                reduced_x1 = true;
        }
        else
                ++number_of_mispredictions;
        }


        //--------------------------------3-red-------------------------------

        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
                BlockDivision[NumOfBlocks] = assumed_norm_of_the_current_block;

                //cout << "enter the next block" << endl;

                //cout << "filtered_list is:" << endl;
                //for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                //    cout << it1->get_sc_prod() << endl;

                //cout << "max_length_of_the_current_block = " << max_length_of_the_current_block << endl;


                //in case non of the list-elements from the new block will end up in filtered_list, the last_element of this block is the last_element of the prev. block
                //last_elements[NumOfBlocks] = last_elements[NumOfBlocks-1];


        }

        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);

        // scaling factor for <px1>: ||p|| * ||x_1||
        //float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        float scale = sqrt ( (float)pApprox.get_approx_norm2() * (float)(it->get_approx_norm2()) );


        //if true_inner_product_px1 is large enough, put it in filtered_list
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {
            //cout  << "true_inner_product_px1 " << true_inner_product_px1 << endl;
            typename FilteredListType::iterator it_filter =  filtered_list.begin();
            typename FilteredListType::iterator it_comparison_flip_filter = filtered_list.end();

            int block_counter = 0;
            bool position_is_found = false;

            //------for debug------

            /*
            cout << "2: filtered_list is:" << endl;
            for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                    cout << it1->get_sc_prod() << endl;

            cout << "2: Number of elements in filtered_list per block: " << endl;
            int i=0;
            while (i<=NumOfBlocks)
            {
                cout << num_of_elements_in_filter[i] << " ";
                i++;
            }
            cout <<endl;
            */
            //----------------------


            while (!reduced_x1 && it_filter != filtered_list.cend())
            {

                LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = (it_filter->getApproxVector()).get_approx_norm2();



                //in res_upper, we store the upper bound on the inner-product of px2 for x2 from filtered_list
                float res_upper = 6000;

                //Compute_one_third_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, res_upper);

                //cout << "res_upper = " << res_upper << endl;


                //check if we even need to iterate over the current block of filtered_list
                //if (it_filter->get_sc_prod() > res_upper)
                //{


                    ApproxTypeNorm2 true_inner_product_x1x2;

                    typename FilteredListType::iterator it_block =  it_filter;

                    int inblock_counter = 0;


                    //while(it_block!=filtered_list.end() && it_block->get_sc_prod() > res_upper && inblock_counter < num_of_elements_in_filter[block_counter])
                    //{

                            //------------------check for reduction-------------------

                            predict = LatticeApproximations::Compare_Sc_Prod_3red(it_block->getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);

                            //cout << "REDUCTION CHECK" << endl;

                            //sgn1 and sgn2 store the signs s.t. || p + sgn1*x1 + sgn*x2 || < || p ||
                            // these values are computed in check3red_signs
                            int sgn1 = 0;
                            int sgn2 = 0;


                            //in fistered_list we store abs(<p, x1>); retrieve the true sign of x1 stored in filtered_list
                            int sgn_px1 = 1;
                            int sgn_px2 = 1 ;
                            int sgn_x1x2 = 1;


                            if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                            if (!it_block->get_sign()) sgn_px2 = -1;
                            if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;

                            //check if reduction
                            if (GaussSieve::check3red_signs(*(it.access_details()), p, (it_block->getApproxVector()).get_details(), sgn_px1, sgn_px2, sgn_x1x2, sgn1, sgn2))
                            {

                                //perfrom actual reduction
                                LatticePoint<ET> reduced = GaussSieve::perform3red(*(it.access_details()), p, (it_block->getApproxVector()).get_details(), sgn1, sgn2);

                                if(reduced.norm2!=0)
                                        main_queue.push(reduced);

                                it = main_list.erase(it); //also makes ++it
                                --current_list_size;

                                //cout << "reduced x1" << endl;
                                //assert(false);

                                //break both while-loops and set the flag not to insert *it into filtered_list
                                reduced_x1 = true;
                                break;
                            }

                            //else
                                //cout << "check3red is false" << endl;


                            //--------------------- end of reduction's check ----------------------


                            //inblock_counter++;
                            //++it_block;

                    //} //end of iteration over the block


                    //insert from the appendix if it_block now points to a smaller element
                    /*
                    if ( ( !appendices[block_counter].empty() && it_block!=filtered_list.end() && (appendices[block_counter].top()).get_sc_prod() < it_block->get_sc_prod() ) )
                    {
                        FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> poped_from_queue =  appendices[block_counter].top();
                        appendices[block_counter].pop();


                        cout << "about to append from appendix "  << endl;

                        //BUG HERE: it_block does not always point to the end of the block
                        filtered_list.insert(it_block,poped_from_queue);
                        num_of_elements_in_filter[block_counter]++;
                        last_elements[block_counter] = --it_block; //TODO: TO CHECK.
                        cout << "num_of_elements_in_filter[block_counter] " << num_of_elements_in_filter[block_counter] << endl;
                        cout << "last_elements[block_counter] is upd to " << last_elements[block_counter]->get_sc_prod() << endl;
                        ++it_block;

                        cout << "2: filtered_list is:" << endl;
                        for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
                            cout << it1->get_sc_prod() << endl;
                    }
                     */

                //} // end of if (it_filter->get_sc_prod() > res_upper)

                //cout << "last_elements[block_counter]->get_sc_prod() " << last_elements[block_counter]->get_sc_prod() << endl;
                //it_filter = last_elements[block_counter];

                //if(it_filter!=filtered_list.end())
                //   cout << "change it_filter to " << last_elements[block_counter]->get_sc_prod() << endl;

                //block_counter++;

                //cout << "block_counter = " << block_counter << endl;
                //if(it_filter!=filtered_list.end())
                    it_filter++;

            } //end of the loop over filtered_list

            if (!reduced_x1)
            {
                bool px1_sign;
                if (true_inner_product_px1>0)
                    px1_sign = true;
                else
                    px1_sign = false;

                FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
                filtered_list.insert(filtered_list.end(),new_filtered_point);
                num_of_elements_in_filter[NumOfBlocks]++;
            }


        } //if (px1 is larger enough)

        if (!reduced_x1)
        ++it;

    }

    /*
    if (filtered_list.size()>60)
    {

         cout << "filtered_list is:" << endl;
         for (auto it1 = filtered_list.cbegin(); it1!=filtered_list.cend(); ++it1)
         cout << it1->get_sc_prod() << endl;

         cout << "Number of elements in filtered_list per block: " << endl;
         int i=0;
         while (i<=NumOfBlocks)
         {
         cout << num_of_elements_in_filter[i] << " ";
         i++;
         }

        assert(false);

    }
    */


     filtered_list.clear();
}



template<class ET>
void Sieve<ET,false>::sieve_3_iteration_new (LatticePoint<ET> &p)
{

    if (p.norm2==0) return;
    ApproxLatticePoint<ET,false> pApprox (p);

    //cout << "------------" << endl;
    //cout << "Run iteration on p of approx norm " << pApprox.get_approx_norm2() << endl;

    int const n = get_ambient_dimension();

    ET scalar;

    typename MainListType::Iterator it_comparison_flip=main_list.cend();

    //if abs( <p, x1> / (|p||x1|) ) < px1bound, do not put it in the filtered_list
    float  px1bound = 0.29; // TO ADJUST

    //'promissingness'
    float x1x2=.33;

    // length_factor determines the difference (mult) between the legnth of the i-th and (i+1)-st blocks
    float length_factor = 2.0; //TO ADJUST
    //float length_factor =10.0; //to debug assume we have only 1 block


    // store target inner-products
    //float px1=.0;
    //float px2=.0;

    //number of blocks
    int NumOfBlocks = 1;

    auto it = main_list.cbegin();

    typename MainListType::Iterator first_element_of_the_block = main_list.cbegin();
    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_the_current_block =  first_element_of_the_block->get_approx_norm2();
    LatticeApproximations::ApproxTypeNorm2 max_length_of_the_current_block = floor(length_factor * assumed_norm_of_the_current_block + 1);


    // main loop over main_list
    // this while-loop runs until the p is longer



    while (it!=main_list.cend())
    {

        //cout << "1: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //check if p is still of max. length
        if (p.norm2 < it.get_true_norm2())
        {
            //cout << "reached the position to insert p" << endl;
            it_comparison_flip = it;
            break;
        }


         //--------------------------------2-red-------------------------------

        //check 2-red
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );


        if(predict){


        ++number_of_exact_scprods;

         // preform 2-reduction
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
            {
                p = GaussSieve::perform2red(p, *(it.access_details()), scalar);
                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);

                //put p back into the queue and break
                if (p.norm2!=0)
                    main_queue.push(p);
                //cout << "put p of size " << p.norm2 << " into the queue " << endl;
                //cout << "pApprox has norm2 " << pApprox.get_approx_norm2() <<  endl;
                //cout << "2-reduced p, break all the loops" << endl;
                return;
                //inner_loop = false;
                //break;
            }
            else
                ++number_of_mispredictions;
        }


        //--------------------------------3-red-------------------------------



        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;

                //cout << "enter the next block" << endl;
                //cout << "assumed_norm_of_current_block = " << assumed_norm_of_the_current_block << endl;
                //cout <<"max_length_of_current_block = " << max_length_of_the_current_block << endl;
                //cout << "filtered_list is:" << endl;

                //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;

                //cout << endl;


        }

        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);

        // scaling factor for <px1>: ||p|| * ||x_1||
        //float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));
        float scale = sqrt ( (float)pApprox.get_approx_norm2() * (float)(it->get_approx_norm2()) );

        //cout  << "true_inner_product_px1 / scale " << (float)true_inner_product_px1 / scale << endl;

        //if true_inner_product_px1 is large enough, put it in filtered_list
        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {

            //loop over all elements x2 from the filtered list to 3-reduce (p, x1, x2)
            auto it_filter = filtered_list2.cbegin();
            auto next_block = filtered_list2.cend();

            //bool filter_loop = true;
            while ( it_filter != filtered_list2.cend())
            {
                LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = get<0>(it_filter->first);


                pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair_bound = make_pair(assumed_norm_of_x1+1, 0);

                //typename FilteredListType2:: iterator next_block = filtered_list2.lower_bound(new_key_pair_bound );
                next_block = filtered_list2.lower_bound(new_key_pair_bound );

                //returns garbadge if there is only one block in the filtered list
                //LatticeApproximations::ApproxTypeNorm2 norm_of_the_next_block = get<0>(next_block->first);
                //cout << "in filter: " << "block-len = " <<  assumed_norm_of_x1 << " next_block = " << norm_of_the_next_block << endl;


                //in res_upper, we store the upper bound on the inner-product of px2 for x2 from filtered_list
                float res_upper = 0.0;

                // as the second argument pass the value assumed_norm_of_the_current_block
                //Compute_px1_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, x1x2, res_upper)
                //the result is stored in res_upper
                Compute_px1_bound(assumed_norm_of_x1, it->get_approx_norm2(), true_inner_product_px1, x1x2, res_upper);

                //cout <<  " res_upper = " << (LatticeApproximations::ApproxTypeNorm2)res_upper << endl;

                typename FilteredListType2:: iterator itup;
                auto it_tmp = it_filter;

                // if the bounds are too large, consider the next length-block
                // otherwise find itup s.t. all inner-products after itup are smaller (i.e. worse) than res_upper. Iterate up until itup;
                // it_filter should point to the element with the smallest inner-product within this block


                new_key_pair_bound = make_pair(assumed_norm_of_x1, (LatticeApproximations::ApproxTypeNorm2)res_upper);

                //find the iterator s.t. all elements < it have <px1> <= res_upper
                itup = filtered_list2.upper_bound(new_key_pair_bound );


                //if (it_tmp->first == itup->first )
                //{
                    //cout << "itup = " << get<1>(itup->first) << endl;

                    //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                    //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;

                    //cout << "itup is equal to it_filter, no iteration" << endl;

                //}
                //else
                //{
                    //debugging

                    //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                    //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;

                    //cout << "itup = " << get<1>(itup->first) << endl;
                    //cout << "it_filter = " << get<1>(it_filter->first) << endl;
                    ApproxTypeNorm2 true_inner_product_x1x2;

                    //while (it_filter != itup && it_filter!=filtered_list2.cend())
                    while(it_filter !=itup && it_filter !=next_block)
                    {

                        //cout << "it_filter = " << get<1>(it_filter->first) << endl;
                        //retrieve x2 from it_filters and compare x1x2
                        //predict = LatticeApproximations::Compare_Sc_Prod_3red((it_filter->second).getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);
                        predict = LatticeApproximations::Compare_Sc_Prod_3red((it_filter->second).getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);

                       // cout << " true_inner_product_x1x2 = " << true_inner_product_x1x2 << endl;

                        if (predict)
                        {
                            //cout << "REDUCTION" << endl;

                            //sgn1 and sgn2 store the signs s.t. || p + sgn1*x1 + sgn*x2 || < || p ||
                            // these values are computed in check3red_signs
                            int sgn1 = 0;
                            int sgn2 = 0;


                            //in fistered_list we store abs(<p, x1>); retrieve the true sign of x1 stored in filtered_list
                            int sgn_px1 = 1;
                            int sgn_px2 = 1 ;
                            int sgn_x1x2 = 1;


                            if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                            if (!(it_filter->second).get_sign()) sgn_px2 = -1;
                            if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;

                            // check
                            //cout << "true_inner_product_px1 " << true_inner_product_px1 << " sgn_px1 " << sgn_px1 << endl;
                            //cout << "true_inner_product_px2 " << (it_filter->second).get_sign() << " sgn_px2 " << sgn_px2 << endl;
                            //cout << "true_inner_product_x1x2 " << true_inner_product_x1x2 << " sgn_x1x2 " << sgn_x1x2 << endl;




                            //check if reduction
                            if (GaussSieve::check3red_signs(p, *(it.access_details()), ((it_filter->second).getApproxVector()).get_details(), sgn_px1, sgn_px2, sgn_x1x2, sgn1, sgn2))
                            {

                                //perfrom actual reduction
                                p = GaussSieve::perform3red(p, *(it.access_details()), ((it_filter->second).getApproxVector()).get_details(), sgn1, sgn2);

                                if (p.norm2!=0) pApprox = static_cast< ApproxLatticePoint<ET,false> >(p);

                                if(p.norm2!=0)
                                    main_queue.push(p);

                                //if(get_best_length2() == 3135573)
                                //{
                                    //cout << "reduced p: RE-START with p on norm2 = " << p.norm2 << endl;

                                //}
                                //cout << "reduced p: RE-START with p on norm2 = " << p.norm2 << endl;
                                //assert(false);


                                filtered_list2.clear();
                                return;

                            }

                        }

                        ++it_filter;
                    }   //while -loop over one length-block

                    //cout << "reached it_up" << endl;

                //} //else-condition

                //go the next length_block
                it_filter = next_block;

                //cout << "go to length " << get<0>(it_filter->first) << endl;

            } // end of while it_filter - loop


            //cout << "filter while-loop is finished "<< endl;

            //insert *it into filtered_list2

            bool px1_sign;
            if (true_inner_product_px1>0)
                px1_sign = true;
            else
                px1_sign = false;

            FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
            //ApproxLatticePoint<ET> point = *it;
            //FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(&point, abs(true_inner_product_px1), px1_sign);

            pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair;
            new_key_pair = make_pair(assumed_norm_of_the_current_block, abs(true_inner_product_px1));

            //cout << "about to insert into filtered_list" << endl;
            //TODO: provide hint to emplace
            filtered_list2.emplace_hint(next_block, new_key_pair, new_filtered_point);


            //cout << "input x with px = " << true_inner_product_px1  << " and px1_sign = " << px1_sign << endl;

            //if(get_best_length2() == 3135573)
            //{
            //    cout <<"filtered_list.size = " << filtered_list2.size() << endl;
            //}
            //if (filtered_list2.size() > 22)
            //    assert(false);

        }

        else
        {
            //cout << "1: scaled px1 = " <<  (float)true_inner_product_px1 / scale << " is smaller than " << px1bound << endl;
        }



        ++it;
    }

    //if(get_best_length2() == 3135573)
   //{
      //  cout << "INSERT p of norm " << p.norm2 << endl;
    //}
    //cout << "INSERT p of norm " << p.norm2 << endl;
    //cout << &it_comparison_flip << endl;

    //INSERT p
    main_list.insert_before(it_comparison_flip,p);
    ++current_list_size;
    //cout << "list_size = " <<current_list_size << endl;
    if(update_shortest_vector_found(p))
    {
        if(verbosity>=2)
        {
            cout << "New shortest vector found. Norm2 = " << get_best_length2() << endl;
        }
    }



    //filtered_list2.clear();

    //keep the old filtered_list2
    //start iteration from it_comparison_flip

    it =it_comparison_flip; //points to the next after p list-element


    //cout << "start the lower part of the list" << endl;

    assumed_norm_of_the_current_block = it->get_approx_norm2();
    max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);

    bool reduced_x1 = false; //if true, do not increase ++it and do not include into filtered_list

    //now we are reducing *it
    while (it!=main_list.cend())

    {
        reduced_x1 = false;
        //cout << "2: consider a list element of approx norm = " << it->get_approx_norm2() << endl;
        //assert(false);


        //--------------------------------2-red-------------------------------
        ++number_of_scprods;
        bool predict = LatticeApproximations::Compare_Sc_Prod(pApprox,*it,it->get_approx_norm2(),2* it->get_length_exponent()-1,n   );

        // preform 2-reduction
        if(predict){


        ++number_of_exact_scprods;
        if ( GaussSieve::check2red_new(p, *(it.access_details()), scalar) )
        {
                LatticePoint<ET> reduced = GaussSieve::perform2red(p, *(it.access_details()), scalar);

                //put reduced back into the queue
                if (reduced.norm2!=0)
                    main_queue.push(reduced);


                //cout << "put x1 of size " << reduced.norm2 << " into the queue " << endl;

                it = main_list.erase(it); //also makes ++it
                --current_list_size;
                reduced_x1 = true;
        }
        else
                ++number_of_mispredictions;
        }


        //--------------------------------3-red-------------------------------

        //check if we reached the next block; if yes, re-compute max_length_of_the_current_block
        if (it->get_approx_norm2() > max_length_of_the_current_block)
        {
                assumed_norm_of_the_current_block = it->get_approx_norm2();
                max_length_of_the_current_block =floor(length_factor * assumed_norm_of_the_current_block + 1);
                NumOfBlocks++;
                //cout << "enter the next block" << endl;
                //cout << "assumed_norm_of_current_block = " << assumed_norm_of_the_current_block << endl;
                //cout <<"max_length_of_current_block = " << max_length_of_the_current_block << endl;
                //cout << "filtered_list is:" << endl;

                //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;

                //cout << endl;


        }

        ApproxTypeNorm2 true_inner_product_px1 = compute_sc_prod(pApprox.get_approx(), it->get_approx(), n);
        float scale = (float)(pow(pApprox.get_approx_norm2(), 0.5)) * (float)(pow (it->get_approx_norm2(), 0.5));


        if (abs((float)true_inner_product_px1 / scale)>px1bound)
        {
                auto it_filter = filtered_list2.cbegin();

                auto next_block = filtered_list2.cend();

                while ( it_filter != filtered_list2.cend() && !reduced_x1 )
                {
                    LatticeApproximations::ApproxTypeNorm2 assumed_norm_of_x1 = get<0>(it_filter->first);
                    pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair_bound = make_pair(assumed_norm_of_x1+1, 0);

                    next_block = filtered_list2.lower_bound(new_key_pair_bound );

                    //returns garbadge if there is only one block in the filtered list
                    //LatticeApproximations::ApproxTypeNorm2 norm_of_the_next_block = get<0>(next_block->first);
                    //cout << "in filter: " << "block-len = " <<  assumed_norm_of_x1 << " next_block = " << norm_of_the_next_block << endl;


                    float res_upper = 0.0;

                    // as the second argument pass the value assumed_norm_of_the_current_block
                    //Compute_px1_bound(assumed_norm_of_x1, assumed_norm_of_the_current_block, true_inner_product_px1, x1x2, res_upper)
                    Compute_px1_bound(pApprox.get_approx_norm2(), assumed_norm_of_x1, true_inner_product_px1, x1x2, res_upper);
                    //cout <<  " res_upper = " << (LatticeApproximations::ApproxTypeNorm2)res_upper << endl;

                    typename FilteredListType2:: iterator itup;
                    auto it_tmp = it_filter;

                    new_key_pair_bound = make_pair(assumed_norm_of_x1, (LatticeApproximations::ApproxTypeNorm2)res_upper);
                    itup = filtered_list2.upper_bound(new_key_pair_bound );

                    if (it_tmp->first == itup->first )
                    {
                        //cout << "itup = " << get<1>(itup->first) << endl;

                        //for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                        //    cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;

                        //cout << "itup is equal to it_filter, no iteration" << endl;

                    }
                    else
                    {
                        //if (filtered_list2.size()>72)
                        //{
                        //    for (auto it1 = filtered_list2.cbegin(); it1!=filtered_list2.cend(); ++it1)
                        //        cout << get<0>(it1->first) << " " << get<1>(it1->first) << endl;

                        //    assert(false);
                        //}
                        //cout << "itup = " << get<1>(itup->first) << endl;
                        //cout << "it_filter = " << get<1>(it_filter->first) << endl;
                        ApproxTypeNorm2 true_inner_product_x1x2;

                        while (it_filter != itup && it_filter!=next_block)
                        {
                            predict = LatticeApproximations::Compare_Sc_Prod_3red((it_filter->second).getApproxVector(), *it, n, x1x2, true_inner_product_x1x2);

                            if (predict)
                            {
                                //cout << "REDUCTION FOR x1" << endl;

                                int sgn1 = 0;
                                int sgn2 = 0;

                                int sgn_px1 = 1;
                                int sgn_px2 = 1 ;
                                int sgn_x1x2 = 1;


                                if (true_inner_product_px1 < 0 ) sgn_px1 = -1;
                                if (!(it_filter->second).get_sign()) sgn_px2 = -1;
                                if (true_inner_product_x1x2 < 0 ) sgn_x1x2 = -1;

                                // check
                                //cout << "true_inner_product_px1 " << true_inner_product_px1 << " sgn_px1 " << sgn_px1 << endl;
                                //cout << "true_inner_product_px2 " << (it_filter->second).get_sign() << " sgn_px2 " << sgn_px2 << endl;
                                //cout << "true_inner_product_x1x2 " << true_inner_product_x1x2 << " sgn_x1x2 " << sgn_x1x2 << endl;


                                //check if reduction
                                if (GaussSieve::check3red_signs(*(it.access_details()), p, ((it_filter->second).getApproxVector()).get_details(), sgn_px1, sgn_x1x2, sgn_px2, sgn1, sgn2))
                                {

                                    //perfrom actual reduction
                                    LatticePoint<ET> reduced = GaussSieve::perform3red(*(it.access_details()), p, ((it_filter->second).getApproxVector()).get_details(), sgn1, sgn2);

                                    if(reduced.norm2!=0)
                                        main_queue.push(reduced);

                                    it = main_list.erase(it); //also makes ++it
                                    --current_list_size;


                                    //break both while-loops and set the flag not to insert *it into filtered_list
                                    reduced_x1 = true;
                                    break;

                                    //cout << "reduced x1:  now  x1 is if norm2 = " << p.norm2 << endl;
                                    //assert(false);

                                }

                            } //if(predict)

                            ++it_filter;
                        } //while -loop over one length-block
                    }//else

                    it_filter = next_block;
                    //cout << "go to length " << get<0>(it_filter->first) << endl;

                }

                //cout << "filter while-loop is finished "<< endl;
                if (!reduced_x1)
                {
                    bool px1_sign;
                    if (true_inner_product_px1>0)
                        px1_sign = true;
                    else
                        px1_sign = false;

                    FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(*it, abs(true_inner_product_px1), px1_sign);
                    //ApproxLatticePoint<ET> point = *it;
                    //FilteredPoint<ET, LatticeApproximations::ApproxTypeNorm2> new_filtered_point(&point, abs(true_inner_product_px1), px1_sign);
                    pair <LatticeApproximations::ApproxTypeNorm2, LatticeApproximations::ApproxTypeNorm2> new_key_pair;
                    new_key_pair = make_pair(assumed_norm_of_the_current_block, abs(true_inner_product_px1));

                    //if(get_best_length2() == 3135573)
                    //{
                    //    cout << get<0>(new_key_pair) << " " << get<1>(new_key_pair) << endl;
                        //cout << "about to insert into filtered_list from lower-part" << endl;

                    //}

                    filtered_list2.emplace_hint(next_block, new_key_pair, new_filtered_point);
                    //cout << "filtered_list2.size" << filtered_list2.size() << endl;


                }

        } //if

        else
        {
            //if(get_best_length2() == 3135573)
            //{
                //cout << "scaled px1 = " <<  (float)true_inner_product_px1 / scale << " is smaller than " << px1bound << endl;
            //}
        }

        //otherwise ++it was already performed during main_list.erase();
        if (!reduced_x1)
            ++it;
    }
    //cout << "finish the main while-loop, clean filtered_list" << endl;
    filtered_list2.clear();

}
