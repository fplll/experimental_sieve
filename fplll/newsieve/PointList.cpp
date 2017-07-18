#pragma deprecated

//continuation from PointList.h
//no include guards or anything needed.

template<class ET>
GaussList<ET,true,-1>::GaussList() :
    mutex_currently_writing(),
    start_sentinel_node (new Node),
    end_sentinel_node   (new Node)
{
        start_sentinel_node->next_node=end_sentinel_node;
        end_sentinel_node->prev_node=start_sentinel_node;
        start_sentinel_node->nodestatus=static_cast<int>(Node::NodeStatus::is_first_node);
        end_sentinel_node->nodestatus  =static_cast<int>(Node::NodeStatus::is_last_node);
}



template<class ET>
GaussList<ET,true,-1>::~GaussList() //called when only one thread is running
{
    //atomic_thread_fence(std::memory_order_acquire);   begin() already does that.
    Iterator next=cbegin();
    Iterator cur (start_sentinel_node);
    while(!next.is_end())
    {
        delete cur.p;
        cur=next;
        ++next;
    }
    delete cur.p;
    delete next.p;
}


template<class ET>
typename GaussList<ET,true,-1>::Iterator GaussList<ET,true,-1>::insert_before(Iterator pos, ExactType const & val) //TODO, type of val might change. Uses emplacement.
{
    Node* newnode = new Node;
    newnode->datum_ptr = new DataType (val);
    Node* nextgood =pos.p; //we work directly with the underlying pointer, not using the iterator, since we do not need atomic loads to traverse here.
    assert(nextgood!=nullptr);
    mutex_currently_writing.lock();
    while(nextgood->is_marked_for_deletion())
    {
        nextgood=nextgood->next_node;
    } //.load(memory_order_relaxed);}
    Node* preced=nextgood->prev_node;
    newnode->next_node=nextgood;
    newnode->prev_node=preced;
    nextgood->prev_node=newnode;
//until here, no other thread can observe our writes, since we did not publish the pointer.
    preced->next_node.store(newnode,memory_order_release);
//this one can be observed (even prior to releasing the lock), and we have to make sure that other non-mutex-protected threads that see this write also see the (non-atomic) writes to newnode.
    mutex_currently_writing.unlock();
    return static_cast<Iterator> (newnode);
}

template <class ET>
void GaussList<ET,true,-1>::unlink(Iterator const & pos, GarbageBin<DataType> &gb)
{
    assert(pos.p!=nullptr);
    {//locked part
        assert(!(pos.p->is_sentinel_node()));
        //std::lock_guard<std::mutex> writelock(mutex_currently_writing);
        mutex_currently_writing.lock();
        if(pos.p->is_plain_node()) //otherwise, already deleted.
        {
            pos.p->nodestatus=static_cast<int>(Node::NodeStatus::is_to_be_deleted);
            NodePointer nextpos=pos.p->next_node; //.load(memory_order_relaxed); //All writes are within locks anyway.
            nextpos->prev_node=pos.p->prev_node;
            pos.p->prev_node->next_node.store(nextpos,memory_order_relaxed); //relaxed should be fine!!!
            mutex_currently_writing.unlock();
            //Put in garbage bin. //TODO: More clever garbage bin, requires changing structs and global counters.
            gb.push(pos.p);
        }
        else
        {
            mutex_currently_writing.unlock();
        }
    }//end of locked part
    return;
}

template<class ET>
void GaussList<ET,true,-1>::sort()
{
    //We perform an in-place merge-sort, using only the prev-links at first.
    //The forward-links are set in the last pass. This is better compatible with multi-threading.
    //We directly work with NodePointers rather than iterators.
    //unsigned long long int current_block_size = 0;
    std::stack< NodePointer > trace_back;
    unsigned long long int pos = 0;
    mutex_currently_writing.lock();
    NodePointer process_next = end_sentinel_node -> prev_node;
    if (process_next == start_sentinel_node)
    {
        cerr << "Warning: sorting empty list.";
        mutex_currently_writing.unlock();
        return;
    }

    //We maintain the following status invariant at the beginning of the while loop.
    //The last pos (non-sentinel) nodes originally reachable by back-traversing from the end are partially ordered in the following sense:
    //Let pos = 2^e_1 + 2^e_2 + ... + 2^e_n  , e_1 > e_2 > ...,
    //Subdivide the last pos Nodes into disjoint blocks B_i of sizes 2^e_i.
    //trace_back contains n NodePointer's r_1, ... r_n, with r_i pointing to the largest element of the B_i. r_i's for smaller B_i's are on top of the stack.
    //The prev_nodes within each B_i point to the next smaller element.
    //There is no guarantee for prev_node pointers traversing block boundaries.
    //process_next contains the pos+1'th element reachable by back-traversing from the end (in the original status).


    while(process_next!=start_sentinel_node)
    {
        assert(!process_next->is_marked_for_deletion());
        NodePointer current_batch = process_next; //to be appended to our stack of blocks.
        process_next  = current_batch->prev_node; //for the next iteration.
        ++pos;
        //We now append the next block, merging blocks of equal size on top of the stack.
        for(int i = 1; pos % (1ULL << i  ) == 0;++i )
        {
            unsigned long long cur_size   = 1ULL<<(i-1); //remaining size of the current block
            unsigned long long stack_size = 1ULL<<(i-1); //remaining size of the block on top of the stack
            NodePointer new_batch = current_batch; //entry point for merged batch.
            if ( *(current_batch->datum_ptr) < *(trace_back.top()->datum_ptr)  )
            {
                new_batch = trace_back.top();
                --stack_size;
                trace_back.top() = new_batch->prev_node;
            }
            else
            {
                new_batch = current_batch;
                --cur_size;
                current_batch = current_batch -> prev_node;
            }

            //correctly link merged block:

            for(NodePointer merge_batch = new_batch; cur_size+stack_size>0;)
            {
                if(cur_size==0)
                {
                    merge_batch->prev_node = trace_back.top();
                    merge_batch = trace_back.top();
                    trace_back.top() = trace_back.top()->prev_node;
                    --stack_size;
                    continue;
                }
                if(stack_size==0)
                {
                    merge_batch->prev_node = current_batch;
                    merge_batch = current_batch;
                    current_batch = current_batch -> prev_node;
                    --cur_size;
                    continue;
                }
                if( *(current_batch->datum_ptr) < *(trace_back.top()->datum_ptr) )
                {
                    merge_batch->prev_node = trace_back.top();
                    merge_batch = trace_back.top();
                    trace_back.top() = trace_back.top()->prev_node;
                    --stack_size;
                }else
                {
                    merge_batch->prev_node = current_batch;
                    merge_batch = current_batch;
                    current_batch = current_batch -> prev_node;
                    --cur_size;
                }
            }
            //merged block has correct backlinks now.
            current_batch = new_batch;
            trace_back.pop();
        }
        trace_back.push(current_batch);
    }
//    mutex_currently_writing.unlock(); return;
    //now the stack contains pointers to blocks, making up all of our list. pos is equal to list-length
    unsigned long long int top_block_size = 1;
    for( ; ( pos%(2*top_block_size) ) ==0 ;  )
    {
        top_block_size *= 2;
    } //top_block_size now gives the size of the last block
    unsigned long long int size_remaining = pos -top_block_size; //rest of the blocks.
    assert(!( trace_back.empty()));
    while(trace_back.size()!=1)
    {
        unsigned long long int next_block_size=1;
        for( ; (size_remaining % (2*next_block_size))==0; )
        {
            next_block_size *= 2;
        } //next_block_size is size of next block.

        //copied from above
        NodePointer current_batch = trace_back.top();
        trace_back.pop();
        unsigned long long cur_size   = top_block_size;  //remaining size of the current block
        unsigned long long other_size = next_block_size; //remaining size of the block on top of the stack
        NodePointer new_batch = nullptr; //entry point for merged batch.
        if ( *(current_batch->datum_ptr) < *(trace_back.top()->datum_ptr)  )
        {
            new_batch = trace_back.top();
            --other_size;
            trace_back.top() = new_batch->prev_node;
        }
        else
        {
            new_batch = current_batch;
            --cur_size;
            current_batch = current_batch -> prev_node;
        }

        for(NodePointer merge_batch = new_batch; cur_size+other_size>0;)
        {
            if(cur_size==0)
            {
                merge_batch->prev_node = trace_back.top();
                merge_batch = trace_back.top();
                trace_back.top() = trace_back.top()->prev_node;
                --other_size;
                continue;
            }
            if(other_size==0)
            {
                merge_batch->prev_node = current_batch;
                merge_batch = current_batch;
                current_batch = current_batch -> prev_node;
                --cur_size;
                continue;
            }
            if( *(current_batch->datum_ptr) < *(trace_back.top()->datum_ptr) )
            {
                merge_batch->prev_node = trace_back.top();
                merge_batch = trace_back.top();
                trace_back.top() = trace_back.top()->prev_node;
                --other_size;
            }
            else
            {
                merge_batch->prev_node = current_batch;
                merge_batch = current_batch;
                current_batch = current_batch -> prev_node;
                --cur_size;
            }
        }
        trace_back.pop();
        trace_back.push(new_batch);
        top_block_size +=next_block_size;
        size_remaining -=next_block_size;
    }
    assert(size_remaining == 0);
    //We have a single block now.
    NodePointer target_next = end_sentinel_node;
    NodePointer iterate = trace_back.top();
    end_sentinel_node->prev_node = iterate;
    for(unsigned long long int i=0;i<pos;++i)
    {
        iterate->next_node.store(target_next,std::memory_order_release);
        target_next=iterate;
        iterate=iterate->prev_node;
    }
    target_next->prev_node = start_sentinel_node;
    start_sentinel_node->next_node.store(target_next,std::memory_order_release);

    mutex_currently_writing.unlock();
}


//iterators:

template<class ET>
GaussIterator<ET,true,-1> GaussIterator<ET,true,-1>::operator++(int) //postfix version
{
    auto tmp=p;
    ++(*this);
    return tmp;
};

template<class ET>
GaussIterator<ET,true,-1>& GaussIterator<ET,true,-1>::operator++() //prefix version
{
    p=p->next_node.load(memory_order_acquire);
    assert(p!=nullptr);
    return *this;
};

template<class ET>
GaussIterator<ET,true,-1>& GaussIterator<ET,true,-1>::operator--() //prefix version
{
    p=p->prev_node.load(memory_order_acquire);
    assert(p!=nullptr);
    return *this;
};










































//old implementation:

//
//template<class DT>
//MTListIterator<DT> MTListIterator<DT>::operator++(int) //postfix version
//{
//    auto tmp=p;
//    ++(*this);
//    return tmp;
//};
//
//template<class DT>
//MTListIterator<DT>& MTListIterator<DT>::operator++()
//{
//    p=p->next_node.load(memory_order_acquire);
//    assert(p!=nullptr);
//    return *this;
//}; //prefix version

/*
template<class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::cbefore_begin() const
{
    return start_sentinel_node; //no atomic load, because
}*/

/*
template<class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::cbegin() const
{
    return start_sentinel_node->next_node.load(std::memory_order_acquire);
};
*/

/*
template<class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::cend() const
{
    return end_sentinel_node;
};
*/

/*
template<class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::insert_before(Iterator const &pos, DT const &val)
{
    DT * const tmp = new DT(val);
    return enlist(pos,tmp);
};
*/

/*
template<class DT>
ListMultiThreaded<DT>::~ListMultiThreaded() //called when only one thread is running
{
    //atomic_thread_fence(std::memory_order_acquire);   begin() already does that.
    Iterator next=cbegin();
    Iterator cur (start_sentinel_node);
    while(!next.is_end())
    {
        delete cur.p;
        cur=next;
        ++next;
    }
    delete cur.p;
    delete next.p;
}
*/

/*
template<class DT>
ListMultiThreaded<DT>::ListMultiThreaded() : //called when only one thread is running
        mutex_currently_writing(),
        start_sentinel_node (new Node),
        end_sentinel_node   (new Node)
    {
        start_sentinel_node->next_node=end_sentinel_node;
        end_sentinel_node->prev_node=start_sentinel_node;
        start_sentinel_node->nodestatus=static_cast<int>(Node::NodeStatus::is_first_node);
        end_sentinel_node->nodestatus  =static_cast<int>(Node::NodeStatus::is_last_node);
    };
*/

/*
template <class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::enlist_before(MTListIterator<DT> const &pos, DT * const valref)
{
    Node* newnode = new Node;
    newnode->datum_ptr = valref;
    Node* nextgood =pos.p; //we work directly with the underlying pointer, not using the iterator, since we do not need atomic loads to traverse here.
    assert(nextgood!=nullptr);
    mutex_currently_writing.lock();
    while(nextgood->is_marked_for_deletion())
    {
        nextgood=nextgood->next_node;
    } //.load(memory_order_relaxed);}
    Node* preced=nextgood->prev_node;
    newnode->next_node=nextgood;
    newnode->prev_node=preced;
    nextgood->prev_node=newnode;
//until here, no other thread can observe our writes, since we did not publish the pointer.
    preced->next_node.store(newnode,memory_order_release);
//this one can be observed (even prior to releasing the lock), and we have to make sure that other non-mutex-protected threads that see this write also see the (non-atomic) writes to newnode.
    mutex_currently_writing.unlock();
    return static_cast<Iterator> (newnode);
}
*/

/*
template <class DT>
void ListMultiThreaded<DT>::unlink(Iterator const & pos, GarbageBin<DT> &gb)
{
    assert(pos.p!=nullptr);
    {//locked part
        assert(!(pos.p->is_sentinel_node()));
        //std::lock_guard<std::mutex> writelock(mutex_currently_writing);
        mutex_currently_writing.lock();
        if(pos.p->is_plain_node()) //otherwise, already deleted.
        {
            pos.p->nodestatus=static_cast<int>(Node::NodeStatus::is_to_be_deleted);
            NodePointer nextpos=pos.p->next_node; //.load(memory_order_relaxed); //All writes are within locks anyway.
            nextpos->prev_node=pos.p->prev_node;
            pos.p->prev_node->next_node.store(nextpos,memory_order_relaxed); //relaxed should be fine!!!
            mutex_currently_writing.unlock();
            //Put in garbage bin. //TODO: More clever garbage bin, requires changing structs and global counters.
            gb.push(pos.p);
        }
        else
        {
            mutex_currently_writing.unlock();
        }
    }//end of locked part
    return;
}
*/

template<class DT>
void GarbageBin<DT>::empty_trash()
{
    while(!this->empty())
    {
        delete this->front();
        this->pop();
    }
};

