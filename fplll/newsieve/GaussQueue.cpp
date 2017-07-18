#ifndef GAUSS_QUEUE_CPP
#define GAUSS_QUEUE_CPP

#include "GaussQueue.h"
#include "SieveGauss.cpp"
#include "ShiSampler.cpp"
#include "Sampler.cpp"

template<class ET,int nfixed> GaussQueue<ET,false,nfixed>::GaussQueue( Sieve<ET,false,nfixed> *caller_sieve)  //constructor
:
main_queue(),
gauss_sieve(caller_sieve),
sampler(nullptr)
{
    assert(caller_sieve!=nullptr);
    std::seed_seq seed{1,2,4}; //just some arbitrary numbers
    //sampler = new EllipticSampler<ET,false, std::mt19937_64, std::seed_seq> (seed);
    cout << "Initializing Sampler" << endl << flush;
    sampler = new ShiSampler<ET,false, std::mt19937_64, std::seed_seq,nfixed> (seed);
    cout << "Finished Initializing Sampler" << endl << flush;
    assert(sampler!=nullptr);
}

//template<class ET>
//GaussQueue<ET,true>::GaussQueue( Sieve<ET,true> *caller_sieve)
//:
//main_queue(),
//gauss_sieve(caller_sieve),
//queue_mutex()
////sampler(nullptr)
//{
//    assert(caller_sieve!=nullptr);
//}


template<class ET,int nfixed> typename GaussQueue<ET,false,nfixed>::RetType GaussQueue<ET,false,nfixed>::true_pop()
{
    if(main_queue.empty())
    {
        ++ (gauss_sieve->number_of_points_sampled);
        ++ (gauss_sieve->number_of_points_constructed);
//        return gauss_sieve->sampler->sample();
        assert(sampler!=nullptr);
        typename GaussQueue<ET,false,nfixed>::RetType ret = sampler->sample();
        return ret;
//        return sampler->sample();
    }
    else
    {
        #ifndef USE_REGULAR_QUEUE
        DataType next_point = std::move( *(main_queue.top()));
        delete main_queue.top();
        #else
        DataType next_point = std::move( *(main_queue.front()));
        delete main_queue.front();
        #endif // USE_REGULAR_QUEUE
        main_queue.pop(); //just removes the pointer (of which we have a copy)
        return next_point;
    }
}

//template<class ET>
//typename GaussQueue<ET,true>::RetType GaussQueue<ET,true>::true_pop()
//{
//    mutex_guard lock(queue_mutex); //global lock. TODO : Enable concurrent sampling.
//    if(main_queue.empty())
//    {
//        ++ (gauss_sieve->number_of_points_sampled); //atomic
//        ++ (gauss_sieve->number_of_points_constructed); //atomic
//        return gauss_sieve->sampler->sample();
//    }
//    else
//    {
//        #ifndef USE_REGULAR_QUEUE
//        LPType next_point = *(main_queue.top());
//        delete main_queue.top();
//        #else
//        LPType next_point = * (main_queue.front());
//        delete main_queue.front();
//        #endif // USE_REGULAR_QUEUE
//        main_queue.pop();
//        return next_point;
//    }
//}


//template<class ET>
//typename GaussQueue<ET,false>::RetType* GaussQueue<ET,false>::pop_take_ownership()
//{
//    if(main_queue.empty())
//    {
//    ++ (gauss_sieve->number_of_points_sampled);
//    ++ (gauss_sieve->number_of_points_constructed);
//    LPType *next_point_ptr = new LPType (gauss_sieve->sampler->sample());
//    return next_point_ptr;
//    }
//    else
//    {
//        #ifndef USE_REGULAR_QUEUE
//        LPType* next_point_ptr = main_queue.top();
//        #else
//        LPType* next_point_ptr = main_queue.front();
//        #endif // USE_REGULAR_QUEUE
//        main_queue.pop(); //remove pointer from queue.
//        return next_point_ptr;
//    }
//}
//
//template<class ET>
//typename GaussQueue<ET,true>::RetType* GaussQueue<ET,true>::pop_take_ownership()
//{
//    assert(false); //currently disabled
//}
//
//template<class ET>
//void GaussQueue<ET,false>::push(LPType const & val)
//{
//    LPType * new_lp = new LPType (val);
//    main_queue.push(new_lp);
//}
//
//template<class ET>
//void GaussQueue<ET,true>::push(LPType const & val)
//{
//    mutex_guard lock(queue_mutex);
//    LPType * new_lp = new LPType (val);
//    main_queue.push(new_lp);
//}


template<class ET, int nfixed>
void GaussQueue<ET,false,nfixed>::push(DataType && val)
{
    DataType * new_lp = new DataType (std::move(val) );
    main_queue.push(new_lp);
}

/*
Commented out, we avoid pointers at this interface level.
(or name the function differently)

template<class ET, int nfixed>
void GaussQueue<ET,false,nfixed>::push(DataType * &val)
{
    main_queue.push(val);
    val = nullptr;
}

*/

//template<class ET>
//void GaussQueue<ET,true>::push(LPType && val)
//{
//    mutex_guard lock(queue_mutex);
//    LPType * new_lp = new LPType (std::move(val) );
//    main_queue.push(new_lp);
//}


//template<class ET>
//void GaussQueue<ET,false>::give_ownership(LPType * const valptr)
//{
//    main_queue.push(valptr);
//}


//template<class ET>
//void GaussQueue<ET,true>::give_ownership(LPType * const valptr)
//{
//assert(false);
//}

template<class ET, int nfixed> GaussQueue<ET,false,nfixed>::~GaussQueue()
{
//TODO: Delete sampler if owned.
    while(! main_queue.empty() )
    {
        #ifndef USE_REGULAR_QUEUE
        delete main_queue.top();
        #else
        delete main_queue.front();
        #endif // USE_REGULAR_QUEUE
        main_queue.pop();
    }
    delete sampler;
}

//template<class ET> //making a second bool template argument does not work. You can not partially specialize member functions. (Workaround is possible, but its syntax is ridiculous).
//GaussQueue<ET,true>::~GaussQueue()
//{
////TODO: Delete sampler if owned.
//    while(! main_queue.empty() )
//    {
//        #ifndef USE_REGULAR_QUEUE
//        delete main_queue.top();
//        #else
//        delete main_queue.front();
//        #endif // USE_REGULAR_QUEUE
//        main_queue.pop();
//    }
//}

template class GaussQueue<Z_NR<long>,false,-1>;

#endif
