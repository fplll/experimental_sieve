#ifndef GAUSS_QUEUE_CPP
#define GAUSS_QUEUE_CPP

#include "GaussQueue.h"
//#include "SieveGauss.cpp"
#include "ShiSampler.h"
#include "Sampler.h"
#include "DebugAll.h"

namespace GaussSieve
{

//constructor
template<class SieveTraits>
GaussQueue<SieveTraits,false>::GaussQueue( Sieve<SieveTraits,false> * const caller_sieve,
  GlobalStaticDataInitializer const &static_data)
:
init_data_type(static_data),
init_ret_type(static_data),
main_queue(),
gauss_sieve(caller_sieve),
sampler(nullptr)
{
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
  assert(caller_sieve==nullptr);
#else
  assert(caller_sieve!=nullptr);
#endif

  std::seed_seq seed{1,2,4}; //just some arbitrary numbers
    //sampler = new EllipticSampler<ET,false, std::mt19937_64, std::seed_seq> (seed);
    std::cout << "Initializing Sampler" << std::endl << std::flush;
  sampler = new ShiSampler<SieveTraits,false, std::mt19937_64, std::seed_seq> (seed);
//    std::cout << "Finished Initializing Sampler" << std::endl << std::flush;
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

template<class SieveTraits>
auto GaussQueue<SieveTraits,false>::true_pop() -> RetType
{
  if(main_queue.empty()) // Queue is empty, sample a new element.
  {
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
    assert(gauss_sieve==nullptr);
#else
    assert(gauss_sieve!=nullptr);
    ++ (gauss_sieve->number_of_points_sampled);
    ++ (gauss_sieve->number_of_points_constructed);
#endif
//        return gauss_sieve->sampler->sample();
    assert(sampler!=nullptr);
    
    //I'VE COMMENTED OUT THIS STATIC ASSERT: we convert from what sampler returns to what queue stores here -- EK
    
    //static_assert(std::is_same<typename SieveTraits::GaussSampler_ReturnType, RetType>::value,
    //  "Sampler must currently return the same type as the queue.");
//    typename SieveTraits::GaussSampler_ReturnType const ret = sampler->sample();
//    return ret;
//        return sampler->sample();
    typename SieveTraits::GaussList_StoredPoint ret = static_cast<typename SieveTraits::GaussList_StoredPoint>(sampler->sample());
    return ret;
    //return static_cast<typename SieveTraits::GaussList_StoredPoint>(sampler->sample());
  }
  else // Queue is not empty, just return stored element.
  {
#ifndef USE_REGULAR_QUEUE
    DataType next_point = std::move( *(main_queue.top()));
    delete main_queue.top();
#else
    DataType next_point = std::move( *(main_queue.front()));
    delete main_queue.front();
#endif // USE_REGULAR_QUEUE
    main_queue.pop(); //This just removes the pointer.
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


template<class SieveTraits>
void GaussQueue<SieveTraits,false>::push(DataType && val)
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

template<class SieveTraits> GaussQueue<SieveTraits,false>::~GaussQueue()
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

//template class GaussQueue<fplll::Z_NR<long>,false,-1>;

}

#endif
