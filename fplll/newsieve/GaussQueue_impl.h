#ifndef GAUSS_QUEUE_CPP
#define GAUSS_QUEUE_CPP

#include "GaussQueue.h"
//#include "SieveGauss.cpp"
#include "GPVSampler.h"
#include "UniformSampler.h"
#include "GPVSamplerExtended.h"
#include "Sampler.h"
#include "DebugAll.h"

namespace GaussSieve
{

// constructor, single-threaded version
template<class SieveTraits>
GaussQueue<SieveTraits,false>::GaussQueue(Sieve<SieveTraits,false> * const caller_sieve,
                                          GlobalStaticDataInitializer const &static_data,
                                          Sampler<SieveTraits,false,std::mt19937_64,std::seed_seq> *user_sampler)
    : init_data_type(static_data),
      init_ret_type(static_data),
      main_queue(),
      gauss_sieve(caller_sieve),
      sampler(user_sampler),
      sampler_owned(user_sampler == nullptr)
{
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
  assert(caller_sieve==nullptr);
#else
  assert(caller_sieve!=nullptr);
#endif
  if (sampler == nullptr)
  {
    std::seed_seq seed{1,2,4}; //just some arbitrary numbers
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing our own sampler. Using GPVSampler")
    sampler = new GPVSampler<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed);
  }
  assert(sampler != nullptr);
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
    gauss_sieve->statistics.increment_number_of_points_sampled();
    gauss_sieve->statistics.increment_number_of_points_constructed();
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




template<class SieveTraits>
void GaussQueue<SieveTraits,false>::push(DataType && val)
{
  DataType * new_lp = new DataType (std::move(val) );
  main_queue.push(new_lp);
}



template<class SieveTraits> GaussQueue<SieveTraits,false>::~GaussQueue()
{
  while(! main_queue.empty() )
  {
#ifndef USE_REGULAR_QUEUE
  delete main_queue.top();
#else
  delete main_queue.front();
#endif // USE_REGULAR_QUEUE
  main_queue.pop();
  }
  if (sampler_owned)
  {
    delete sampler;
  }
}

}  // end namespace GaussSieve

#endif  // include guard
