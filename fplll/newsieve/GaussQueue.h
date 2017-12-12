#ifndef GAUSS_QUEUE_H
#define GAUSS_QUEUE_H
/* defines the classes used for the main Queues in the Gauss Sieve */

#include "DefaultIncludes.h"
#include "SieveUtility.h"
//#include "LatticePoint.h"
//#include "LatticePointsNew.h"
#include "Typedefs.h"
//#include "EllipticSampler.h"
#include "Sampler.h"
#include "GlobalStaticData.h"

namespace GaussSieve{

// declares GaussQueue class, to be specialized.
template <class SieveTraits, bool MT> class GaussQueue;


/**
TODO: Documentation, enable sorted queue.
*/

//template<class ET,int nfixed> class IsLongerVector_class //TODO : Move to GaussQueue.h
//{
//    public: bool operator() (LatticePoint<ET> const &A, LatticePoint<ET> const & B)
//    {
//     return (A.get_norm2() > B.get_norm2() );
//    }
//};
//
//template<class ET> class IsLongerVector_classPtr
//{
//    public: bool operator() (LatticePoint<ET>* const &A, LatticePoint<ET>* const & B)
//    {
//        return (A->get_norm2() > B->get_norm2() );
//    }
//};

//forward-declarations:
template<class SieveTraits, bool MT> class Sieve; // TODO: Change To ST

template<class SieveTraits, bool MT, class Engine, class Sseq> class Sampler;

/*
template<class ET,bool MT, int nfixed> class IsLongerVector_ExactPtr
{
    public: bool operator() (CompressedPoint<ET,MT,nfixed>* const &A, CompressedPoint<ET,MT,nfixed>* const & B)
    {
        return (A->get_exact_norm2() > B->get_exact_norm2() );
    }
};
*/

template<class SieveTraits> //single-threaded version:
class GaussQueue<SieveTraits,false>
{
public:
    using DataType = typename SieveTraits::GaussQueue_DataType;    //Type of Data internally stored
    using RetType=   typename SieveTraits::GaussQueue_ReturnType;    //Type of Data returned
    using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;

    static_assert(std::is_same<DataType,RetType>::value, "Currently, DataType and RetType must be identical.");
#ifndef USE_REGULAR_QUEUE
    // TODO: Make this one work (and actually a template argument)
    static_assert(false, "Only regular queue might work at the moment");
    using QueueType = DOES NOT WORK : std::priority_queue< DataType* , std::vector<DataType* >, IsLongerVector_ExactPtr<ET,false,nfixed> >;
#else
    using QueueType = std::queue<DataType*>;
#endif
    GaussQueue()=delete;
    explicit inline GaussQueue(Sieve<SieveTraits,false>* const caller_sieve,
      GlobalStaticDataInitializer const &static_data); //only constructor
    GaussQueue(GaussQueue const &) = delete;
    GaussQueue(GaussQueue &&) = delete;
    GaussQueue& operator= (GaussQueue const &)=delete;
    GaussQueue& operator= (GaussQueue &&) = delete;
    inline ~GaussQueue();

    // we might as well always return false (or make this private)!
    // if the internal queue is empty, we just sample a new vector.
    [[deprecated("The queue is never empty from the callers POV.")]]
    inline bool empty() const  {return main_queue.empty();};

    //returns size of queue (used for diagnostics and statistics only)
    inline long long size() const {return main_queue.size();}; //TODO: MAY BE MORE THAN LONGLONG
    void push(DataType const &val) = delete; //puts a copy of val in the queue : deleted
    inline void push(DataType && val);     //uses move semantics for that.

    // allow pushing via pointer?

    inline auto true_pop() -> RetType; //removes front element from queue *and returns it*.

private:
    StaticInitializer<DataType> const init_data_type;
    StaticInitializer<RetType>  const init_ret_type;
    QueueType main_queue;           //actual queue of lattice points to be processed.
    Sieve<SieveTraits,false>* const gauss_sieve;   //pointer to caller object.
public:
    Sampler<SieveTraits,false,std::mt19937_64, std::seed_seq> * sampler; //or a type derived from it.
};


// clang-format off


/*
template<class ET> //multi-threaded version:
class GaussQueue<ET,true>
{
public:
    using LengthType = ET; //entries of lattice points
    using LPType = LatticePoint<ET>; //Type of Data internally stored
    using RetType= LatticePoint<ET>; //Type of Data returned
    using mutex_guard = std::lock_guard<std::mutex>;
    #ifndef USE_REGULAR_QUEUE
    using QueueType =      std::priority_queue< LPType* , std::vector<LPType* >, IsLongerVector_classPtr<ET> >;
    #else
    using QueueType =      std::queue<LPType*>;
    #endif
    using size_type = typename QueueType::size_type;
    //using SamplerType =    KleinSampler<typename ET::underlying_data_type, FP_NR<double> > ;

    GaussQueue()=delete;
    GaussQueue(Sieve<ET,true> *caller_sieve); //only constructor, not thread-safe
    GaussQueue(GaussQueue const &old) = delete;
    GaussQueue(GaussQueue &&old) = delete;
    GaussQueue& operator= (GaussQueue const &old)=delete;
    GaussQueue& operator= (GaussQueue &&old) = delete;
    ~GaussQueue(); //not thread-safe

    //TODO: Fix const - correctness

    bool empty()                      {mutex_guard lock(queue_mutex); return main_queue.empty();}; //checks whether the queue is currently empty. Blocks
    size_type size()                  {mutex_guard lock(queue_mutex); return main_queue.size();};  //returns size of queue (used for diagnostics and statistics only). Blocks
    void push(LPType const &val); //puts a copy of val in the queue
    void push(LPType && val);     //uses move semantics for that.
    [[deprecated("Ownership transfer clashes with compressed storage.")]]
    void give_ownership(LPType * const valptr); //takes a pointer to a list point and puts the point into the queue, moves ownership (avoids copying)
    RetType  true_pop(); //removes front element from queue *and returns it*.
    [[deprecated("Use copy elison rather than ownership transfer.")]]
    RetType* pop_take_ownership() ; //removes front elements from queue and returns handle to it.
                                   //Transfers ownership to the caller. Return type might change, but should be dereferencable, deleteable.
                                   //might become deprecated

private:
    QueueType main_queue;
    Sieve<ET,true>* gauss_sieve; //caller object.
    std::mutex queue_mutex; //global lock. We do not differentiate reads and writes.
    //SamplerType *sampler; //controlled by the GaussSieve currently. TODO: Change that
};

*/

}

#endif // GAUSS_QUEUE_H
//clang-format on
