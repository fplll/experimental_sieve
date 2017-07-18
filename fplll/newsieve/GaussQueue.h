#ifndef GAUSS_QUEUE_H
#define GAUSS_QUEUE_H
/* defines the classes used for the main Queues in the Gauss Sieve */

template <class ET, bool MT, int nfixed> class GaussQueue;

//template <class ET> class IsLongerVector_class; //class wrapper to compare vectors by length. Needs to be wrapped in a class to work seamlessly with some STL containers.
//template <class ET> class IsLongerVector_classPtr;


//template <class ET>
//class GaussQueue<ET,true>;
//
//template <class ET>
//class GaussQueue<ET,false>;
//
//template <class ET>
//using GaussQueueST=GaussQueue<ET,false>;
//
//template <class ET>
//using GaussQueueMT=GaussQueue<ET,true>;

#include <mutex>
#include <atomic>
#include <queue>
#include <utility>
#include <random>
//#include "LatticePoint.h"
#include "LatticePointsNew.h"
#include "assert.h"
#include "Typedefs.h"
//#include "EllipticSampler.h"


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
template<class ET,bool MT, int nfixed> class Sieve;
template<class ET,bool MT, class Engine, class Sseq, int nfixed> class Sampler;

/*
template<class ET,bool MT, int nfixed> class IsLongerVector_ExactPtr
{
    public: bool operator() (CompressedPoint<ET,MT,nfixed>* const &A, CompressedPoint<ET,MT,nfixed>* const & B)
    {
        return (A->get_exact_norm2() > B->get_exact_norm2() );
    }
};
*/

template<class ET,int nfixed> //single-threaded version:
class GaussQueue<ET,false,nfixed>
{
public:
    using DataType = typename GaussSieve::GaussQueue_DataType<ET,false,nfixed>;    //Type of Data internally stored
    using RetType=   typename GaussSieve::GaussQueue_ReturnType<ET,false,nfixed>;    //Type of Data returned
    #ifndef USE_REGULAR_QUEUE
    using QueueType =      std::priority_queue< DataType* , std::vector<DataType* >, IsLongerVector_ExactPtr<ET,false,nfixed> >;
    #else
    using QueueType =      std::queue<DataType*>;
    #endif
    using size_type = typename QueueType::size_type;
    //using SamplerType =    KleinSampler<typename ET::underlying_data_type, FP_NR<double> > ;
    GaussQueue()=delete;
    GaussQueue(Sieve<ET,false,nfixed> *caller_sieve); //only constructor
    GaussQueue(GaussQueue const &old) = delete;
    GaussQueue(GaussQueue &&old) = delete;
    GaussQueue& operator= (GaussQueue const &old)=delete;
    GaussQueue& operator= (GaussQueue &&old) = delete;
    ~GaussQueue();

    [[deprecated("The queue is never empty from the callers POV")]]
    bool empty() const              {return main_queue.empty();};  //we might as well always return false (or make this private)!
    size_type size() const          {return main_queue.size();};   //returns size of queue (used for diagnostics and statistics only)
    void push(DataType const &val) = delete; //puts a copy of val in the queue : deleted
    void push(DataType && val);     //uses move semantics for that.
    //void push(DataType * &val); //uses move semantics! val is changed to nullptr

//    [[deprecated("Ownership transfer clashes with compressed storage.")]]
//    void give_ownership(LPType * const valptr); //takes a pointer to a list point and puts the point into the queue, moves ownership (avoids copying)

    RetType true_pop(); //removes front element from queue *and returns it*.
    //RetType true_pop() = delete;

//    [[deprecated("Use copy elison rather than ownership transfer.")]]
//    RetType* pop_take_ownership() ; //removes front elements from queue and returns handle to it.
                                   //Transfers ownership to the caller. Return type might change, but should be dereferencable, deleteable.
                                   //might become deprecated

private:
    QueueType main_queue;           //actual queue of lattice points to be processed.
    Sieve<ET,false,nfixed>* gauss_sieve;   //pointer to caller object.
    //SamplerType *sampler; //controlled by the GaussSieve currently. TODO: Change that
public:
    Sampler<ET,false,std::mt19937_64, std::seed_seq,nfixed> * sampler; //or a type derived from it.
};





/*
template<class ET> //multi-threaded version:
class GaussQueue<ET,true>
{
public:
    using EntryType = ET; //entries of lattice points
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



#endif // GAUSS_QUEUE_H
