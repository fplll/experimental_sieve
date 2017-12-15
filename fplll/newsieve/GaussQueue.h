#ifndef GAUSS_QUEUE_H
#define GAUSS_QUEUE_H
/**
 This file declares the classes used for the main Queues in the Gauss Sieve
*/

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"
#include "Sampler.h"
#include "SieveUtility.h"
#include "Typedefs.h"

/**
  This file defines a class for the queue that stores the yet-unprocessed points of the Sieve.
  It essentially wraps around a std::(priority?) queue, with the added feature that
  whenever the queue is empty, we call the sampler to get a fresh point.

  Note that for compatibility with multi-threading, we provide a function true_pop that obtains the
  next element and removes it from the queue in a single operation (for std::queue, these are
  separate operations, but for multi-threaded data structures, it is better to make it one)

  TODO: Currently, we use an (unsorted) queue if USE_REGULAR_QUEUE is set, a (sorted) priority queue
        otherwise (in the single-threaded case). The unsorted queue is faster overall.
  TODO: Make the decision a template argument?

  NOTE: The Queue's constructor takes a pointer to the caller sieve, in order to access some of its
        data (notably, it needs to forward that pointer to the sampler, which might be
        user-provided, so we do not know which data it needs)
        For unit testing, there is the debug symbol DEBUG_SIEVE_STANDALONE_QUEUE that allows
        (forces, actually) to use nullptr for that pointer. It deactivates some features, naturally.
  TODO: Enable user-provided samplers
        (The infrastructure for this feature exists, we just need to add parameter to constructors
        and ways to set the sampler)
*/

/**
  Note: Due to circular dependencies between the Queue and the main sieve object, most functionality
        is implemented in GaussQueue_impl.h
        ( Since everything is header-only, there is some uglyness involved here:
          Notably, GaussQueue_impl must be included *only* at the root of the inclusion tree,
          *after* all non-impl_h files. )
  TODO: Revise / merge the above DOC into a global structure readme.
*/

namespace GaussSieve
{

// forward-declarations:
template <class SieveTraits, bool MT> class Sieve;
template <class SieveTraits, bool MT, class Engine, class Sseq> class Sampler;

// declares GaussQueue class, to be specialized.
template <class SieveTraits, bool MT> class GaussQueue;

template <class SieveTraits>  // single-threaded version:
class GaussQueue<SieveTraits, false>
{
private:  // aliases to avoid typing long names
  using DataType                    = typename SieveTraits::GaussQueue_DataType;
  using RetType                     = typename SieveTraits::GaussQueue_ReturnType;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;
#ifndef USE_REGULAR_QUEUE
  // If this symbol is not set, we use a priority queue to always process the smallest element from
  // the queue.
  // However, there is the issue that std::priortiy_queue does not allow non-const access to its
  // elements (presumably, because modifying would destroy the sorting)
  // This is true even if we erase that element immediately afterwards.

  // In particular, we can not use move semantics to get the next element, which is bad because
  // lattice points are not copyable. (This is a known issue with std::priority_queue)
  // As a solution, we store pointers to lattice points in the priority queue (which we can copy)
  // (Alternatively, use const_cast; while it is extremely likely to work, this depends on the
  // implementation details of priority_queue and might be UB. Also, we would need to overwrite
  // norm2 in the lattice to its old value, because the pop() operation (i.e. delete the top
  // element) of priority_queue might perform comparisons with the to-be-deleted element to update
  // its data structure, depending on the data structure it keeps)

  // Comparator is a local nested class that encapsulate the function by which to order.
  // Note that std::priority_queue defaults to < and allows access to the largest element,
  // so we have to reverse the ordering.
  struct Comparator
  {
    bool operator()(DataType *const &lhs, DataType *const &rhs) { return *lhs > *rhs; }
  };
  using QueueType = std::priority_queue<DataType *, std::vector<DataType *>, Comparator>;
#else
  using QueueType = std::queue<DataType>;
#endif

public:
  // We only ever want 1 queue. No copying / moving / assigning.
  // clang-format off
  GaussQueue()                              = delete;
  GaussQueue(GaussQueue const &)            = delete;
  GaussQueue(GaussQueue &&)                 = delete;
  GaussQueue &operator=(GaussQueue const &) = delete;
  GaussQueue &operator=(GaussQueue &&)      = delete;
  // clang-format on

  // Sole constructor:
  // we take take a pointer to the caller sieve as an argument to initialize the sampler and to
  // access statistics. static_data is needed to initialize the used lattice point classes.

  // user_sampler is a pointer to a user-provided samplers. If == nullptr, we create our own
  // seed_sampler is only used if we create our own sampler.
  // TODO: Consider storing a reference to the statistics class instead.
  // clang-format off
  explicit inline GaussQueue(Sieve<SieveTraits, false> *const caller_sieve,
                             GlobalStaticDataInitializer const &static_data,
                             int seed_sampler,
                             Sampler<SieveTraits, false, std::mt19937_64, std::seed_seq> *user_sampler = nullptr);
  // clang-format on
  inline ~GaussQueue();

  // we might as well always return false (or make this private)!
  // if the internal queue is empty, we just sample a new vector.
  // clang-format off
  [[deprecated("The queue is never empty from the users point of view.")]]
  NODISCARD inline bool empty() const
  // clang-format on
  {
    return main_queue.empty();
  }

  // returns size of queue (used for diagnostics and statistics only)
  typename QueueType::size_type size() const { return main_queue.size(); }

  // since we cannot / do not copy point, users have to use (possibly explicit) move semantics
  // (i.e. push(std::move(point)), which makes point unusable for the caller.
  // Note that there is a potential conversion at the call site.
  // TODO: Allow to push other lattice point types and explicitly convert.
  // (this one only sees implicit conversions)
  // clang-format off
  inline void push(DataType const &val) = delete;  // puts a copy of val in the queue : deleted
  inline void push(DataType      &&val);           // uses move semantics for that.
  // clang-format on

  // removes front element from queue *and returns it*.
  inline auto true_pop() -> RetType;

private:
  // clang-format off
  StaticInitializer<DataType> const init_data_type;
  StaticInitializer<RetType>  const init_ret_type;
  QueueType main_queue;                          // actual queue of lattice points to be processed.
  Sieve<SieveTraits, false> *const gauss_sieve;  // pointer to caller object.
  // clang-format on

  // NOTE: the sampler is public, because some modules of the sieve need to communicate with the
  // sampler (e.g. when using progressive rank).
  // TODO: Make private, add friends.
public:
  // a pointer to the Sampler. Note that we store a pointer, because we want to allow
  // user-provided samplers (which we do not own). Note that the Sampler class has
  // virtual members and *sampler might be a type derived from it.
  Sampler<SieveTraits, false, std::mt19937_64, std::seed_seq> *sampler;
  bool sampler_owned;  // Is the above pointer owning. Required to correctly delete it.
                       // (because we should not do that for user-provided ones).
                       // TODO: Use smart pointers instead.
};

}  // end namespace GaussSieve
#endif  // GAUSS_QUEUE_H
