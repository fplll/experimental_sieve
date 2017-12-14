// clang-format off

#ifndef POINT_LIST_NEW_H
#define POINT_LIST_NEW_H

#include "DefaultIncludes.h"


#include "Typedefs.h"
#include <list>
#include "SieveUtility.h"

//Class for (weakly?) sorted list of lattice points.
//includes thread-safe variant(s). May need experiments which implementation is best. (mutex on whole structure on every write, lock-free,...)
//Note that the application is rather lenient on requirement:
// -- keeping the list sorted is not strictly required and can be worked around. It should be _somewhat_ sorted for efficiency reasons.
//Reading and actually using vectors that are marked-for-deletion does not hurt too much.
//The structure of the algorithm allows for relatively simple garbage collection.
//Note that we do NOT assume a sequentially constistent memory model here. Relaxing from that should gain a bit of performance:
//Be aware that reading time from the lattice point list even asymptotically leading order.

namespace GaussSieve{

//forward declarations
template <class SieveTraits, bool MT> class GaussListNew;
template <class SieveTraits, bool MT> class GaussIteratorNew;


template <class SieveTraits>
class GaussListNew<SieveTraits, false>
{
public:
  friend GaussIteratorNew<SieveTraits,false>;
  using StoredPoint = typename SieveTraits::GaussList_StoredPoint;
  using ReturnType  = typename SieveTraits::GaussList_ReturnType;

  using UnderlyingContainer = std::list<StoredPoint>;
  using Iterator = GaussIteratorNew<SieveTraits,false>;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;

  explicit GaussListNew(GlobalStaticDataInitializer const &static_data)
      : init_stored_point(static_data), init_return_type(static_data), actual_list() {}
  GaussListNew(GaussListNew const & old) = delete;
  GaussListNew(GaussListNew && old) = delete;
  GaussListNew & operator= (GaussListNew const &other) = delete;
  GaussListNew & operator= (GaussListNew &&other) = delete;
  ~GaussListNew() = default; //FIXME: CHECK!!!

  Iterator cbegin() { return static_cast<Iterator>(actual_list.begin()); };
  Iterator cend()   { return static_cast<Iterator>(actual_list.end()); };

    //These functions insert (possibly a copy of val) into the list.
    //TODO: include ownership transfer semantics to avoid some copying, possibly include refcounts in LatticePoints.
    //This becomes really tricky if overwrite is allowed in multithreaded case.

    //no automatic copying: You can use insert_before(point.make_copy() );
    //The issue is that copying should be explicit.
  Iterator insert_before(Iterator pos, StoredPoint const & val) = delete;

  template<class LatticePoint, TEMPL_RESTRICT_DECL(IsALatticePoint<LatticePoint>::value)>
  Iterator insert_before(Iterator pos, LatticePoint && val)
  { return static_cast<Iterator> (actual_list.emplace(pos.it, static_cast<StoredPoint>( std::move(val)))); }

    /*Iterator insert_before_give_ownership(Iterator pos, DetailType * const val) = delete;  //TODO
    Iterator insert_before(Iterator pos, DataType const & val) =  delete; //TODO
    Iterator insert_before_give_ownership(Iterator pos, DataType const & val) = delete;                          //TODO
    */
    //Iterator insert_before(Iterator pos, DetailType && val)           {return actual_list.insert(pos,std::move(val));};


  Iterator erase(Iterator pos) {return static_cast<Iterator> (actual_list.erase(pos.it));} //only for single-threaded
    //void unlink(Iterator pos)=delete;     //MT only                              //{actual_list.erase(pos);};

    //TODO: Use aux_data, sort by calling comparison function
  void sort() {actual_list.sort();}  //only for single-threaded (for now).

private:
  StaticInitializer<StoredPoint> const init_stored_point;
  StaticInitializer<ReturnType>  const init_return_type;
  UnderlyingContainer actual_list;
};

//behaves like a const_iterator to approximate points.

template <class SieveTraits>
class GaussIteratorNew<SieveTraits,false>
{
    friend GaussListNew<SieveTraits,false>;
    public:

    /*
    What this class wraps around:
    */
    using UnderlyingIterator = typename GaussListNew<SieveTraits,false>::UnderlyingContainer::iterator; //non-const
    using CUnderlyingIterator= typename GaussListNew<SieveTraits,false>::UnderlyingContainer::const_iterator; //const-version

    /*
    This is what dereferencing will give us (modulo cv - qualifiers)
    */
    using DerefType  = typename SieveTraits::GaussList_ReturnType; //without cv - spec.

    GaussIteratorNew () = delete; // ???
    GaussIteratorNew (GaussIteratorNew const & other) = default;
    GaussIteratorNew (GaussIteratorNew && other) = default;
    GaussIteratorNew & operator= (GaussIteratorNew const & other) = default;
    GaussIteratorNew & operator= (GaussIteratorNew && other) = default;
    ~GaussIteratorNew() = default;

    //ExactLatticePoint<ET,nfixed> const & dereference_exactly_r() {return it -> access_exact_point_r();};
    //LatticePoint<ET> * access_details() { assert(  (*this)->get_details_ptr_rw()!=nullptr );  return (*this)-> get_details_ptr_rw() ;};
    explicit GaussIteratorNew(UnderlyingIterator const & other) : it(other) {}; //make iterator from pointer

    GaussIteratorNew&  operator++() {++it; return *this;}; //prefix version
    GaussIteratorNew  operator++(int){return it++;}; //postfix version

    [[deprecated("Decrement operator is deprecated") ]] GaussIteratorNew&  operator--() {--it; return *this;}; //For tests only
    bool operator==( GaussIteratorNew const & other) const {return it==(other.it);};
    bool operator!=( GaussIteratorNew const & other) const {return it!=(other.it);};

    bool is_end() const = delete; //not implemented for single-threaded yet.
    //intrinsic check for end, validity etc?
    //operator DerefType const *();
    DerefType const &  operator*() const    {return *it;};
    DerefType const *  operator->() const   {return &(*it);};
    //CUnderlyingIterator const operator->() const {return static_cast<CUnderlyingIterator>(it);};

    //Note : access_details may become deprecated at some point, if details only store difference to approximation or coefficients wrt other vectors.
    //(i.e. we might store the details in a compressed fashion to reduce memory complexity)

    //DetailType * access_details() {return it->get_details_ptr_rw();}; //Note: In multithreaded environment, there is no write access.
    //get_exact_point() const {return *(it->get_details_ptr());}; //retrieves a copy of the exact value. preferably use this one.
    //ET get_true_norm2() const  {return it->get_exact_norm2();}; // Use this function to get the true norm2 value.
    //DerefType & deref_rw() {return *it;};
    private:
    UnderlyingIterator it;
    //may require handle to container object.
};


//Multithreaded:

/*

template <class ET>
class GaussList<ET, true, -1>
{
public:
    //friend GaussIterator<ET,true,-1>;
    using LengthType= ET;
    using DataType = ApproxLatticePoint<ET,false,-1>; //Type of data stored in the list (from caller's POV) //TODO: distinguish MT case in data type.
    using DataPointer=DataType *;
    using AtomicDataPointer = std::atomic<DataType *>;
    using Node=ListMTNode<DataType>; //will change
    using NodePointer = Node *;
    using AtomicNodePointer = std::atomic<ListMTNode<DataType> *>;
    using UnderlyingContainer = std::list<DataType>;
    //using Iterator = GaussIterator<ET,true,-1>;
    using Iterator=GaussIterator<ET,true,-1>;
    using DetailType = typename DataType::DetailType;
    using ExactType = LatticePoint<ET>;
    explicit GaussList(); //creates empty list. Details in cpp.
    GaussList(GaussList const & old) = delete; //no moving / copying due to mutex.
    GaussList(GaussList && old) = delete;
    GaussList & operator= (GaussList const &other) = delete;
    GaussList & operator= (GaussList &&other) = delete;
    //TODO: Create from single-threaded list.
    ~GaussList();   //Deletes the underlying container. Details in cpp. Not thread-safe


    Iterator cbegin()                                                   {return start_sentinel_node->next_node.load(std::memory_order_acquire);};
    Iterator cend()                                                     {return end_sentinel_node;};

    Iterator insert_before(Iterator pos, ExactType const & val); //TODO, type of val might change. Uses emplacement.

    //inserts a copy of val before pos. The iterator pos remains valid.
    //NOTE: If *pos is marked-for-deletion, val will be inserted before
    //the next non-marked position pos' that is reachable by increasing pos.
    //Even in this case, incrementing pos will eventually reach pos' *WITHOUT* seeing the newly inserted value.

    Iterator erase(Iterator pos) =delete; //single-threaded only.
    void unlink(Iterator const &pos, GarbageBin<DataType> &gb); //marks the element pointed at by pos as deleted. Such elements will (eventually)
                                                          //become unreachable by traversing the list. Any Iterators to *pos (including pos)
                                                          //remain valid. *pos is put onto gb, whose job is eventually freeing memory.
                                                          //TODO : Include mechanism to forwarding arguments to gb.

    void sort();



    //TODO: Constructor from SingleThreaded variant and vice versa?

private:
    std::mutex mutex_currently_writing;    //global lock for writing threads.
    NodePointer const start_sentinel_node; //node before the start of the list. This is never modified, so no atomic here.
    NodePointer const end_sentinel_node;   //node after the end of the list, could probably do with a single sentinel.
};

template<class ET>
class GaussIterator<ET,true,-1>
{
public:
    friend GaussList<ET,true,-1>;
    friend void swap(GaussIterator &A, GaussIterator &B)              {std::swap(A.p,B.p);};
    using DataType    = typename GaussList<ET,true,-1>::DataType;
    using Node=ListMTNode< DataType >; //data representation
    using NodePointer = Node *;
    using DataPointer = DataType*;
    using CDataPointer= DataType const *;
    using DerefType  = DataType; //without cv - spec.
    using DetailType = LatticePoint<ET>;
    using ExactType  = LatticePoint<ET>; //need not be the same!

    GaussIterator() = delete; //should always init with valid object.
    GaussIterator(NodePointer const & _p)                                   : p(_p)         {assert(_p!=nullptr);};
    GaussIterator(GaussIterator<ET,true,-1> const &old)                 : p(old.p)      {};
    GaussIterator(GaussIterator<ET,true,-1> && old)=default;
    GaussIterator& operator=(GaussIterator const &other) = default;
    GaussIterator& operator=(GaussIterator &&other) = default;
    ~GaussIterator() {};

    GaussIterator& operator++();    //prefix version
    GaussIterator  operator++(int); //postfix version

    GaussIterator& operator--();    ////For tests. Not needed in MT. Tests are so far performed on ST. To be deleted

    DataPointer const operator->() const                                {return p->datum_ptr;}; //Note weird semantics of -> overload cause datum_ptr (which is a pointer) to get dereferenced as well.
    DerefType const & operator*()   const                               {return *(p->datum_ptr);};
    //operator DataPointer() const                                        {return p->datum_ptr;}; //converts from Iterator to pointer to lattice point.
    bool operator==(GaussIterator<ET,true,-1> const & other) const      {return p==other.p;};  //test for equality of iterators.
    bool operator!=(GaussIterator<ET,true,-1> const & other) const      {return p!=other.p;};  //test for inequality of iterators. Note: Even if it1 != it2, it may still be the case that *it1 == *it2.
    bool is_end() const                                                 {return p->check_for_end_node();};
    bool is_good() const                                                {return !(p->is_marked_for_deletion());}; //Not: Note reliable if caller did not lock mutex. May be made private.

//  DetailType*  access_details() {return it->get_details_ptr_rw();}; //Note: In multithreaded environment, there is no write access.
    ExactType get_exact_point() const {return p->datum_ptr->get_details();}; //retrieves a copy of the exact value. preferably use this one.
    ET get_true_norm2() const {return p->datum_ptr->get_details_ptr_r()->get_norm2();}; // Use this function to get the true norm2 value.


private:
    NodePointer p; //does not own. Need not be atomic. Is NOT a const - pointer!
};

*/

/*

template <class DT> //need to redo. //TODO: Thou shalt not inherit from STL containers (no virtual destructors).
class GarbageBin : public std::queue< ListMTNode<DT> * >
{
public:
    GarbageBin() = default;
    GarbageBin(GarbageBin const &old) = delete;
    GarbageBin(GarbageBin &&old) = default;
    GarbageBin& operator=(GarbageBin const &old) = delete;
    GarbageBin& operator=(GarbageBin &&old) = default;
    ~GarbageBin()   {empty_trash();}
    void empty_trash();
};

*/

/*


template <class DT>
class ListMultiThreaded
{
//friend PointListMTNode<ET>;
public:
    using Node=ListMTNode<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    using AtomicDataPointer = std::atomic<DT *>;
    using NodePointer = ListMTNode<DT> *;
    using AtomicNodePointer = std::atomic<ListMTNode<DT> *>;
    using Iterator=MTListIterator<DT>;
    explicit ListMultiThreaded(); //called when only one thread is running. Creates empty list.
    ListMultiThreaded(ListMultiThreaded const &old)=delete; //No copying via constructor! (due to mutex)
    ListMultiThreaded(ListMultiThreaded && old)=delete; //No moving (mutex)
    ListMultiThreaded & operator=(ListMultiThreaded const & old) = delete; //No copy assignment. (due to mutex)
    ListMultiThreaded & operator=(ListMultiThreaded && old)=delete; //dito
    //TODO: Constructor from SingleThreaded variant and vice versa?
    ~ListMultiThreaded();
    Iterator cbegin() const; //returns iterator to first element. On empty list, returns valid past-the-end iterator.
    Iterator cend() const;   //returns past-the-end iterator. Must not be dereferenced. (same behaviour as stl::list)
    //Iterator cbefore_begin() const; //returns iterator before the first element. May not be dereferenced. Included for completeness.
    void unlink(Iterator const &pos, GarbageBin<DT> &gb); //marks the element pointed at by pos as deleted. Such elements will (eventually)
                                                          //become unreachable by traversing the list. Any Iterators to *pos (including pos)
                                                          //remain valid. *pos is put onto gb, whose job is eventually freeing memory.
                                                          //TODO : Allow forwarding additional arguments to gb.
        //  For each of the following insertion routines, the
        //  return value is an iterator to the newly inserted element.
    Iterator insert_before(Iterator const &pos, DT const &val);     //inserts a copy of val before pos. The iterator pos remains valid.
                                                                    //NOTE: If *pos is marked-for-deletion, val will be inserted before
                                                                    //the next non-marked position pos' that is reachable by increasing pos.
                                                                    //Even in this case, incrementing pos will eventually reach pos' without seeing the newly inserted value.
    Iterator enlist_before(Iterator const &pos, DT * const valref); //moves *valref just before pos (as above), transfering ownership to the list. Avoids copying of *valref.
//same as above, but inserts after the corresponding element. Return values is an iterator to the newly inserted element.
//In case *pos is marked-for-deletion, inserts after the next non-marked position pos' that is reachable by increasing pos.
//Due to these reasons, it is not guaranteed that pos+1 == retval.
//Furthermore, insert_before(pos+1,val) may differ from insert_after(pos,val).
    Iterator insert_after(Iterator  const &pos, DT const &val) = delete;     //not implemented yet.
    Iterator enlist_after(Iterator  const &pos, DT * const valref) = delete; //not implemented yet.

private:
    std::mutex mutex_currently_writing;
    NodePointer const start_sentinel_node; //node before the start of the list. This is never modified, so no atomic here.
    NodePointer const end_sentinel_node; //node after the end of the list, could probably do with a single sentinel.
};

*/

//restriction: Iterator itself is thread-local.


//template<class DT>
//class MTListIterator
//{
//public:
//    friend GaussList<typename DT::LengthType,true,-1>;
//    friend ListMultiThreaded<DT>;
//    friend void swap(MTListIterator &A, MTListIterator &B)      {std::swap(A.p,B.p);};
//    using Node=ListMTNode<DT>;
//    using DataType    = DT;
//    using DataPointer = DT*;
//    using CDataPointer= DT const *;
//    MTListIterator() = delete; //should always init with valid object.
//    MTListIterator(Node * const & _p): p(_p)                    {assert(_p!=nullptr);};
//    MTListIterator(MTListIterator<DT> const &old): p(old.p)     {};
//    MTListIterator(MTListIterator<DT> && old)=default;
//    ~MTListIterator() {};
//    MTListIterator<DT>& operator=(MTListIterator<DT> const &other) = default;
//    MTListIterator<DT>& operator=(MTListIterator<DT> &&other) = default;
//    MTListIterator<DT>& operator++(); //prefix version
//    MTListIterator<DT>  operator++(int); //postfix version
//    DataPointer operator->() const                              {return p->datum_ptr;}; //Note weird semantics of -> overload cause datum_ptr (which is a pointer) to get dereferenced as well.
//    DataType & operator*()   const                              {return *(p->datum_ptr);};
//    operator DataPointer() const                                {return p->datum_ptr;}; //converts from Iterator to pointer to lattice point.
//    bool operator==(MTListIterator<DT> const & other) const     {return p==other.p;};  //test for equality of iterators.
//    bool operator!=(MTListIterator<DT> const & other) const     {return p!=other.p;};  //test for inequality of iterators. Note: Even if it1 != it2, it may still be the case that *it1 == *it2.
//    bool is_end() const                                         {return p->check_for_end_node();};
//    bool is_good() const                                        {return !(p->is_marked_for_deletion() ); //Not: Note reliable if caller did not lock mutex. May be made private.
//    };
//private:
//    Node * p; //does not own. Need not be atomic. Is NOT a const - pointer!
//};
////TODO: define iterator_traits specialisation (or derive from std::iterator)
//

/*

//custom linked list node. For internal use only. Owns an object of type DT
template <class DT>
class ListMTNode
{
//    friend GaussList<typename DT::LengthType,true,-1>;
//    friend ListMultiThreaded<DT>;
//    friend MTListIterator<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    using AtomicDataPointer = std::atomic<DT *>;
    using NodePointer = ListMTNode<DT> *;
    using AtomicNodePointer = std::atomic<ListMTNode<DT> *>;
public:
    ListMTNode() : next_node(nullptr),prev_node(nullptr), datum_ptr(nullptr),nodestatus(0)  {} ;
    ListMTNode(ListMTNode const &old) = delete;
    ListMTNode(ListMTNode &&old) = delete;
    ListMTNode & operator=(ListMTNode const &old) =delete;
    ListMTNode & operator=(ListMTNode &&old) = delete;
    ~ListMTNode()                                               {delete datum_ptr;}; //destructor
    bool check_for_end_node() const                             {return nodestatus == static_cast<int>(NodeStatus::is_last_node);};
    bool is_marked_for_deletion() const                         {return nodestatus == static_cast<int>(NodeStatus::is_to_be_deleted);};
    bool is_sentinel_node() const                               {return (nodestatus == static_cast<int>(NodeStatus::is_first_node))||(nodestatus == static_cast<int>(NodeStatus::is_last_node) );};
    bool is_plain_node() const                                  {return nodestatus == 0;};
//private: --not intended for outside use, anyway.
public:
    AtomicNodePointer next_node;
    NodePointer       prev_node;
    DataPointer       datum_ptr; //actual data. We use a pointer here for potential atomicity. This is hidden from the user.
    int nodestatus; //We use int rather than an enum-type to be on the safe side regarding that atomic-operations be possible.
public:
    enum class NodeStatus
    {
        is_to_be_deleted=1,
        is_first_node=2,
        is_last_node=4
    }; //meaning of values for nodestatus. Comparison requires typecast.
};

*/

//template <class ET>
//using PointListMultiThreaded= ListMultiThreaded<LatticePoint<ET>>;
//template <class ET> class PointListIterator;


//--------------------------------- OLD CPP FILE

//continuation from PointList.h
//no include guards or anything needed.

/*

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

*/

/*

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

*/


//iterators:

/*

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

*/

/*

template<class DT>
void GarbageBin<DT>::empty_trash()
{
    while(!this->empty())
    {
        delete this->front();
        this->pop();
    }
};

*/

}

#endif

//clang-format on
