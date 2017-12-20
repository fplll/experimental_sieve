#ifndef ARRAY_LIST_H
#define ARRAY_LIST_H

#include "DefaultIncludes.h"

namespace GaussSieve
{

namespace ArrayListDetails
{

template<class T, unsigned int blocksize>
struct ArrayListBlock
{
  std::array<T,blocksize> block;
  unsigned int usedsize;

  ArrayListBlock() : block(), usedsize(0) {}
  ArrayListBlock(T const &value): block(), usedsize(1)
  {
    block[0] = value;
  }
  ArrayListBlock(T &&value): block(), usedsize(1)
  {
    block[0] = std::move(value);
  }

  ArrayListBlock(ArrayListBlock const &) = default;
  ArrayListBlock(ArrayListBlock &&) = default;
};

}  // end namespace ArrayListDetails

template<class T, unsigned int blocksize> class ArrayList;
template<class T, unsigned int blocksize> class ArrayListConstIterator;

template<class T, unsigned int blocksize>
class ArrayList
{
  static_assert(blocksize > 0, "");
  static_assert(std::is_default_constructible<T>::value, "" );

private:
  using Block = ArrayListDetails::ArrayListBlock<T,blocksize>;
  std::list<Block> blocks;
  std::size_t total_size;
public:
  using const_iterator = ArrayListConstIterator<T,blocksize>;
  using size_type = std::size_t;
  using value_type = T;
  friend const_iterator;

  ArrayList() : blocks(), total_size(0) {}
  ArrayList(ArrayList const &) = default;
  ArrayList(ArrayList &&) = default;

  size_type size() const noexcept { return total_size; }
  NODISCARD bool empty() const noexcept { return total_size==0; }
  //void clear() noexcept { blocks.clear(); total_size = 0; }
  const_iterator cbegin() const noexcept { return const_iterator{ blocks.cbegin(),0 }; }
  const_iterator cend() const noexcept { return const_iterator{blocks.cend(),0}; }
  void push_back(T const &value)
  {
    if(total_size==0)  // special case
    {
      blocks.emplace_back(T);
    }
    else
    {
      // TODO: Can improve?
      auto it = blocks.cend();
      --it;
      if( it.index == it.current_block_size)
      {

      }
    }

  }

};

template<class T, unsigned int blocksize>
class ArrayListConstIterator
{
  friend ArrayList<T, blocksize>;
private:
  using Block = ArrayListDetails::ArrayListBlock<T,blocksize>;
  using Array = std::array<T,blocksize>;
  using ArrayIt = typename std::array<T,blocksize>::iterator;
  using ListIt  = typename std::list<Block>::iterator;

  ListIt list_it;
  unsigned int index;
  unsigned int current_block_size;

  void update_max_index() { current_block_size = list_it->used_size; }

  // Default constructors and destructors
public:

  bool operator==(ArrayListConstIterator const &other) const { return (list_it==other.list_it) && (index==other.index); }
  bool operator!=(ArrayListConstIterator const &other) const { return (list_it!=other.list_it) || (index!=other.index); }

  ArrayListConstIterator(ListIt const &new_list_it, unsigned int const &new_index)
      : list_it(new_list_it), index(new_index), current_block_size(new_list_it->used_size) {}

  ArrayListConstIterator &operator++()  // prefix only for now.
  {
    ++index;
    if(index == current_block_size)
    {
      ++list_it;
      index = 0;
      update_max_index();
    }
    return *this;
  }

  T const & operator*() const { return list_it->block[index]; }
  T const * operator->() const { return list_it->block.data()+index; }

};



}  // end namespace


#endif
