#ifndef GAUSS_SIEVE_BITAPPROX_H
#define GAUSS_SIEVE_BITAPPROX_H

#include <random>
#include "GlobalStaticData.h"
#include "DefaultIncludes.h"
#include <bitset>
#include <math.h> //for log2

#include <boost/dynamic_bitset.hpp> // maybe remove at some point.

namespace GaussSieve{ namespace SimHash {

/*****************************************************************************
is_a_power_of_two(n) checks whether n is a power of 2.
Slow, but designed for use in static asserts
(This is why it's a recursive one-line function, to be C++11 - constexpr)
******************************************************************************/

template<class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
constexpr bool is_a_power_of_two(Integer const n)
{
  // one-line function to be C++11 - constexpr. Slow, but only used in static_asserts anyway.
  return (n>0) &&  ( (n==1) || ( (n%2==0) && is_a_power_of_two(n/2)));
}

/**********************************************************
Global constants, accessible as SimHash::sim_hash_len, etc.
Consider moving (some) of them to SieveTraits
**********************************************************/

unsigned int constexpr sim_hash_len = 64;
//unsigned int constexpr sim_hash2_len = 64;
//unsigned int constexpr sim_hash_number_of_coos = 4; // probably unused.
unsigned int constexpr num_of_levels = 3; //number of approximation level (bitapprox is of size sim_hash_len*num_of_levels
unsigned int constexpr num_of_transforms = 2; //per 1 approximation (to get dim-dimensional binary vector)
constexpr std::array<unsigned int, num_of_levels> threshold_lvls = {2, 2, 2};

// Default Parameters are JUST FOR TESTING. Change these!
template<class SieveTraits = void, bool MT = false> class CoordinateSelection; // forward declare

}}



namespace GaussSieve{ namespace SimHash {

/*
 stores (num_of_levels*num_of_transforms) permutation matrices of dim=ambient_dim
  must be instantiated once at the start of the sieve
 */
class PMatrix
{
  private:
  // vector is of size ambient_dimension. The individual permutation[i]'s are different.
  std::vector<uint_fast16_t> permutation;
  public:
  PMatrix() = default;
  PMatrix(unsigned int dim) //initialize with a random dim-dimensional permutation.
  {
    //TODO: TO CHECK!
//    PMatrix::total_num_of_matrices = static_cast<int>(( SimHash::sim_hash_len/ ambient_dim + 1)*
//                                                          SimHash::num_of_levels*SimHash::num_of_transforms );

  DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to initialize P matrix")
  permutation.resize(dim);
  for(uint_fast16_t i=0; i<dim; ++i) permutation[i] = i;
  std::random_shuffle(permutation.begin(), permutation.end());
#ifdef DEBUG_SIEVE_TRACEINITIATLIZATIONS
  print();
#endif
}
  // apply the stored permutation to the given vector. Changes the argument.
  template<class T> void apply(std::vector<T> &vec) const
  {
    auto const dim = vec.size();
    assert(dim == permutation.size());
    std::vector<T> out(dim);
    for(unsigned int i=0;i<dim;++i)
    {
      out[i] = vec[permutation[i]];
    }
    using std::swap;
    swap(out,vec);
    return;
  }

  inline void print(std::ostream & os = std::cout) const
  {
    // os << "P  [" << i <<"] is: " << std::endl;
    for (uint_fast16_t k=0; k<permutation.size(); ++k)
    {
      os << permutation[k] << " ";
    }
    os << std::endl;
  }
};

}} // end namespace GaussSieve::SimHash

namespace GaussSieve { namespace SimHash {

/**
  A diagonal matrix with +/-1 on the diagonal.
*/

class DMatrix
{
  //friend StaticInitializer<DMatrix>;
  private:
  std::vector<uint_fast8_t> diagonal;
  public:
  DMatrix() = default;
  DMatrix(unsigned int const dim)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to fill-up the D matrix")

    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_int_distribution<uint_fast8_t> distr(0, 1); //0 means no sign-flip, 1 means sign-flip.
    diagonal.resize(dim);
    for(uint_fast16_t i=0;i<dim;++i)
    {
      diagonal[i] = distr(rng);
    }
#ifdef DEBUG_SIEVE_TRACEINITIATLIZATIONS
    print();
#endif
  }
  template<class T> void apply(std::vector<T> &vec) const
  {
    auto const dim = vec.size();
    assert(dim == diagonal.size());
    for(unsigned int i = 0;i<dim;++i)
    {
      if(diagonal[i]!=0) vec[i] = -vec[i]; //otherwise, keep vec[i] as it is.
    }
  }

  inline void print(std::ostream & os = std::cout) const
  {
    for (uint_fast16_t k=0;k<diagonal.size(); ++k )
    {
      os << diagonal[k] << " ";
    }
    os << std::endl;
  }
};
}} // end namespace GaussSieve::SimHash

namespace GaussSieve{ namespace SimHash {

/********************************************************************
SimHash::CoordinateSelection<SieveTraits,IsMultithreaded>
is a class template holding the data related to bitapproximation.
This class only contains static data and static member functions.

Note: We may template this by the user lattice point class as well or
store these data inside some other class which uses the data.
Initializiation is performed by StaticInitializer<CoordinateSelection>.
This is called by ExactLatticePoint at the moment, but this is subject
to change.

Note2: Template arguments are subject to change.

TODO: Consider writing a wrapper around either a std::bitset or a std::dynamic_bitset
that actually has these as static data and initialize it via the main sieve.
********************************************************************/

template<class SieveTraits, bool MT>
class CoordinateSelection
{
  public:
  friend StaticInitializer<CoordinateSelection<SieveTraits,MT>>;

  CoordinateSelection() = delete; // cannot be instantiated.

  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  inline static auto transform_and_bitapprox(LatP const &point)
    -> std::array<std::bitset<sim_hash_len>,num_of_levels>;
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  inline static auto transform_and_bitapprox_2nd_layer(LatP const &point)
    -> std::array<std::bitset<sim_hash_len>,num_of_levels>;
  private:
  static unsigned int number_of_blocks;
  static unsigned int fast_walsh_hadamard_len;
  static std::vector<std::array<PMatrix,num_of_transforms> > pmatrices;
  static std::vector<std::array<DMatrix,num_of_transforms> > dmatrices;
};

template<class SieveTraits,bool MT>
unsigned int CoordinateSelection<SieveTraits,MT>::number_of_blocks;
template<class SieveTraits,bool MT>
std::vector<std::array<PMatrix,num_of_transforms> > CoordinateSelection<SieveTraits,MT>::pmatrices;
template<class SieveTraits,bool MT>
std::vector<std::array<DMatrix,num_of_transforms> > CoordinateSelection<SieveTraits,MT>::dmatrices;
template<class SieveTraits,bool MT>
unsigned int CoordinateSelection<SieveTraits,MT>::fast_walsh_hadamard_len;
}} // end namespace GaussSieve::BitApprox


namespace GaussSieve{

/**************************************************************************************
Initialization of CoordinateSelection is done here.

Note: This is done in namespace GaussSieve rather than GaussSieve::SimHash, because
      the StaticInitializer template that is specialized is in the GaussSieve namespace.
***************************************************************************************/

template<class SieveTraits,bool MT>
class StaticInitializer<SimHash::CoordinateSelection<SieveTraits,MT> >
: public DefaultStaticInitializer< SimHash::CoordinateSelection<SieveTraits,MT> >
{
  using Parent = DefaultStaticInitializer< SimHash::CoordinateSelection<SieveTraits,MT> >;
  using Data = SimHash::CoordinateSelection<SieveTraits,MT>;
  public:

  template<class T,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const & initializer) : StaticInitializer(initializer.dim) {}

  StaticInitializer(int ambient_dimension)
  {
    assert(Parent::user_count > 0);
    if(Parent::user_count>1)
    {
      // check that the value is what is should be (see below)
      // TODO: Store dimension and verify it was the same as in the last invocation?
      assert(Data::number_of_blocks == (SimHash::num_of_levels * SimHash::sim_hash_len -1) / ambient_dimension + 1);
    }
    else
    {
      Data::number_of_blocks =  (SimHash::num_of_levels * SimHash::sim_hash_len -1) / ambient_dimension + 1;
      // num_of_levels * sim_hash_len is the total number of required bits in the output simhashes.
      // We need ceil(total number of required bits / ambient_dimension) blocks (i.e. repetitions of the input
      // vector) to get at least that many output bits.

      Data::fast_walsh_hadamard_len = static_cast<unsigned int>( pow(2, floor(log2(ambient_dimension)) ) );

      Data::pmatrices.resize(Data::number_of_blocks);
      Data::dmatrices.resize(Data::number_of_blocks);
      for(unsigned int i=0;i<Data::number_of_blocks;++i)
      {
        for(unsigned int j=0; j< SimHash::num_of_transforms; ++j)
        {
          Data::pmatrices[i][j] = static_cast<SimHash::PMatrix>(ambient_dimension);
          Data::dmatrices[i][j] = static_cast<SimHash::DMatrix>(ambient_dimension);
        }
      }
      // TODO
      //Data::print();
    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing CoordinateSelection; Counter is " << Parent::user_count )
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing CoordinateSelection; Counter is " << Parent::user_count )
  }
};
} // end namespace GaussSieve


namespace GaussSieve{ namespace SimHash{

/**
 fast_partial_walsh_hadamard<len>(input) performs (fast) Walsh-Hadamard Transform on the first
 len coordinates of its input. len is enforced to be a power of two and input must be either
 a std::vector or a std::array. The entries of the vector / array must support arithmentic (+,-) and
 be default-constructible and swappable.
*/
template<class T>
inline auto  fast_partial_walsh_hadamard(std::vector<T> input, unsigned int len) // Note: Pass by value is intentional. We modify the local copy.
-> std::vector<T>
{
  assert( std::bitset<std::numeric_limits<unsigned long>::digits>{len}.count() == 1 ); // len must be a power of two.
  //static_assert(is_a_power_of_two(len), "len must be a power of two");

  assert(len <= input.size() ); //maybe static
  std::vector<T> output(input.size() );

  for (uint_fast16_t i =len; i<output.size(); ++i)
  {
    output[i] = input[i];
  }

  //on [0...len-1] coordinates perform WH

  for (uint_fast16_t i = len >> 1; i> 0; i>>=1 )
  {
    for(uint_fast16_t j = 0; j < len ; j++)
    {
      // Note: The static_cast<T> is required to be able to use mpz_class, because mpz_class internally
      // used expression templates (so the types of a+b and a-b don't match for mpz_classes a,b)
      output[j]= ((j/i)%2!=0) ? static_cast<T>(input[j-i] - input[j]) : static_cast<T>(input[j] + input[i+j]);
    }
    std::swap(input,output);
  }
  return output;
}

template<class T, std::size_t arraylen>
inline auto fast_partial_walsh_hadamard(std::array<T,arraylen> input, unsigned int len) // Note: Pass by value is intentional. We modify the local copy.
-> std::array<T,arraylen>
{
  assert( std::bitset<std::numeric_limits<unsigned long>::digits>{len}.count() == 1 ); // len must be a power of two.
  assert( len <= arraylen);
  std::array<T,arraylen> output(); // assumes that entries are default-constructible.
  for (uint_fast16_t i = len >> 1; i> 0; i>>=1 )
  {
    for(uint_fast16_t j = 0; j < len ; j++)
    {
      // Note: The static_cast<T> is required to be able to use mpz_class, because mpz_class internally
      // used expression templates (so the types of a+b and a-b don't match for mpz_classes a,b)
      output[j]= ((j/i)%2!=0) ? static_cast<T>(input[j-i] - input[j]) : static_cast<T>(input[j] + input[i+j]);
    }
    std::swap(input,output);
  }
  return output;
}


template<class SieveTraits, bool MT>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto CoordinateSelection<SieveTraits,MT>::transform_and_bitapprox(LatP const &point) //, uint_fast16_t iteration)
-> std::array<std::bitset<sim_hash_len>,num_of_levels>
{
  unsigned int dim = static_cast<unsigned int>(point.get_dim());
  using ET = Get_AbsoluteCooType<LatP>;
  using Block  = std::vector<ET>; // will be of length dim.
  Block copy_of_point(dim);
  for(unsigned int i = 0;i<dim;++i)
  {
    copy_of_point[i] = point.get_absolute_coo(i);
  }
  std::vector<Block> blocks(number_of_blocks,std::move(copy_of_point)); //holds number_of_blocks many copies of the input point.

  // each block gets transformed by applying WH * D * P num_transform many times.
  for(uint_fast16_t i=0; i<number_of_blocks;++i ) // for each block:
  {
    for(uint_fast8_t j=0; j<num_of_transforms;++j) // repeat num_of_transforms many times:
    {
      pmatrices[i][j].apply(blocks[i]);
      dmatrices[i][j].apply(blocks[i]);
      blocks[i] = fast_partial_walsh_hadamard(blocks[i],fast_walsh_hadamard_len);
    }
  }
  // put together the blocks into an array of bitsets.
  std::array<std::bitset<sim_hash_len>,num_of_levels> ret;
  for(unsigned int n=0;n<num_of_levels;++n)
  {
    for(unsigned int m=0;m<sim_hash_len;++m)
    {
      // index of the bit currently considered
      // if we use only one level of indexing.
      unsigned int const flat_bit_count = n*sim_hash_len + m;

      ret[n][m] = (blocks[flat_bit_count / dim][flat_bit_count % dim] > 0);
    }
  }
  return ret;
}

// TODO: Merge this function with the one above to avoid recomputations.

template<class SieveTraits, bool MT>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto CoordinateSelection<SieveTraits,MT>::transform_and_bitapprox_2nd_layer(LatP const &point) //, uint_fast16_t iteration)
-> std::array<std::bitset<sim_hash_len>,num_of_levels>
{
  using std::abs;
  using std::max;
  unsigned int dim = static_cast<unsigned int>(point.get_dim());
  using ET = Get_AbsoluteCooType<LatP>;
  using Block  = std::vector<ET>; // will be of length dim.
  Block copy_of_point(dim);
  for(unsigned int i = 0;i<dim;++i)
  {
    copy_of_point[i] = point.get_absolute_coo(i);
  }
  std::vector<Block> blocks(number_of_blocks,std::move(copy_of_point)); //holds number_of_blocks many copies of the input point.

  // each block gets transformed by applying WH * D * P num_transform many times.
  for(uint_fast16_t i=0; i<number_of_blocks;++i ) // for each block:
  {
    for(uint_fast8_t j=0; j<num_of_transforms;++j) // repeat num_of_transforms many times:
    {
      pmatrices[i][j].apply(blocks[i]);
      dmatrices[i][j].apply(blocks[i]);
      blocks[i] = fast_partial_walsh_hadamard(blocks[i],fast_walsh_hadamard_len);
    }
  }
  // put together the blocks into an array of bitsets.
  std::array<std::bitset<sim_hash_len>,num_of_levels> ret;
  for(unsigned int n=0;n<num_of_levels;++n)
  {
  // maximal entry of the outbut block.
  // Note: The walsh-hadamard transform used above is not normalized, so the elements are larger
  // than the original vector. That's why we do not compare to norm2() currently.
  // TODO : Reconsider the point where we change the sign. Currently, it's at
  // half the maximal entry. This may not be ideal.
    ET maxentry = 0;
    for(unsigned int m=0;m<sim_hash_len;++m)
    {
      unsigned int const flat_bit_count = n*sim_hash_len + m;
      maxentry = max(maxentry,abs(blocks[flat_bit_count/dim][flat_bit_count % dim]));
    }

    for(unsigned int m=0;m<sim_hash_len;++m)
    {
      // index of the bit currently considered
      // if we use only one level of indexing.
      unsigned int const flat_bit_count = n*sim_hash_len + m;

      ret[n][m] = (abs(blocks[flat_bit_count / dim][flat_bit_count % dim])*2 > maxentry);
    }
  }
  return ret;
}


/**
  This class stores the result of computing a scalar product of bitwise
  approximations.
  Essentially, just wraps around an int.
  Note: We might want to wrap an approximate scalar product of t as value = #bits - t. -- Gotti
*/

//#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX
class BitApproxScalarProduct
{

  public:
  using BitApproxScalarProduct_WrappedType       = uint_fast32_t;


  BitApproxScalarProduct(BitApproxScalarProduct const &old) = delete; // why not copy ints?
  BitApproxScalarProduct(BitApproxScalarProduct &&old)      = default;

  explicit constexpr BitApproxScalarProduct(BitApproxScalarProduct_WrappedType const rhs):value(rhs) {}
  explicit operator BitApproxScalarProduct_WrappedType() { return value; }

  BitApproxScalarProduct &operator=(BitApproxScalarProduct const &other) = delete; // Why?
  BitApproxScalarProduct &operator=(BitApproxScalarProduct &&other) = default;

  //TODO: operator >=, <=

  inline bool operator<=(BitApproxScalarProduct_WrappedType && rhs)
  {
    return  this->value <= rhs;
  }

  /*
  friend std::ostream& operator<<(std::ostream &os, BitApproxScalarProduct const &value)
  {
    os<< value;
    return os;
  }
   */
  //member
  BitApproxScalarProduct_WrappedType value;
};
//#endif



}} //end namespace GaussSieve::SimHash



// OLD UNUSED CODE


/*
class RelevantCoordinates;

uint_fast16_t constexpr  sim_hash_len = 64;
uint_fast16_t constexpr  sim_hash2_len = 64;

uint_fast16_t constexpr  num_of_coord = 4;
*/

/*
//TODO: REPLACE rand() by a proper rand
//template<int nfixed>
class RelevantCoordinates
{
  friend StaticInitializer<RelevantCoordinates>;

  public:
  RelevantCoordinates() = delete;

  RelevantCoordinates(RelevantCoordinates const &) = delete;
  RelevantCoordinates(RelevantCoordinates &&obj)   = delete;


  RelevantCoordinates  &operator=(RelevantCoordinates const &obj) = delete;
  RelevantCoordinates  &operator=(RelevantCoordinates &obj)       = delete;


  // for 0<=i<sim_hash_len; 0<=j<num_of_coord
  static inline uint_fast16_t get_ij_value(uint_fast16_t i, uint_fast16_t j)
  {
    return (rel_coo[i][j]);
  }

  public:
  static inline void print()
  {
    std::cout << "relevant matrix is: " << std::endl;
    for (uint_fast16_t i=0; i<sim_hash_len; i++)
    {
      for (uint_fast16_t j=0; j<num_of_coord; ++j)
      {
        std::cout << rel_coo[i][j] << ", ";
      }
      std::cout << std::endl;
      //std::cout << rel_coo[i][0] << "," << rel_coo[i][1] << "," << rel_coo[i][2] << "," << rel_coo[i][3] <<std::endl;
    }
  }

  //member
  private:
  static std::array<uint_fast16_t, num_of_coord> rel_coo[sim_hash_len];
};
std::array<uint_fast16_t,num_of_coord> RelevantCoordinates::rel_coo[sim_hash_len] = {};
*/


/*
// Static Initializer:
// template<int nxfixed>
template<>
class StaticInitializer<class RelevantCoordinates>
: public DefaultStaticInitializer<RelevantCoordinates>
{
  using Parent = DefaultStaticInitializer<RelevantCoordinates>;
  public:

  //template<class T,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  //StaticInitializer(T const & initializer) : StaticInitializer(initializer.dim) {} //<-WHAT IS IT FOR?

  StaticInitializer(int ambient_dimension)
  {
    assert(Parent::user_count > 0);
    if(Parent::user_count>1)
    {
      std::cout << "user_cout for  RelevantCoordinates > 1" << std::endl;
    }
    else
    {

      std::mt19937 rng;
      rng.seed(std::random_device()());
      std::uniform_int_distribution<std::mt19937::result_type> distr(0,ambient_dimension-1);

      for (uint_fast16_t i=0; i<sim_hash_len; ++i)
      {
        RelevantCoordinates::rel_coo[i][0] = distr(rng);
        RelevantCoordinates::rel_coo[i][1] = distr(rng);
        RelevantCoordinates::rel_coo[i][2] = distr(rng);
        RelevantCoordinates::rel_coo[i][3] = distr(rng);
      }
      RelevantCoordinates::print();
    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing RelevantCoordinates; Counter is " << Parent::user_count )
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing RelevantCoordinates; Counter is " << Parent::user_count )
  }
};
*/


/*
 stores (num_of_levels*num_of_transforms) diagonal matrices of dim=ambient_dim
  with {-1, 0, 1} elements on the main diagonal
  must be instantiated once at the start of the sieve
 */

/*
template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
inline std::bitset<sim_hash_len> compute_fixed_bitapproximation(LatP const &point)
{
  //using RelevantCoords = GaussSieve::RelevantCoordinates;
  using ET = Get_CoordinateType<LatP>;

  //std::cout << "rel_coo_matrix used: " <<std::endl;
  //RelevantCoordinates::print();
  //assert(false);

  std::bitset<sim_hash_len> ret;

  uint_fast16_t bound = mystd::constexpr_min(static_cast<uint_fast16_t>(point.get_dim()), static_cast<uint_fast16_t>(sim_hash_len));

  for(uint_fast16_t i=0;i<bound;++i)
  {
    ret[i] = (point.get_absolute_coo(i)>=0) ? 1 : 0;
  }

  //for(uint_fast16_t i=0;i<sim_hash_len;++i)
  for(uint_fast16_t i=point.get_dim();i<sim_hash_len;++i)
  {
  // Gotti: 3x + and 1x - seems strange to me.
    ET res =  point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,0)) +
              point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,1)) +
              point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,2)) -
              point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,3));
      //std::cout << "i = " << i << " res =" << res << " ";
      ret[i] = (res>=0) ? 1: 0;
    }

    //assert(false);
    return ret;
}
*/

/*
template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
static inline std::bitset<sim_hash2_len> compute_2order_fixed_bitapproximation(LatP const &point)
{
  std::bitset<sim_hash2_len> ret;

  using ET = Get_CoordinateType<LatP>;
  std::array<ET, sim_hash2_len> input_vector;

  //assume sim_hash_len=sim_hash2_len
  for(uint_fast16_t i=0;i<sim_hash_len;++i)
  {
    input_vector[i] = point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,0)) -
                      point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,1)) +
                      point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,2)) -
                      point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,3));
  }
  assert(false); // TODO: This looks wrong. Please check! - -Gotti
    // (applying WH-Trafo *after* selecting relevant coos seems an error)

//    std::array<ET, sim_hash2_len> hadamard = fast_walsh_hadamard<ET>(input_vector);
//    for(uint_fast16_t i=0;i<sim_hash_len;++i)
//    {
//      ret[i] = (hadamard[i]>=0) ? 1: 0;
//    }
//
//    return ret;
//
}
*/

/*
template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
static inline boost::dynamic_bitset<> compute_2nd_order_bitapproximation(LatP const &point)
{
  using ET = Get_CoordinateType<LatP>;
  using std::abs;
  ET max_coo = 0;
  auto dim = point.get_dim();
  boost::dynamic_bitset<> ret{static_cast<size_t>(dim) };
  //find the max fist
  for(uint_fast16_t i=0;i<dim;++i)
  {
    // equivalent, but has the problem of not working well for all coo types (in particular, mpz_class)
//         max_coo = max_coo ^ ((max_coo ^ abs(point.get_absolute_coo(i))) & -(max_coo < abs(point.get_absolute_coo(i)) )); //<-works for positive coeffs

      // Note : The static_cast is needed to deactivate lazy evaluation inside mpz_class.
    using std::max;
    max_coo = max(max_coo, static_cast<ET>(abs(point.get_absolute_coo(i))));
  }
    //std::cout << "max_coo = " << max_coo << std::endl;
    //compute the 2nd-order approximation
  for(uint_fast16_t i=0;i<dim;++i)
  {
    ret[i] = (abs(2*point.get_absolute_coo(i))  >= max_coo) ? 1 : 0;
  }
  return ret;
}
*/

// helpers for the function below, do not use.
namespace GaussSieve{ namespace SimHash { namespace Helpers{
template<int SizeOfBitSet> struct MakeBitApprox_Helper
{
  static_assert(SizeOfBitSet>=0, "Only for fixed-size bit-sets.");
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline std::bitset<SizeOfBitSet> compute_bitapproximation(LatP const &point)
  {
    auto const dim = point.get_dim();
    static_assert(dim == SizeOfBitSet, "Only usable if size of bitset equals fixed ambient dim of vector");
    std::bitset<SizeOfBitSet> ret;
    for(uint_fast16_t i=0; i<dim; ++i )
    {
      ret[i] = (point.get_absolute_coo(i)>=0) ? 1 :0;
    }
    return ret;
  }
};

// variable-dim version:
template<> struct MakeBitApprox_Helper<-1>
{
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline boost::dynamic_bitset<> compute_bitapproximation(LatP const &point)
  {
    auto const dim = point.get_dim();
    boost::dynamic_bitset<> ret{ static_cast<size_t>(dim) };
//    ret.resize(dim);
    for(uint_fast16_t i=0;i<dim;++i)
    {
      ret[i] = (point.get_absolute_coo(i)>=0) ? 1 : 0;
    }
    return ret;
  }
};
}}} // end namespace GaussSieve::SimHash::Helpers

namespace GaussSieve{ namespace SimHash{

/**
  compute_coordinate_wise_bitapproximation<length>(point)
  computes the coordinate wise 1-bit-approximation of the lattice point point.
  It returns a std::bitset<length> for lenght>= 0 or boost::dynamic_bitset for length == -1.
  Currently, the length parameter (if >=0) has to match the (fixed!) dimension of the lattice point.
  For length==1, the dynamic_bitset that is returned has length equal to the dimension of point.
*/
template<int SizeOfBitSet,class LatP>
auto compute_coordinate_wise_bitapproximation(LatP const &point)
-> decltype(Helpers::MakeBitApprox_Helper<SizeOfBitSet>::compute_bitapproximation(std::declval<LatP>() ))
{
  static_assert(IsALatticePoint<mystd::decay_t<LatP>>::value,"Not a lattice point.");
  return Helpers::MakeBitApprox_Helper<SizeOfBitSet>::compute_bitapproximation(point);
}
}}// end namespace GaussSieve::SimHash




namespace GaussSieve{ namespace SimHash{

//assume sim_hash2_len is a power-of-two.

// deprecated. Use versions below (faster, more flexible)
// not for now
template<class T>
inline std::vector<T> fast_walsh_hadamard_ext(std::vector<T> input, unsigned int len)
{

  std::vector<T> output (input.size());
  std::vector<T> tmp = input;
  uint_fast16_t i, j, s;

  for (i = len>>1; i > 0; i>>=1)
  {
    for (j = 0; j < len; j++)
    {
      s = (j/i)%2;
      output[j]=input[(s?-i:0)+j]+(s?-1:1)*input[(s?0:i)+j];
    }
    //tmp = inp; inp = out; out = tmp;
    std::swap(input, output);
  }

  //lower matrix with 1's on the main diag
  for (i =len; i<output.size(); ++i)
  {
    output[i] = tmp[i];
  }

  return output;
}

}} // end namespace GaussSieve::SimHash


#endif // include guards
