#ifndef GAUSS_SIEVE_BITAPPROX_H
#define GAUSS_SIEVE_BITAPPROX_H

#include <random>
#include "GlobalStaticData.h"
#include "DefaultIncludes.h"
#include <bitset>

#include <boost/dynamic_bitset.hpp> // maybe remove at some point.

namespace GaussSieve::SimHash
{

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
unsigned int constexpr sim_hash2_len = 64;
unsigned int constexpr sim_hash_number_of_coos = 4; // probably unused.

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
  static inline void print(std::ostream &os = std::cout)
  {
    os << "relevant matrix is: " << std::endl;
    for (uint_fast16_t i=0; i<sim_hash_len; i++)
    {
      for (uint_fast16_t j=0; j< sim_hash_number_of_coos; ++j)
      {
        os << rel_coo[i][j] << ", ";
      }
      os << std::endl;
    }
  }
  static inline uint_fast16_t get_ij_value(uint_fast16_t i, uint_fast16_t j)
  {
    return rel_coo[i][j];
  }

  private:
  static std::array<uint_fast16_t, sim_hash_number_of_coos> rel_coo[sim_hash_len];
};

template<class SieveTraits,bool MT>
std::array<uint_fast16_t,sim_hash_number_of_coos> CoordinateSelection<SieveTraits,MT>::rel_coo[sim_hash_len] = {};

} // end namespace GaussSieve::BitApprox


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
    }
    else
    {
      std::mt19937 rng;
      rng.seed(std::random_device()()); // TODO: Make this deterministic?
      std::uniform_int_distribution<std::mt19937::result_type> distr(0,ambient_dimension-1);

      for (uint_fast16_t i=0; i< SimHash::sim_hash_len; ++i)
      {
        for(uint_fast16_t j=0; i<SimHash::sim_hash_number_of_coos; ++j ) Data::rel_coo[i][j] = distr(rng);
      }
      Data::print();
    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing CoordinateSelection; Counter is " << Parent::user_count )
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing CoordinateSelection; Counter is " << Parent::user_count )
  }
};
}


namespace GaussSieve
{

//template<int nfixed>
class RelevantCoordinates;

uint_fast16_t constexpr  sim_hash_len = 64;
uint_fast16_t constexpr  sim_hash2_len = 64;

uint_fast16_t constexpr  num_of_coord = 4;

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
        std::cout << rel_coo[i][j] << ", ";

      std::cout << std::endl;
      //std::cout << rel_coo[i][0] << "," << rel_coo[i][1] << "," << rel_coo[i][2] << "," << rel_coo[i][3] <<std::endl;
    }
  }

  //member
  private:
  static std::array<uint_fast16_t, num_of_coord> rel_coo[sim_hash_len];
};
std::array<uint_fast16_t,num_of_coord> RelevantCoordinates::rel_coo[sim_hash_len] = {};


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

} // end namespace GaussSieve


// helpers for the function below, do not use.
namespace GaussSieve::SimHash::Helpers{
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
} // end namespace GaussSieve::SimHash::Helpers

namespace GaussSieve::SimHash
{

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
}// end namespace GaussSieve::SimHash

namespace GaussSieve::SimHash
{

//assume sim_hash2_len is a power-of-two.

// deprecated. Use versions below (faster, more flexible)
template<class ET>
[[deprecated]] inline std::array<ET,sim_hash2_len> fast_walsh_hadamard(std::array<ET,sim_hash2_len> const &input)
{
  std::array<ET,sim_hash2_len> inp = input;
  std::array<ET,sim_hash2_len> out;
  std::array<ET,sim_hash2_len> tmp;
  uint_fast16_t i, j, s;

  for (i = sim_hash2_len>>1; i > 0; i>>=1)
  {
    for (j = 0; j < sim_hash2_len; j++)
    {
      s = (j/i)%2;
      out[j]=inp[(s?-i:0)+j]+(s?-1:1)*inp[(s?0:i)+j];
    }
    tmp = inp; inp = out; out = tmp;
  }

  return out;
}


/**
 fast_partial_walsh_hadamard<len>(input) performs (fast) Walsh-Hadamard Transform on the first
 len coordinates of its input. len is enforced to be a power of two and input must be either
 a std::vector or a std::array. The entries of the vector / array must support arithmentic (+,-) and
 be default-constructible and swappable.
*/
template< unsigned int len, class T>
inline std::vector<T>  fast_partial_walsh_hadamard(std::vector<T> input) // Note: Pass by value is intentional. We modify the local copy.
{
  static_assert(is_a_power_of_two(len), "len must be a power of two");
  assert(len >= input.size() ); //maybe static
  std::vector<T> output(input.size() );
  for (uint_fast16_t i = len >> 1; i> 0; i>>=1 )
  {
    for(uint_fast16_t j = 0; j < len ; j++)
    {
      output[j]= ((j/i)%2!=0) ? input[j-i] - input[j] : input[j] + input[i+j];
    }
    std::swap(input,output);
  }
  return output;
}

template< unsigned int len, class T, std::size_t arraylen>
inline std::array<T,arraylen>  fast_partial_walsh_hadamard(std::array<T,arraylen> input) // Note: Pass by value is intentional. We modify the local copy.
{
  static_assert(is_a_power_of_two(len), "len must be a power of two");
  static_assert(arraylen >= len, "Cannot perfor WH-Transform on array of smaller size");
  std::array<T,arraylen> output(); // assumes that entries are default-constructible.
  for (uint_fast16_t i = len >> 1; i> 0; i>>=1 )
  {
    for(uint_fast16_t j = 0; j < len ; j++)
    {
      output[j]= ((j/i)%2!=0) ? input[j-i] - input[j] : input[j] + input[i+j];
    }
    std::swap(input,output);
  }
  return output;
}


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
    std::array<ET, sim_hash2_len> hadamard = fast_walsh_hadamard<ET>(input_vector);
    for(uint_fast16_t i=0;i<sim_hash_len;++i)
    {
      ret[i] = (hadamard[i]>=0) ? 1: 0;
    }

    return ret;
  }

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



}


#endif // include guards
