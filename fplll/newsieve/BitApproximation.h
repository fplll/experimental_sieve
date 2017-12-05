#ifndef GAUSS_SIEVE_BITAPPROX_H
#define GAUSS_SIEVE_BITAPPROX_H

#include <bitset>
#include <boost/dynamic_bitset.hpp> // maybe remove at some point.
#include <random>
#include <math.h> //for log2, sqrt

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"

// clang-format off

namespace GaussSieve
{
namespace SimHash
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
  return (n > 0) && ( (n == 1) || ( (n%2 == 0) && is_a_power_of_two(n/2) ) );
}

/**********************************************************
Global constants, accessible as SimHash::sim_hash_len, etc.
Consider moving (some) of them to SieveTraits
**********************************************************/

unsigned int constexpr sim_hash_len = 64;
// unsigned int constexpr sim_hash2_len = 64;
// unsigned int constexpr sim_hash_number_of_coos = 4; // probably unused.
unsigned int constexpr num_of_levels = 2;  // number of approximation level
                                           // (bitapprox is of size sim_hash_len*num_of_levels)
// per level, we perform a random orthogonal transformation.
// For improved speed, we approximate a uniformly random matrix by applying
// WH * D * P num_of_transforms many times. (WH = Walsh-Hadamard, D = diagonal, P = permutation)
unsigned int constexpr num_of_transforms = 2;
constexpr std::array<unsigned int, num_of_levels> threshold_lvls_2sieve = {{4,8}};
constexpr std::array<unsigned int, num_of_levels> threshold_lvls_3sieve = {{0,0}};

//constexpr int num_of_coos = 4;

// Default Parameters are JUST FOR TESTING. Change these!

template<class SieveTraits = void, bool MT = false> class CoordinateSelection;  // forward declare

/**
 PMatrix creates and stores a random permutation (matrix) in dim dimensions.
 Can be apply()'d to a vector.
 */
class PMatrix
{
private:
  std::vector<uint_fast16_t> permutation;
public:
  PMatrix() = default;
  PMatrix(unsigned int dim)  // initialize with a random dim-dimensional permutation.
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to initialize P matrix")
    permutation.resize(dim);
    for (uint_fast16_t i = 0; i < dim; ++i)
    {
      permutation[i] = i;
    }
    std::random_shuffle(permutation.begin(), permutation.end());
#ifdef DEBUG_SIEVE_TRACEINITIATLIZATIONS
    print();
#endif
  }
  // apply the stored permutation to the given vector. Changes the argument.
  template<class T> inline void apply(std::vector<T> &vec) const
  {
    auto const dim = vec.size();
    assert(dim == permutation.size());
    std::vector<T> out(dim);
    for (unsigned int i = 0; i < dim; ++i)
    {
      out[i] = vec[permutation[i]];
    }
    using std::swap;
    swap(out, vec);
    return;
  }

  inline void print(std::ostream &os = std::cout) const
  {
    // os << "P  [" << i <<"] is: " << std::endl;
    for (uint_fast16_t k = 0; k < permutation.size(); ++k)
    {
      os << permutation[k] << " ";
    }
    os << std::endl;
  }
};

/**
  Creates and stores a random diagonal matrix with +/-1 on the diagonal.
  Can be apply()'ed to a vector.
*/

class DMatrix
{
private:
  std::vector<uint_fast8_t> diagonal;  // 0 means no sign-flip, 1 means sign-flip:
public:
  DMatrix() = default;
  DMatrix(unsigned int const dim)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to fill-up the D matrix")

    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_int_distribution<uint_fast8_t> distr(0, 1);
    diagonal.resize(dim);
    for (uint_fast16_t i = 0; i < dim; ++i)
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
    for (unsigned int i = 0; i < dim; ++i)
    {
      if (diagonal[i] != 0)  // otherwise, keep vec[i] as it is.
      {
        vec[i] = -vec[i];
      }
    }
  }

  inline void print(std::ostream &os = std::cout) const
  {
    for (uint_fast16_t k = 0; k < diagonal.size(); ++k)
    {
      os << static_cast<int>(diagonal[k]) << " ";
    }
    os << std::endl;
  }
};

/************************************
 * OLD RELEVANT COORDITATES (FOR COMPARISON)
 ************************************/
/*
class RMatrix
{
  private:
  std::array<uint_fast16_t,SimHash::num_of_coos> coos [SimHash::sim_hash_len];
  public:
  RMatrix() = default;
  RMatrix(unsigned int const dim)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to fill-up the R matrix")

    std::mt19937 rng;
    rng.seed(std::random_device()());
    std::uniform_int_distribution<std::mt19937::result_type> distr(0,dim-1);
    for (uint_fast16_t i=0; i<SimHash::sim_hash_len; ++i)
    {
      for (uint_fast16_t j=0; j<SimHash::num_of_coos; ++j)
        coos[i][j] = distr(rng);
    }

  }

  uint_fast16_t get_value(uint_fast16_t i, uint_fast16_t j) const
  {
       return coos[i][j];
  }

};
*/
/********************************************************************
SimHash::CoordinateSelection<SieveTraits,IsMultithreaded>
is a class template holding the data related to bitapproximation.
This class only contains static data and static member functions.

Notably, it determines the random directions we are using for SimHash.
This is essentially encoded in a set of orthogonal matrices (up to scaling), each of which is
(for efficiency reason) defined and stored as a sequence of special matrices (permutations,
diagonal matrices and the Walsh-Hadamard matrix) that can be applied to a vector fast.

Note: We may template this by the user lattice point class as well or
store these data inside some other class which uses the data.
Initializiation is performed by StaticInitializer<CoordinateSelection>.
This is called by ExactLatticePoint at the moment, but this is subject
to change.

Note2: Template arguments are subject to change.

Note3: Currently, we only use CoordinateSelection<void,false>
       This will change once control over bitapproximations is moved into the main list.

TODO: Consider writing a wrapper around either a std::bitset or a std::dynamic_bitset
      that actually has these as static data and initialize it via the main sieve.
********************************************************************/

template<class SieveTraits, bool MT>
class CoordinateSelection
{
public:
  friend StaticInitializer< CoordinateSelection<SieveTraits,MT> >;

  CoordinateSelection() = delete;  // cannot be instantiated.

  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  inline static auto transform_and_bitapprox(LatP const &point)
      -> std::array< std::bitset<sim_hash_len>, num_of_levels >;

  //template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  //inline static auto transform_and_bitapprox_simple(LatP const &point)
  //    -> std::array< std::bitset<sim_hash_len>, num_of_levels >;

  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  inline static auto transform_and_bitapprox_2nd_layer(LatP const &point)
      -> std::array< std::bitset<sim_hash_len>, num_of_levels >;

private:
  static unsigned int number_of_blocks;
  static unsigned int fast_walsh_hadamard_len;
  static std::vector< std::array<PMatrix,num_of_transforms> > pmatrices;
  static std::vector< std::array<DMatrix,num_of_transforms> > dmatrices;
  static unsigned int ambient_dimension;  // TODO: Remove this in favor of a static scratch space

  //static std::array<RMatrix,SimHash::num_of_levels> rmatrices;
};

template<class SieveTraits, bool MT>
unsigned int CoordinateSelection<SieveTraits,MT>::number_of_blocks;
template<class SieveTraits, bool MT>
std::vector< std::array<PMatrix,num_of_transforms> > CoordinateSelection<SieveTraits,MT>::pmatrices;
template<class SieveTraits, bool MT>
std::vector< std::array<DMatrix,num_of_transforms> > CoordinateSelection<SieveTraits,MT>::dmatrices;
template<class SieveTraits, bool MT>
unsigned int CoordinateSelection<SieveTraits,MT>::fast_walsh_hadamard_len;
template<class SieveTraits, bool MT>
unsigned int CoordinateSelection<SieveTraits,MT>::ambient_dimension;

//template<class SieveTraits, bool MT>
//std::array<RMatrix,SimHash::num_of_levels> CoordinateSelection<SieveTraits,MT>::rmatrices={};

}  // end namespace (GaussSieve::)BitApprox
   // now in namespace GaussSieve

/**************************************************************************************
Initialization of CoordinateSelection is done here.

Note: This is done in namespace GaussSieve rather than GaussSieve::SimHash, because
      the StaticInitializer template that is specialized is in the GaussSieve namespace.
***************************************************************************************/

template<class SieveTraits, bool MT>
class StaticInitializer< SimHash::CoordinateSelection<SieveTraits,MT> > final
: public DefaultStaticInitializer< SimHash::CoordinateSelection<SieveTraits,MT> >
{
  using Parent = DefaultStaticInitializer< SimHash::CoordinateSelection<SieveTraits,MT> >;
  using Data   = SimHash::CoordinateSelection<SieveTraits,MT>;
public:

  template<class T, TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const &initializer) : StaticInitializer(initializer.dim) {}

  StaticInitializer(int const ambient_dimension)
  {
    assert(Parent::user_count > 0);
    if (Parent::user_count > 1)
    {
      // check that the value is what is should be (see below)
      // TODO: Store dimension and verify it was the same as in the last invocation?
      // assert(Data::ambient_dimension == ambient_dimension);
      assert(Data::number_of_blocks ==
             (SimHash::num_of_levels * SimHash::sim_hash_len - 1) / ambient_dimension + 1);
    }
    else
    {
      // num_of_levels * sim_hash_len is the total number of required bits in the output simhashes.
      // We need ceil(total number of required bits / ambient_dimension) blocks
      // (i.e. repetitions of the input vector) to get at least that many output bits.
      // (X-1)/ Y + 1 is the same as ceil(X/Y)
      Data::ambient_dimension = ambient_dimension;
      Data::number_of_blocks =
          ((SimHash::num_of_levels * SimHash::sim_hash_len - 1) / ambient_dimension) + 1;

      Data::fast_walsh_hadamard_len =
          static_cast<unsigned int>(    pow(  2,  floor(log2(ambient_dimension))  )    );

      Data::pmatrices.resize(Data::number_of_blocks);
      Data::dmatrices.resize(Data::number_of_blocks);



      for (unsigned int i = 0; i < Data::number_of_blocks; ++i)
      {
        for (unsigned int j = 0; j < SimHash::num_of_transforms; ++j)
        {
          Data::pmatrices[i][j] = static_cast<SimHash::PMatrix>(ambient_dimension);
          Data::dmatrices[i][j] = static_cast<SimHash::DMatrix>(ambient_dimension);
        }
      }
      // TODO: Print if DEBUG symbol is set.

      //for (unsigned int j = 0; j < SimHash::num_of_levels; ++j)
      //  Data::rmatrices[j] = static_cast<SimHash::RMatrix>(ambient_dimension);
    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing CoordinateSelection; Counter is " << Parent::user_count )
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing CoordinateSelection; Counter is " << Parent::user_count )
  }
};


namespace SimHash  // entering GaussSieve::SimHash
{

/**
 fast_partial_walsh_hadamard<len>(input) performs (fast) Walsh-Hadamard Transform on the first
 len coordinates of its input. len is enforced to be a power of two and input must be either
 a std::vector or a std::array. The entries of the vector / array must support arithmentic (+,-) and
 be default-constructible and swappable.

 Note that the transformation is an orthogonal transformation, combined with a scaling by
 sqrt(len).
 TODO: Change the scaling?
*/

// TODO: Make part of the CoordinateSelection class

template<class T>
inline auto fast_partial_walsh_hadamard(std::vector<T> input, unsigned int const len)
    -> std::vector<T>
{
  //static // thread_local
  //std::vector<T> input=input_arg;  // TODO: Move to other static data.
  //input = input_arg;

  using std::swap;
  // len must be a power of two:
  assert(std::bitset< std::numeric_limits<unsigned long>::digits >{len}.count() == 1);
  // static_assert(is_a_power_of_two(len), "len must be a power of two");
  assert(len <= input .size());  // maybe static
  std::vector<T> output(input.size());

//  const double lengthfactor = std::sqrt(len);  // we have to properly rescale the unmodified coos.
  for (uint_fast16_t i = len; i < output.size(); ++i)
  {
//    input [i] = static_cast<T>(input[i] * lengthfactor);
    output[i] = input[i];
  }

  // on [0...len-1] coordinates perform WH

  for (uint_fast16_t i = len >> 1; i > 0; i >>= 1)
  {
    for (uint_fast16_t j = 0; j < len; j++)
    {
      // Note: static_cast<T> is required to be able to use mpz_class, because mpz_class internally
      // uses expression templates (so the types of a+b and a-b don't match for mpz_classes a,b)
      output[j] = ( (j/i)%2 != 0) ? static_cast<T>(input[j-i] - input[j]  )   // for j/i odd
                                  : static_cast<T>(input[j]   + input[i+j]);  // for j/i even
    }
    swap(input, output);
  }
  const double lengthfactor = std::sqrt(len);  // we have to properly rescale the modified coos.
  for(uint_fast16_t i = 0; i < len; ++i)
  {
    output[i] /= lengthfactor;
  }
  return output;
}

/*
// unused, did not incorporate recent changes.
// (unused) variant of the above function, for arrays instead of vectors.
template<class T, std::size_t arraylen>
inline auto fast_partial_walsh_hadamard(std::array<T,arraylen> input, unsigned int const len)
    -> std::array<T,arraylen>
{
  using std::swap;
  assert(std::bitset< std::numeric_limits<unsigned long>::digits>{len}.count() == 1);
  assert(len <= arraylen);
  std::array<T,arraylen> output{}; // assumes that entries are default-constructible.

  const double lengthfactor = std::sqrt(len);  // we have to properly rescale the unmodified coos.
  for (uint_fast16_t i = len; i < arraylen; ++i)
  {
    input [i] = static_cast<T>(input[i] * lengthfactor);
    output[i] = input[i];
  }

  for (uint_fast16_t i = len >> 1; i > 0; i >>= 1)
  {
    for (uint_fast16_t j = 0; j < len; j++)
    {
      // Note: static_cast<T> is required to be able to use mpz_class, because mpz_class internally
      // uses expression templates (so the types of a+b and a-b don't match for mpz_classes a,b)
      output[j] = ( (j/i)%2 != 0) ? static_cast<T>(input[j-i] - input[j]  )   // for j/i odd
                                  : static_cast<T>(input[j]   + input[i+j]);  // for j/i even
    }
    swap(input, output);
  }
  return output;
}
*/

template<class SieveTraits, bool MT>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto CoordinateSelection<SieveTraits,MT>::transform_and_bitapprox(LatP const &point)
    -> std::array< std::bitset<sim_hash_len>,num_of_levels >
{
  unsigned int const dim = static_cast<unsigned int>(point.get_dim());

  using ET    = Get_AbsoluteCooType<LatP>;
  using Block = std::vector<ET>;  // will be of length dim.
  Block copy_of_point(dim);
  for (unsigned int i = 0; i < dim; ++i)
  {
    copy_of_point[i] = point.get_absolute_coo(i);
  }

  // initialize blocks by number_of_blocks many copies of the input point.
  std::vector<Block> blocks(number_of_blocks, std::move(copy_of_point));

  // each block gets transformed by applying WH * D * P num_transform many times.
  for (uint_fast16_t i = 0; i < number_of_blocks; ++i)  // for each block:
  {
    for (uint_fast8_t j = 0; j < num_of_transforms; ++j)  // repeat num_of_transforms many times:
    {
      pmatrices[i][j].apply(blocks[i]);
      dmatrices[i][j].apply(blocks[i]);
      blocks[i] = fast_partial_walsh_hadamard(blocks[i], fast_walsh_hadamard_len);
    }
  }
  // put together the blocks into an array of bitsets.
  std::array< std::bitset<sim_hash_len>,num_of_levels > ret;
  for (unsigned int n = 0; n < num_of_levels; ++n)
  {
    for (unsigned int m = 0; m < sim_hash_len; ++m)
    {
      // index of the bit currently considered
      // if we use only one level of indexing.
      unsigned int const flat_bit_count = n * sim_hash_len + m;

      ret[n][m] = (blocks[flat_bit_count / dim][flat_bit_count % dim] > 0);
    }
  }
  return ret;
}

/*
template<class SieveTraits, bool MT>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto CoordinateSelection<SieveTraits,MT>::transform_and_bitapprox_simple(LatP const &point)
    -> std::array< std::bitset<sim_hash_len>,num_of_levels >
{
  using ET    = Get_AbsoluteCooType<LatP>;
  ET res;

  std::array< std::bitset<sim_hash_len>, num_of_levels > ret;
  for (unsigned int i = 0; i < num_of_levels; ++i)
  {
    for (unsigned int j = 0; j < sim_hash_len; ++j)
    {
      for (unsigned int k=0; k < SimHash::num_of_coos; ++k)
      {
        uint_fast16_t coo = CoordinateSelection<SieveTraits,MT>::rmatrices[i].get_value(i, k);
        res+=point[coo];
      }
      ret[i][j] = res>0;

    }

  }
  return ret;

}
 */
// TODO: Merge this function with the one above to avoid recomputations.
// NEW TODO: to delete, not sensitive

#if 0  // deactivated
template<class SieveTraits, bool MT>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto CoordinateSelection<SieveTraits,MT>::transform_and_bitapprox_2nd_layer(LatP const &point)
    -> std::array< std::bitset<sim_hash_len>,num_of_levels >
{
  using std::abs;
  using std::max;
  unsigned int const dim = static_cast<unsigned int>(point.get_dim());

  using ET    = Get_AbsoluteCooType<LatP>;
  using Block = std::vector<ET>;  // will be of length dim.
  Block copy_of_point(dim);
  for (unsigned int i = 0; i < dim; ++i)
  {
    copy_of_point[i] = point.get_absolute_coo(i);
  }
  std::vector<Block> blocks(number_of_blocks, std::move(copy_of_point));

  for (uint_fast16_t i = 0; i < number_of_blocks; ++i)  // for each block:
  {
    for (uint_fast8_t j = 0; j < num_of_transforms; ++j)  // repeat num_of_transforms many times:
    {
      pmatrices[i][j].apply(blocks[i]);
      dmatrices[i][j].apply(blocks[i]);
      blocks[i] = fast_partial_walsh_hadamard(blocks[i], fast_walsh_hadamard_len);
    }
  }
  // put together the blocks into an array of bitsets.
  std::array< std::bitset<sim_hash_len>,num_of_levels > ret;
  for (unsigned int n = 0; n < num_of_levels; ++n)
  {
    // maximal entry of the outbut block.
    // Note: The walsh-hadamard transform used above is not normalized, so the elements are larger
    // than the original vector. That's why we do not compare to norm2() currently.
    // TODO : Reconsider the point where we change the sign. Currently, it's at
    // half the maximal entry. This may not be ideal.
    /*
    ET maxentry = 0;
    for (unsigned int m = 0; m < sim_hash_len; ++m)
    {
      unsigned int const flat_bit_count = n * sim_hash_len + m;
      maxentry = max(maxentry, static_cast<ET>(abs(blocks[flat_bit_count / dim][flat_bit_count % dim])));
    }
    */
    ET const norm2 = point.get_norm2();

    for (unsigned int m = 0; m < sim_hash_len; ++m)
    {
      // index of the bit currently considered
      // if we use only one level of indexing.
      unsigned int const flat_bit_count = n * sim_hash_len + m;
      ret[n][m] = ( blocks[flat_bit_count / dim][flat_bit_count % dim]
                  * blocks[flat_bit_count / dim][flat_bit_count % dim] * dim * 6 > norm2);
    }
//    std::cout << ret[0].count() << std::endl << std::flush;
  }
  return ret;
}
#endif // 0

/**
  This class stores the result of computing a scalar product of bitwise
  approximations.
  Essentially, just wraps around an int.
  Note: We might want to wrap an approximate scalar product of t as value = #bits - t. -- Gotti
*/

class BitApproxScalarProduct
{

public:
  using BitApproxScalarProduct_WrappedType = uint_fast32_t;


  BitApproxScalarProduct(BitApproxScalarProduct const &) = delete;  // why not copy ints? -- Gotti
  BitApproxScalarProduct(BitApproxScalarProduct &&)      = default;

  explicit constexpr BitApproxScalarProduct(BitApproxScalarProduct_WrappedType const rhs):value(rhs) {}
  explicit operator  BitApproxScalarProduct_WrappedType() { return value; }

  BitApproxScalarProduct &operator=(BitApproxScalarProduct const &other) = delete;  // Why? -- Gotti
  BitApproxScalarProduct &operator=(BitApproxScalarProduct &&other) = default;

  // TODO: operator >=, <=

  inline bool operator<=(BitApproxScalarProduct_WrappedType &&rhs) const
  {
    return (this->value <= rhs);
  }

//  friend std::ostream& operator<<(std::ostream &os, BitApproxScalarProduct const &value)
//  {
//    os << value;
//    return os;
//  }

  BitApproxScalarProduct_WrappedType value;
};

template<class SieveTraits, bool MT>
using SimHashBlock = std::bitset<SieveTraits::sim_hash_len>;



}  // end namespace (GaussSieve::)SimHash
}  // end namespace GaussSieve


// OLD UNUSED CODE


// helpers for the function below, do not use.
namespace GaussSieve{ namespace SimHash { namespace Helpers{
template<int SizeOfBitSet> struct MakeBitApprox_Helper
{
  static_assert(SizeOfBitSet >= 0, "Only for fixed-size bit-sets.");
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline std::bitset<SizeOfBitSet> compute_bitapproximation(LatP const &point)
  {
    auto const dim = point.get_dim();
    static_assert(dim == SizeOfBitSet, "size of bitset must equal fixed ambient dim of vector");
    std::bitset<SizeOfBitSet> ret;
    for (uint_fast16_t i = 0; i < dim; ++i)
    {
      ret[i] = (point.get_absolute_coo(i) >= 0) ? 1 : 0;
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
    for (uint_fast16_t i = 0; i < dim; ++i)
    {
      ret[i] = (point.get_absolute_coo(i) >= 0) ? 1 : 0;
    }
    return ret;
  }
};
}}} // end namespaces GaussSieve::SimHash::Helpers

namespace GaussSieve{ namespace SimHash{

/**
  compute_coordinate_wise_bitapproximation<length>(point)
  computes the coordinate wise 1-bit-approximation of the lattice point point.
  It returns a std::bitset<length> for lenght>= 0 or boost::dynamic_bitset for length == -1.
  Currently, the length parameter (if >=0) has to match the (fixed!) dimension of the lattice point.
  For length==1, the dynamic_bitset that is returned has length equal to the dimension of point.
*/
template<int SizeOfBitSet, class LatP>
auto compute_coordinate_wise_bitapproximation(LatP const &point)
    -> decltype(  Helpers::MakeBitApprox_Helper<SizeOfBitSet>::compute_bitapproximation( std::declval<LatP>() )  )
{
  static_assert(IsALatticePoint<mystd::decay_t<LatP>>::value,"Not a lattice point.");
  return Helpers::MakeBitApprox_Helper<SizeOfBitSet>::compute_bitapproximation(point);
}
}}// end namespace GaussSieve::SimHash




namespace GaussSieve{ namespace SimHash{

// assume sim_hash2_len is a power-of-two.

// deprecated. Use fast_partial_walsh_hadamard instead.
template<class T>
[[deprecated]] inline std::vector<T> fast_walsh_hadamard_ext(std::vector<T> input, unsigned int len)
{
  std::vector<T> output (input.size());
  std::vector<T> tmp = input;
  uint_fast16_t i, j, s;

  for (i = (len >> 1); i > 0; i >>= 1)
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
  for (i = len; i < output.size(); ++i)
  {
    output[i] = tmp[i];
  }

  return output;
}

}}  // end namespace GaussSieve::SimHash


#endif // include guards

// clang-format on
