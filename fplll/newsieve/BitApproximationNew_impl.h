#pragma once
#ifndef GAUSS_SIEVE_BITAPPROX_NEW_H
#error wrong usage
#endif

namespace GaussSieve
{
namespace SimHashNew
{

template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
CoordinateSelection<sim_hash_len,sim_hash_num,MT,DimensionType_arg>::
    CoordinateSelection(DimensionType const dim, unsigned int random_seed)
{
  std::mt19937 rng;
  rng.seed(random_seed);
// sim_hash_num * sim_hash_len is the total number of required bits in the output simhashes.
// We need ceil(total number of required bits / ambient_dimension) blocks
// (i.e. repetitions of the input vector) to get at least that many output bits.
// (X-1)/ Y + 1 is the same as ceil(X/Y)
  unsigned int const ambient_dimension = dim;
  assert(ambient_dimension!=0);
  number_of_orthogonal_blocks =((sim_hash_num * sim_hash_len - 1) / ambient_dimension) + 1;
  fast_walsh_hadamard_loglen = static_cast<unsigned int>(floor(log2(ambient_dimension)));
  fast_walsh_hadamard_len =
      static_cast<unsigned int>(    pow(  2,  floor(log2(ambient_dimension))  )    );

  pmatrices.resize(number_of_orthogonal_blocks);
  dmatrices.resize(number_of_orthogonal_blocks);
  for (unsigned int i = 0; i < number_of_orthogonal_blocks; ++i)
  {
    for (unsigned int j = 0; j < num_of_transforms; ++j)
    {
      pmatrices[i][j] = PMatrix(ambient_dimension,rng);
      dmatrices[i][j] = DMatrix(ambient_dimension,rng);
    }
  }
}


/**
 fast_partial_walsh_hadamard<len>(input) performs (fast) Walsh-Hadamard Transform on the first
 len coordinates of its input. len is enforced to be a power of two and input must be either
 a std::vector or a std::array. The entries of the vector / array must support arithmentic (+,-) and
 be default-constructible and swappable.

 Note that the transformation is an orthogonal transformation, combined with a scaling by
 sqrt(len).
 TODO: Change the scaling?
*/

template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
template<class T>
inline auto CoordinateSelection<sim_hash_len,sim_hash_num,MT,DimensionType_arg>::
    fast_partial_walsh_hadamard(std::vector<T> input) const
    -> std::vector<T>
{
  unsigned int const len = fast_walsh_hadamard_len;
  //static // thread_local
  //std::vector<T> input=input_arg;  // TODO: Move to other static data.
  //input = input_arg;

  using std::swap;
  // len must be a power of two:
  assert(is_a_power_of_two(len));
  assert(len <= input.size());  // maybe static
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

template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
template<class T>
void inline CoordinateSelection<sim_hash_len, sim_hash_num, MT, DimensionType_arg>::
    faster_almost_partial_walsh_hadamard_inplace(std::vector<T> & input) const
{
  unsigned int const len = fast_walsh_hadamard_len;
  assert(is_a_power_of_two(len));
  assert(len <= input.size());  // maybe static


  // on [0...len-1] coordinates perform WH
  T tmp;

  // The following 3 versions are equivalent, up to permutation and signs
  // experimentally, the one selected was the fastest on my machine...

  // Version 1:

  // writing the indices in binary, for the WH-Transform, we perform
  // output[****0*****] = input[*****0*****] + input[*****1*****]
  // output[****1*****] = input[*****0*****] - input[*****1*****]
  // and then set input to output, where the affected "special" bit ranges over all possibilites,
  // starting with the msb.

  // Changing the processing order for the "special bit" and where the -sign appears above does
  // not matter for our application:
  // It is equivalent to a fixed permutation of coos and multiplying some coos by +/-1 before/after
  // the WH transformation; choose a variant that requires slightly less computation.

  // (Note that where the - sign appears here does not matter for our application)

  // i            = 0000010000 in the above notation
  // higher_bits is *****00000
  // lower_bits  is 000000****


  for (uint_fast16_t i = len >> 1; i > 0; i >>= 1)
//  for (uint_fast16_t i = 1; i < len; i<<=1) // loop order does not matter...
  {
    // this order of the 2 inner loops is 5-10% faster. No idea why.
    for(uint_fast16_t lower_bits = 0; lower_bits < i; ++lower_bits)
    {
      for(uint_fast16_t higher_bits =0; higher_bits < len; higher_bits += (2*i) )
      {
        tmp = input[higher_bits+0+lower_bits];
        input[higher_bits+0+lower_bits]+=input[higher_bits+i+lower_bits];
//        input[higher_bits+i+lower_bits] =tmp - input[higher_bits+i+lower_bits];  // <- This is WH trafo
        input[higher_bits+i+lower_bits]-=tmp;// - input[higher_bits+i+lower_bits]; // <- This is just as good, but 5-10% faster.
      }
    }
  }


  /*
  // Version 2 :
  unsigned int const log_len = fast_walsh_hadamard_loglen;
  for(uint_fast16_t i=0; i < log_len;++i)
  {
    unsigned int const ip = (1<<i);

    for(uint_fast16_t j=0; j < (1 << i); ++j)
    {
      for(uint_fast16_t k=0; k < (1 << (log_len-1-i)); ++k)
      {
        tmp = input[j+ (k<<(i+1))];
        input[j+(k<<(i+1))] += input[j+(k<<(i+1))+ip];
        input[j+(k<<(i+1))+ip]-=tmp;
      }
    }
  }
  */

  // Version 3:
  /*
  for(uint_fast16_t j = 0; j<len; ++j)
  {
    if ((j & i) !=0) continue;
    tmp = input[j];
    input[j]+=input[j^i];
    input[j^i]-=tmp;
  }
  */

  const double lengthfactor = std::sqrt(len);  // we have to properly rescale the modified coos.
  for(uint_fast16_t i = 0; i < len; ++i)
  {
    input[i] /= lengthfactor;
  }
  return;
}


template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto CoordinateSelection<sim_hash_len,sim_hash_num,MT,DimensionType_arg>::
    compute_all_bitapproximations(LatP const &point) const -> SimHashes
{
  unsigned int const dim = static_cast<unsigned int>(point.get_dim());

  using ET    = Get_AbsoluteCooType<LatP>;
  using Block = std::vector<ET>;  // will be of length dim.
  Block copy_of_point(dim);
  for (unsigned int i = 0; i < dim; ++i)
  {
    copy_of_point[i] = point.get_absolute_coo(i);
  }

  // initialize blocks by number_of_orthogonal_blocks many copies of the input point.
  std::vector<Block> blocks(number_of_orthogonal_blocks, std::move(copy_of_point));

  // each block gets transformed by applying WH * D * P num_transform many times.
  for (uint_fast16_t i = 0; i < number_of_orthogonal_blocks; ++i)  // for each block:
  {
    for (uint_fast8_t j = 0; j < num_of_transforms; ++j)  // repeat num_of_transforms many times:
    {
      pmatrices[i][j].apply(blocks[i]);
      dmatrices[i][j].apply(blocks[i]);
//      blocks[i] = fast_partial_walsh_hadamard(blocks[i]);
      faster_almost_partial_walsh_hadamard_inplace(blocks[i]);
    }
  }
  // put together the blocks into an array of bitsets.
  SimHashes ret;
  for (unsigned int n = 0; n < sim_hash_num; ++n)
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




PMatrix::PMatrix(unsigned int dim, std::mt19937 &rng)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to initialize P matrix")
  permutation.resize(dim);
  for (uint_fast16_t i = 0; i < dim; ++i)
  {
    permutation[i] = i;
  }
  std::shuffle(permutation.begin(), permutation.end(), rng);
#ifdef DEBUG_SIEVE_TRACEINITIATLIZATIONS
//  print();
#endif
}

template<class T>
inline void PMatrix::apply(std::vector<T> &vec) const
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

inline void PMatrix::print(std::ostream &os) const
{
  // os << "P  [" << i <<"] is: " << std::endl;
  for (uint_fast16_t k = 0; k < permutation.size(); ++k)
  {
    os << permutation[k] << " ";
  }
  os << std::endl;
}


DMatrix::DMatrix(unsigned int const dim, std::mt19937 &rng)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to fill-up the D matrix")
  std::uniform_int_distribution<uint_fast8_t> distr(0, 1);
  diagonal.resize(dim);
  for (uint_fast16_t i = 0; i < dim; ++i)
  {
    diagonal[i] = distr(rng);
  }
#ifdef DEBUG_SIEVE_TRACEINITIATLIZATIONS
//  print();
#endif
}

template<class T> inline void DMatrix::apply(std::vector<T> &vec) const
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

inline void DMatrix::print(std::ostream &os) const
{
  for (uint_fast16_t k = 0; k < diagonal.size(); ++k)
  {
    os << static_cast<int>(diagonal[k]) << " ";
  }
  os << std::endl;
}




}  // end namespace SimHash (within GaussSieve::)
}  // end namespace GaussSieve


/**
  UNUSED CODE STARTS HERE:
*/


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

#if 0  // deactivated
template<class SieveTraits, bool MT>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto CoordinateSelection<SieveTraits,MT>::transform_and_bitapprox_2nd_layer(LatP const &point)
    -> SimHashes<SieveTraits,MT>
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
  std::vector<Block> blocks(number_of_orthogonal_blocks, std::move(copy_of_point));

  for (uint_fast16_t i = 0; i < number_of_orthogonal_blocks; ++i)  // for each block:
  {
    for (uint_fast8_t j = 0; j < num_of_transforms; ++j)  // repeat num_of_transforms many times:
    {
      pmatrices[i][j].apply(blocks[i]);
      dmatrices[i][j].apply(blocks[i]);
      blocks[i] = fast_partial_walsh_hadamard(blocks[i], fast_walsh_hadamard_len);
    }
  }
  // put together the blocks into an array of bitsets.
  SimHashes<SieveTraits,MT> ret;
  for (unsigned int n = 0; n < sim_hash_num ; ++n)
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


/** OLD UNUSED CODE
*/

/*
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
*/


/**
  compute_coordinate_wise_bitapproximation<length>(point)
  computes the coordinate wise 1-bit-approximation of the lattice point point.
  It returns a std::bitset<length> for lenght>= 0 or boost::dynamic_bitset for length == -1.
  Currently, the length parameter (if >=0) has to match the (fixed!) dimension of the lattice point.
  For length==1, the dynamic_bitset that is returned has length equal to the dimension of point.
*/

/*
namespace GaussSieve{ namespace SimHash{
template<int SizeOfBitSet, class LatP>
auto compute_coordinate_wise_bitapproximation(LatP const &point)
    -> decltype(  Helpers::MakeBitApprox_Helper<SizeOfBitSet>::compute_bitapproximation( std::declval<LatP>() )  )
{
  static_assert(IsALatticePoint<mystd::decay_t<LatP>>::value,"Not a lattice point.");
  return Helpers::MakeBitApprox_Helper<SizeOfBitSet>::compute_bitapproximation(point);
}
}}// end namespace GaussSieve::SimHash
*/


/*
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
*/
