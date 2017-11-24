#ifndef GAUSS_SIEVE_BITAPPROX_H
#define GAUSS_SIEVE_BITAPPROX_H

#include <random>
#include "GlobalStaticData.h"
#include "DefaultIncludes.h"

namespace GaussSieve{

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
//template<int nxfixed>
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
      //assert(false);
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


// Fixed-dim version:
template<int SizeOfBitSet> struct BitApproximation
{
  static_assert(SizeOfBitSet>=0, "Only for fixed-size bit-sets.");
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline std::bitset<SizeOfBitSet> compute_bitapproximation(LatP const &point)
  {
    auto dim = point.get_dim();
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
template<> struct BitApproximation<-1>
{

  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline boost::dynamic_bitset<> compute_bitapproximation(LatP const &point)
  {
    auto dim = point.get_dim();
    boost::dynamic_bitset<> ret{static_cast<size_t>(dim) };
//    ret.resize(dim);
    for(uint_fast16_t i=0;i<dim;++i)
    {
      ret[i] = (point.get_absolute_coo(i)>=0) ? 1 : 0;
    }
    return ret;
  }

  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline std::bitset<sim_hash_len> compute_fixed_bitapproximation(LatP const &point)
  {
    //using RelevantCoords = GaussSieve::RelevantCoordinates;
    using ET = Get_CoordinateType<LatP>;




    //std::cout << "rel_coo_matrix used: " <<std::endl;
    //RelevantCoordinates::print();
    //assert(false);

    std::bitset<sim_hash_len> ret;

    uint_fast16_t bound = std::min(static_cast<uint_fast16_t>( point.get_dim()), sim_hash_len);


    for(uint_fast16_t i=0;i<bound;++i)
    {
        ret[i] = (point.get_absolute_coo(i)>=0) ? 1 : 0;
    }

    //for(uint_fast16_t i=0;i<sim_hash_len;++i)
    for(uint_fast16_t i=point.get_dim();i<sim_hash_len;++i)
    {
      ET res = point.get_absolute_coo(RelevantCoordinates::get_ij_value(i,0)) +
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

    std::array<ET, sim_hash2_len> hadamard = fast_walsh_hadamard<ET>(input_vector);
    for(uint_fast16_t i=0;i<sim_hash_len;++i)
    {
      ret[i] = (hadamard[i]>=0) ? 1: 0;
    }

    return ret;
  }


  //assume sim_hash2_len is a power-of-two
  template<class ET>
  static std::array<ET,sim_hash2_len> fast_walsh_hadamard(std::array<ET,sim_hash2_len> const &input)
  {
    std::array<ET,sim_hash2_len> inp = input;
    std::array<ET,sim_hash2_len> out;
    std::array<ET,sim_hash2_len> tmp;

    uint_fast16_t i, j, s;

    for (i = sim_hash2_len>>1; i > 0; i>>=1) {
        for (j = 0; j < sim_hash2_len; j++) {
            s = j/i%2;
            out[j]=inp[(s?-i:0)+j]+(s?-1:1)*inp[(s?0:i)+j];
        }
        tmp = inp; inp = out; out = tmp;
    }

    return out;
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

  /*
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline boost::dynamic_bitset<> compute_2nd_order_bitapproximation(LatP const &point)
  {
    boost::dynamic_bitset<> ret{64};
    using ET = Get_CoordinateType<LatP>;

    return ret;
  }
   */
};




}


#endif // include guards
