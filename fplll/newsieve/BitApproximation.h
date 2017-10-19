#ifndef BITSET_APPROXIMATION_H
#define BITSET_APPROXIMATION_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include <cstdint>
#include <limits>
#include <cmath>
#include <vector>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include "GlobalStaticData.h"


/**
  This defines a lattice point approximation, where the approximation is a binary vector of sgns:
   * approx[i]=0 if true_vector[i]>=0
*/

namespace GaussSieve{

template<int nfixed> class BitApproximation;
  
  
#define FOR_FIXED_DIM template <int X = nfixed, typename std::enable_if<X >= 0, int>::type = 0>
#define FOR_VARIABLE_DIM template <int X = nfixed, typename std::enable_if<X == -1, int>::type = 0>


template<int nfixed>
class BitApproximation
{
  friend
  StaticInitializer<BitApproximation<nfixed>>;
  public:

  //TODO
  using ScalarProductType = int_fast32_t;

  private:
  using ApproxEntryType = bool;
  using ApproxNorm2Type = int_fast32_t;

  using Container = typename std::conditional<nfixed >= 0,
                      std::bitset<nfixed >=0 ? nfixed:0>,  // if nfixed >= 0
                      boost::dynamic_bitset<>  >                       // if nfixed <  0
                      ::type;

  public:
  template<class LatticePoint> 
  explicit BitApproximation(LatticePoint const &exact_point);
  
  BitApproximation() = delete;
  BitApproximation(BitApproximation const &old) = delete;
  BitApproximation(BitApproximation && old) = default;
  BitApproximation& operator=(BitApproximation const &other) = delete;
  BitApproximation& operator=(BitApproximation && other) = default;
  ~BitApproximation() {};
  
  ApproxEntryType &operator[](uint_fast16_t idx) { return data[idx]; };
  ApproxEntryType const &operator[](uint_fast16_t idx) const { return data[idx]; };
  
  FOR_FIXED_DIM
  static constexpr MaybeFixed<nfixed> get_dim() { return MaybeFixed<nfixed>(nfixed); }
  
  FOR_VARIABLE_DIM
  static MaybeFixed<-1> get_dim() { return dim; }
  
  FOR_FIXED_DIM
  void reserve_size() {}
  
  FOR_VARIABLE_DIM
  void reserve_size() {data.reserve(static_cast<typename Container::size_type>(dim));}
  
public:
  
  
  Container data; // array of bits
private:
  static MaybeFixed<nfixed> dim; // dimension of data
  // TODO: Remove (DefaultStaticInit's job)
#ifdef DEBUG_SIEVE_LP_INIT
  static bool class_initialized;
#endif // DEBUG_SIEVE_LP_INIT


};


}
#endif

