#ifndef BITSET_APPROXIMATION_H
#define BITSET_APPROXIMATION_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include <cstdint>
#include <limits>
#include <cmath>
#include <vector>
#include <bitset>
#include "GlobalStaticData.h"


/**
  This defines a lattice point approximation, where the approximation is a binary vector of sgns:
   * approx[i]=0 if true_vector[i]>=0
*/

namespace GaussSieve{

template<int nfixed> class BitApproximation;

template<int nfixed>
class BitApproximation
{
  friend
  StaticInitializer<BitApproximation<nfixed>>;
  public:

  //TODO
  using ScalarProductType = int_fast32_t;

  private:
  //using ApproxEntryType = typename EMVApproximationTraits::ApproxEntryType;
  //using ApproxNorm2Type = typename EMVApproximationTraits::ApproxNorm2Type;

  using Container = typename std::conditional<nfixed >= 0,
                      std::bitset<nfixed >=0 ? nfixed:0>,  // if nfixed >= 0
                      std::vector<bool> >                       // if nfixed <  0
                      ::type;

  public:
  template<class LatticePoint> 
  explicit BitApproximation(LatticePoint const &exact_point);


}


}
#endif

