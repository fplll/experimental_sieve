/**
  Very simple lattice point class that stores precomputed values for its length.
  Otherwise, just a wrapper around a vector, really.
*/

#ifndef EXACT_LATTICE_POINT_H
#define EXACT_LATTICE_POINT_H

#include "DebugAll.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include "assert.h"
#include <array>
#include <vector>
#include <iostream>
#include "PlainLatticePoint.h" // for conversions

#define FOR_FIXED_DIM template <int X = nfixed, typename std::enable_if<X >= 0, int>::type = 0>
#define FOR_VARIABLE_DIM template <int X = nfixed, typename std::enable_if<X == -1, int>::type = 0>

namespace GaussSieve
{

// most simple Lattice Point that just wraps around a vector / array of ET's.
// dimension is static
template <class ET, int nfixed> class ExactLatticePoint;

template <class ET, int nfixed> class LatticePointTraits<ExactLatticePoint<ET, nfixed>>
{
public:
  using AuxDataType             = MaybeFixed<nfixed>;
  using ScalarProductReturnType = ET;
  using CoordinateVector        = std::true_type;
  using CoordinateAccess        = std::true_type;
  using AbsoluteCoos            = std::true_type;
  using CoordinateType          = ET;
  using CheapNorm2              = std::true_type;
  using CheapNegate             = std::true_type;
};

template <class ET, int nfixed>
class ExactLatticePoint : public GeneralLatticePoint<ExactLatticePoint<ET, nfixed>>
{
public:
  friend StaticInitializer<ExactLatticePoint<ET,nfixed>>;
  using LatticePointTag         = std::true_type;
  using AuxDataType             = typename GetAuxDataType<ExactLatticePoint>::type;
  using ScalarProductReturnType = ET;
  using Container =
    typename std::conditional<nfixed >= 0, std::array<ET, nfixed >=0 ? nfixed:0>,  // if nfixed >= 0
                                std::vector<ET>>::type;               // if nfixed <0

/*
  static bool class_init(AuxDataType const aux_data)
  {
    if(user_counter>0)
    {
      if(dim!=aux_data)
      {
#ifdef DEBUG_SIEVE_LP_INIT
        assert(class_initialized);
        std::cerr << "Warning: Class initialization failed for " << class_name()
                  << std::endl << std::flush;
#endif
        return false;
      }
    }
#ifdef DEBUG_SIEVE_LP_INIT
    class_initialized = true;
#endif
    ++user_counter;
    dim = aux_data;
    return true;
  };

  static bool class_uninit()
  {
#ifdef DEBUG_SIEVE_LP_INIT
    assert(class_initialized);
#endif
    assert(user_counter>0);
    --user_counter;
#ifdef DEBUG_SIEVE_LP_INIT
    if(user_counter==0)
    {
      class_initialized=false;
    }
#endif
    return(user_counter==0);
  }
  */

  FOR_FIXED_DIM
  static constexpr MaybeFixed<nfixed> get_dim()
  {
    static_assert(X == nfixed, "");
    return MaybeFixed<nfixed>(nfixed);
  }

  FOR_VARIABLE_DIM
  static MaybeFixed<-1> get_dim()
  {
    static_assert(nfixed == -1, "");
    return dim;
  }

#ifdef DEBUG_SIEVE_LP_INIT
  FOR_FIXED_DIM
  explicit ExactLatticePoint()
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((ExactLatticePoint<ET, nfixed>::class_initialized));
  }

  FOR_VARIABLE_DIM
  explicit ExactLatticePoint() : data(static_cast<unsigned int>(get_dim()))
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((ExactLatticePoint<ET, nfixed>::class_initialized));
  }
#else
  FOR_FIXED_DIM
  explicit ExactLatticePoint() = default;

  FOR_VARIABLE_DIM
  explicit ExactLatticePoint() : data(static_cast<unsigned int>(get_dim())) {}
#endif

  FOR_FIXED_DIM
  explicit ExactLatticePoint(MaybeFixed<nfixed>)
  {
#ifdef DEBUG_SIEVE_LP_INIT
    assert((ExactLatticePoint<ET, nfixed>::class_initialized));
#endif
  };

  FOR_VARIABLE_DIM
  explicit ExactLatticePoint(MaybeFixed<nfixed> dim) : data(static_cast<unsigned int>(dim))
  {
    static_assert(nfixed == -1, "");
#ifdef DEBUG_SIEVE_LP_INIT
    assert((ExactLatticePoint<ET, -1>::class_initialized));
#endif
#ifdef DEBUG_SIEVE_LP_MATCHDIM
    assert(dim == get_dim());
#endif
  }

// TODO: Debug output and validation.

  explicit ExactLatticePoint(PlainLatticePoint<ET,nfixed> &&plain_point)
  : data(std::move(plain_point.data))
  {
    sanitize();
  }

  ExactLatticePoint(ExactLatticePoint const &old) = delete;
  ExactLatticePoint(ExactLatticePoint &&old)      = default;
  ExactLatticePoint &operator=(ExactLatticePoint const &other) = delete;
  ExactLatticePoint &operator=(ExactLatticePoint &&other) = default;
  ET &operator[](uint_fast16_t idx) { return data[idx]; };
  ET const &operator[](uint_fast16_t idx) const { return data[idx]; };
  static std::string class_name() { return "Exact Lattice Point"; };

  void sanitize() { norm2 = compute_sc_product(*this, *this); }
  void sanitize( ScalarProductReturnType const & new_norm2 ) { norm2 = new_norm2; }

  ET get_norm2() const { return norm2; }


  ET compute_sc_product(ExactLatticePoint const &lp1, ExactLatticePoint const &lp2)
  {
  ET res1 = 0;
  ET res2 = 0;
  ET res3 = 0;
  ET res4 = 0;
  uint_fast16_t dim = get_dim();
  for(uint_fast16_t i=0; i < (dim/4) * 4;i+=4)
  {
    res1 += lp1[i+0] * lp2[i+0];
    res2 += lp1[i+1] * lp2[i+1];
    res3 += lp1[i+2] * lp2[i+2];
    res4 += lp1[i+3] * lp2[i+3];
  }
  for(uint_fast8_t i= (dim/4) * 4; i < dim; ++i)
  {
    res1+= lp1[i] * lp2[i];
  }
  res1+=res2;
  res1+=res3;
  res1+=res4;
  return res1;
  }


private:
  static MaybeFixed<nfixed> dim;  // note that for nfixed != -1, this variable is actually unused.
#ifdef DEBUG_SIEVE_LP_INIT
  static bool class_initialized;
#endif  // DEBUG_SIEVE_LP_INIT
//  static unsigned int user_counter;
  Container data;
  ET norm2;
};

// Static Initializer:
template<class ET, int nfixed> class StaticInitializer<ExactLatticePoint<ET,nfixed>>
{
  public:
  StaticInitializer(MaybeFixed<nfixed> const new_dim)
  {
    if(user_counter>0)
    {
      assert((new_dim == ExactLatticePoint<ET,nfixed>::dim));
      // TODO: Throw exception!
    }
    else
    {
      ExactLatticePoint<ET,nfixed>::dim = new_dim;
    }
    ++user_counter;
#ifdef DEBUG_SIEVE_LP_INIT
    ExactLatticePoint<ET,nfixed>::class_initialized = true;
#endif
  }

  ~StaticInitializer()
  {
    assert(user_counter > 0);
    --user_counter;
#ifdef DEBUG_SIEVE_LP_INIT
    ExactLatticePoint<ET,nfixed>::class_initialized = (user_counter > 0);
#endif
  }

  static bool is_initialized(){ return user_counter > 0; }; // Does an object exist?

  private:
  static unsigned int user_counter; // counts number of StaticInitializer instances.
};


// initialize static data:

template <class ET, int nfixed>
MaybeFixed<nfixed> ExactLatticePoint<ET, nfixed>::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);

template <class ET, int nfixed>
unsigned int StaticInitializer<ExactLatticePoint<ET,nfixed>>::user_counter = 0;


#ifdef DEBUG_SIEVE_LP_INIT
template <class ET, int nfixed> bool ExactLatticePoint<ET, nfixed>::class_initialized = false;
#endif  // DEBUG_SIEVE_LP_INIT


}  // end of namespace

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM

#endif
