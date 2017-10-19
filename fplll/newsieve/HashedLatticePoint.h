/**
  Very simple lattice point class that stores precomputed values for its length.
  Otherwise, just a wrapper around a vector, really.
*/

#ifndef HASHED_LATTICE_POINT_H
#define HASHED_LATTICE_POINT_H

#include "DefaultIncludes.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include <array>
#include <vector>
#include "PlainLatticePoint.h" // for conversions
#include "GlobalStaticData.h"

#define FOR_FIXED_DIM template <int X = nfixed, typename std::enable_if<X >= 0, int>::type = 0>
#define FOR_VARIABLE_DIM template <int X = nfixed, typename std::enable_if<X == -1, int>::type = 0>


namespace GaussSieve
{

// most simple Lattice Point that just wraps around a vector / array of ET's.
// dimension is static
template <class ET, int nfixed> class HashedLatticePoint;

template <class ET, int nfixed> class LatticePointTraits<HashedLatticePoint<ET, nfixed>>
{
public:
  using Trait_ScalarProductStorageType= ET;
  using Trait_CoordinateVector        = std::true_type;
  using Trait_CoordinateAccess        = std::true_type;
  using Trait_AbsoluteCoos            = std::true_type;
  using Trait_CoordinateType          = ET;
  using Trait_CheapNorm2              = std::true_type;
  using Trait_CheapNegate             = std::true_type;
};

template <class ET, int nfixed>
class HashedLatticePoint : public GeneralLatticePoint<HashedLatticePoint<ET, nfixed>>
{
public:
  friend StaticInitializer<HashedLatticePoint<ET,nfixed>>;
  using LatticePointTag         = std::true_type;
  using ScalarProductStorageType = ET;
  using Container = typename std::conditional<nfixed >= 0,
        std::array<ET, nfixed >=0 ? nfixed:0>,  // if nfixed >= 0
        std::vector<ET>  >::type;               // if nfixed <0
        // Note : The nfixed >=0 ? nfixed:0 is always nfixed;
        // The ?: expression is only needed to silence compiler errors/warnings.


  //TODO: MAKE PROPERLY
  static constexpr unsigned short number_of_hash_tables = 102;

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

  FOR_FIXED_DIM
  explicit HashedLatticePoint()
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((StaticInitializer<HashedLatticePoint<ET,nfixed>>::is_initialized));
  }

  FOR_VARIABLE_DIM
  explicit HashedLatticePoint() : data(static_cast<unsigned int>(get_dim()))
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((StaticInitializer<HashedLatticePoint<ET,nfixed>>::is_initialized));
  }

  FOR_FIXED_DIM
  explicit HashedLatticePoint(MaybeFixed<nfixed>) : HashedLatticePoint() {}

  FOR_VARIABLE_DIM
  explicit HashedLatticePoint(MaybeFixed<nfixed> dim) : HashedLatticePoint()
  {
#ifdef DEBUG_SIEVE_LP_MATCHDIM
    assert(dim == get_dim());
#endif
  }

// TODO: Debug output and validation.

  explicit HashedLatticePoint(PlainLatticePoint<ET,nfixed> &&plain_point)
  : data(std::move(plain_point.data))
  {
    sanitize();
  }

  HashedLatticePoint(HashedLatticePoint const &old) = delete;
  HashedLatticePoint(HashedLatticePoint &&old)      = default;
  HashedLatticePoint &operator=(HashedLatticePoint const &other) = delete;
  HashedLatticePoint &operator=(HashedLatticePoint &&other) = default;
  ET &operator[](uint_fast16_t idx) { return data[idx]; };
  ET const &operator[](uint_fast16_t idx) const { return data[idx]; };
  static std::string class_name() { return "Exact Lattice Point"; };

  void sanitize() { norm2 = compute_sc_product(*this, *this); }
  void sanitize( ScalarProductStorageType const & new_norm2 ) { norm2 = new_norm2; }

  ET get_norm2() const { return norm2; }


  inline ET compute_sc_product(HashedLatticePoint const &lp1, HashedLatticePoint const &lp2)
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

  public:
    std::array<int, number_of_hash_tables> positions;

private:
  static MaybeFixed<nfixed> dim;  // note that for nfixed != -1, this variable is actually unused.
  Container data;
  ET norm2;
};

// initialize static data:
template <class ET, int nfixed>
MaybeFixed<nfixed> HashedLatticePoint<ET, nfixed>::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);


// Static Initializer:
template<class ET, int nfixed> class StaticInitializer<HashedLatticePoint<ET,nfixed>>
  : public DefaultStaticInitializer<HashedLatticePoint<ET,nfixed>>
{
  using Parent = DefaultStaticInitializer<HashedLatticePoint<ET,nfixed>>;
  public:

  template<class T,TEMPL_RESTRICT_DECL(IsArgForStaticInitializer<T>::value)>
  StaticInitializer(T const & initializer) : StaticInitializer(initializer.dim) {}

  StaticInitializer(MaybeFixed<nfixed> const new_dim)
  {
    assert(Parent::user_count > 0);
    if(Parent::user_count>1)
    {
      assert((new_dim == HashedLatticePoint<ET,nfixed>::dim));
      // TODO: Throw exception!
    }
    else
    {
      HashedLatticePoint<ET,nfixed>::dim = new_dim;
    }
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing HashedLatticePoint with nfixed =" << nfixed << "REALDIM = " << new_dim << "Counter is" << Parent::user_count )
  }
  ~StaticInitializer()
  {
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing HashedLatticePoint with nfixed =" << nfixed << "Counter is" << Parent::user_count )
  }
};

}  // end of namespace

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM

#endif
