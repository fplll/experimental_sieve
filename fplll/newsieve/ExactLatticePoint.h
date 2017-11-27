/**
  Very simple lattice point class that stores precomputed values for its length.
  Otherwise, just a wrapper around a vector, really.
*/

#ifndef EXACT_LATTICE_POINT_H
#define EXACT_LATTICE_POINT_H


//#define EXACT_LATTICE_POINT_HAS_BITAPPROX
//#define EXACT_LATTICE_POINT_HAS_BITAPPROX_2ND_ORDER
#define EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED

#include "DefaultIncludes.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include <array>
#include <vector>
#include "PlainLatticePoint.h" // for conversions
#include "GlobalStaticData.h"

#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX
  #include <bitset> //for approximation
  #include <boost/dynamic_bitset.hpp> //for approximation
#endif

//#include "RelevantCoords.h"
#include "BitApproximation.h"

#define FOR_FIXED_DIM template <int X = nfixed, typename std::enable_if<X >= 0, int>::type = 0>
#define FOR_VARIABLE_DIM template <int X = nfixed, typename std::enable_if<X == -1, int>::type = 0>

namespace GaussSieve
{

// most simple Lattice Point that just wraps around a vector / array of ET's.
// dimension is static
template <class ET, int nfixed> class ExactLatticePoint;

// result of bitwise approximate scalar product
// wraps around int_fast32_t

template <class ET, int nfixed> class LatticePointTraits<ExactLatticePoint<ET, nfixed>>
{
public:
  using Trait_ScalarProductStorageType= ET;
  using Trait_AbsoluteCoos            = std::true_type;
  using Trait_CoordinateType          = ET;
  using Trait_CheapNorm2              = std::true_type;
  using Trait_CheapNegate             = std::true_type;

  using Trait_InternalRepIsAbsolute   = std::true_type;
  using Trait_InternalRepByCoos       = std::true_type;
  using Trait_InternalRepLinear       = std::true_type;
  using Trait_InternalRep_RW          = std::true_type;
  using Trait_AccessNorm2             = std::true_type;
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX
  using Trait_BitApprox               = std::true_type;
#endif
};

// Note: I renamed Approx to BitApprox here. We have two types of Approximations and everything related to Bits should have
// BitApprox in its name to reduce confusion. -- Gotti

template <class ET, int nfixed>
class ExactLatticePoint : public GeneralLatticePoint<ExactLatticePoint<ET, nfixed>>
{
public:
  friend StaticInitializer<ExactLatticePoint<ET,nfixed>>;
  using LatticePointTag         = std::true_type;
  using Container = mystd::conditional_t<nfixed >= 0,
        std::array<ET, nfixed >=0 ? nfixed:0>,  // if nfixed >= 0
        std::vector<ET>  >;                     // if nfixed <0
        // Note : The nfixed >=0 ? nfixed:0 is always nfixed;
        // The ?: expression is only needed to silence compiler errors/warnings.


#ifdef  EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
  using BitApproxContainerFixed   = std::bitset<SimHash::sim_hash_len>;
#endif

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
  explicit ExactLatticePoint()
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((StaticInitializer<ExactLatticePoint<ET,nfixed>>::is_initialized));
  }

  FOR_VARIABLE_DIM
  explicit ExactLatticePoint() :
  data(static_cast<unsigned int>(get_dim()))
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((StaticInitializer<ExactLatticePoint<ET,nfixed>>::is_initialized));
  }

  // TODO: Remove these constructors?
  FOR_FIXED_DIM
  explicit ExactLatticePoint(MaybeFixed<nfixed>) : ExactLatticePoint() {}

  FOR_VARIABLE_DIM
  explicit ExactLatticePoint(MaybeFixed<nfixed> dim) : ExactLatticePoint()
  {
#ifdef DEBUG_SIEVE_LP_MATCHDIM
    assert(dim == get_dim());
#endif
  }

// TODO: Debug output and validation.

  explicit ExactLatticePoint(PlainLatticePoint<ET,nfixed> &&plain_point)
  :
  data(std::move(plain_point.data))
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

  void sanitize()
  {
    norm2 = compute_sc_product(*this, *this);

 #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
    //fixed_bitapprox_data = SimHash::compute_fixed_bitapproximation(*this);
    fixed_bitapprox_data_level = SimHash::compute_fixed_bitapprox_level(*this);
#endif
  }
  void sanitize( ET const & new_norm2 )
  {
    norm2 = new_norm2;

#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
    //fixed_bitapprox_data = SimHash::compute_fixed_bitapproximation(*this);
    
    fixed_bitapprox_data_level = SimHash::compute_fixed_bitapprox_level(*this);
#endif

  }

  ET get_norm2() const { return norm2; }

  std::ostream& write_lp_to_stream (std::ostream &os, bool const include_norm2=true, bool const include_approx=true) const;


#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
  inline SimHash::BitApproxScalarProduct do_compute_sc_product_bitapprox_fixed(ExactLatticePoint const & another) const;
  
  inline SimHash::BitApproxScalarProduct do_compute_sc_product_bitapprox_level(ExactLatticePoint const & another, int lvl) const;
#endif

  inline ET do_compute_sc_product(ExactLatticePoint const &lp2) const
  {
  ET res1 = 0;
  ET res2 = 0;
  ET res3 = 0;
  ET res4 = 0;
  uint_fast16_t dim = get_dim();
  for(uint_fast16_t i=0; i < (dim/4) * 4;i+=4)
  {
    res1 += (*this)[i+0] * lp2[i+0];
    res2 += (*this)[i+1] * lp2[i+1];
    res3 += (*this)[i+2] * lp2[i+2];
    res4 += (*this)[i+3] * lp2[i+3];
  }
  for(uint_fast8_t i= (dim/4) * 4; i < dim; ++i)
  {
    res1+= (*this)[i] * lp2[i];
  }
  res1+=res2;
  res1+=res3;
  res1+=res4;
  return res1;
  }

  // moved to LatticePointConcept.h
  // The reason is that we want a compute_bitapproximation function
  // that also works with other lattice point classes:
  // In the constructors, we want to be able to set bitapprox_data before data
  // in the initializer list. (I was getting tons of compiler warnings...) This only works well with a non-member function.
  // Alternatively, move it back to ExactLatticePoint.h, but make it a static function.

private:
  static MaybeFixed<nfixed> dim;  // note that for nfixed != -1, this variable is actually unused.

/*
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX
  BitApproxContainer bitapprox_data;
  #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_2ND_ORDER
  BitApproxContainer bitapprox2_data;
  #endif
#endif
*/


#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
  BitApproxContainerFixed fixed_bitapprox_data;

  std::array<BitApproxContainerFixed, SimHash::num_of_levels> fixed_bitapprox_data_level;
#endif

  Container data;
  ET norm2;
};

/*
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX
// bit-approximate scalar product
template <class ET, int nfixed>
inline SimHash::BitApproxScalarProduct ExactLatticePoint<ET, nfixed>::do_compute_sc_product_bitapprox(ExactLatticePoint const &another) const
{
  auto const dim = this->get_dim();
//  BitApproxScalarProduct result(0);
//  result = dim - (this.bitapprox_data ^ another.bitapprox_data).count();
//  return result;
  return SimHash::BitApproxScalarProduct{ static_cast<size_t>(dim - (this->bitapprox_data ^ another.bitapprox_data).count()) };
}

// 2nd order bit-approximate scalar product
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_2ND_ORDER
template <class ET, int nfixed>
inline SimHash::BitApproxScalarProduct ExactLatticePoint<ET, nfixed>::do_compute_sc_product_bitapprox_2nd_order(ExactLatticePoint const &another) const
{
  BitApproxContainer not_xor_res = (bitapprox_data ^ another.bitapprox_data).flip();
  //std::cout << not_xor_res << std::endl;
  //BitApproxContainer not_xor_res = (bitapprox_data ^ another.bitapprox_data);
  return SimHash::BitApproxScalarProduct{ static_cast<size_t>( (not_xor_res & another.bitapprox2_data).count() +
                                                      (not_xor_res & this->bitapprox2_data).count() ) +
                                                      not_xor_res.count() };
}
#endif
#endif
*/
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
template <class ET, int nfixed>
inline SimHash::BitApproxScalarProduct ExactLatticePoint<ET, nfixed>::do_compute_sc_product_bitapprox_fixed(ExactLatticePoint const &another) const
{
  //std::cout << "sim-hash1 = " << this->fixed_bitapprox_data << std::endl;
  //std::cout << "sim-hash1 = " << another.fixed_bitapprox_data << std::endl;
  return SimHash::BitApproxScalarProduct {static_cast<size_t>(sim_hash_len - (this->fixed_bitapprox_data ^ another.fixed_bitapprox_data).count()) };
}

template <class ET, int nfixed>
inline SimHash::BitApproxScalarProduct ExactLatticePoint<ET, nfixed>::do_compute_sc_product_bitapprox_level(ExactLatticePoint const & another, int lvl) const
{
  return SimHash::BitApproxScalarProduct {static_cast<size_t>(SimHash::sim_hash_len - (
                                            this->fixed_bitapprox_data_level[lvl] ^ another.fixed_bitapprox_data_level[lvl]).count())};
}


#endif

template <class ET, int nfixed>
inline std::ostream& ExactLatticePoint<ET, nfixed>::write_lp_to_stream(std::ostream &os, bool const include_norm2, bool const include_approx) const
{

  GeneralLatticePoint<ExactLatticePoint<ET,nfixed>>::write_lp_to_stream(os,include_norm2,include_approx);

#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX
  if (include_approx)
  {
    os << " bit-approx = [ ";
    for (uint_fast16_t i =0; i<dim; ++i)
    {
      os << bitapprox_data[i] << " ";
    }
    os << "]";
    os << std::endl;

    #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_2ND_ORDER
    os << " bit-approx2 = [ ";
    for (uint_fast16_t i =0; i<dim; ++i)
    {
      os << bitapprox2_data[i] << " ";
    }
    os << "]";
    #endif
  }
#endif
#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
  
  if (include_approx)
  {
    os << std::endl;
    for (uint_fast16_t lvl =0; lvl<SimHash::num_of_levels; ++lvl)
    {
    os<< "lvl =  " << lvl << " [ ";
    for (uint_fast16_t i=0; i<SimHash::sim_hash_len; i++) {
      os << fixed_bitapprox_data_level[lvl][i] << " ";
    }
    os << "] "<<std::endl;
    }
    
  }
   
#endif
  os << std::endl;
  return os;
}


template <class ET, int nfixed>
std::ostream& operator<<(std::ostream &os, ExactLatticePoint<ET, nfixed> const &LatP)
{
  return LatP.write_lp_to_stream(os, true, true);
}

// initialize static data:
template <class ET, int nfixed>
MaybeFixed<nfixed> ExactLatticePoint<ET, nfixed>::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);

// Static Initializer:
template<class ET, int nfixed> class StaticInitializer<ExactLatticePoint<ET,nfixed>>
  : public DefaultStaticInitializer<ExactLatticePoint<ET,nfixed>>
{
  using Parent = DefaultStaticInitializer<ExactLatticePoint<ET,nfixed>>;
  public:

  template<class T,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const & initializer) : StaticInitializer(initializer.dim) {}

  StaticInitializer(MaybeFixed<nfixed> const new_dim):   init_D_matrices(new_dim), init_P_matrices(new_dim)
  {

    assert(Parent::user_count > 0);
    if(Parent::user_count>1)
    {
      assert((new_dim == ExactLatticePoint<ET,nfixed>::dim));
      // TODO: Throw exception!
    }
    else
    {
      ExactLatticePoint<ET,nfixed>::dim = new_dim;

    }
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing ExactLatticePoint with nfixed = " << nfixed << " REALDIM = " << new_dim << " Counter is" << Parent::user_count )
  }
  ~StaticInitializer()
  {
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing ExactLatticePoint with nfixed = " << nfixed << " Counter is " << Parent::user_count )
  }
  
  
  GaussSieve::StaticInitializer<DMatrix> init_D_matrices;
  GaussSieve::StaticInitializer<PMatrix> init_P_matrices;
  //GaussSieve::StaticInitializer<RelevantCoordinates> init_relevant_coo_matrix;
};

}  // end of namespace

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM

#endif
