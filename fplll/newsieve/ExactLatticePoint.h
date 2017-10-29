/**
  Very simple lattice point class that stores precomputed values for its length.
  Otherwise, just a wrapper around a vector, really.
*/

#ifndef EXACT_LATTICE_POINT_H
#define EXACT_LATTICE_POINT_H

#include "DefaultIncludes.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include <array>
#include <vector>
#include "PlainLatticePoint.h" // for conversions
#include "GlobalStaticData.h"

#include <bitset> //for approximation
#include <boost/dynamic_bitset.hpp> //for approximation


#define FOR_FIXED_DIM template <int X = nfixed, typename std::enable_if<X >= 0, int>::type = 0>
#define FOR_VARIABLE_DIM template <int X = nfixed, typename std::enable_if<X == -1, int>::type = 0>

namespace GaussSieve
{

// most simple Lattice Point that just wraps around a vector / array of ET's.
// dimension is static
template <class ET, int nfixed> class ExactLatticePoint;

// result of approximate scalar product
// wraps around int_fast32_t
class ApproxScProduct;

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
  using Trait_BitApprox               = std::true_type;
};

class ApproxScProduct
{
  
  public:
    using ApproxScProductReturn       = int_fast32_t;
    
  
  ApproxScProduct(ApproxScProduct const &old) = delete;
  ApproxScProduct(ApproxScProduct &&old)      = default;
  
  ApproxScProduct(ApproxScProductReturn rhs) { this->value = rhs;};
  
  ApproxScProduct &operator=(ApproxScProduct const &other) = delete;
  ApproxScProduct &operator=(ApproxScProduct &&other) = default;
  
  //TODO: operator >=, <=
  
  inline bool operator<=(ApproxScProductReturn && rhs)
  {
    return  this->value <= rhs;
  } 
  
  //member
  ApproxScProductReturn value;
  
    
};

template <class ET, int nfixed>
class ExactLatticePoint : public GeneralLatticePoint<ExactLatticePoint<ET, nfixed>>
{
public:
  friend StaticInitializer<ExactLatticePoint<ET,nfixed>>;
  using LatticePointTag         = std::true_type;
  using ScalarProductStorageType = ET;
  using Container = typename std::conditional<nfixed >= 0,
        std::array<ET, nfixed >=0 ? nfixed:0>,  // if nfixed >= 0
        std::vector<ET>  >::type;               // if nfixed <0
        // Note : The nfixed >=0 ? nfixed:0 is always nfixed;
        // The ?: expression is only needed to silence compiler errors/warnings.
  
  /*
  using ApproxContainer = typename std::conditional<nfixed >= 0,
                          std::bitset<nfixed >=0 ? nfixed:0>,  // if nfixed >= 0
                          boost::dynamic_bitset<>  >::type;                   // if nfixed <  0
  */
  using ApproxContainer = boost::dynamic_bitset<>;

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
  explicit ExactLatticePoint() : data(static_cast<unsigned int>(get_dim()))
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((StaticInitializer<ExactLatticePoint<ET,nfixed>>::is_initialized));
  }

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
  : data(std::move(plain_point.data)), bitapprox_data()
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

  void sanitize() {
      norm2 = compute_sc_product(*this, *this);
      compute_approximation(*this);
  }
  void sanitize( ScalarProductStorageType const & new_norm2 ) 
  { 
      norm2 = new_norm2; 
      compute_approximation(*this);
  }

  ET get_norm2() const { return norm2; }
  
  std::ostream& write_lp_to_stream (std::ostream &os, bool const include_norm2=true, bool const include_approx=true) const;

  inline ApproxScProduct compute_sc_product_bitapprox(ExactLatticePoint const & another);

  inline ET compute_sc_product(ExactLatticePoint const &lp1, ExactLatticePoint const &lp2)
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
  
  
  inline void compute_approximation(ExactLatticePoint &point)
  {
      
      uint_fast16_t dim = get_dim();
      bitapprox_data.resize(dim);
      for(uint_fast16_t i=0;i<dim;++i)
      {
        bitapprox_data[i] = (point[i]>=0) ? 1 : 0;
      }
      
  }
  

private:
  static MaybeFixed<nfixed> dim;  // note that for nfixed != -1, this variable is actually unused.

  ApproxContainer bitapprox_data;

  Container data;
  ET norm2;
};


//approximate scalar product
template <class ET, int nfixed>
inline ApproxScProduct ExactLatticePoint<ET, nfixed>::compute_sc_product_bitapprox(ExactLatticePoint const &another)
  {
    ApproxScProduct result(0);
    auto const dim = this.get_dim();
    result = dim - (this.bitapprox_data ^ another.bitapprox_data).count();
    return result;
}

template <class ET, int nfixed>
inline std::ostream& ExactLatticePoint<ET, nfixed>::write_lp_to_stream(std::ostream &os, bool const include_norm2, bool const include_approx) const
{
// Note: include_approx is ignored, because classes that have an approximation overload this anyway.
  DEBUG_TRACEGENERIC("Using generic writer (absolute) for " << LatP::class_name() )
  //std::cout << "call from ExactLatticePoint" << std::endl;
  auto const dim = get_dim();
  //os << "dim = " << dim << std::endl;
  os << "[ "; // makes spaces symmetric
  for (uint_fast16_t i =0; i<dim; ++i)
  {
    os << data[i] << " ";
  }
  os << "]";
  
  if(include_norm2)
  {
    os <<", norm2= " << norm2 << " ";
  }
  
  if (include_approx)
  {
    os << " approx = [ ";
    for (uint_fast16_t i =0; i<dim; ++i)
    {
      os << bitapprox_data[i] << " ";
    }
    os << "]";
  }
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

  StaticInitializer(MaybeFixed<nfixed> const new_dim)
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
};

}  // end of namespace

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM

#endif
