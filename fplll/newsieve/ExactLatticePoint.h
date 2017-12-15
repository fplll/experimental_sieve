/**
  Very simple lattice point class that stores precomputed values for its length.
  Otherwise, just a wrapper around a vector, really.
  Dimension is static (i.e. the same for all objects)
  It is a lattice point, so has to follow the interface specified in
  LatticePointConcept.h
*/

// clang-format checked and adjusted manually -- Gotti

#ifndef EXACT_LATTICE_POINT_H
#define EXACT_LATTICE_POINT_H

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"
#include "LatticePointConcept.h"
#include "PlainLatticePoint.h"  // for conversions
#include "SieveUtility.h"

// local defines, undef at the end of file, used to enable functions only for (non-)fixed dimension.
// clang-format off
#define FOR_FIXED_DIM    template<int nfixed_copy = nfixed, TEMPL_RESTRICT_DECL(nfixed_copy >= 0)>
#define FOR_VARIABLE_DIM template<int nfixed_copy = nfixed, TEMPL_RESTRICT_DECL(nfixed_copy == -1)>
// clang-format on

// to make some functions C++14-constexpr (if supported)
// depends on DEBUG_SYMBOL
// TODO: remove this macro and write it out everywhere, disambiguate
// C++11 and C++14 constexpr
// clang-format off
#ifdef DEBUG_SIEVE_LP_INIT
  #define CONSTEXPR_IN_NON_DEBUG_LP_INIT CPP14CONSTEXPR
#else
  #define CONSTEXPR_IN_NON_DEBUG_LP_INIT
#endif
// clang-format on

namespace GaussSieve
{

// forward declaration
template <class ET, int nfixed> class ExactLatticePoint;

/**
  lattice point traits for ExactLatticePoint:
*/
template <class ET, int nfixed> struct LatticePointTraits<ExactLatticePoint<ET, nfixed>>
{
public:
  // set appropriate traits:
  // clang-format off
  using Trait_ScalarProductStorageType = ET;
  using Trait_CoordinateType           = ET;
  using Trait_CheapNorm2               = std::true_type;  // promises that Norm2 is precomputed
  using Trait_CheapNegate              = std::true_type;  // no update of Norm2 if we flip the sign

  // These traits just mean that we have an RW-operator[] that outputs absolute coordinate
  // (i.e. with respects to the ambient space) and these determine the point uniquely.
  using Trait_AbsoluteCoos             = std::true_type;
  using Trait_InternalRepIsAbsolute    = std::true_type;
  using Trait_InternalRepByCoos        = std::true_type;
  using Trait_InternalRepLinear        = std::true_type;
  using Trait_InternalRep_RW           = std::true_type;
  // clang-format on
};

/**
  ExactLatticePoint<ET,nfixed> is a lattice point class that wraps around a vector / array of ET's.
  The parameter ET determines the type of the entries.
  The parameter nfixed determines whether the dimension is fixed at compile-time or run-time.
  If nfixed >= 0, it is fixed to that value at compile time.
  If nfixed == -1,it is fixed at runtime.
  Even in the latter case, the dimension of all ExactLatticePoint<ET,-1> objects (for given ET) is
  the same at any given point in time, because dimension is static
  (it is set by the static initializer)
  ExactlatticePoints store and maintain a precomputed norm^2.
  For most of the functionality, we use the default implementation from the parent class
  GeneralLatticePoint<ET,nfixed> (cf. LatticePointGeneric.h)
*/

template <class ET, int nfixed>
class ExactLatticePoint final : public GeneralLatticePoint<ExactLatticePoint<ET, nfixed>>
{
  friend StaticInitializer<ExactLatticePoint<ET, nfixed>>;

public:
  using LatticePointTag = std::true_type;

private:
  // Container type used to store the actual point
  using Container =
      mystd::conditional_t<(nfixed >= 0),                               // depends on nfixed:
                           std::array<ET, (nfixed >= 0 ? nfixed : 0)>,  // if nfixed >= 0
                           std::vector<ET>>;                            // if nfixed <0
  // Note : The nfixed >=0 ? nfixed : 0 is always nfixed, of course.
  // The ?: expression is only needed to silence compiler errors/warnings.

public:
  // get dimension
  FOR_FIXED_DIM
  static constexpr MaybeFixed<nfixed> get_dim()
  {
    static_assert(nfixed_copy == nfixed, "");  // nfixed_copy defined in FOR_FIXED_DIM
    return MaybeFixed<nfixed>(nfixed);
  }

  FOR_VARIABLE_DIM
  static constexpr MaybeFixed<-1> const &get_dim()
  {
    static_assert(nfixed == -1, "");
    return dim;
  }

  // The macro confuses clang-format completely...
  // clang-format off
  FOR_FIXED_DIM  // This introduces template params, so we do not =default
  constexpr explicit ExactLatticePoint() {}
  // clang-format on

  // Gotti: If this gives problems on clang, remove the
  // CONSTEXPR_IN_NON_DEBUG_LP_INIT
  FOR_VARIABLE_DIM
  CONSTEXPR_IN_NON_DEBUG_LP_INIT explicit ExactLatticePoint()  // make noexcept?
      : data(static_cast<unsigned int>(get_dim()))
  {
#ifdef DEBUG_SIEVE_LP_INIT
    // For some reason, clang-format messes up the indentations of its own clang-format commands...
    // clang-format off
    // double (( )) because assert is a macro and it mis-parses the ","
    assert((StaticInitializer< ExactLatticePoint<ET,nfixed> >::is_initialized));
// clang-format on
#endif
  }

  // TODO: Remove these constructors?
  FOR_FIXED_DIM
  constexpr explicit ExactLatticePoint(MaybeFixed<nfixed>) : ExactLatticePoint() {}

  FOR_VARIABLE_DIM
  CONSTEXPR_IN_NON_DEBUG_LP_INIT explicit ExactLatticePoint(MaybeFixed<nfixed> dim)
      : ExactLatticePoint()
  {
#ifdef DEBUG_SIEVE_LP_MATCHDIM
    assert(dim == get_dim());
#endif
  }

  // TODO: Debug output and validation.

  explicit ExactLatticePoint(PlainLatticePoint<ET, nfixed> &&plain_point) noexcept
      : data(std::move(plain_point.data))
  {
    sanitize();  // compute norm2
  }

  // lattice points are moveable but not copyable (by design, to catch unwanted copying)
  // (note that these declarations just make that explicit, if we did not write these,
  // the auto-generated copy constructor would be deleted, because the parent's is.)
  // clang-format off
  ExactLatticePoint(ExactLatticePoint const &)            = delete;
  ExactLatticePoint(ExactLatticePoint &&)       	        = default;
  ExactLatticePoint &operator=(ExactLatticePoint const &) = delete;
  ExactLatticePoint &operator=(ExactLatticePoint &&)      = default;
  // clang-format on

  // behaves just like a normal container and can be used as such
  // clang-format off
  ET       &operator[](uint_fast16_t idx)       { return data[idx]; }
  ET const &operator[](uint_fast16_t idx) const { return data[idx]; }
  // clang-format on
  static std::string class_name() { return "Exact Lattice Point"; }

  // brings the internal data into the correct state, i.e. we recompute the norm.
  // clang-format off
  void sanitize()
  {
    sanitize(compute_sc_product(*this, *this));
  }
  // clang-format on

  // version of the above where we already know norm2
  void sanitize(ET const &new_norm2) { norm2 = new_norm2; }

  ET get_norm2() const { return norm2; }

  // Compute scalar product with another point.
  // We implicitly assert the other point has a "compatible" type
  template <class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline ET do_compute_sc_product(LatP2 const &lp2) const
  {
    // Compute the sum of (*this)[i]  * lp2[i] over i.
    // Naively, we would start with res = 0 and increment
    // res+= (*this)[i] * lp2[i] in a loop over i.

    // Instead, we split the computation of res into a sum of 4 sub-terms, corresponding to the
    // contribution from a subset of the coordinates (depending on i mod 4).
    // The reason is that inside the for-loop,
    // the res1+= ,... res4+= statements can be computed in parallel by SIMD.
    // with only a single res, this is not the case, because each +=operation has to wait
    // for the previous one to finish.
    // whether this helps or not depends on compiler / architecture.
    // Note that to get the "naive" implementation, just comment out this function altogether,
    // the default inherited from the parent will do exactly that.
    ET res1 = 0;
    ET res2 = 0;
    ET res3 = 0;
    ET res4 = 0;

    uint_fast16_t dim = get_dim();
    for (uint_fast16_t i = 0; i < (dim / 4) * 4; i += 4)  // perform 4 += in one go
    {
      res1 += (*this)[i + 0] * lp2.get_absolute_coo(i + 0);
      res2 += (*this)[i + 1] * lp2.get_absolute_coo(i + 1);
      res3 += (*this)[i + 2] * lp2.get_absolute_coo(i + 2);
      res4 += (*this)[i + 3] * lp2.get_absolute_coo(i + 3);
    }
    for (uint_fast16_t i = (dim / 4) * 4; i < dim;
         ++i)  // rest of the loop if dim is not divisible by 4
    {
      res1 += (*this)[i] * lp2.get_absolute_coo(i);
    }
    res1 += res2;  // add up the sub-sums.
    res1 += res3;
    res1 += res4;
    return res1;
  }

private:
  static MaybeFixed<nfixed> dim;  // note that for nfixed != -1, this variable is actually unused.
  Container data;                 // either a vector or an array that stores the actual point.
  ET norm2;                       // stores the precomputed norm2
};

// define (and compile-time initialize) static data of template classes:
template <class ET, int nfixed>
MaybeFixed<nfixed> ExactLatticePoint<ET, nfixed>::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);

/************************************************************
    Static Initializer: This class sets the dimension of the point
************************************************************/

// clang-format off
template <class ET, int nfixed>
class StaticInitializer< ExactLatticePoint<ET,nfixed> > final
    : public DefaultStaticInitializer<ExactLatticePoint<ET, nfixed>>
{
  using Parent = DefaultStaticInitializer< ExactLatticePoint<ET,nfixed> >;

  // clang-format on
public:
  template <class T, TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const &initializer) : StaticInitializer(initializer.dim)
  {
  }

  StaticInitializer(MaybeFixed<nfixed> const new_dim)
  {
    assert(Parent::user_count > 0);
    if (Parent::user_count > 1)
    {
      assert((new_dim == ExactLatticePoint<ET, nfixed>::dim));
      // TODO: Throw exception!
    }
    else
    {
      ExactLatticePoint<ET, nfixed>::dim = new_dim;
    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing ExactLatticePoint with nfixed = "
                                      << nfixed << " REALDIM = " << new_dim << " Counter is"
                                      << Parent::user_count)
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing ExactLatticePoint with nfixed = "
                                      << nfixed << " Counter is " << Parent::user_count)
  }
};  // end of static initializer

}  // end of namespace

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM
#undef CONSTEXPR_IN_NON_DEBUG_LP_INIT

#endif  // include guard
