#ifndef GSO_LATTICE_POINT_H
#define GSO_LATTICE_POINT_H

#include "DefaultIncludes.h"
#include "LatticePointConcept.h"
#include "LatticePointGeneric.h"

namespace GaussSieve
{

template <class AbsoluteEntries, class RelativeCoos, class Basis, int nfixed>
struct LatticePointTraits<GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>>
{
public:
  using Trait_ScalarProductStorageType      = AbsoluteEntries;
  using Trait_ScalarProductStorageType_Full = AbsoluteEntries;
  using Trait_ExposesCoos                   = std::false_type;
  //using Trait_CoordinateType                = void;
  //using Trait_Coos_RW                       = std::false_type;
  using Trait_AbsoluteCooType               = AbsoluteEntries;
  using Trait_RepCooType                    = AbsoluteEntries;
  using Trait_ExposesInternalRep            = std::true_type;
  using Trait_InternalRepLinear             = std::true_type;
  using Trait_InternalRep_RW                = std::true_type;
  using Trait_InternalRepByCoos             = std::false_type;
  using Trait_InternalRepIsAbsolute         = std::false_type;
  //using Trait_AbsoluteCoos                  = std::false_type;
  using Trait_CheapNorm2                    = std::true_type;
  using Trait_CheapNegate                   = std::true_type;
  using Trait_BitApprox                     = std::false_type;
  using Trait_Leveled                       = std::false_type;
};


// local defines, undef at the end of file, used to enable functions only for (non-)fixed dimension.
// clang-format off
#define FOR_FIXED_DIM    template<int nfixed_copy = nfixed, TEMPL_RESTRICT_DECL(nfixed_copy >= 0)>
#define FOR_VARIABLE_DIM template<int nfixed_copy = nfixed, TEMPL_RESTRICT_DECL(nfixed_copy == -1)>
// clang-format on


template <class AbsoluteEntries, class RelativeCoos, class Basis, int nfixed>
    final : public GeneralLatticePoint<GSOLatticePoint< AbsoluteEntries, RelativeCoos, Basis, nfixed>>
class GSOLatticePoint
{
  // constructors are defaulted.
  // Note that copy constructor is automatically deleted, because the parent's is.

  friend StaticInitializer<GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>>;

public:

  using LatticePointTag = std::true_type;

private:

  // Container type used to store the coefficients of the point wrt an orthonormal basis
  using AbsoluteContainer =
      mystd::conditional_t<(nfixed >= 0),                               // depends on nfixed:
                           std::array<AbsoluteEntries, (nfixed >= 0 ? nfixed : 0)>,  // if nfixed >= 0
                           std::vector<AbsoluteEntries>>;                            // if nfixed <0

  using CoefficientContainer =
      mystd::conditional_t<(nfixed >=0),
              std::array<RelativeCoos, (nfixed >= 0 ? nfixed : 0)>,
              std::vector<RelativeCoos>>;

  static Basis const *basisptr;
  static MaybeFixed<nfixed> dim;

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

  /*
    Constructors:
  */

  // clang-format off
  FOR_FIXED_DIM  // This introduces template params, so we do not =default
  constexpr explicit GSOLatticePoint() {}
  // clang-format on

  FOR_FIXED_DIM
  constexpr explicit GSOLatticePoint(MaybeFixed<nfixed> dim) {}

  FOR_VARIABLE_DIM
  explicit GSOLatticePoint(MaybeFixed<nfixed> dim)  // TODO: constexpr
      : GSOLatticePoint()
  {
#ifdef DEBUG_SIEVE_LP_MATCHDIM
    assert(dim == get_dim());
#endif
  }

  FOR_VARIABLE_DIM
  explicit GSOLatticePoint()  // make noexcept?
      : onb_coos(get_dim()), relative_coos(get_dim())
  {
    // TODO : assert initialization
  }

  static std::string class_name() { return "Lattice point with coeffs wrt. supplied basis"; }

  AbsoluteEntries get_norm2() const { return norm2; }

  void sanitize() &
  {
    sanitize(compute_sc_product(*this, *this));
  }
  // clang-format on

  // version of the above where we already know norm2. Note that sanitize() just forwards to this.
  void sanitize(AbsoluteEntries const &new_norm2)
  {
    norm2 = new_norm2;
  }

  ET get_norm2() const { return norm2; }


private:
  AbsoluteContainer onb_coos;  // coordinates wrt to an orthonormal basis
  CoefficientContainer relative_coos;  // coordinates relative to the basis indicated by the static
                                       // pointer
  AbsoluteEntries norm2;  // precomputed norm^2
};

/**
  static member intitalization, compile-time (or startup):
*/
template<class AbsoluteEntries, class RelativeCoos, class Basis, int nfixed>
Basis const * GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>::basisptr = nullptr;
template<class AbsoluteEntries, class RelativeCoos, class Basis, int nfixed>
MaybeFixed<nfixed> GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);

/**
  Run-time initalization, via RAII wrapper:
*/

template <class ET, int nfixed>
class StaticInitializer<GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>> final
    : public DefaultStaticInitializer<GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>>
{
  using Parent = DefaultStaticInitializer<GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>>;

public:
  template <class T, TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const &initializer) : StaticInitializer(initializer.dim)
  {
  }

  StaticInitializer(MaybeFixed<nfixed> const new_dim, )
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
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing GSOLatticePoint with nfixed = "
                                      << nfixed << " REALDIM = " << new_dim << " Counter is"
                                      << Parent::user_count)
  }

  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing GSOLatticePoint with nfixed = "
                                      << nfixed << " Counter is " << Parent::user_count)
  }
};  // end of static initializer



}  // end namespace GaussSieve

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM

#endif  // include guard

}  // end of namespace

