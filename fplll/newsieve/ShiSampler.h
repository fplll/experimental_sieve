/**
 An adaption of the sampler of Shi's old implemenation to our classes and conventions.

 Has a small bug wrt. usage of what GSO returns, notably the scaling by the length of the
 GS vectors is probably not done as intented. This means that the Gaussian is not spherical.
 This was present in the old version and we try to match it for comparisons.

 Use EllipticSampler if you want a spherical Gaussian.
*/

#ifndef SHI_SAMPLER_H  // equivalent to Shi's old sampler, using new framework.
#define SHI_SAMPLER_H

#include "DefaultIncludes.h"
#include "Sampler.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"
#include <random>
#include <vector>
#include "LatticeBases.h"


namespace GaussSieve
{

// forward declaration
template <class SieveTraits, bool MT> class Sieve;


template <class SieveTraits, bool MT, class Engine, class Sseq> class ShiSampler;

template <class SieveTraits, bool MT, class Engine, class Sseq>
class ShiSampler final : public Sampler<SieveTraits, MT, Engine, Sseq>
{
public:
  using DimensionType = typename SieveTraits::DimensionType;
  using EntryType     = typename SieveTraits::EntryType;
  using RetType       = typename SieveTraits::GaussSampler_ReturnType;

  // Note: Sampler::sieveptr is only initialized during Sampler::init.
  // Consequently, some member fields will only be set during custom_init.
  explicit ShiSampler(Sseq &seq, double const _cutoff = 2.0)
      : Sampler<SieveTraits, MT, Engine, Sseq>(seq), cutoff(_cutoff), initialized(false),
      static_init_rettype(nullptr), static_init_plainpoint(nullptr)
      {
        DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing ShiSampler.")
      };
  virtual SamplerType sampler_type() const override { return SamplerType::shi_sampler; };
  virtual ~ShiSampler()
  {
    if(initialized)
    {
      delete static_init_plainpoint;
      delete static_init_rettype;
//      RetType::class_uninit();
//      SieveTraits::PlainPoint::class_uninit();
    }
  };
  virtual inline RetType sample(int thread = 0) override;

private:
  inline virtual void custom_init(SieveLatticeBasis<SieveTraits,MT> const & input_basis) override;
//  fplll::ZZ_mat<typename ET::underlying_data_type> current_basis;
//  std::vector<MyLatticePoint<ET, nfixed>> helper_current_basis;  // TODO: Use different type
  std::vector<std::vector<double>> mu_matrix; // copied from basis.
//  fplll::Matrix<fplll::FP_NR<double>> mu;
  // stores standard dev. for each dimension, already squared and divided by pi.
  std::vector<double> s2pi;
  std::vector<double> maxdeviations;  // stores s*cutoff for each dimension.
  DimensionType dim;
  uint_fast16_t lattice_rank;
  double cutoff;
  bool initialized;

protected:
  // bring into scope from parent
  using Sampler<SieveTraits, MT, Engine, Sseq>::sieveptr;
  using Sampler<SieveTraits, MT, Engine, Sseq>::engine;
  std::vector<typename SieveTraits::PlainPoint> basis;

  StaticInitializer<RetType> *static_init_rettype;
  StaticInitializer<typename SieveTraits::PlainPoint> *static_init_plainpoint;

  // TODO: Store basis in a different way.
  //    vector<MyLatticePoint> basis;
};

}  // end namespace

#endif
