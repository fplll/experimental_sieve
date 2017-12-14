// clang-format status: OK

/**
 Gaussian GPV sampler from Gentry, Peikert, Vaikuntanathan 'How to Use a Short Basis:
 Trapdoors for Hard Lattices and New Cryptographic Constructions'.
*/

#ifndef GPV_SAMPLER_H
#define GPV_SAMPLER_H

#include "DefaultIncludes.h"
#include "Sampler.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"

#include "LatticeBases.h"

namespace GaussSieve
{

// forward declaration
template <class SieveTraits, bool MT> class Sieve;

template <class SieveTraits, bool MT, class Engine, class Sseq> class ShiSampler;

template <class SieveTraits, bool MT, class Engine, class Sseq>
class GPVSampler final : public Sampler<SieveTraits, MT, Engine, Sseq>
{
public:
  using DimensionType = typename SieveTraits::DimensionType;
  using LengthType    = typename SieveTraits::LengthType;
  using RetType       = typename SieveTraits::GaussSampler_ReturnType;

  // Note: Sampler::sieveptr is only initialized during Sampler::init.
  // Consequently, some member fields will only be set during custom_init.
  // The constructor takes a seed of type Sseq, e.g. std::seed_seq seed{3,2,1},
  // and a cut_off parameter for maxdeviations, namely
  // maxdeviations = sqrt(s2pi)*_cutoff
  // clang-format off
  explicit GPVSampler(Sseq &seq, double const _cutoff = 2.0)
      : Sampler<SieveTraits, MT, Engine, Sseq>(seq),
        cutoff(_cutoff),
        initialized(false),
        static_init_rettype(nullptr),
        static_init_plainpoint(nullptr)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing GPVSampler.")
  }
  // clang-format on

  virtual SamplerType sampler_type() const override { return SamplerType::GPV_sampler; };
  virtual ~GPVSampler()
  {
    if (initialized)
    {
      delete static_init_plainpoint;
      delete static_init_rettype;
    }
  };
  virtual inline RetType sample(int const thread = 0) override;

private:
  inline virtual void custom_init(SieveLatticeBasis<SieveTraits, MT> const &input_basis) override;

  std::vector<std::vector<double>> mu_matrix;  // copied from basis.
  std::vector<double> s2pi;
  std::vector<double> maxdeviations;  // stores s*cutoff for each dimension.
  DimensionType dim;                  // ambient dimension
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
};

}  // end namespace

#endif
