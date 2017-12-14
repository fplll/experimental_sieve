// clang-format status: OK

#ifndef GPV_SAMPLER_EXTENDED_H
#define GPV_SAMPLER_EXTENDED_H

#include "DefaultIncludes.h"
#include "LatticeBases.h"
#include "Sampler.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"

/* GPVExtended sampler centered at 0 with st. dev. parameter s = ||b*_1|| * sqrt(n);
   Same as GPSSampler but with additional member uint_fast16_t start_babai;
   The sampling routine of GPVExtended runs GPV on dimensions {n,..., n-start_babai},
   (creates a CVP insance as the result) and on {start_babai, ..., 1} runs Babai's
   algorithm.

   TODO: Extend to smth better (qualitatively) than Babai, e.g. enumeration, Random Sampling...

 *  */

namespace GaussSieve
{

template <class SieveTraits, bool MT> class Sieve;

// forward declaration:
template <class SieveTraits, bool MT, class Engine, class Sseq> class GPVSamplerExtended;

template <class SieveTraits, bool MT, class Engine, class Sseq>
class GPVSamplerExtended final : public Sampler<SieveTraits, MT, Engine, Sseq>
{
public:
  using DimensionType = typename SieveTraits::DimensionType;
  using LengthType    = typename SieveTraits::LengthType;
  using RetType       = typename SieveTraits::GaussSampler_ReturnType;

  // clang-format off
  explicit GPVSamplerExtended(Sseq &seq, uint_fast16_t babai, double const _cutoff = 2.0)
      : Sampler<SieveTraits,MT,Engine,Sseq>(seq),
        start_babai(babai),
        cutoff(_cutoff),
        initialized(false),
        static_init_rettype(nullptr),
        static_init_plainpoint(nullptr)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing GPVSampler.")
  }
  // clang-format on

  virtual SamplerType sampler_type() const override { return SamplerType::GPVExtended_sampler; };
  virtual ~GPVSamplerExtended()
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

  // mu_i,j = r_i,j / ||b*_j||^2.. lower triangular matrix
  std::vector<std::vector<double>> mu_matrix;

  std::vector<double> s2pi;           // st. dev. per dimension, squared and divided by pi.
  std::vector<double> maxdeviations;  // [s2pi * cutoff] - for rejection sampling
  DimensionType dim;

  uint_fast16_t lattice_rank;
  uint_fast16_t start_babai;  // start_babai < lattice_rank; use Babai on { b_{start_babai},...b_1}

  double cutoff;
  bool initialized;

protected:
  using Sampler<SieveTraits, MT, Engine, Sseq>::sieveptr;
  using Sampler<SieveTraits, MT, Engine, Sseq>::engine;
  std::vector<typename SieveTraits::PlainPoint> basis;

  StaticInitializer<RetType> *static_init_rettype;
  StaticInitializer<typename SieveTraits::PlainPoint> *static_init_plainpoint;
};

}  // namespace GaussSieve

#endif
