#ifndef UNIFORM_SAMPLER_H
#define UNIFORM_SAMPLER_H

#include "DefaultIncludes.h"
#include "LatticeBases.h"
#include "Sampler.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"

/* Sparce uniform sampler
 * Chooses sparcity-many '1'-coefficients at random from {0,..,dim-1}
 * and sparcity-many '-1'-coefficients from {0,..., dim-1}
 *  */

namespace GaussSieve
{
// forward declarations:
template <class SieveTraits, bool MT> class Sieve;
template <class SieveTraits, bool MT, class Engine, class Sseq> class UniformSampler;

template <class SieveTraits, bool MT, class Engine, class Sseq>
class UniformSampler final : public Sampler<SieveTraits, MT, Engine, Sseq>
{
public:
  using DimensionType = typename SieveTraits::DimensionType;
  using EntryType     = typename SieveTraits::EntryType;
  using RetType       = typename SieveTraits::GaussSampler_ReturnType;
  // clang-format off
  explicit UniformSampler(Sseq &seq, unsigned int inp_sparcity)
      : Sampler<SieveTraits,MT,Engine,Sseq>(seq),
        // dim, lattice_rank un-/ or default-initialized.
        sparcity(inp_sparcity),
        initialized(false),
        static_init_rettype(nullptr),
        static_init_plainpoint(nullptr)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing Uniform sampler.");
  }
  // clang-format on

  virtual SamplerType sampler_type() const override { return SamplerType::uniform_sampler; }

  // TODO: destructor
  virtual ~UniformSampler()
  {
    if (initialized)
    {
      delete static_init_plainpoint;
      delete static_init_rettype;
    }
  }

  virtual inline RetType sample(int const thread = 0) override;

private:
  inline virtual void custom_init(SieveLatticeBasis<SieveTraits, MT> const &input_basis) override;

  DimensionType dim;
  uint_fast16_t lattice_rank;

  unsigned int sparcity;
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
