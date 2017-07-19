/**
 An adaption of the sampler of Shi's old implemenation to our classes and conventions.

 Has a small bug wrt. usage of what GSO returns, notably the scaling by the length of the
 GS vectors is probably not done as intented. This means that the Gaussian is not spherical.
 This was present in the old version and we try to match it for comparisons.

 Use EllipticSampler if you want a spherical Gaussian.
*/

#ifndef SHI_SAMPLER_H  // equivalent to Shi's old sampler, using new framework.
#define SHI_SAMPLER_H

#include "Sampler.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include <random>
#include <vector>

template <class ET, bool MT, class Engine, class Sseq, int nfixed> class ShiSampler;

template <class ET, bool MT, class Engine, class Sseq, int nfixed>
class ShiSampler final : public Sampler<ET, MT, Engine, Sseq, nfixed>
{
public:
  explicit ShiSampler(Sseq &seq, double const _cutoff = 2.0)
      : Sampler<ET, MT, Engine, Sseq, nfixed>(seq), dim(nfixed < 0 ? 0 : nfixed), cutoff(_cutoff){};
  virtual SamplerType sampler_type() const override { return SamplerType::shi_sampler; };
  virtual ~ShiSampler();
  virtual typename GaussSieve::GaussSampler_ReturnType<ET, MT, nfixed>
  sample(int thread = 0) override;

private:
  inline virtual void custom_init() override;
  fplll::ZZ_mat<typename ET::underlying_data_type> current_basis;
  std::vector<MyLatticePoint<ET, nfixed>> helper_current_basis;  // TODO: Use different type
  fplll::Matrix<fplll::FP_NR<double>> mu;

  // stores standard dev. for each dimension, already squared and divided by pi.
  std::vector<double> s2pi;
  std::vector<double> maxdeviations;  // stores s*cutoff for each dimension.
  Dimension<nfixed> dim;
  unsigned int rank;
  double cutoff;

protected:
  using Sampler<ET, MT, Engine, Sseq, nfixed>::sieveptr;
  using Sampler<ET, MT, Engine, Sseq, nfixed>::engine;
  // TODO: Store basis in a different way.
  //    vector<MyLatticePoint> basis;
};

#endif
