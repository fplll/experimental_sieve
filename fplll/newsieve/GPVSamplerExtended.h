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
#include <random>
#include <vector>
#include "LatticeBases.h"

/* Flexible GPV sampler centered at 0 with st. dev. parameter s = ||b*_1|| * sqrt(n);
    GPV sampler does the following for B = Q * R (R - upper-triag), r_{i,i} =||b*_1||^2

    for i from dim to 1
        shift[i] = - \sum_{j>i} x_j r_{i,j}
        x_i = SampleZ (s / r_{i,i}), c_i / r_{i,i}
        b = b + x_i & b_i

    return b

 *  */


namespace GaussSieve
{

template <class SieveTraits, bool MT> class Sieve;


template <class SieveTraits, bool MT, class Engine, class Sseq> class GPVSampler; //declaration

template <class SieveTraits, bool MT, class Engine, class Sseq>
class GPVSampler final : public Sampler<SieveTraits, MT, Engine, Sseq>
{
public:
  using DimensionType = typename SieveTraits::DimensionType;
  using LengthType     = typename SieveTraits::LengthType;
  using RetType       = typename SieveTraits::GaussSampler_ReturnType;

  explicit GPVSampler(Sseq &seq, uint_fast16_t babai, double const _cutoff = 2.0)
      : Sampler<SieveTraits,MT,Engine,Sseq>(seq),
        start_babai(babai), cutoff(_cutoff), initialized(false),
        static_init_rettype(nullptr), static_init_plainpoint(nullptr)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing GPVSampler.")
  };
  virtual SamplerType sampler_type() const override { return SamplerType::GPV_sampler; };
  virtual ~GPVSampler()
  {
    if(initialized)
    {
      delete static_init_plainpoint;
      delete static_init_rettype;
    }
  };
  virtual inline RetType sample(int const thread = 0) override;

private:
  inline virtual void custom_init(SieveLatticeBasis<SieveTraits,MT> const & input_basis) override;


  //std::vector<std::vector<double>> q_matrix;
  //std::vector<std::vector<double>> r_matrix; //B = Q*R, r_ij = mu_ij * | b*_j|^2

  std::vector<std::vector<double>> mu_matrix; // mu_i,j = r_i,j / ||b*_j||^2.. lower triangular matrix


  std::vector<double> s2pi; // stores standard dev. for each dimension, already squared and divided by pi.
  std::vector<double> maxdeviations; // [s2pi * cutoff] - for rejection sampling
  DimensionType dim;

  uint_fast16_t lattice_rank;
  uint_fast16_t start_babai; // start_babai < lattice_rank; use Babai on { b_{start_babai},...b_1}

  double cutoff;
  bool initialized;

protected:

  using Sampler<SieveTraits, MT, Engine, Sseq>::sieveptr;
  using Sampler<SieveTraits, MT, Engine, Sseq>::engine;
  std::vector<typename SieveTraits::PlainPoint> basis;

  StaticInitializer<RetType> *static_init_rettype;
  StaticInitializer<typename SieveTraits::PlainPoint> *static_init_plainpoint;


};

} //namespace

#endif
