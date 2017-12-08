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

/* Flexible GPV sampler
 TODO: describe
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
    using EntryType     = typename SieveTraits::EntryType;
    using RetType       = typename SieveTraits::GaussSampler_ReturnType;
    
    explicit GPVSampler(Sseq &seq, double const _cutoff = 2.0)
    : Sampler<SieveTraits, MT, Engine, Sseq>(seq), cutoff(_cutoff), initialized(false),
    static_init_rettype(nullptr), static_init_plainpoint(nullptr)
    {
      DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing GPVSampler.")
    };
    virtual SamplerType sampler_type() const override { return SamplerType::GPV_sampler; };
    virtual GPVSampler()
    {
      if(initialized)
      {
        delete static_init_plainpoint;
        delete static_init_rettype;
        //      RetType::class_uninit();
        //      SieveTraits::PlainPoint::class_uninit();
      }
    };
    virtual inline RetType sample(int const thread = 0) override;
    
  private:
    inline virtual void custom_init(SieveLatticeBasis<SieveTraits,MT> const & input_basis) override;
    
    
    std::vector<std::vector<double>> Q_matrix; // copied from basis.
    std::vector<std::vector<double>> R_matrix;
    //  fplll::Matrix<fplll::FP_NR<double>> mu;
    // stores standard dev. for each dimension, already squared and divided by pi.
    std::vector<double> s2pi;
    std::vector<double> maxdeviations;
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

    
  };
  
} //namespace

#endif