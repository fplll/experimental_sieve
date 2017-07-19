// clang-format off

#ifndef SHI_SAMPLER_H //equivalent to Shi's old sampler, using new framework.
#define SHI_SAMPLER_H

template<class ET, bool MT, class Engine, class Sseq, int nfixed> class ShiSampler;

#include <random>
#include "Sampler.h"
#include "Typedefs.h"
#include "SieveUtility.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include <vector>
//only defined for nfixed==-1 for now.

template<class ET, bool MT, class Engine, class Sseq,int nfixed>
class ShiSampler: public Sampler<ET,MT, Engine, Sseq,nfixed>
{
    public:
    ShiSampler(Sseq & seq, double const _cutoff = 2.0): Sampler<ET,MT,Engine,Sseq,nfixed>(seq), dim(nfixed<0?0:nfixed), cutoff(_cutoff) {};
    virtual SamplerType  sampler_type() const override                          {return SamplerType::shi_sampler;};
    virtual ~ShiSampler();
    virtual typename GaussSieve::GaussSampler_ReturnType<ET,MT,nfixed> sample(int thread=0) override;
    private:
    virtual void custom_init() override;
    fplll::ZZ_mat<typename ET::underlying_data_type> current_basis;
    std::vector<MyLatticePoint <ET,nfixed> > helper_current_basis;
    fplll::Matrix<fplll::FP_NR<double > > mu;
    std::vector<double> s2pi; //stores standard dev. for each dimension, already squared and multiplied by pi.
    std::vector<double> maxdeviations; //stores s*cutoff for each dimension.
    Dimension<nfixed> dim;
    unsigned int rank;
    double cutoff;
    protected:
    using Sampler<ET,MT,Engine,Sseq,nfixed>::sieveptr;
    using Sampler<ET,MT,Engine,Sseq,nfixed>::engine;
//    vector<MyLatticePoint> basis;
};



#endif

//clang-format on
