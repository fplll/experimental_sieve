#ifndef SHI_SAMPLER_H //equivalent to Shi's old sampler, using new framework.
#define SHI_SAMPLER_H

template<class ET, bool MT, class Engine, class Sseq, int nfixed> class ShiSampler;

#include <random>
#include "Sampler.h"
#include "Typedefs.h"
#include "Utility.h"


//only defined for nfixed==-1 for now.

template<class ET, bool MT, class Engine, class Sseq,int nfixed>
class ShiSampler: public Sampler<ET,MT, Engine, Sseq,nfixed>
{
    public:
    ShiSampler(Sseq & seq, double const _cutoff = 2.0)
        :   Sampler<ET,MT,Engine,Sseq>(seq), dim(0), cutoff(_cutoff) {};
    virtual SamplerType  sampler_type() const override                          {return SamplerType::shi_sampler;};
    virtual ~ShiSampler();
    virtual typename GaussSieve::GaussSampler_ReturnType<ET,MT,nfixed> sample(int thread=0) override;
    private:
    virtual void custom_init() override;
    ZZ_mat<typename ET::underlying_data_type> current_basis;
    vector<MyLatticePoint <ET,nfixed> > helper_current_basis;
    Matrix<FP_NR<double > > mu;
    vector<double> s2pi; //stores standard dev. for each dimension, already squared and multiplied by pi.
    vector<double> maxdeviations; //stores s*cutoff for each dimension.
    Dimension<nfixed> dim;
    unsigned int rank;
    double cutoff;
    protected:
    using Sampler<ET,MT,Engine,Sseq>::sieveptr;
    using Sampler<ET,MT,Engine,Sseq>::engine;
//    vector<MyLatticePoint> basis;
};



#endif
