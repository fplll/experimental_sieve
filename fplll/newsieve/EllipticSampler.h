//NOT WOKING ATM, does not yet used new lattic point classes.

#ifndef ELLIPTIC_SAMPLER_H
#define ELLIPTIC_SAMPLER_H

#include "Sampler.h"

#include <random>

template<class ET, bool MT, class Engine, class Sseq>
class EllipticSampler;

template<class ET, bool MT, class Engine, class Sseq>
class EllipticSampler: public Sampler<ET,MT, Engine, Sseq>
{
    public:
    EllipticSampler(Sseq & seq, double const s=0.8, double const cutoff = 10.0)
        :   Sampler<ET,MT,Engine,Sseq>(seq),s2pi(s*s/GaussSieve::pi),maxdeviations(s*cutoff) {};
    virtual SamplerType  sampler_type() const override                          {return SamplerType::elliptic_sampler;};
    virtual ~EllipticSampler();
    virtual LatticePoint<ET> sample(int thread=0) override;
    private:
    virtual void custom_init() override;
    ZZ_mat<typename ET::underlying_data_type> current_basis;
    Matrix<FP_NR<double > > mu;
    double s2pi; //stores standard dev. for each dimension, already squared and multiplied by pi.
    double maxdeviations; //stores s*cutoff for each dimension.
    protected:
    using Sampler<ET,MT,Engine,Sseq>::sieveptr;
    using Sampler<ET,MT,Engine,Sseq>::engine;
};

template<class ET,bool MT, class Engine, class Sseq>
void EllipticSampler<ET,MT,Engine,Sseq>::custom_init()
{
    current_basis = sieveptr->get_original_basis();
    Matrix<ET> u, u_inv; //intentionally uninitialized.
    MatGSO<ET, FP_NR<double> > GSO(current_basis, u, u_inv, MatGSOFlags::GSO_INT_GRAM);
    GSO.update_gso(); //todo: raise exception in case of error.
    mu = GSO.get_mu_matrix();
//    g  = pGSO.get_g_matrix();
//    F maxbistar2 = pGSO.get_max_bstar();
//    F tmp;
//    variances.resize(rank); variances.shrink_to_fit();
//    maxdeviations.resize(rank); maxdeviations.shrink_to_fit();
//

}
template<class ET,bool MT, class Engine, class Sseq>
EllipticSampler<ET,MT,Engine, Sseq>::~EllipticSampler()
{

}

template<class ET,bool MT, class Engine, class Sseq>
LatticePoint<ET> EllipticSampler<ET,MT,Engine, Sseq>::sample(int thread)
{
    assert(sieveptr!=nullptr);
    unsigned int const dim = sieveptr->get_ambient_dimension();
    unsigned int const rank = sieveptr->get_lattice_rank();
    NumVect<ET> vec(dim); vec.fill(0); //current vector built up so far.
    vector<double> shifts(rank, 0.0); //shift, expressed in coordinates wrt the Gram-Schmidt basis.
    for(int j=rank-1; j>=0; --j)
    {
        long const newcoeff = GaussSieve::sample_z_gaussian_VMD<long,Engine>(s2pi,shifts[j],engine.rnd(),maxdeviations); //coefficient of b_j in vec.
        //vec+= current_basis[j].get_underlying_row(); //build up vector
        vec.addmul_si(current_basis[j].get_underlying_row(), newcoeff);
        for(int i=0;i<j;++i) //adjust shifts
        {
            shifts[i]-=newcoeff* (mu[j][i].get_d() );
        }
    }
    return vec; //converts to LatticePoint<ET>
}



#endif
