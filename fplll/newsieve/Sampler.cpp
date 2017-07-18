#ifndef SAMPLER_CPP
#define SAMPLER_CPP

#include "Sampler.h"

template<class Engine, class Sseq> void MTPRNG<Engine,true,Sseq>::reseed(Sseq & _seq)
{
    seeder.seed(_seq);
    unsigned int old_threads=num_threads;
    num_threads=0;
    init(old_threads); //will restart all engines, because num_threads = 0;
};

template<class Engine, class Sseq> void MTPRNG<Engine,true,Sseq>::init(unsigned int const _num_threads)
{
    if(_num_threads<=num_threads) //no need to initalize.
    {
        return;
    }
    engines.resize(_num_threads);
    engines.shrink_to_fit();
    uint32_t per_engine_seed[seed_length];
    //else initialize remaining threads
    for(unsigned int i=num_threads;i<_num_threads;++i)
    {

        for(unsigned int j=0;j<seed_length;++j)
        {
            per_engine_seed[j] = seeder();
        }
        std::seed_seq per_engine_see_seq(per_engine_seed, per_engine_seed+seed_length);
        engines[i].seed(per_engine_see_seq);
    }
    num_threads = _num_threads;
}

template<class Engine, class Sseq> void MTPRNG<Engine,false,Sseq>::reseed(Sseq & _seq)
{
    std::mt19937_64 seeder(_seq);
    uint32_t per_engine_seed[seed_length];
    for(unsigned int j=0;j<seed_length;++j)
    {
        per_engine_seed[j] = seeder();
    }
    std::seed_seq derived_seed_seq(per_engine_seed, per_engine_seed+seed_length);
    engine.seed(derived_seed_seq);
};

/**
 * sampling integral Gaussians by rejection sampling
 */

 //Z must be an integral POD class

template<class Z, class Engine> Z GaussSieve::sample_z_gaussian(double s, double const center, Engine & engine, double const cutoff)
{
//Note : The following allows to access / modify floating point exceptions and modes.
//#pragma STDC FENV_ACCESS on
//This is too compiler/implementation-specific and does not work most of the time...

    static_assert(is_integral<Z>::value,"Return type for sample_z_gaussian must be POD integral type.");
    Z maxdev = static_cast<Z>(std::ceil(s * cutoff)); //maximum deviation of the Gaussian from the center. Note that maxdev may be 1.
    std::uniform_int_distribution<Z> uniform_in_range (std::floor(center-maxdev),std::ceil(center+maxdev));
    std::uniform_real_distribution<double> rejection_test(0.0,1.0); //defaults to value from [0,1), used in rejection sampling.
    Z closest_int = std::round(center); //closest int to center, i.e. most likely value.
    double adj = -(center-closest_int)*(center-closest_int); //negative squared distance to most likely value. Used to scale up the Gaussian weight function s.t. it is 1 at the most likely value.
    s = s*s/pi; //overwriting s.
    //std::fenv_t env;
    //feholdexcept( &env); //This disables all floating-point exceptions.

//use rejection sampling
    while(true)
    {
        Z result = uniform_in_range(engine); //sample uniform result.
        double dist = result - center;
    //compute Gaussian weight. std::fma(dist,dist,adj) computes dist^2 + adj = (result-center)^2  - MIN{(result-center)^2 | result integral}.
    //s was overwritten to be s^2/pi.

    //Note that the argument of the exp-function might be a tiny positive value due to numeric error
    //(Even if result==closest_int, adj = ROUND((closest_int-center)^2), the computation of std::fma(dist,dist,adj) does not round the intermediate dist^2, leading to a non-zero argument)
    //In particular, it is conceivable that floating point underruns occur in the std::fma - call.
    //Furthermore, if cutoff is large or if s<<1 (in this case, the issue is the rounding when we determined the range), the argument to exp can be extremely small, leading to further potential underruns.
    //We do not care about this for now...

        if( rejection_test(engine) <  std::exp(-std::fma(dist,dist,adj)/s))
        {
            //std::feclearexcept(FE_UNDERFLOW);
            //std::feupdateenv(&env);
            return result;
        }
    }
}

template<class Z, class Engine> Z GaussSieve::sample_z_gaussian_VMD(double const s2pi, double const center, Engine & engine, double const maxdeviation)
{
//Note : The following allows to access / modify floating point exceptions and modes.
//#pragma STDC FENV_ACCESS on
//This is too compiler/implementation-specific and does not work most of the time...

    static_assert(is_integral<Z>::value,"Return type for sample_z_gaussian must be POD integral type.");
    std::uniform_int_distribution<Z> uniform_in_range (std::floor(center-maxdeviation),std::ceil(center+maxdeviation));
    std::uniform_real_distribution<double> rejection_test(0.0,1.0); //defaults to value from [0,1), used in rejection sampling.
    Z closest_int = std::round(center); //closest int to center, i.e. most likely value.
    double adj = -(center-closest_int)*(center-closest_int); //negative squared distance to most likely value. Used to scale up the Gaussian weight function s.t. it is 1 at the most likely value.

    //std::fenv_t env;
    //feholdexcept( &env); //This disables all floating-point exceptions.

//use rejection sampling
    while(true)
    {
        Z result = uniform_in_range(engine); //sample uniform result.
        double dist = result - center;
    //compute Gaussian weight. std::fma(dist,dist,adj) computes dist^2 + adj = (result-center)^2  - MIN{(result-center)^2 | result integral}.
    //s was overwritten to be s^2/pi.

    //Note that the argument of the exp-function might be a tiny positive value due to numeric error
    //(Even if result==closest_int, adj = ROUND((closest_int-center)^2), the computation of std::fma(dist,dist,adj) does not round the intermediate dist^2, leading to a non-zero argument)
    //In particular, it is conceivable that floating point underruns occur in the std::fma - call.
    //Furthermore, if cutoff is large or if s<<1 (in this case, the issue is the rounding when we determined the range), the argument to exp can be extremely small, leading to further potential underruns.
    //We do not care about this for now...

        if( rejection_test(engine) <  std::exp(-std::fma(dist,dist,adj)/s2pi))
        {
            //std::feclearexcept(FE_UNDERFLOW);
            //std::feupdateenv(&env);
            return result;
        }
    }
}

template class MTPRNG<std::mt19937_64,false, std::seed_seq>;
//template class MTPRNG<std::mt19937,true,  std::seed_seq>;
template class Sampler<Z_NR<long>, false, std::mt19937_64,std::seed_seq>;
//template class Sampler<Z_NR<long>, true,  std::mt19937,std::seed_seq>;


#endif
