// clang-format off

#ifndef SAMPLER_CPP
#define SAMPLER_CPP

#include "Sampler.h"
#include <type_traits>

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


template class MTPRNG<std::mt19937_64,false, std::seed_seq>;
//template class MTPRNG<std::mt19937,true,  std::seed_seq>;
template class Sampler<fplll::Z_NR<long>, false, std::mt19937_64,std::seed_seq, -1>;
//template class Sampler<Z_NR<long>, true,  std::mt19937,std::seed_seq>;


#endif

//clang-format on
