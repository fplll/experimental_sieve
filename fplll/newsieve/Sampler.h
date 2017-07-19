// This file provided the interface for the lattice point sampler.
// To this end, this file provides a virtual interface template class Sampler, from which the actual
// samplers are derived.

// Dependencies:
// The sampler instance contains a pointer back to the sieve main object, creating circular
// dependencies.
// For this reason, parts of the implementation are in Sampler_impl.h or Sampler.cpp (TODO: Decide
// on whether we want header-only)

#ifndef SAMPLER_H
#define SAMPLER_H

// forward declarations
#include "Typedefs.h"
#include <iostream>
#include <random>
#include <type_traits>
//#include <cfenv>

// forward declarations

template <class ET, bool MT, int nfixed> class Sieve;

template <class ET, bool MT, class Engine, class Sseq, int nfixed> class Sampler;

// printing
template <class ET, bool MT, class Engine, class Sseq, int nfixed>
inline std::ostream &operator<<(std::ostream &os,
                                Sampler<ET, MT, Engine, Sseq, nfixed> *const samplerptr);

// reading (may also be used by constructor from istream)
template <class ET, bool MT, class Engine, class Sseq, int nfixed>
inline std::istream &operator>>(std::istream &is,
                                Sampler<ET, MT, Engine, Sseq, nfixed> *const samplerptr);

enum class SamplerType
{
  user_defined     = 0,
  elliptic_sampler = 1,
  shi_sampler      = 2,
  gauss_sampler    = 3
};

template <class Engine, bool MT, class Sseq>
class MTPRNG;  // wrapper around (a vector of) random number engines of type Engine
               // This is used to unify the single and multi-threaded case

namespace GaussSieve
{

/**

 These functions sample from a discrete Gaussian distribution with parameter s and center c
 on (the 1-dimensional lattice given by) the integers Z.

 We cutoff the Gaussian at s*cutoff. This means that the distribution is discrete on a subset of
 Z with  output probability for x being proportional to exp(-pi(x-c)^2/s^2). Note the scaling by pi
 in the exponent.

 For reasons of numerical stability, center should not be very large in absolute value (it is
 possible to reduce to |center|<1 anyway), such that  center +/- cutoff * s does not overflow.
 Z must be an integral POD type (e.g. short, int, long, long long).
 center needs to be representable as an (exact) double.
 We do NOT support mpz_t here! The output will take the role of coefficients wrt a given basis.

 We only support double for the floating point numbers. For sieving algorithms, there is no really
 good reason for now to support different precisions, as sampling does not dominate anyway.

 Note: if one ever wants to have higher precision, one also needs to adjust the PRNGs to actually
 output high precision.

 engine is supposed to be a random number engine (as defined by the STL).

 The variant sample_z_gaussian_VMD takes
 s2pi = s^2 * pi and maxdeviation = cutoff * s as parameters.

 **/

template <class Z, class Engine>
Z sample_z_gaussian(double s, double const center, Engine &engine, double const cutoff);

template <class Z, class Engine>
Z sample_z_gaussian_VMD(double const s2pi, double const center, Engine &engine,
                        double const maxdeviation);
}

// clang-format off

//The class MTPRNG is just a wrapper around a PRNG engine to facilitate switching to multi-threaded.
//Due to the fact that we support multi-threading, MTPRNG<Engine,true,.> is a wrapper around
//a vector of Engines, whereas, MTPRNG<Engine,false,.> is a wrapper around Engine.
//reseed seeds *all* engines. Use rnd(thread-id) to obtain the underlying Engine.
//Thread-safety: Init and reseed are not thread-safe. Concurrent calls to rnd are fine.
//You may concurrently call rnd and use init to increase (but not decrease) the number of threads.

//Randomness path:
//The global (master) seed that is input to the constructor resp. to reseed is used to create 20x32 bit per-thread-seeds for each actual Engine.
//The output from theses engine(s) is then accessed using rnd(thread-number).
//This is done even in the single-threaded case to ensure consistency.
//Note that for obtaining the per-thread seeds from the master seeds, we use a fixed Mersenne twister engine and not the engine given as template parameter.

template<class Engine, class Sseq>  class MTPRNG<Engine,true, Sseq> //multithreaded case of MTPRNG
{
    public:
    MTPRNG(Sseq & _seq = {}) : seeder(_seq), engines(0), num_threads(0)       {}; //constructs an uninitialized MTPRNG
    void reseed(Sseq & _seq);
    void init(unsigned int const _num_threads); //will make sure at least _num_threads engines are actually running, starting new ones as desired.
                                                //Will never reseed/restart already running engines.
                                                //Reducing the number of threads and increasing it back saves the random state (unless we reseed).
    Engine & rnd(unsigned int const thread)                                 {return engines[thread];};
    private:
    std::mt19937_64 seeder; //seeded with initial seq and consecutively used to seed the children PRNGs.
    std::vector<Engine> engines;
    unsigned int num_threads; //number of initialized engines. May differ from size of the vector. In particular, num_threads = 0 means uninitialized.
    Sseq seq;
    static unsigned int constexpr seed_length = 20; //number of 32bit values to use as seeds for the underlying engine(s).
                                                    //Technically, we could use state_size if the engine provides it, but not all default engines do.
};


template<class Engine, class Sseq>  class MTPRNG<Engine, false, Sseq>           //singlethreaded case of MTPRNG: just wrapper around Engine
{
    public:
    MTPRNG(Sseq & _seq ={}) : engine()                      {reseed(_seq);};
    void reseed(Sseq & _seq);
    void init(unsigned int const = 1)                       {} //does nothing.
    Engine & rnd(unsigned int const)                        {return engine;};   //Argument is number of thread. It is ignored.
    Engine & rnd()                                          {return engine;};   //Version without thread-id
    private:
    Engine engine;
    static unsigned int constexpr seed_length = 20; //number of 32bit values to use as (per-thread) seed for the underlying engine.
                                                    //Technically, we could use state_size if the engine provides it, but not all default engines do.
};      //End of MTPRNG





//generic Sampler. All other sampler are derived from it.
template<class ET,bool MT, class Engine, class Sseq, int nfixed> //Sseq is supposed to satisfy the C++ concept "SeedSequence". The standard library has std::seed_seq as a canonical example.
                                                     //Engine is supposed to satisfy the C++ concept of a "Random Number Engine". <random> provides several of those, e.g. std::mt19937_64.
class Sampler
//Note :    In multi-threaded environment, we only have 1 sampler object. thread-number is given to sample();
{
    public:
    using GaussSampler_ReturnType = typename GaussSieve::GaussSampler_ReturnType<ET,MT,nfixed>;
    friend std::ostream & operator<< <ET,MT>(std::ostream &os, Sampler<ET,MT,Engine,Sseq,nfixed>* const samplerptr);
    friend std::istream & operator>> <ET,MT>(std::istream &is, Sampler<ET,MT,Engine,Sseq,nfixed>* const samplerptr);

    Sampler<ET,MT,Engine,Sseq,nfixed> (Sseq & initial_seed): engine(initial_seed), sieveptr(nullptr)                      {}
    //We call init first, then custom_init (via init).
    void init(Sieve<ET,MT,nfixed> * const sieve);
    virtual ~Sampler()=0; //needs to be virtual
    virtual SamplerType  sampler_type() const {return SamplerType::user_defined;};    //run-time type information.
                                                                    //This may be used to determine how to interpret a dump file.
                                                                    //defaults to user-defined.
                                                                    //Other values mean that the GaussSieve dumping routine is aware of the type, simplifying the syntax for dumping / reading.
    //TODO : Write some explanation how to do that.
    virtual GaussSampler_ReturnType sample(int thread=0)=0; //thread is the index of the calling thread (we need to keep separate PRNGs for each thread)
    //TODO : Allow sampling in subspaces, updating basis.

    private:
    virtual void custom_init()                                                                  {}         //called before any points are sampled;
    virtual std::ostream & dump_to_stream(std::ostream &os)  {return os;};    //dummy implementation of << operator.
    virtual std::istream & read_from_stream(std::istream &is){return is;};    //dummy implementation of >> operator.
    protected:
    MTPRNG<Engine, MT, Sseq> engine; //or engines
    Sieve<ET,MT,nfixed> * sieveptr; //pointer to parent sieve. Set in init();
};


#include "SieveGauss.h"


template <class ET,bool MT, class Engine, class Sseq, int nfixed> Sampler<ET,MT, Engine,Sseq,nfixed>::~Sampler() {} //actually needed, even though destructor is pure virtual as the base class destructor is eventually called implicitly.

template <class ET,bool MT, class Engine, class Sseq, int nfixed> void Sampler<ET,MT,Engine,Sseq,nfixed>::init(Sieve<ET,MT,nfixed> * const sieve)
{
    sieveptr = sieve;
//    std::cout << "Initializing RNGS engines" << std::endl << flush;
    engine.init(sieve->get_num_threads());
//    cout << "Done. Starting custom initialization of specific sampler" << endl << flush;
    custom_init();
//    cout << "Finished custom initialization" << endl << flush;
}

template<class ET,bool MT, class Engine, class Sseq, int nfixed> inline std::ostream & operator<<(std::ostream &os,Sampler<ET,MT,Engine,Sseq,nfixed>* const samplerptr){return samplerptr->dump_to_stream(os);};
template<class ET,bool MT, class Engine, class Sseq, int nfixed> inline std::istream & operator>>(std::istream &is,Sampler<ET,MT,Engine,Sseq,nfixed>* const samplerptr){return samplerptr->read_from_stream(is);};



template<class Z, class Engine> Z GaussSieve::sample_z_gaussian(double s, double const center, Engine & engine, double const cutoff)
{
//Note : The following allows to access / modify floating point exceptions and modes.
//#pragma STDC FENV_ACCESS on
//This is too compiler/implementation-specific and does not work most of the time...

    static_assert(std::is_integral<Z>::value,"Return type for sample_z_gaussian must be POD integral type.");
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

    static_assert(std::is_integral<Z>::value,"Return type for sample_z_gaussian must be POD integral type.");
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

        if( rejection_test(engine) < std::exp(-std::fma(dist,dist,adj)/s2pi))
        {
            //std::feclearexcept(FE_UNDERFLOW);
            //std::feupdateenv(&env);
            return result;
        }
    }
}


#endif

//clang-format on
