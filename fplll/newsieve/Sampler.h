#ifndef SAMPLER_H
#define SAMPLER_H

//forward declarations
#include <random>
#include <iostream>
#include <cfenv>
#include <type_traits>
#include "Typedefs.h"

//forward declarations

template<class ET,bool MT, int nfixed> class Sieve;

template<class ET,bool MT, class Engine, class Sseq, int nfixed=-1> class Sampler;
template<class ET,bool MT, class Engine, class Sseq, int nfixed> std::ostream & operator<<(std::ostream &os, Sampler<ET,MT, Engine, Sseq, nfixed>* const samplerptr); //printing
template<class ET,bool MT, class Engine, class Sseq, int nfixed> std::istream & operator>>(std::istream &is, Sampler<ET,MT, Engine, Sseq, nfixed>* const samplerptr); //reading (may also be used by constructor from istream)
enum class SamplerType
{
    user_defined = 0,
    elliptic_sampler= 1,
    shi_sampler = 2,
    gauss_sampler =3
};

template<class Engine, bool MT, class Sseq> class MTPRNG;       //wrapper around (a vector of) random number engines of type Engine
                                                                //This is used to unify the single and multi-threaded case
namespace GaussSieve
{
template<class Z, class Engine>     Z sample_z_gaussian(double s, double const center, Engine & engine, double const cutoff);
    template<class Z, class Engine>     Z sample_z_gaussian_VMD(double const s2pi, double const center, Engine & engine, double const maxdeviation);
    //samples from a discrete Gaussian distribution with parameter s and center c. We cutoff the Gaussian at s*cutoff.
    //i.e. the distribution is discrete on Z with output probability for x being proportional to exp(- pi(x-c)^2/s^2). Note the scaling by pi in the exponent.
    //For reasons of numerical stability, center should not be very large in absolute value (it is possible to reduce to |center|<1 anyway), s.t.
    //center +/- cutoff * s does not overflow.
    //Z must be among one of short, int, long, long long. center needs to be representable as an (exact) double.
    //We do NOT support mpz_t here! Note that the output takes the role of coefficients wrt a given basis.
    //We only support double. For sieving algorithms, there is no really good reason to support higher precision.
    //Note: if one wants to have higher precision, one also needs to adjust the PRNGs to actually output high precision.

    //The variant sample_z_gaussian_VMD takes s2pi = s^2 * pi and cutoff * s as parameters.

}

//includes
#include "LatticePointsNew.h"

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
    cout << "Initializing RNGS engines" << endl << flush;
    engine.init(sieve->get_num_threads());
    cout << "Done. Starting custom initialization of specific sampler" << endl << flush;
    custom_init();
    cout << "Finished custom initialization" << endl << flush;
}

template<class ET,bool MT, class Engine, class Sseq> std::ostream & operator<<(std::ostream &os,Sampler<ET,MT,Engine,Sseq>* const samplerptr){return samplerptr->dump_to_stream(os);};
template<class ET,bool MT, class Engine, class Sseq> std::istream & operator>>(std::istream &is,Sampler<ET,MT,Engine,Sseq>* const samplerptr){return samplerptr->read_from_stream(is);};





#endif
