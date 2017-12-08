#ifndef GPV_SAMPLER_IMPL_H
#define GPV_SAMPLER_IMPL_H

#include "DefaultIncludes.h"
#include "Sampler.h"
#include "GPVSampler.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include <random>
#include <vector>
#include <math.h>
#include "LatticeBases.h"

namespace GaussSieve
{
  
  template <class SieveTraits, bool MT, class Engine, class Sseq>
  void GPVSampler<SieveTraits, MT, Engine, Sseq>::custom_init(SieveLatticeBasis<SieveTraits,MT> const & input_basis)
  {
    assert(!initialized);
#ifndef DEBUG_SIEVE_STANDALONE_SAMPLER
    assert(sieveptr!=nullptr);
#else
    assert(sieveptr==nullptr);
#endif
    
    dim           = input_basis.ambient_dimension;
    lattice_rank  = input_basis.lattice_rank;
    
    mu_matrix     = input_basis.get_mu_matrix();
    
    
    s2pi.resize(lattice_rank);
    maxdeviations.resize(lattice_rank);
    basis.resize(lattice_rank);
    
    auto const maxbistar2 = input_basis.get_maxbistar2();
    
    double st_dev = maxbistar2*2.2; // square of the st.dev guaranteed by GPV
    
    //std::cout << st_dev << std::endl;
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      
      double maxdev_nonsc = st_dev / convert_to_double(input_basis.get_g(i,i)); // g_(i,i) = ||b^*_i||^2
      
      s2pi[i] = maxdev_nonsc / GaussSieve::pi;
      maxdeviations[i] = sqrt(maxdev_nonsc) * cutoff;
      
      //std::cout <<"maxdeviations[i] = " << maxdeviations[i] << " ";
      
      basis[i] = input_basis.get_basis_vector(i).make_copy();
    }
    
    using RetType = typename SieveTraits::GaussSampler_ReturnType;
    
    if(static_init_plainpoint!=nullptr)
    {
      assert(false);
    }
    if(static_init_rettype!=nullptr)
    {
      assert(false);
    }
    
    static_init_rettype   = new StaticInitializer<RetType>(MaybeFixed<SieveTraits::get_nfixed>{dim});
    static_init_plainpoint= new StaticInitializer<typename SieveTraits::PlainPoint>(MaybeFixed<SieveTraits::get_nfixed>{dim});
    
    initialized = true;


  }
  
  
  //TODO: progressive rank + babai_start
  template <class SieveTraits, bool MT, class Engine, class Sseq>
  typename SieveTraits::GaussSampler_ReturnType
  GPVSampler<SieveTraits, MT, Engine, Sseq>::sample(int const thread)
  {
    assert(initialized);
#ifdef DEBUG_SIEVE_STANDALONE_SAMPLER
    assert(sieveptr==nullptr);
#else
    assert(sieveptr!=nullptr);
#endif
    
    typename SieveTraits::PlainPoint vec;
    vec.fill_with_zero();
    
    
    
    uint_fast16_t i = lattice_rank;
    
    std::vector<double> shifts(lattice_rank, 0.0);
    
    while ( vec.is_zero())
    {
      while(i>0)
      {
        --i;
        long const newcoeff = sample_z_gaussian_VMD<long, Engine>(
        s2pi[i], shifts[i], engine.rnd(), maxdeviations[i]);  // coefficient of b_j in vec.
        
        //std::cout<< i << " newcoeff = " << newcoeff << std::endl;
        
        vec += basis[i] * newcoeff;
        
        for (uint_fast16_t j = lattice_rank-1; j >i; --j)  // adjust shifts
        {
          shifts[i] -= newcoeff * (mu_matrix[j][i]);
          //std::cout << mu_matrix[j][i] <<std::endl;
        }
      }
    }
   
    typename SieveTraits::GaussSampler_ReturnType ret;
    ret = make_from_any_vector<typename SieveTraits::GaussSampler_ReturnType>(vec, dim);
    
    return ret;
  }
  
}

#endif
