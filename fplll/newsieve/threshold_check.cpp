#define USE_REGULAR_QUEUE //comment out if you use priority-queue

//#define USE_LSH
#define PROGRESSIVE // progressive-sieve on

//#define USE_APPROXPOINT //<-does not work yet

#include "fplll.h"
#include <thread>
#include <chrono>
#include "SieveGauss_main.h"
#include <iostream>
#include <fstream>
#include <random>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include <unistd.h>



using namespace fplll;
using namespace GaussSieve;

template <class ZT> void test_run_sieve(int dim, std::ofstream &ofs)
{
  
}

int main(int argc, char **argv)
{
  
  
  const unsigned int bases_per_dim  = 4;
  const unsigned int dim_min = 45;
  const unsigned int dim_max = 59;
  
  std::array<double, dim_max-dim_min+1> av_time;
  std::array<int, dim_max-dim_min+1> dims;
  std::array<unsigned long, dim_max-dim_min+1> av_list_size;
  
  bool constexpr multithreaded = false;
  int k = 2;
  
  using Traits = GaussSieve::DefaultSieveTraits<int32_t, false, -1, ZZ_mat< mpz_t > >;
  TerminationCondition<Traits,multithreaded> * termcond;
  termcond = new MinkowskiTerminationCondition<Traits, multithreaded>;
  
  ZZ_mat<mpz_t> B;
 
  
  for (unsigned dim = dim_min; dim<=dim_max;++dim)
  {
    B.resize(dim, dim);
    
    double time_cnt = 0;
    unsigned long list_size_cnt = 0;
    
    for (unsigned int seed = 1; seed <=bases_per_dim; ++seed)
    {
      string input_file_name = "Inputs/dim"+std::to_string(dim)+"_seed"+std::to_string(seed);
      std::ifstream input_file(input_file_name);
      std::cout << "reading B from file ..." << std::endl;
            
      input_file >> B;
      input_file.close();
      
      lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);
      
      Sieve<Traits, multithreaded> Test_2Sieve (B, k, 0);
      Test_2Sieve.set_termination_condition(termcond);
      
      auto start = std::chrono::high_resolution_clock::now();
      Test_2Sieve.run();
      
      auto finish = std::chrono::high_resolution_clock::now();
      auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
      
      time_cnt+=microseconds.count()/1000000.0;
      
      list_size_cnt+=Test_2Sieve.statistics.get_current_list_size();
      //list_size_cnt+=Test_2Sieve.get_final_list_size();
    }
    
    dims[dim-dim_min] = dim;
    av_time[dim-dim_min] = time_cnt/ bases_per_dim;
    av_list_size[dim-dim_min] = list_size_cnt / bases_per_dim;
  }
  
  for (unsigned int i=0; i<av_time.size(); ++i)
  {
    std::cout << dims[i] <<" " << av_time[i] << " " <<av_list_size[i] << std::endl;
  }

}