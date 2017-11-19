#ifndef RELEVANT_COORDS_H
#define RELEVANT_COORDS_H

#include <random>
#include "GlobalStaticData.h"

namespace GaussSieve{

//template<int nfixed>
class RelevantCoordinates;

// SIMPLE SINGLETON
// Thread-safe version according to http://www.modernescpp.com/index.php/thread-safe-initialization-of-a-singleton

//class RelevantCoordinates
//{
//private:
//  //TODO: pass a seed to rand()
//  RelevantCoordinates(int ambient_dim)
//  {
//    //std::cout << "inside constructor " << std::endl;
//    for (uint_fast16_t i=0; i<64; ++i)
//    {
//      rel_coo[i][0] = rand() % ambient_dim;
//      rel_coo[i][1] = rand() % ambient_dim;
//      rel_coo[i][2] = rand() % ambient_dim;
//      rel_coo[i][3] = rand() % ambient_dim;
//    }
//    //std::cout << "finish initializing rel_coo " << std::endl;
//  }
//
//  RelevantCoordinates(RelevantCoordinates const &) = delete;
//  void operator=(RelevantCoordinates const &)      = delete;
//
//  ~RelevantCoordinates()
//  {
//    std::cout << "Dtor" << std::endl;
//  }
//
//public:
//  static RelevantCoordinates& get_instance(int ambient_dim)
//  {
//    //std::cout << "inside get_instance " << std::endl;
//    static RelevantCoordinates single_instance(ambient_dim);
//    return single_instance;
//  }
//
//  // for 0<=i<=63; 0<=j<=4
//  uint_fast16_t get_ij_value(uint_fast16_t i, uint_fast16_t j)
//  {
//    return (rel_coo[i][j]);
//  }
//
//  //member
//  std::array<uint_fast16_t, 4> rel_coo[64];
//
//};

  
  uint_fast16_t constexpr  sim_hash_len = 64;
  uint_fast16_t constexpr  sim_hash2_len = 64;
  
  uint_fast16_t constexpr  num_of_coord = 4;

//TODO: REPLACE rand() by a proper rand
//template<int nfixed>
class RelevantCoordinates
{
  friend StaticInitializer<RelevantCoordinates>;
  
public:
  RelevantCoordinates() = delete; 

  RelevantCoordinates(RelevantCoordinates const &) = delete;
  RelevantCoordinates(RelevantCoordinates &&obj)   = delete;


  RelevantCoordinates  &operator=(RelevantCoordinates const &obj) = delete;
  RelevantCoordinates  &operator=(RelevantCoordinates &obj)       = delete;


  // for 0<=i<sim_hash_len; 0<=j<num_of_coord
  static uint_fast16_t get_ij_value(uint_fast16_t i, uint_fast16_t j)
  {
    return (rel_coo[i][j]);
  }
  
public:
  static void print()
  {
    std::cout << "relevant matrix is: " << std::endl;
    for (uint_fast16_t i=0; i<sim_hash_len; i++)
    {
      for (uint_fast16_t j=0; j<num_of_coord; ++j)
        std::cout << rel_coo[i][j] << ", ";
        
      std::cout << std::endl;
      //std::cout << rel_coo[i][0] << "," << rel_coo[i][1] << "," << rel_coo[i][2] << "," << rel_coo[i][3] <<std::endl;
    }
  }

  //member
  private:
  static std::array<uint_fast16_t, num_of_coord> rel_coo[sim_hash_len];

};

std::array<uint_fast16_t,num_of_coord> RelevantCoordinates::rel_coo[sim_hash_len] = {};


// Static Initializer:
//template<int nxfixed>
template<>
class StaticInitializer<class RelevantCoordinates>
: public DefaultStaticInitializer<RelevantCoordinates>
{
  using Parent = DefaultStaticInitializer<RelevantCoordinates>;
public:

  //template<class T,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  //StaticInitializer(T const & initializer) : StaticInitializer(initializer.dim) {} //<-WHAT IS IT FOR?

  StaticInitializer(int abmient_dimension)
  {
    assert(Parent::user_count > 0);
    if(Parent::user_count>1)
    {
      //assert(false);
    }
    else
    {
      std::mt19937 rng;
      rng.seed(std::random_device()());
      std::uniform_int_distribution<std::mt19937::result_type> distr(0,abmient_dimension-1);
      
      for (uint_fast16_t i=0; i<sim_hash_len; ++i)
      {
        RelevantCoordinates::rel_coo[i][0] = distr(rng);
        RelevantCoordinates::rel_coo[i][1] = distr(rng);
        RelevantCoordinates::rel_coo[i][2] = distr(rng);
        RelevantCoordinates::rel_coo[i][3] = distr(rng);
        /*
        RelevantCoordinates::rel_coo[i][0] = rand() % abmient_dimension;
        RelevantCoordinates::rel_coo[i][1] = rand() % abmient_dimension;
        RelevantCoordinates::rel_coo[i][2] = rand() % abmient_dimension;
        RelevantCoordinates::rel_coo[i][3] = rand() % abmient_dimension;
         */
      }
      RelevantCoordinates::print();

    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing RelevantCoordinates; Counter is " << Parent::user_count )
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing RelevantCoordinates; Counter is " << Parent::user_count )
  }
};


}


#endif
