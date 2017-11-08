#ifndef RELEVANT_COORDS_H
#define RELEVANT_COORDS_H


#include "GlobalStaticData.h"

#define NUM_OF_HASH_VALUES 64



namespace GaussSieve{

template<int nfixed>
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


//TODO: REPLACE rand() by a proper rand
template<int nfixed>
class RelevantCoordinates
{
  friend StaticInitializer<RelevantCoordinates<nfixed>>;

  RelevantCoordinates(RelevantCoordinates const &) = delete;
  RelevantCoordinates(RelevantCoordinates &&obj)   = delete;


  RelevantCoordinates  &operator=(RelevantCoordinates const &obj) = delete;
  RelevantCoordinates  &operator=(RelevantCoordinates &obj)       = delete;


  // for 0<=i<=63; 0<=j<=4
  uint_fast16_t get_ij_value(uint_fast16_t i, uint_fast16_t j)
  {
    return (rel_coo[i][j]);
  }

  //member
  private:
  static std::array<uint_fast16_t, 4> rel_coo[NUM_OF_HASH_VALUES];
  static int test;

};

template<int nfixed> int RelevantCoordinates<nfixed>::test = 10;
template<int nfixed> std::array<uint_fast16_t,4> RelevantCoordinates<nfixed>::rel_coo[NUM_OF_HASH_VALUES] = {};


// Static Initializer:
template<int nfixed> class StaticInitializer<RelevantCoordinates<nfixed>>
: public DefaultStaticInitializer<RelevantCoordinates<nfixed>>
{
  using Parent = DefaultStaticInitializer<RelevantCoordinates<nfixed>>;
public:

  template<class T,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const & initializer) : StaticInitializer(initializer.dim) {} //<-WHAT IS IT FOR?

  StaticInitializer()
  {
    assert(Parent::user_count > 0);
    if(Parent::user_count>1)
    {
      assert(false);
    }
    else
    {
      RelevantCoordinates<nfixed>::test = 10;

      for (uint_fast16_t i=0; i<64; ++i)
      {
        RelevantCoordinates<nfixed>::rel_coo[i][0] = rand() % nfixed;
        RelevantCoordinates<nfixed>::rel_coo[i][1] = rand() % nfixed;
        RelevantCoordinates<nfixed>::rel_coo[i][2] = rand() % nfixed;
        RelevantCoordinates<nfixed>::rel_coo[i][3] = rand() % nfixed;
      }

    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing RelevantCoordinates with nfixed = " << nfixed  << " Counter is" << Parent::user_count )
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing RelevantCoordinates with nfixed = " << nfixed << " Counter is " << Parent::user_count )
  }
};


}


#endif
