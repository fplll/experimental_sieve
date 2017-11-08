#ifndef RELEVANT_COORDS_H
#define RELEVANT_COORDS_H


//Singleton

namespace GaussSieve{
  
// Thread-safe version according to http://www.modernescpp.com/index.php/thread-safe-initialization-of-a-singleton
class RelevantCoordinates;

class RelevantCoordinates
{
private:
  //TODO: pass a seed to rand()
  RelevantCoordinates(int ambient_dim)
  {
    //std::cout << "inside constructor " << std::endl;
    for (uint_fast16_t i=0; i<64; ++i)
    {
      rel_coo[i][0] = rand() % ambient_dim;
      rel_coo[i][1] = rand() % ambient_dim;
      rel_coo[i][2] = rand() % ambient_dim;
      rel_coo[i][3] = rand() % ambient_dim;
    }
    //std::cout << "finish initializing rel_coo " << std::endl;
  }
  
  RelevantCoordinates(RelevantCoordinates const &) = delete;
  void operator=(RelevantCoordinates const &)      = delete;
  
  ~RelevantCoordinates()
  {
    std::cout << "Dtor" << std::endl;
  }
  
public:
  static RelevantCoordinates& get_instance(int ambient_dim)
  {
    //std::cout << "inside get_instance " << std::endl;
    static RelevantCoordinates single_instance(ambient_dim);
    return single_instance;
  }
  
  // for 0<=i<=63; 0<=j<=4
  uint_fast16_t get_ij_value(uint_fast16_t i, uint_fast16_t j)
  {
    return (rel_coo[i][j]);
  }
  
  //member
  std::array<uint_fast16_t, 4> rel_coo[64];
  
};
  
}


#endif
