#ifndef GAUSS_SIEVE_STATISTICS_H
#define GAUSS_SIEVE_STATISTICS_H

#include "DefaultIncludes.h"

/**
 This file defines the SieveStatistics class that is responsible for collecting various types of
 statistics.
*/

namespace GaussSieve{

template<class SieveTraits, bool MT> struct GaussSieveStatistics;
template<class SieveTraits>
struct GaussSieveStatistics<SieveTraits,false>
{

  GaussSieveStatistics()
  : number_of_collisions(0),
    number_of_points_sampled(0),
    number_of_points_constructed(0),
    current_list_size(0),
    filtered_list_size(0),
    number_of_scprods_level1(0),
    number_of_scprods_level2(0),
    number_of_scprods_level3(0),

    number_of_exact_scprods(0),
    number_of_mispredictions(0)
    {}
// TODO: Move parts of these statistics into the actual object they relate to.
// e.g. there is no reason to maintain list sizes...
// Only the getters should remain.

  inline unsigned long int get_number_of_collisions() const {return number_of_collisions;};
  unsigned long int number_of_collisions;
  inline void increment_number_of_collisions() {++number_of_collisions; }

  unsigned long int number_of_points_sampled;
  inline unsigned long int get_number_of_points_sampled() const {return number_of_points_sampled;};
  inline void increment_number_of_points_sampled() { ++number_of_points_sampled; }

  unsigned long long int number_of_points_constructed; //sampling  + succesful pairs
  inline unsigned long long get_number_of_points_constructed() const {return number_of_points_constructed;};
  inline void increment_number_of_points_constructed() { ++number_of_points_constructed; }

  unsigned long int get_current_list_size() const {return current_list_size;};
  unsigned long int current_list_size;
  inline void increment_current_list_size() { ++current_list_size; }
  inline void increment_current_list_size_by(long int const amount) {current_list_size+=amount;}
  inline void decrement_current_list_size() { --current_list_size; }

  inline unsigned long int get_filtered_list_size() const    {return filtered_list_size;};
  unsigned long int filtered_list_size; //size of filtered_list
  inline void set_filtered_list_size(unsigned long int const to) { filtered_list_size = to;}
  inline void increment_filtered_list_size() {++filtered_list_size;}



  unsigned long long int number_of_scprods_level1; //for k=2 case
  unsigned long long get_number_of_scprods_level1() const  {return number_of_scprods_level1;}
  inline void increment_number_of_scprods_level1() {++number_of_scprods_level1;}

  unsigned long long int number_of_scprods_level2; //for k=2,3 cases
  inline unsigned long long get_number_of_scprods_level2() const {return number_of_scprods_level2;}
  inline void increment_number_of_scprods_level2() {++number_of_scprods_level2;}

  unsigned long long int number_of_scprods_level3; //for k=2,3,4 cases
  inline unsigned long long get_number_of_scprods_level3() const {return number_of_scprods_level3;}




  unsigned long long int number_of_exact_scprods;
  unsigned long long int number_of_mispredictions; //could not reduce in spite of approximation saying so

//TODO: total time spent?

/* Old MT Variant
//note: we might collect statistics per-thread and merge occasionally. This means these statistics might be inaccurate.
    atomic_ulong number_of_collisions;
    atomic_ulong number_of_points_sampled;
    atomic_ullong number_of_points_constructed;
    atomic_ulong current_list_size;
    atomic_ullong number_of_scprods;
    atomic_ullong number_of_total_scprods_level1;
    atomic_ullong number_of_exact_scprods;
    atomic_ullong number_of_mispredictions;
*/




};












} // end namespace GaussSieve

#endif // include guard
