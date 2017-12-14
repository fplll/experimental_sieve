#ifndef GAUSS_SIEVE_STATISTICS_H
#define GAUSS_SIEVE_STATISTICS_H

#include "DefaultIncludes.h"

/**
 This file defines the SieveStatistics class that is responsible for collecting various types of
 statistics.
*/

namespace GaussSieve
{

template <class SieveTraits, bool MT> class Sieve;

template <class SieveTraits, bool MT> struct GaussSieveStatistics;

template <class SieveTraits>
struct GaussSieveStatistics<SieveTraits, false>
{

  GaussSieveStatistics(Sieve<SieveTraits,false> *backptr)
      : sieveptr(backptr),
        number_of_collisions(0),
        number_of_points_sampled(0),
        number_of_points_constructed(0),
        // current_list_size(0),
        filtered_list_size(0),
        number_of_scprods_level1(0),
        number_of_scprods_level2(0),
        number_of_scprods_level3(0),
        number_of_2reds(0),
        number_of_3reds(0),
        number_of_exact_scprods(0),
        number_of_mispredictions(0)
  {
  }
// TODO: Move parts of these statistics into the actual object they relate to.
// e.g. there is no reason to maintain list sizes outside of the actual list objects...
// Only the getters should remain.

  Sieve<SieveTraits, false> *const sieveptr;

  inline unsigned long int get_number_of_collisions() const { return number_of_collisions; }
  unsigned long int number_of_collisions;
  inline void increment_number_of_collisions() { ++number_of_collisions; }

  unsigned long int number_of_points_sampled;
  inline unsigned long int get_number_of_points_sampled() const { return number_of_points_sampled; }
  inline void increment_number_of_points_sampled() { ++number_of_points_sampled; }

  unsigned long long int number_of_points_constructed;  // sampling  + succesful pairs
  inline unsigned long long get_number_of_points_constructed() const { return number_of_points_constructed; }
  inline void increment_number_of_points_constructed() { ++number_of_points_constructed; }

  unsigned long int get_current_list_size() const { return sieveptr->main_list.size(); }

  inline unsigned long int get_filtered_list_size() const    {return filtered_list_size;};
  unsigned long int filtered_list_size;  // size of filtered_list
  inline void set_filtered_list_size(unsigned long int const to) { filtered_list_size = to; }
  inline void increment_filtered_list_size() { ++filtered_list_size; }

  unsigned long int get_current_queue_size()                  { return sieveptr->main_queue.size(); }  //TODO : fix const-correctness


  /*Stats for exact scalar products*/
  unsigned long long int number_of_scprods_level1;  // for k=2 case
  unsigned long long get_number_of_scprods_level1() const  { return number_of_scprods_level1; }
  inline void increment_number_of_scprods_level1() { ++number_of_scprods_level1; }

  unsigned long long int number_of_scprods_level2;  // for k=2,3 cases
  inline unsigned long long get_number_of_scprods_level2() const { return number_of_scprods_level2; }
  inline void increment_number_of_scprods_level2() { ++number_of_scprods_level2; }

  unsigned long long int number_of_scprods_level3; // for k=2,3,4 cases
  inline unsigned long long get_number_of_scprods_level3() const { return number_of_scprods_level3; }

  /* Stats for approximate scalar products */
  unsigned long long int number_of_approx_scprods_level1; //for k=2 case
  unsigned long long get_number_of_approx_scprods_level1() const  {return number_of_approx_scprods_level1;}
  inline void increment_number_of_approx_scprods_level1() {++number_of_approx_scprods_level1;}

  unsigned long long int number_of_approx_scprods_level2; //for k=2,3 cases
  inline unsigned long long get_number_of_approx_scprods_level2() const {return number_of_approx_scprods_level2;}
  inline void increment_number_of_approx_scprods_level2() {++number_of_approx_scprods_level2;}

  /* for 3-sieve, to compare number of 2-reds vs. 3-reds */

  unsigned int number_of_2reds;
  unsigned int get_number_of_2reds() const  {return number_of_2reds;}
  inline void increment_number_of_2reds() {++number_of_2reds;}


  unsigned int number_of_3reds;
  unsigned int get_number_of_3reds() const  {return number_of_3reds;}
  inline void increment_number_of_3reds() {++number_of_3reds;}


  unsigned long long int number_of_exact_scprods;
  unsigned long long int number_of_mispredictions; //could not reduce in spite of approximation saying so

  inline void dump_status_to_stream(std::ostream &of, int howverb);

  // Temporary code, to be removed

  #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
  inline void compute_and_print_statistics_lvl(std::ostream &of, int lvl, bool do_print);
  inline void compute_and_print_statistics_all(std::ostream &of);
  inline void compute_and_print_statistics_all_innloop(std::ostream &of);


    // THIS IS ONLY TO GET STATISTICS FOR BITAPPROX. to be deleted

  //#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
  std::array< std::vector<int>, SimHash::num_of_levels> no_red_stat;
  std::array< std::vector<int>, SimHash::num_of_levels> red_stat;

  std::array< std::vector<int>, SimHash::num_of_levels> no_red_stat_innloop;
  std::array< std::vector<int>, SimHash::num_of_levels> red_stat_innloop;

  //std::array< std::vector<int>, SimHash::num_of_levels> no_red_stat_layer;
  //std::array< std::vector<int>, SimHash::num_of_levels> red_stat_layer;
  #endif



//TODO: total time spent?

/* Old MT Variant : Do not delete
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

/**
 To be deleted
*/

/*
 Computes statistics for sim-hash
 */

#ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
template<class SieveTraits>
inline void GaussSieveStatistics<SieveTraits,false>::compute_and_print_statistics_lvl(std::ostream &of, int lvl, bool do_print)
  {
    using std::endl;


    unsigned long long sum_no_red = 0;
    unsigned long long sum_red = 0;

    for (unsigned int i=0; i<no_red_stat[lvl].size(); ++i)
    {
      sum_no_red+=no_red_stat[lvl][i];
      sum_red+=red_stat[lvl][i];
    }

    float pdf_no_red[no_red_stat[lvl].size()];
    float pdf_red[red_stat[lvl].size()];

    float cdf_no_red[no_red_stat[lvl].size()];
    float cdf_red[red_stat[lvl].size()];

    float accum_cdf_no_red = 0;
    float accum_cdf_red = 0;

    for (unsigned int i=0; i<no_red_stat[lvl].size(); ++i)
    {
      pdf_no_red[i] = (float)no_red_stat[lvl][i]/ (float)sum_no_red;
      pdf_red[i] = (float)red_stat[lvl][i] / (float) sum_red;

      accum_cdf_no_red+=no_red_stat[lvl][i];
      cdf_no_red[i] = accum_cdf_no_red/sum_no_red;

      accum_cdf_red+=red_stat[lvl][i];
      cdf_red[i] = accum_cdf_red/sum_red;

    }


    if (do_print)
    {
      std::ofstream myfile;
      myfile.open ("newsieve/Statistics/Statistics.txt");


      myfile << "Statistics for dim = " << sieveptr->lattice_rank << endl;
      myfile <<  std::setw(40) << " NO REDUCTION "<< std::setw(30) << " REDUCTION " << endl;

      for (unsigned int i=0; i<no_red_stat.size(); ++i)
      {
        myfile << " | " <<std::setw(3) << i <<"  | " << std::setw(10) << no_red_stat[lvl][i] << " | " <<
        std::setw(16) << pdf_no_red[i]  << " | " << std::setw(16) << cdf_no_red[i]  << " ||" <<
        std::setw(7) << red_stat[lvl][i] << " | " <<
        std::setw(16) << pdf_red[i]  << " | " << std::setw(16) << cdf_red[i]  << " |" << endl;
      }

      myfile.close();
    }
}


template<class SieveTraits>
inline void GaussSieveStatistics<SieveTraits,false>::compute_and_print_statistics_all(std::ostream &of)
  {
    using std::endl;
    std::ofstream myfile;
    myfile.open ("newsieve/Statistics/Statistics for k="+std::to_string(sieveptr->sieve_k) + "_dim="+std::to_string(sieveptr->lattice_rank));
    myfile << "Statistics for dim = " << sieveptr->lattice_rank << endl;

    /*for levels*/
    std::array<std::vector<float>, SimHash::num_of_levels> cdf_no_red;
    std::array<std::vector<float>, SimHash::num_of_levels> cdf_red;

    std::array<std::vector<float>, SimHash::num_of_levels> pdf_no_red;
    std::array<std::vector<float>, SimHash::num_of_levels> pdf_red;


    for (unsigned int lvl=0; lvl<SimHash::num_of_levels; ++lvl)
    {

      unsigned long long sum_no_red = 0;
      unsigned long long sum_red = 0;



      for (unsigned int i=0; i<no_red_stat[lvl].size(); ++i)
      {
        sum_no_red+=no_red_stat[lvl][i];
        sum_red+=red_stat[lvl][i];

      }

      /*for levels*/
      pdf_no_red[lvl].resize(no_red_stat[lvl].size());
      pdf_red[lvl].resize(no_red_stat[lvl].size());

      cdf_no_red[lvl].resize(no_red_stat[lvl].size());
      cdf_red[lvl].resize(no_red_stat[lvl].size());


      float accum_cdf_no_red = 0;
      float accum_cdf_red = 0;



      for (unsigned int i=0; i<no_red_stat[lvl].size(); ++i)
      {
        /*for levels*/
        pdf_no_red[lvl][i] = (float)no_red_stat[lvl][i]/ (float)sum_no_red;
        pdf_red[lvl][i] = (float)red_stat[lvl][i] / (float) sum_red;

        accum_cdf_no_red+=no_red_stat[lvl][i];
        cdf_no_red[lvl][i] = accum_cdf_no_red/sum_no_red;

        accum_cdf_red+=red_stat[lvl][i];
        cdf_red[lvl][i] = accum_cdf_red/sum_red;

      }

    }

    for (unsigned int lvl=0; lvl<SimHash::num_of_levels; ++lvl)
    {
      myfile << std::setw(45) << "lvl = " << lvl;
    }
    myfile << endl;

    for (unsigned int i=0; i<no_red_stat[0].size(); ++i)
    {
        myfile << " | " <<std::setw(3) << i <<"  | ";
        for (unsigned int lvl=0; lvl<SimHash::num_of_levels; ++lvl)
        {
          myfile<< std::setw(10) << no_red_stat[lvl][i] << " | " << std::setw(13) <<
          pdf_no_red[lvl][i]  << " ||" <<std::setw(7) <<
          red_stat[lvl][i] << " | " << std::setw(13) <<
          pdf_red[lvl][i]  << " |||";
        }
        myfile << endl;
    }

    myfile.close();

}

/* only for k>2 */
template<class SieveTraits>
inline void GaussSieveStatistics<SieveTraits,false>::compute_and_print_statistics_all_innloop(std::ostream &of)
{
  using std::endl;
  std::ofstream myfile;
  myfile.open ("newsieve/Statistics/Inner loop Statistics for k="+std::to_string(sieveptr->sieve_k) + "_dim="+std::to_string(sieveptr->lattice_rank));
  myfile << "Inner loop Statistics for dim = " << sieveptr->lattice_rank << endl;

  std::array<std::vector<float>, SimHash::num_of_levels> pdf_no_red_innloop;
  std::array<std::vector<float>, SimHash::num_of_levels> pdf_red_innloop;


  for (unsigned int lvl=0; lvl<SimHash::num_of_levels; ++lvl)
  {

    unsigned long long sum_no_red = 0;
    unsigned long long sum_red = 0;



    for (unsigned int i=0; i<no_red_stat_innloop[lvl].size(); ++i)
    {
      sum_no_red+=no_red_stat_innloop[lvl][i];
      sum_red+=red_stat_innloop[lvl][i];

    }
    pdf_no_red_innloop[lvl].resize(no_red_stat_innloop[lvl].size());
    pdf_red_innloop[lvl].resize(red_stat_innloop[lvl].size());

    for (unsigned int i=0; i<no_red_stat[lvl].size(); ++i)
    {

      pdf_no_red_innloop[lvl][i] = (float)no_red_stat_innloop[lvl][i]/ (float)sum_no_red;
      pdf_red_innloop[lvl][i] = (float)red_stat_innloop[lvl][i] / (float) sum_red;

    }

  }

  for (unsigned int lvl=0; lvl<SimHash::num_of_levels; ++lvl)
  {
    myfile << std::setw(45) << "lvl = " << lvl;
  }
  myfile << endl;

  for (unsigned int i=0; i<no_red_stat_innloop[0].size(); ++i)
  {
    myfile << " | " <<std::setw(3) << i <<"  | ";
    for (unsigned int lvl=0; lvl<SimHash::num_of_levels; ++lvl)
    {
      myfile<< std::setw(10) << no_red_stat_innloop[lvl][i] << " | " << std::setw(13) <<
      pdf_no_red_innloop[lvl][i]  << " ||" <<std::setw(7) <<
      red_stat_innloop[lvl][i] << " | " << std::setw(13) <<
      pdf_red_innloop[lvl][i]  << " |||";
    }
    myfile << endl;
  }

  myfile.close();

}


#endif

template<class SieveTraits>
inline void GaussSieveStatistics<SieveTraits,false>::dump_status_to_stream(std::ostream &of, int howverb)
{
  using std::endl;
  if(howverb >=1 ) of << "Number of collisions=" << number_of_collisions << endl;
  if(howverb >=1 ) of << "Number of points Sampled=" << number_of_points_sampled << endl;
  if(howverb >=1 ) of << "Number of points Constructed=" << number_of_points_constructed << endl;
  // if(howverb >=1 ) of << "Number of approx. scalar products=" << number_of_scprods << endl;

  if(howverb >= 1) of << "Number of exact scalar products=" << number_of_exact_scprods << endl;
  if(howverb >= 1) of << "Number of scalar products level 1=" << number_of_scprods_level1 << endl;
  if(howverb >= 1) of << "Number of scalar products level 2=" << number_of_scprods_level2 << endl;
  if(howverb >= 1) of << "Number of scalar products level 3=" << number_of_scprods_level3 << endl;

  if(howverb >= 1) of << "Number of 2-reductions inside 3-sieve=" << number_of_2reds << endl;
  if(howverb >= 1) of << "Number of 3-reductions inside 3-sieve=" << number_of_3reds << endl;

  if(howverb >= 1) of << "Number of mispredictions=" << number_of_mispredictions << endl;
  if(howverb >= 1) of << "Final List Size=" << get_current_list_size() << endl;
  if(howverb >= 1) of << "Final Queue Size="<< get_current_queue_size()<< endl;

// clang-format off
  //ONLY TO TEST BITAPPROX. TO BE DELETED
    #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX
    if(howverb >= 1)
    {
                   of << "No reduction: ";
                   for (unsigned int i=0; i!=no_red_stat.size(); ++i) of <<no_red_stat[i] << " ";
                   of << endl;
                   #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_2ND_ORDER
                   of << "2nd order: ";
                   for (unsigned int i=0; i!=no_red_stat2.size(); ++i) of <<no_red_stat2[i] << " ";
                   #endif
                   of << endl;
    }
    if(howverb>=1)
    {
                   of << "Reduction: ";
                   for (unsigned int i=0; i!=red_stat.size(); ++i) of <<red_stat[i] << " ";
                   of << endl;
                   #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_2ND_ORDER
                   of << "2nd order: ";
                   for (unsigned int i=0; i!=red_stat2.size(); ++i) of <<red_stat2[i] << " ";
                   #endif
                   of << endl;
    }
    #endif
    #ifdef EXACT_LATTICE_POINT_HAS_BITAPPROX_FIXED
    /*
    if(howverb>=1)
    {
      of << "SIM-HASH No reduction: ";
      for (unsigned int i=0; i!=no_red_stat_sim_hash.size(); ++i) of <<no_red_stat_sim_hash[i] << " ";
      of << endl;
      of << "SIM-HASH Reduction: ";
      for (unsigned int i=0; i!=red_stat_sim_hash.size(); ++i) of <<red_stat_sim_hash[i] << " ";
      of << endl;
    }
    compute_statistics(of);
    compute_statistics_2nd_order(of);
    */
    compute_and_print_statistics_all(of);
    compute_and_print_statistics_all_innloop(of);
    #endif
// clang-format on
}

} // end namespace GaussSieve

#endif // include guard
