#define USE_REGULAR_QUEUE //comment out to not use priority-queue
/*
  This provides an implementation of Gauss sieving, including using
  tuples of vectors in fplll framework. The Gauss Sieve code is
  based on Panagiotis Voulgaris's implementation of the Gauss sieve.
*/

#include "GaussQueue.h"
#include "Sampler.h"
#include "sieve_main.h"
#include "fplll.h"
//#include "LatticePoint.h"
//#include "PointList.h"
#include <thread>
#include <chrono>
#include "SieveGauss.cpp"
//#include "LatticePoint2.h"
#include "TermCondNew.h"
//#include "FilteredPoint.h"

//Don't do this:
using namespace fplll;
//using namespace LatticeApproximations;

/**
 * help function
 */
//static void main_usage(char *myself)
//{
//  cout << "Usage: " << myself << " [options]\n"
//       << "List of options:\n"
//       << "  -a [2|3|4]\n"
//       << "     2- or 3- or 4-sieve;\n"
//       << "  -f filename\n"
//       << "     Input filename\n"
//       << "  -r nnn\n"
//       << "     Generate a random instance of dimension nnn\n"
//       << "  -t nnn\n"
//       << "     Targeted norm^2=nnn\n"
//       << "  -s nnn\n"
//       << "     Using seed=nnn\n"
//       << "  -b nnn\n"
//       << "     BKZ preprocessing of blocksize=nnn\n"
//       << "  -v\n"
//       << "     Verbose mode\n";
//  exit(0);
//}

/**
 * run sieve
 */

#if 0
template <class ZT> int main_run_sieve(ZZ_mat<ZT> B, Z_NR<ZT> goal_norm, int alg, int ver, int seed)
{
  /*
   GaussSieve<ZT, FP_NR<double>> gsieve(B, alg, ver, seed);
  gsieve.set_goal_norm2(goal_norm);
  if (gsieve.alg == 3)
    gsieve.run_3sieve();
  else if (gsieve.alg == 4)
    gsieve.run_4sieve();
  else
    gsieve.run_2sieve();
   */
  return 0;
}
#endif

/**
 * main function
 */


/*
template<class DT>
void ListTester(ListMultiThreaded<DT> * const Z, GarbageBin<DT> * const gb, int id,int verbose)
{
for (int i=0; i <25; ++i)
{
    int count =0;
    int insertions=0;
    int deletions =0;
    for(MTListIterator<DT> it = Z->begin(); !it.is_end(); ++it)
    {
        ++count;
        DT tmp = *it;
        if(( tmp * tmp + i + 33*id) % 39 == 11)
        {
        Z->unlink(it,*gb);
        if (verbose>=2) cout << "Deleted";
        ++deletions;
        }
        if(( tmp * tmp + i + 18*id) %57 ==13  )
        {
        Z->insert(it, 1000000*id + 100* i + ((tmp * tmp )%93) );
        if (verbose>=2) cout << "Inserted";
        ++insertions;
        }
    }
    if (verbose>=1)
    {
    cout << "Thread " << id <<", iteration " << i << "count: " <<count <<" deletes: "<<deletions << "inserts: "<< insertions<<endl;
    }
    MTListIterator<DT> it = Z->end();
    Z->insert(it, 1000000*id + 100*i + 1);
}
return;
}

*/



int main(int argc, char **argv)
{
    char *target_norm_string = NULL;
    int opt, dim = 10;
    Z_NR<mpz_t> target_norm;
    target_norm = 0;
    int k=2;
    int bs=2;

    if (argc == 1)
    {
        cout << "Please provide the dimension." << endl;
        return -1;
    }

    while ((opt = getopt(argc, argv, "d:t:k:")) != -1) {
        switch (opt) {
            case 'd':
                dim = atoi(optarg);
                break;
            case 't':
                target_norm_string = optarg;
                break;
            case 'k':
                k=atoi(optarg);
            case 'b':
                bs = atoi(optarg);
            break;
        }
    }

    if (target_norm_string!=NULL)
    {
        target_norm.set_str(target_norm_string);
    }
    if(target_norm > 0)
    {
        cout << "target norm set: " << target_norm << endl;
    }

    // ZZ_mat is an integer row-oriented matrix. See /nr/matrix.h
    ZZ_mat<mpz_t> B;

    B.resize(dim, dim);

    //generates a lower-triangular matrix B; the argument determines (in a complicated way) the bit-size of entries
    //B.gen_trg(1.1);

    srand (1);
    //generates GM lattice
    B.gen_qary_prime(1, 10*dim);

    //KleinSampler<ZT, F> is templated by two classes; returns NumVect<Z_NR<ZT> of dim = B.NumCols()
    //    KleinSampler<mpz_t, FP_NR<double>> *Sampler = new KleinSampler<mpz_t, FP_NR<double>>(B, 0, 234234);
    //    NumVect<Z_NR<mpz_t> > sample(dim);
    //    sample = Sampler->sample();
    //    cout << sample << endl;

    /* preprocessing of basis */
    clock_t stime = clock();
    if (bs > 0)
        bkz_reduction(B, bs, BKZ_DEFAULT, FT_DEFAULT, 0);
    else
        lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);

    clock_t etime = clock();
    double secs   = (etime - stime) / (double)CLOCKS_PER_SEC;

    if (bs > 0)
        cout << "# [info] BKZ took time " << secs << " s" << endl;
    else
        cout << "# [info] LLL took time " << secs << " s" << endl;

    //ApproxLatticePoint<Z_NR<mpz_t>,false,-1> X ( conv_matrixrow_to_lattice_point(B[0]) );

    //ApproxLatticePoint<Z_NR<mpz_t>,false,-1> p ( conv_matrixrow_to_lattice_point(B[1]) );
    //int32_t inner_prod = LatticeApproximations::compute_sc_prod(X.get_approx(), p.get_approx(), dim);

    //cout << inner_prod << endl;

    //cout << "A point from a filtered list: " << l.getApproxVector() << " sc_prod: " << l.get_sc_prod() << endl;

    //cout << X;
    bool constexpr multithreaded = false;

    #ifndef USE_REGULAR_QUEUE
        cout << "Use Priority Queue" << endl;
    #else
        cout << "Use Standard Queue" << endl;
    #endif
    //cout << "run sieve on B[0] = " << B[0] << endl;
    //cout << "B[1] = " << B[1] << endl;
    auto start = std::chrono::high_resolution_clock::now();


	Sieve<Z_NR< mpz_t > , multithreaded> Test_2Sieve (B);
	TerminationCondition<Z_NR<mpz_t>,multithreaded> * termcond = new MinkowskiTerminationCondition<Z_NR<mpz_t>, multithreaded>;
	Test_2Sieve.set_termination_condition(termcond);
	if(target_norm!=0)
    {
        delete termcond;
        termcond = new LengthTerminationCondition<Z_NR<mpz_t>,multithreaded>(target_norm);
        Test_2Sieve.set_termination_condition(termcond);
        cout << "Setting target norm2 to" << target_norm << endl << flush;
    }


    //Test_2Sieve.run();
	//cout << "sv is " << endl;
    //Test_2Sieve.get_shortest_vector_found().printLatticePoint();
    //Test_2Sieve.print_status();




    Sieve<Z_NR<mpz_t>, false> Test_3Sieve(B, k);
    TerminationCondition<Z_NR<mpz_t>,false> * termcond2 = new MinkowskiTerminationCondition<Z_NR<mpz_t>, false>;
    Test_3Sieve.set_termination_condition(termcond2); //if this is commented out, the next line segfaults
    Test_3Sieve.run();


    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
    cout << " Time taken: " << microseconds.count()/1000000.0 << "sec" << endl;
    delete termcond;
#if 0
#if 0
  dot_time = 0;
  dot_num = 0;
  count_bad = 0;
#endif

  alg = 2;

  /* parse */
  if (argc == 1)
  {
    main_usage(argv[0]);
    return -1;
  }
  while ((option = getopt(argc, argv, "a:f:r:t:s:b:v")) != -1)
  {
    switch (option)
    {
    case 'a':
      alg = atoi(optarg);
      if (alg != 2 && alg != 3 && alg != 4)
      {
        cout << " Error, only support 2-, 3- and 4-sieve" << endl;
        exit(1);
      }
    case 'f':
      input_file_name = optarg;
      flag_file       = true;
      break;
    case 'r':
      dim       = atoi(optarg);
      flag_file = false;
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'b':
      bs = atoi(optarg);
      break;
    case 'v':
      flag_verbose = true;
      break;
    case 't':
      // ngoal_norm = atol(optarg);
      cout << optarg << endl;
      goal_norm_s = optarg;
      break;
    case 'h':
      main_usage(argv[0]);
      return -1;
    case '?':
      main_usage(argv[0]);
      return -1;
    case ':':
      main_usage(argv[0]);
      return -1;
    }
  }

  /* set lattice */
  ZZ_mat<mpz_t> B;
  if (flag_file)
  {
    ifstream input_file(input_file_name);
    if (input_file.is_open())
    {
      input_file >> B;
      input_file.close();
    }
    else
    {
      cin >> B;
    }
    if (flag_verbose)
    {
      cout << "# [info] reading lattice of dimension " << B.get_rows() << "x" << B.get_cols()
           << endl;
    }
  }
  else
  {
    if (flag_verbose)
    {
      cout << "# [info] generating random lattice of dimension " << dim << endl;
    }
    srand(time(NULL));
    B.resize(dim, dim);
    //B.gen_trg(1.1);
    BTest.gen_qary_prime(1, 10*dim);
  }

  /* set targeted norm */
  Z_NR<mpz_t> goal_norm, max;
  if (goal_norm_s != NULL)
  {
    goal_norm.set_str(goal_norm_s);
  }
  if (goal_norm < 0)
    goal_norm = 0;
  if (flag_verbose)
    cout << "# [info] goal norm^2 is " << goal_norm << endl;

  /* preprocessing of basis */
  clock_t stime = clock();
  if (bs > 0)
    bkz_reduction(B, bs, BKZ_DEFAULT, FT_DEFAULT, 0);
  else
      lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);

  clock_t etime = clock();
  double secs   = (etime - stime) / (double)CLOCKS_PER_SEC;
  if (flag_verbose)
  {
    if (bs > 0)
      cout << "# [info] BKZ took time " << secs << " s" << endl;
    else
      cout << "# [info] LLL took time " << secs << " s" << endl;
  }
  // cout << B << endl;

  /* decide integer type */
  stime = clock();
  max   = B.get_max();

#if 1
  if (max < std::numeric_limits<int>::max())
  {
    long goal_norm_l = abs(goal_norm.get_si());
    Z_NR<long> goal_norm_lt;
    goal_norm_lt = goal_norm_l;
    ZZ_mat<long> B2(B.get_rows(), B.get_cols());
    for (int i = 0; i < B.get_rows(); i++)
      for (int j = 0; j < B.get_cols(); j++)
        B2(i, j) = B(i, j).get_si();
    main_run_sieve<long>(B2, goal_norm_lt, alg, flag_verbose, seed);
  }
  else
#endif
 //   main_run_sieve<mpz_t>(B, goal_norm, alg, flag_verbose, seed); //commented out for test

  etime = clock();
  secs  = (etime - stime) / (double)CLOCKS_PER_SEC;
  if (flag_verbose)
  {
    cout << "# [info] sieve took time " << secs << " s" << endl;
/* dot product time */
#if 0
    cout << "# [info] dot_time " << dot_time << endl;
    cout << "# [info] dot_num " << dot_num << endl;
    cout << "# [info] dot_time/dot_number " << (double) dot_time/dot_num << endl;
#endif
  }
#endif
  return 1;
}
