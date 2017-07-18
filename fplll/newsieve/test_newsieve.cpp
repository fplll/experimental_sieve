#define USE_REGULAR_QUEUE //comment out if you use priority-queue

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
#include <iostream>
#include <fstream>
#include <random>
#include "MyLatticePointClass.cpp"
#include "LatticePointConcept.h"
#include "PlainLatticePoint.h"



using namespace fplll;

// Z_NR can be either long, double, or mpz_t

template <class ZT> void test_run_sieve(int dim, std::ofstream &ofs)
{
    //TerminationConditions< Z_NR<ZT> > term_cond; //sets target-legnth as Mink. bound


    //lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);
    //Sieve<Z_NR<ZT>, false> Test_Queue (B);
    //Test_Queue.run_2_sieve();

    ZZ_mat<ZT> BTest;
    BTest.resize(dim, dim);
    //BTest.gen_trg(1.1);
    srand (1);
    BTest.gen_qary_prime(1, 10*dim);



    if (dim >= 60)
        bkz_reduction(BTest, 8, BKZ_DEFAULT, FT_DEFAULT, 0);
    else
        lll_reduction(BTest, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);

    Sieve<Z_NR<ZT>, false> Test_Queue (BTest);

    ofs << "sieve is run on B[0]" << BTest[0] << endl;

    auto start = std::chrono::high_resolution_clock::now();
    Test_Queue.run();
    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);

    Test_Queue.print_status(-1,ofs);


    ofs << " Time taken: " << microseconds.count()/1000000.0 << "sec" << endl;
}

template <class Z> void sample_gaussians(int number, double s, double center, double cutoff)
{
    std::mt19937_64 engine;
    for (int i=0; i<number;++i)
    {
        cout << GaussSieve::sample_z_gaussian<Z, std::mt19937_64>(s,center,engine, cutoff) << endl;
    }
}

template <class ET,int nfixed>
void TestMyLatticePointClass()
{
    Dimension<-1> dim (10);
    MyLatticePoint<Z_NR<mpz_t>, -1>::class_init(dim);
    ZZ_mat<mpz_t> B;
    B.resize(10, 10);
    //generates a lower-triangular matrix B; the argument determines (in a complicated way) the bit-size of entries
    //B.gen_trg(1.1);

    srand (1);
    //generates GM lattice

    B.gen_qary_prime(1, 15);


    MyLatticePoint<Z_NR<mpz_t>, -1> test_pointA =  MyLatticePoint<Z_NR<mpz_t>, -1>(B[0]);
    MyLatticePoint<Z_NR<mpz_t>, -1> test_pointB =  MyLatticePoint<Z_NR<mpz_t>, -1>(B[1]);
    test_pointA.write_to_stream(cout);
    test_pointB.write_to_stream(cout);


    //unsigned int dim = 10;

    MyLatticePoint<Z_NR<mpz_t>, -1> test_pointSum = add(test_pointA, test_pointB);

    //MyLatticePoint<Z_NR<mpz_t>, -1> test_pointSum1 = add(test_pointA, test_pointB, 10); --DOES NOT WORK. WHY DOES THE ABOVE WORK?


    test_pointSum.write_to_stream(cout);

    MyLatticePoint<Z_NR<mpz_t>, -1> test_pointSub = sub(test_pointA, test_pointB);
    test_pointSub.write_to_stream(cout);

    MyLatticePoint<Z_NR<mpz_t>, -1> test_pointA_neg = negative_of(test_pointA);
    test_pointA_neg.write_to_stream(cout);

    Z_NR<mpz_t> mult;
    mult = 2;
    MyLatticePoint<Z_NR<mpz_t>, -1> test_pointA_mult = scalar_mult(test_pointA, mult);
    test_pointA_mult.write_to_stream(cout);

    //B.gen_qary_prime(1, 15);

    MyLatticePoint<Z_NR<mpz_t>, -1> test_point =  MyLatticePoint<Z_NR<mpz_t>, -1>(B[0]);


    test_pointA.write_to_stream(cout);
    test_pointB.write_to_stream(cout);

//    static_assert(IsALatticePoint< MyLatticePoint< Z_NR<mpz_t> ,-1> >::value," XXX");
//    cout << test_pointA;

    MyLatticePoint<Z_NR<mpz_t>, -1> test_pointA1;
    MyLatticePoint<Z_NR<mpz_t>, -1> test_pointA2;
    
    Z_NR<mpz_t> sc_prod = compute_sc_product(test_pointA, test_pointB);
    cout << sc_prod << endl;
    
    Z_NR<mpz_t> sc_prod1 = compute_sc_product(test_pointA1, test_pointA2);
    cout << "sc_prod1 " << sc_prod1 << endl;


    Z_NR<mpz_t> sc_prod_target;
    sc_prod_target = 5000;
    cout << compare_sc_product(test_pointA, test_pointB, sc_prod_target) << endl;
    cout << compare_abs_sc_product(test_pointA, test_pointB, sc_prod_target) << endl;

    cout << test_pointA.get_norm2() << endl;

    test_pointA.fill_with_zero();
    test_pointA.write_to_stream(cout);


}

//template <class ZT>
void test_2sieve()
{
    const int dim = 20;
    ZZ_mat<mpz_t> BTest;
    BTest.resize(dim, dim);
    //BTest.gen_trg(1.1);
    srand (1);
    BTest.gen_qary_prime(1, 10*dim);
    
    if (dim >= 60)
        bkz_reduction(BTest, 8, BKZ_DEFAULT, FT_DEFAULT, 0);
    else
        lll_reduction(BTest, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);
    
    Sieve<Z_NR<mpz_t>, false> Test_2Sieve(BTest);
    
    bool constexpr multithreaded = false;
    TerminationCondition<Z_NR<mpz_t>,multithreaded> * termcond = new MinkowskiTerminationCondition<Z_NR<mpz_t>, multithreaded>;
    Test_2Sieve.set_termination_condition(termcond);
    //Test_2Sieve.run();
}


int main(int argc, char **argv)
{
//sample_gaussians<long>(50, 10.0, 0.3, 4.0);
//sample_gaussians<long>(50, 0.0000001, 0.48, 4.0); //should still be fast.

    PlainLatticePoint<Z_NR<mpz_t>,-1> X;

    TestMyLatticePointClass<Z_NR<mpz_t>, -1>();
    
    test_2sieve();


//        //int dim[] = {52, 54, 56, 58, 60, 62, 64};
//        int dim = 62;
//        int length = 7;
//
//    	#ifdef USE_REGULAR_QUEUE
//        std::ofstream ofs("test_sieve_PQ_dim" +to_string(dim) + ".txt");
//		ofs << "WITH PRIORITY QUEUE" << endl;
//	#else
//		std::ofstream ofs("test_sieve_dim"+to_string(dim) + ".txt");
//		ofs << "WITH STANDARD QUEUE" << endl;
//	#endif
//	for (int i=0; i<length; i++) {
//
//
//		ofs << "start sieve on lattice of dim = " << dim[i] << endl;
//        	test_run_sieve<mpz_t>(dim[i], ofs);
//		ofs << "----------------------------------------" << endl;
//	}

//    for (int i=0; i<1; i++)
//    {
//        ofs << "start sieve on lattice of dim =  " << dim << endl;
//        test_run_sieve<mpz_t>(dim, ofs);
//        ofs << "----------------------------------------" << endl;
//    }
//   ofs.close();
  return 1;
}
