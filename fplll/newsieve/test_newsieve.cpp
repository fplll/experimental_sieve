// clang-format off

#define USE_REGULAR_QUEUE //comment out if you use priority-queue

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
    
    char *target_norm_string = NULL;
    char* input_file_name = NULL;
    bool flag_file = false;
    
    int opt, dim = 10;
    int b = 2;
    int k = 3;
    Z_NR<mpz_t> target_norm;
    mpz_class target_norm_conv;
    target_norm = 0;
    
    if (argc == 1)
    {
        cout << "Please provide the dimension." << endl;
        return -1;
    }

    while ((opt = getopt(argc, argv, "d:t:k:f:b:")) != -1) {
        switch (opt) {
            case 'd':
                dim = atoi(optarg);
                break;
            case 't':
                target_norm_string = optarg;
                break;
            case 'f':
                input_file_name = optarg;
                flag_file = true;
                break;
            case 'k':
                k=atoi(optarg);
            case 'b':
                b = atoi(optarg);
            break;
        }
    }
    
    ZZ_mat<mpz_t> B;
    B.resize(dim, dim);
    
    if (flag_file) {
        ifstream input_file(input_file_name);
        if (input_file.is_open()) {
            cout << "reading B from file ..." << endl;
            input_file >> B;
            input_file.close();
        }
    }
    else {
        srand (1);
        //generates GM lattice
        B.gen_qary_prime(1, 10*dim);
    }
    
    if (target_norm_string!=NULL)
    {
        target_norm.set_str(target_norm_string);
        target_norm_conv = mpz_class(target_norm_string);
    }
    
    if(target_norm > 0)
    {
        cout << "target norm set: " << target_norm << endl;
    }
    
    /* preprocessing of basis */
    clock_t stime = clock();
    if (b > 2)
        bkz_reduction(B, b, BKZ_DEFAULT, FT_DEFAULT, 0);
    else
        lll_reduction(B, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);

    clock_t etime = clock();
    double secs   = (etime - stime) / (double)CLOCKS_PER_SEC;

    if (b > 2)
        cout << "# [info] BKZ took time " << secs << " s" << endl;
    else
        cout << "# [info] LLL took time " << secs << " s" << endl;
        
    bool constexpr multithreaded = false;
    using Traits = GaussSieve::DefaultSieveTraits<mpz_class, false, -1>;
    
    
    
    auto start = std::chrono::high_resolution_clock::now();
    

	Sieve<Traits, multithreaded> Test_3Sieve (B, k, 0);
	//TerminationCondition<Traits,multithreaded> * termcond = new MinkowskiTerminationCondition<Traits, multithreaded>;
    TerminationCondition<Traits,multithreaded> * termcond = new LengthTerminationCondition<Traits, multithreaded> (target_norm_conv);
	Test_3Sieve.set_termination_condition(termcond);
    
    
    Test_3Sieve.run();
    cout << "sv is " << endl;
    Test_3Sieve.print_status();
    
    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
    cout << " Time taken: " << microseconds.count()/1000000.0 << "sec" << endl;
    delete termcond;

}

//clang-format on
