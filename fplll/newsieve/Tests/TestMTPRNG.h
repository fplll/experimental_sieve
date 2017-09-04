#ifndef TEST_MTPRNG_H
#define TEST_MTPRNG_H

#include "../MTPRNG.h"
#include <vector>
#include <iostream>
#include <random>

bool test_mtprng()
{

std::seed_seq sseq1({1,2,3});
std::seed_seq sseq2({1,2,3});
std::seed_seq sseq3({1,2,3});
std::seed_seq sseq4({1,2,3});
std::seed_seq sseq5({1,2,3});
std::seed_seq sseq6({1,2,3});


// Note: std::seed_seq is non-copyable and using it to initialize anything actually modifies it.

GaussSieve::MTPRNG< std::mt19937_64, true, std::seed_seq> mtprng(sseq1);
GaussSieve::MTPRNG< std::mt19937_64, false, std::seed_seq> stprng(sseq2);

mtprng.init(1);
stprng.init(1);

std::vector<uint_fast64_t> results1(30);
std::vector<uint_fast64_t> results2(30);

for(unsigned int i=0;i<10;++i)
{
  results1[i] = mtprng.rnd(0)();
  results2[i] = stprng.rnd(0)();
}
mtprng.init(2);

for(unsigned int i=10;i<20;++i)
{
  results1[i] = mtprng.rnd(0)();
  results2[i] = stprng.rnd(0)();
  mtprng.rnd(1)();
  mtprng.rnd(1)();
  mtprng.rnd(1)();
  mtprng.rnd(1)();
}

mtprng.init(0);
mtprng.init(1);

for(unsigned int i=20;i<30;++i)
{
  results1[i] = mtprng.rnd(0)();
  results2[i] = stprng.rnd(0)();
}

assert(results1==results2);



mtprng.reseed(sseq3);
for(unsigned int i=0;i<30;++i)
{
  results1[i]=mtprng.rnd(0)();
}
assert(results1==results2);

stprng.reseed(sseq4);

for(unsigned int i=0;i<30;++i)
{
  results2[i]=stprng.rnd(0)();
}
assert(results1==results2);


// This is just to see that
// a) the values look random
// b) that they actually depend on the seed

/*
std::seed_seq sseqother{};
mtprng.reseed(sseqother);

for(unsigned int i=0;i<30;++i)
{
  results1[i]=mtprng.rnd(0)();
}

*/


/*
for(unsigned int i=0;i<30;++i)
{
std::cout << results1[i] << " ";
}
std::cout << std::endl;
for(unsigned int i=0;i<30;++i)
{
std::cout << results2[i] << " ";
}
std::cout << std::endl;
*/


return true;
}


#endif
