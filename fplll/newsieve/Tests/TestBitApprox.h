#ifndef TEST_BITAPPROX_H
#define TEST_BITAPPROX_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "../BitApproximation.h"
#include <limits>
#include "../ExactLatticePoint.h"
#include "gmpxx.h"
//#include <math.h> 

bool test_bit_approx()
{
  int constexpr dim = 20;
  int constexpr dimfixed=20;
  
  using LP = GaussSieve::ExactLatticePoint<mpz_class, dimfixed>;
  using GaussSieve::MaybeFixed;
  using GaussSieve::BitApproximation;
  
  std::array<mpz_class,dim> arr;

  for(int i=0;i<dim;++i)
  {
    arr[i] = std::pow(-1, i) * i;
    std::cout << arr[i] << " ";
  }
  

  LP latp = GaussSieve::make_from_any_vector<LP>(arr,MaybeFixed<dimfixed>{dim});

  BitApproximation<dim> approx(latp);
  
  
  std::cout << "approximated point = :" << approx; 
  
  return true;
}


#endif