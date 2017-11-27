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
  
  using LP = GaussSieve::ExactLatticePoint<long, dim>;
  using GaussSieve::MaybeFixed;
  
  std::array<long,dim> arr;
  std::array<long,dim> arr2;
  for(int i=0;i<dim;++i)
  {
    arr[i] = std::pow(-1, i+1) * i;
    arr2[i] = std::pow(-1, i) * (i-1);
  }
  
  
  LP latp = GaussSieve::make_from_any_vector<LP>(arr,MaybeFixed<dim>{dim});
  
  //LP latp2 = GaussSieve::make_from_any_vector<LP>(arr2,MaybeFixed<dim>{dim});
  
  std::cout << "exact point with approximations: " <<  latp << std::endl;

  /*
 
  int constexpr dimfixed=20;

 
  
  using GaussSieve::BitApproximation;




  LP latp = GaussSieve::make_from_any_vector<LP>(arr,MaybeFixed<dimfixed>{dim});
  LP latp2 = GaussSieve::make_from_any_vector<LP>(arr2,MaybeFixed<dimfixed>{dim});

  std::cout << "exact point: " <<  latp << std::endl;
  BitApproximation<dim> approx(latp);
  std::cout << std::endl << "approximated point = :" << approx << std::endl;

  std::cout << "exact point: " << latp2 << std::endl;
  BitApproximation<dim> approx2(latp2);
  std::cout << std::endl << "approximated point = :" << approx2 << std::endl;

  int scprod = approximate_scalar_product (approx, approx2);

  std::cout << "approx <,> = " << scprod << std::endl;
  */
  return true;
}


#endif
