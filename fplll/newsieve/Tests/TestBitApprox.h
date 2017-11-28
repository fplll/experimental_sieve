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

  int constexpr dim = 67;
  
  typedef GaussSieve::ExactLatticePoint<long, -1> LP;
  using GaussSieve::MaybeFixed;
  
  std::array<long,dim> arr;
  std::array<long,dim> arr2;
  for(int i=0;i<dim;++i)
  {
    arr[i] = std::pow(-1, i+1) * i;
    arr2[i] = std::pow(-1, i) * (i-1)+13;
  }
  
  GaussSieve::StaticInitializer<LP> init1 (MaybeFixed<-1>{dim});
  LP latp = GaussSieve::make_from_any_vector<LP>(arr,MaybeFixed<-1>{dim});
  LP latp2 = GaussSieve::make_from_any_vector<LP>(arr2,MaybeFixed<-1>{dim});
  
  std::cout << "latp = " << latp << std::endl;
  std::cout << "latp2 = " << latp2 << std::endl;
  
  
  for (int i = 0; i<GaussSieve::SimHash::num_of_levels; ++i)
  {
    GaussSieve::SimHash::BitApproxScalarProduct approx_scprod_res = compute_sc_product_bitapprox_level(latp, latp2, i);
    std::cout << static_cast<uint_fast32_t>(approx_scprod_res) <<std::endl;
  }
  
  int constexpr dim_var = 70;
  std::array<long,dim_var> arr3;
  
  for(int i=0;i<dim_var;++i)
  {
    arr3[i] = std::pow(-1, i+1) * i - (i/2);
  }
  
  //typedef GaussSieve::ExactLatticePoint<long, 70> LPFixed;
  //GaussSieve::StaticInitializer<LPFixed> init2 (MaybeFixed<70>{dim_var});
  
  //LPFixed latp3 = GaussSieve::make_from_any_vector<LPFixed>(arr3,MaybeFixed<70>{dim_var});
  //std::cout << "latp3 = " << latp3 << std::endl;
  
  
  
  
  //LP latp2 = GaussSieve::make_from_any_vector<LP>(arr2,MaybeFixed<dim>{dim});
  
  //std::cout << "exact point with approximations: " <<  latp << std::endl;

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
