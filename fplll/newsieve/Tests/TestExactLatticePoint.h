#ifndef TEST_EXACT_LATTICE_POINT_H
#define TEST_EXACT_LATTICE_POINT_H

#include <type_traits>
#include "../ExactLatticePoint.h"

#include "fplll/defs.h"
#include "fplll/nr/nr.h"
#include "vector"
#include "gmpxx.h"

bool test_exact_LP()
{
  using GaussSieve::ExactLatticePoint;
  using GaussSieve::MaybeFixed;

  typedef GaussSieve::ExactLatticePoint<long, -1> LPvar;
  typedef GaussSieve::ExactLatticePoint<long, 10> LPfix;
  typedef GaussSieve::ExactLatticePoint<mpz_class,10> LPGMP;
  static_assert(GaussSieve::IsALatticePoint<LPvar>::value,"");
  static_assert(GaussSieve::IsCooVector<LPvar>::value,"");
  LPvar::class_init(MaybeFixed<-1>{10});
  LPfix::class_init(MaybeFixed<10>{10});
  LPGMP::class_init(MaybeFixed<10>{10});
  LPvar X1;
  LPvar X2(10);
  LPvar X3(MaybeFixed<-1>(10));

  LPfix Y1;
  LPfix Y2(10);
  LPfix Y3(MaybeFixed<10>(10));

  std::vector<long> vec;
  vec.reserve(10);
  std::array<long,10> arr;
  long carr[10];
  long carr2[10];

  for(unsigned int i=0;i<10;++i)
  {
  vec[i] = i;
  arr[i] = i*i;
  carr[i]=10*i;
  carr2[i]=9;
  }

  X1 = GaussSieve::make_from_any_vector<LPvar>(vec, MaybeFixed<-1>(10));
  X2 = GaussSieve::make_from_any_vector<LPvar>(vec, MaybeFixed<-1>(10));
  X3 = GaussSieve::make_from_any_vector<LPvar>(vec, MaybeFixed<-1>(10));
  assert(X1==X2);
  assert(X1==X3);
  Y1 = GaussSieve::make_from_any_vector<LPfix>(vec, MaybeFixed<10>());
  Y2 = GaussSieve::make_from_any_vector<LPfix>(arr, MaybeFixed<10>());
  Y3 = GaussSieve::make_from_any_vector<LPfix>(carr, MaybeFixed<10>());
  assert(Y1==X1);
  assert(Y1!=Y2);
  assert(Y1!=Y3);
  assert(X1 <= X2);
  assert(X2 >= X1);
  assert(Y1 < Y2);
  assert(Y2 > Y1);
  assert(!(X1<X2));
  assert(X1.get_norm2() == 285);
  assert(Y3.get_norm2() == 28500);
  X2 = X1.make_copy();
  X3 = GaussSieve::make_from_any_vector<LPvar>(carr2, MaybeFixed<-1>(10));
  X2 = X3 - X1; // so X2 is now X1 in reverse.
  assert(X2.get_norm2() == 285);
  assert(X2+ X1 == X3);
  assert( compute_sc_product(X1,X3) == 9*45 );
  assert( X1.get_dim() == 10);
  assert( X1.get_vec_size() == 10);
  assert( Y1.get_dim() == 10);
  assert( Y1.get_vec_size() == 10);
  assert(X1.is_zero() == false);
  X1 = X1 - X1;
  assert(X1.is_zero() == true);
  X2.fill_with_zero();
  assert(X2.is_zero() == true);
  std::array<fplll::Z_NR<mpz_t>,10> mpzvec1, mpzvec2;
  LPGMP Z1;
  LPGMP Z2(10);
  LPGMP Z3( MaybeFixed<10>(10) );


  for(int i=0;i<10;++i)
  {
    mpzvec1[i] = i;
    mpzvec2[i] = 2*i;
  }

  Z1 = GaussSieve::make_from_znr_vector<LPGMP>(mpzvec1, MaybeFixed<10>());
  Z2 = Z1.make_copy();
  Z3 = GaussSieve::make_from_znr_vector<LPGMP>(mpzvec2, MaybeFixed<10>());
//  Z3 = GaussSieve::make_from_znr_vector<LPGMP>(mpzvec2, MaybeFixed<10>());
  assert(Z1 == Z2);

  assert(Z1+Z2 == Z3);

  assert(Z1 != Z3);
  assert(Z1.get_norm2() == 285);
  Z3 = Z1 - Z2;
  assert(Z3.is_zero());
  assert(!Z1.is_zero());
  Z1.fill_with_zero();
  assert(Z1.is_zero());

//  Y2 = make_from_any_vector(vec, MaybeFixed<-1>(10));
//  assert(Y1==Y2);
//  assert(X1==Y1);

  std::cout << X1 << X2 << X3 << std::endl;
  std::cout << Y1 << Y2 << Y3 << std::endl;
  std::cout << Z1 << Z2 << Z3 << std::endl << std::flush;
  return true;
};

#endif
