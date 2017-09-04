#ifndef TEST_MAYBE_FIXED_INT_H
#define TEST_MAYBE_FIXED_INT_H

#include <type_traits>
#include "assert.h"
#include "../SieveUtility.h"

bool test_maybe_fixed_int()
{
  using GaussSieve::MaybeFixed;
  MaybeFixed<-1> n1; n1=2;
  MaybeFixed<-1> n12{2};
  MaybeFixed<-1> n13(2);
  MaybeFixed<-1> n2 = 3;
  n13 = n2;
  MaybeFixed<3> n3;
  MaybeFixed<4> n4;
  MaybeFixed<3> n33;

  assert(n1==n12);
  assert(! (n1!= n12));
  assert(n1 != n2);
  assert(! (n1==n2));
  assert(n1 < n2);
  assert(n1 <=n2);
  assert(n2 > n1);
  assert(n2 >=n1);

  assert(n1 < n3);
  assert(n1 != n3);
  assert(n3 > n1);
  assert(n3 != n1);
  assert(!(n3==n1));
  assert(!(n1==n3));

  assert(n2 == n3);
  assert(n3 == n2);
  assert(!(n2!=n3));
  assert(!(n3!=n2));

// NOTE: n3 etc is not declared constexpr above.
// If we actually overload operator== etc., this causes gcc to raise an error due to that.

  static_assert(n3==n33);
  static_assert(n3< n4);
  static_assert(n3<=n4);
  static_assert(n4>=n3);
  static_assert(n4> n3);
  static_assert(!(n3==n4));
  static_assert(n3!=n4);

  return true;
}

#endif
