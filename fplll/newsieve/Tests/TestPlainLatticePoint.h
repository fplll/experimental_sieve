#include <type_traits>
#include "../PlainLatticePoint.h"

#include "fplll/defs.h"
#include "fplll/nr/nr.h"
#include "vector"

bool test_plain_LP()
{
  using GaussSieve::PlainLatticePoint;
  using GaussSieve::Dimension;

  typedef GaussSieve::PlainLatticePoint<fplll::Z_NR<long>, -1> PLPvar1;
  typedef GaussSieve::PlainLatticePoint<fplll::Z_NR<long>, 10> PLPfix1;
  static_assert(GaussSieve::IsALatticePoint<PLPvar1>::value,"");
  static_assert(GaussSieve::IsCooVector<PLPvar1>::value,"");
  PLPvar1::class_init(Dimension<-1>{10});
  PLPfix1::class_init(Dimension<10>{10});
  PLPvar1 X1;
  PLPvar1 X2(10);
  PLPvar1 X3(Dimension<-1>(10));

  PLPfix1 Y1;
  PLPfix1 Y2(10);
  PLPfix1 Y3(Dimension<10>(10));

  std::cout << X1 << X2 << X3 << std::endl;
  std::cout << Y1 << Y2 << Y3 << std::endl;
  return true;
};
