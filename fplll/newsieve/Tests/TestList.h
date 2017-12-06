#ifndef TEST_LIST_H
#define TEST_LIST_H

//#include "../PointListNew.h"
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"
#include <iostream>
#include "../GaussListBitapprox.h"

bool test_list()
{
  using Traits = GaussSieve::DefaultSieveTraits<mpz_class, false, -1>;
  using List = GaussSieve::GaussListWithBitApprox<Traits,false>;
  using StoredPoint = typename Traits::GaussList_StoredPoint;
  using ReturnType  = typename Traits::GaussList_ReturnType;
  unsigned int constexpr dim = 30;
  typename Traits::GlobalStaticDataInitializer init_arg (dim);

  List main_list(init_arg,1);

  StoredPoint point;
  std::array<mpz_class,dim> v;
  for(unsigned int i = 0; i < dim; ++i) { v[i] = (1-2*(i%2))*i*i*i; }
  point = GaussSieve::make_from_any_vector<StoredPoint>(v,dim);
  auto it = main_list.cend();
  assert(main_list.empty());
  main_list.insert_before(it, std::move(point));
  assert(main_list.size() == 1);
  it = main_list.cbegin();
  std::cout << *it << std::endl;
  std::cout << it.access_bitapproximation(0) << std::endl;
  std::cout << it.get_approx_norm2() << std::endl;

  std::cout << "No tests yet for list." << std::endl;
  return true;
}


// Old Code:

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




#endif
