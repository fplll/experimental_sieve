#ifndef TEST_LIST_H
#define TEST_LIST_H

#include "../PointListNew.h"
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"
#include <iostream>

bool test_list()
{
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
