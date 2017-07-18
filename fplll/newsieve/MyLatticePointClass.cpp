#ifndef MY_LATTICE_POINT_CLASS_CPP
#define MY_LATTICE_POINT_CLASS_CPP


#include "MyLatticePointClass.h"


template <class ET,int nfixed>
MyLatticePoint<ET, nfixed> add (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET,nfixed> const &B)
{
    MyLatticePoint<ET, nfixed> sum;

    //MyLatticePoint<ET, nfixed> sum = make_copy(A, auxdata);

    for (unsigned int i=0; i<MyLatticePoint<ET,nfixed>::dim; ++i)
    {
        //sum.data[i] = sum.data[i]+B.data[i];
        sum.data[i].add(A.data[i], B.data[i]); //IS CORRECT?
    }
    sum.update_norm2();
    return sum;
}

template <class ET, int nfixed>
MyLatticePoint<ET,nfixed> sub (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET, nfixed> const &B)
{
    MyLatticePoint<ET, nfixed> sum;

    for (unsigned int i=0; i<MyLatticePoint<ET,nfixed>::dim; ++i)
    {
        sum.data[i].sub(A.data[i], B.data[i]);
    }
    sum.update_norm2();
    return sum;
}

template <class ET,int nfixed> MyLatticePoint<ET,nfixed> negative_of(MyLatticePoint<ET,nfixed> const &A)
{
    ET zero;
    zero = 0;
    MyLatticePoint<ET, nfixed> neg;
    for (unsigned int i=0; i<MyLatticePoint<ET,nfixed>::dim; ++i)
    {
        neg.data[i].sub(zero, A.data[i]);
    }
    return neg;
}

template <class ET,int nfixed>

MyLatticePoint<ET,nfixed> scalar_mult (MyLatticePoint<ET,nfixed> const &A, ET const & multiple)

{
    MyLatticePoint<ET, nfixed> res;

    for (unsigned int i=0; i<MyLatticePoint<ET,nfixed>::dim; ++i)
    {
        res.data[i].mul(A.data[i], multiple);
    }
    res.update_norm2();
    return res;
}

/**
    Checks whether <A,B>  > target
    In our use-cases, target should be positive.
**/
template <class ET,int nfixed>
bool compare_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET const & target)
{
    ET sc_product = compute_sc_product(A, B);
    return (sc_product > target)? (true) : (false);
}


/**
    Checks whether |<A,B>| > target.
    Only really useful if target>=0
 **/
template <class ET,int nfixed>
bool compare_abs_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, ET const & target)
{
    ET sc_product = compute_sc_product(A, B);
    ET abs_sc_prod;
    abs_sc_prod.abs(sc_product);
    return  (abs_sc_prod>target) ? (true) : (false);
}

template <class ET,int nfixed>
ET compute_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B)
{
    ET res;
    res = 0;
    for (unsigned int i=0; i<MyLatticePoint<ET,nfixed>::dim; ++i)
    {
        res.addmul(A.data[i], B.data[i]);
    }
    return res;
}


template <class ET,int nfixed>
MyLatticePoint<ET, nfixed> MyLatticePoint<ET,nfixed>::make_copy() const
{
    MyLatticePoint tmp_copy;
    for(unsigned int i=0; i<MyLatticePoint<ET,nfixed>::dim; ++i)
    {
        tmp_copy.data[i] = data[i];
    }
    tmp_copy.norm2 = norm2;
    return tmp_copy;
}

/*
template <class ET,int nfixed> MyLatticePoint<ET, nfixed>
void print (std::ostream &os, MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata)
{

}
*/
#endif
