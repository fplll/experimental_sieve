#ifndef EXACT_LATTICE_POINT_CPP
#define EXACT_LATTICE_POINT_CPP

#include "ExactLatticePoint.h"

/*  header file continuation from ExactLatticePoint.h
    This file contains the implementations
*/

/*
template<class ET, int nfixed>
ExactLatticePoint<ET,nfixed>::ExactLatticePoint(int n) //creates vector of length n
    :
    NumVect<ET>(n), //creates all-zero vector of length n
    norm2()
{
    this->data.resize(n);
    this->fill(0);
    norm2 = 0;
}
*/

template<class ET,int nfixed>
ExactLatticePoint<ET,nfixed>::ExactLatticePoint(int n, long fillwith) : NumVect<ET>(n),norm2() // for debugging
{
    this->data.resize(n);
    this->fill(fillwith);
    normalize();
    //norm2 = fillwith*fillwith*n; gives type problems.
}


template <class ET,int nfixed>
ExactLatticePoint<ET> operator+ (ExactLatticePoint<ET,nfixed> const &A, ExactLatticePoint<ET,nfixed> const &B)
{

	ExactLatticePoint<ET,nfixed> C(A);
	//length-check is done in by add in numvect
	C.NumVect<ET>::add(B, A.size());
	C.Normalize();
	return C;

}

template <class ET, int nfixed>
ExactLatticePoint<ET,nfixed> operator- (ExactLatticePoint<ET,nfixed> const &A, ExactLatticePoint<ET> const &B)
{

	ExactLatticePoint<ET,nfixed> C(A);
	C.NumVect<ET>::sub(B, A.size());
	C.Normalize();
	return C;
}

template <class ET,int nfixed> ExactLatticePoint<ET,nfixed> operator- (ExactLatticePoint<ET,nfixed> const &A) //unary-
{
    ExactLatticePoint<ET,nfixed> ans(A);
    unsigned int const n = A.get_dim();
    for(unsigned int i=0;i<n;++i)
    {
    ans[i].neg();
    }
    return ans;
}

template <class ET,int nfixed>
void scalar_mult (ExactLatticePoint<ET,nfixed> &A, ET const & multiple)
{
    A.NumVect<ET>::mul(A, A.size(), multiple);
    A.Normalize(); //we could do smarter.
}


//Simple dot_product
//template <class ET,int nfixed>
//void compute_exact_sc_product (ET &result, const ExactLatticePoint<ET,nfixed> &p1, const ExactLatticePoint<ET,nfixed> &p2)
//{
//   dot_product(result, p1, p2);
//}

template<class ET, int nfixed>
ET exact_scalar_product(ExactLatticePoint<ET,nfixed> const &p1, ExactLatticePoint<ET,nfixed> const & p2)
{
    ET res;
    dot_product(res,p1,p2);
    return res;
}

template <class ET,int nfixed>
void ExactLatticePoint<ET,nfixed>::normalize() //sets norm2 to the correct value. Use after using operations from the NumVect<ET> base class.
{
    this->norm2 = exact_scalar_product(*this,*this);
};

//Convert MatrixRow to LatticePoint
template <class ET,int nfixed> ExactLatticePoint<ET,nfixed> conv_matrixrow_to_lattice_point (MatrixRow<ET> const &row)
{
	ExactLatticePoint<ET,nfixed> res(row.get_underlying_row());
	//NumVect<ET> tmp(row.get_underlying_row());
	return res;
}

#endif
