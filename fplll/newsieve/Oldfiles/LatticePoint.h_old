//
//  LatticePoint.h
//

#ifndef LATTICE_VECTOR_CLASS_H
#define LATTICE_VECTOR_CLASS_H

#include "sieve_common.h" //needed (at least for convertions from MatrixRow (the header has to be revised);

/**
 * Class for list element

    data types:

    class Z_NR: stores integers; defined in type_nr_Z.h
    class NumVect: uses STL's container of type <vector>; defined in type_numvec.h


    members:
    v: instance of NumVect
    norm: l^2-norm of v

**/

template <class ET> class LatticePoint;

template <class ET>
void sc_product (ET &result, const LatticePoint<ET> &p1, const LatticePoint<ET> &p2);



template<class ET> //ET: entries of the individual vectors. Required to be copy-constructible. Use E = Z_NR<mpz_t> rather than E=mpz_t.
class LatticePoint : public NumVect<ET>
{
       using NV = NumVect<ET>;

    /* square L2 norm of the vector */
public:
        ET norm2;

public:
    LatticePoint()=default;
    LatticePoint(const LatticePoint &Point) = default; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    LatticePoint(LatticePoint &&Point) = default ; // : NumVect<ET> (), norm2(0) {swap(Point.NV::data), swap(Point.norm2);}
    LatticePoint(int n); //creates all-zero vector of length n
    LatticePoint(int n, long fillwith);    //TODO: fillwith should be ET, not long.
    LatticePoint(NumVect<ET> vector_) : NumVect<ET>(vector_)
    {
        Normalize();
    }
    void Normalize(); //sets norm2 to the correct value. Use after using operations from the NumVect<ET> base class.

    LatticePoint& operator=(LatticePoint const &that) =default;
    LatticePoint& operator=(LatticePoint && that) =default;
    /*
    LatticePoint& operator=(LatticePoint that)
    {
        swap(*this,that);
        return *this;
    }
    */
    ~LatticePoint() {}

    friend void swap(LatticePoint &A, LatticePoint &B)
    {
        std::swap(A.data, B.data);
        std::swap(A.norm2,B.norm2);
    }

    inline NV& getVector() const {return this->data.get();}
    inline ET get_norm2() const {return norm2;}
    //inline void get_norm2 (ET &norm_to_return) const {norm_to_return = norm2;}
    inline void setNorm2 (ET norm) {this->norm2 = norm;} //should not be required, actually.

    void printLatticePoint(ostream & os = cout) const
    {
        os << * (static_cast<NumVect<ET> const *>(this)) << " of norm: " << this->norm2 << endl;
    }
    bool operator< (LatticePoint<ET> const &other ) const {return (this->norm2 < other.norm2);};
//    void subtr (LatticePoint<ET> const A)
//    {
//        this->this.sub(A, A.size());
//    }
};

template<class ET>
LatticePoint<ET>::LatticePoint(int n)
    :
    NumVect<ET>(n), //creates all-zero vector of length n
    norm2()
{
    this->data.resize(n);
    this->fill(0);
    norm2 = 0; //why? LatticePoint(int n) : NumVect<ET>(n), norm(0) errs
}

template<class ET>
LatticePoint<ET>::LatticePoint(int n, long fillwith) : NumVect<ET>(n),norm2() // for debugging
{
    this->data.resize(n);
    this->fill(fillwith);
    sc_product(this->norm2, *this, *this);
    //norm2 = fillwith*fillwith*n; gives type problems.
}


template <class ET>
LatticePoint<ET> operator+ (LatticePoint<ET> const &A, LatticePoint<ET> const &B)
{

	LatticePoint<ET> C(A);
	//length-check is done in by add in numvect
	C.NumVect<ET>::add(B, A.size());
	C.Normalize();
	return C;

}

template <class ET>
LatticePoint<ET> operator- (LatticePoint<ET> const &A, LatticePoint<ET> const &B)
{

	LatticePoint<ET> C(A);
	C.NumVect<ET>::sub(B, A.size());
	C.Normalize();
	return C;

}

template <class ET>
void scalar_mult (LatticePoint<ET> &A, ET multiple)
{
    A.NumVect<ET>::mul(A, A.size(), multiple);
    A.Normalize(); //we could do smarter.
}


//Simple dot_product
template <class ET>
void sc_product (ET &result, const LatticePoint<ET> &p1, const LatticePoint<ET> &p2)
{
   dot_product(result, p1, p2);
}

//Convert MatrixRow to LatticePoint

template <class ET>
void LatticePoint<ET>::Normalize() //sets norm2 to the correct value. Use after using operations from the NumVect<ET> base class.
{
    sc_product(this->norm2, *this,*this);
};

//template <class ET>
//LatticePoint<ET> conv_matrixrow_to_lattice_point (MatrixRow<ET> const &row)
//{
//	LatticePoint<ET> res(row.get_underlying_row());
//	//NumVect<ET> tmp(row.get_underlying_row());
//	return res;
//}

// Convert sample() result NumVect to LatticePoint


template<class ET>
inline LatticePoint<ET> conv_sample_to_lattice_point (NumVect<ET> const &sample)
{
    LatticePoint<ET> p (sample);
    return p;
}

#endif
