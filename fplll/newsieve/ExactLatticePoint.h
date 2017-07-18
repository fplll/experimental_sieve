#ifndef EXACT_LATTICE_POINT_H
#define EXACT_LATTICE_POINT_H

#include <iostream>

/*
    ExactLatticePoint stores an (exact) lattice point, together with its squared length.
    We also have some elementary arithmetic (addition, substraction etc. ) defined on them.
*/

template <class ET,int nfixed> ExactLatticePoint<ET,nfixed> conv_matrixrow_to_lattice_point (MatrixRow<ET> const &row);
template <class ET,int nfixed> std::istream & operator>> (std::istream & is, ExactLatticePoint<ET,nfixed> & exact_point);
template <class ET,int nfixed> std::ostream & operator<< (std::ostream & os, ExactLatticePoint<ET,nfixed> const & exact_point)
{
exact_point.print_exact_point(os); return os;
}


template<class ET,int nfixed> //ET: entries of the individual vectors. Required to be copy-constructible. Use E = Z_NR<mpz_t> rather than E=mpz_t.
class ExactLatticePoint : public NumVect<ET>
{
       using NV = NumVect<ET>;

    /* square L2 norm of the vector */
public:
        ET norm2;
public:
    ExactLatticePoint()=default;
    ExactLatticePoint(const ExactLatticePoint &Point) = default;
    ExactLatticePoint(ExactLatticePoint &&Point) = default ;
    //ExactLatticePoint(int n); //creates all-zero vector of length n
    ExactLatticePoint(int n, long fillwith);    //TODO: fillwith should be ET, not long.
    explicit ExactLatticePoint(NumVect<ET> const & vector_) : NumVect<ET>(vector_) {normalize();};          //Creates an exact lattice point from a NumVect
    explicit ExactLatticePoint(NumVect<ET> && vector_) : NumVect<ET>(std::move(vector_)) {normalize();};
    void normalize(); //sets norm2 to the correct value. Use after using operations from the NumVect<ET> base class.
    ExactLatticePoint& operator=(ExactLatticePoint const &that) =default;
    ExactLatticePoint& operator=(ExactLatticePoint && that) =default;
    ~ExactLatticePoint() {}

    friend void swap(ExactLatticePoint<ET,nfixed> &A, ExactLatticePoint<ET,nfixed> &B)
    {
        std::swap(A.data, B.data);
        std::swap(A.norm2,B.norm2);
    }
    friend std::ostream & operator<< <ET,nfixed>(std::ostream &os, ExactLatticePoint<ET,nfixed> const & exact_point);     //printing
    friend std::istream & operator>> <ET,nfixed>(std::istream &is, ExactLatticePoint<ET,nfixed> & exact_point) = delete;  //reading (may also be used by constructor from istream)


    NV const & access_vector() const {return *this;} //Probably not needed anyway.
    ET const & access_norm2() const {return norm2;} //Note: Does not copy
    unsigned int get_dim() const {return NumVect<ET>::size();} //returns dimension
    //inline void get_norm2 (ET &norm_to_return) const {norm_to_return = norm2;}
    //inline void setNorm2 (ET norm) {this->norm2 = norm;} //should not be required, actually.

    void print_exact_point(std::ostream & os = cout) const
    {
        os << * (static_cast<NumVect<ET> const *>(this)) << " of norm: " << this->norm2 << endl;
    }
    bool operator< (ExactLatticePoint const &other ) const {return (this->norm2 < other.norm2);};
    bool operator<=(ExactLatticePoint const &other ) const {return (this->norm2 <=other.norm2);};
    bool operator> (ExactLatticePoint const &other ) const {return (this->norm2 > other.norm2);};
    bool operator>=(ExactLatticePoint const &other ) const {return (this->norm2 >=other.norm2);};
};

template <class ET,int nfixed> ExactLatticePoint<ET> operator+ (ExactLatticePoint<ET,nfixed> const &A, ExactLatticePoint<ET,nfixed> const &B);
template <class ET,int nfixed> ExactLatticePoint<ET,nfixed> operator- (ExactLatticePoint<ET,nfixed> const &A, ExactLatticePoint<ET> const &B);
template <class ET,int nfixed> ExactLatticePoint<ET,nfixed> operator- (ExactLatticePoint<ET,nfixed> const &A); //unary-
template <class ET,int nfixed> void scalar_mult (ExactLatticePoint<ET,nfixed> &A, ET const & multiple); // A <- A*multiple
//template <class ET,int nfixed> void compute_exact_sc_product (ET &result, const ExactLatticePoint<ET,nfixed> &p1, const ExactLatticePoint<ET,nfixed> &p2); only use version above, for simplicity. With Return Value Optimization mandated by now, there is no reason for this.
template <class ET,int nfixed> ExactLatticePoint<ET,nfixed> conv_matrixrow_to_lattice_point (MatrixRow<ET> const &row);
template <class ET,int nfixed> ET exact_scalar_product(ExactLatticePoint<ET,nfixed> const &p1, ExactLatticePoint<ET,nfixed> const & p2);

#endif
