#ifndef APPROXIMATE_LATTICE_POINT_H
#define APPROXIMATE_LATTICE_POINT_H


/*
    ApproximateLatticePoint stores an approximation to a lattice point in the form
    2^exponent * (mantissa), where exponent is shared among coordinates and
    (mantissa) is a vector of type LatticeApproximations::ApproxType (should be approx 16-32 bits)
    We also store a copy of the norm^2 of the mantissa (so the real norm is 2^2exponent * norm_of_mantissa
    Note that is does not really depend much on the underlying type ET, which is only used to select the correct conversion routines and to know whether exponent may be negative.
*/

template <class ET> class ApproximateLatticePoint<ET,-1>
{
    public:
    using ApproxTypeNorm2   = LatticeApproximationsNew::ApproximationNorm2Type;
    using ExactEntryType    = ET;
    using ApproxEntryType   = LatticeApproximationsNew::ApproximationEntriesType;

    public: //consider making some constructors private and befriend the list class(es).
    ApproximateLatticePoint() : length_exponent(0),approx(nullptr), approxn2(0){}; //creates an empty approximation
    ApproximateLatticePoint(ApproximateLatticePoint<ET,-1> const & other) = delete; //as long as we don't store the length, we can't directly copy
    ApproximateLatticePoint(ApproximateLatticePoint<ET,-1> && other)                //moving is fine
            : length_exponent(other.length_exponent), approx(other.approx), approxn2(other.approxn2) {other.approx=nullptr;};
    ApproximateLatticePoint& operator= (ApproximateLatticePoint<ET,-1> const &other) =delete;
    ApproximateLatticePoint& operator= (ApproximateLatticePoint<ET,-1> && other)
    {
        length_exponent = other.length_exponent;approxn2=other.approxn2;swap(approx,other.approx);return *this;
    };
    explicit ApproximateLatticePoint(ExactLatticePoint<ET,-1> const & exact_point);
    //ApproximateLatticePoint(ExactLatticePoint const & exact_point); //create approximation from exact point.
    ApproximateLatticePoint(ApproximateLatticePoint<ET,-1> const & other, int const dim);
    ApproximateLatticePoint(ApproximateLatticePoint<ET,-1> && other, int const dim)
            : length_exponent(other.length_exponent), approx(other.approx), approxn2(other.approxn2) {other.approx=nullptr;};
    ~ApproximateLatticePoint(){delete approx;};

    ApproxEntryType*    access_vectors_mantissa()   const               {return approx;};
    ApproxTypeNorm2     get_norm2_mantissa()   const                    {return approxn2;};
    signed int get_vectors_exponent()  const                            {return length_exponent;};
    signed int get_norm2_exponent() const                               {return 2*length_exponent;};
    unsigned int get_dim() const                                        =delete; //Not implemented! Dimension is not stored to save memory (because it's the same for each one and we store lots of Approx. Points)
    void replace_by(ApproximateLatticePoint<ET,-1> const & other, int const dim);
    double get_norm2_d() const;
    ApproximateLatticePoint<ET,-1> make_copy(int const dim) const;

    protected: //internal data
    signed int length_exponent; //note : May be negative
    ApproxEntryType *approx; //array for the mantissa
    ApproxTypeNorm2 approxn2; //mantissa of the norm2
    //The approximation itself is given by 2^length_exponent * approx
    //The approximation to the norm^2 is given by 2^2length_exponent * approxn2
    //Note that we need to care about overflows here by truncating accordingly.
};

template<class ET, int nfixed> ostream & operator<< (ostream &os, ApproximateLatticePoint<ET,nfixed> const & approx_point) = delete; //doesn't work without dim
template<class ET, int nfixed> istream & operator>> (istream &os, ApproximateLatticePoint<ET,nfixed> & approx_point) = delete;


//template<class ET, int nfixed> ostream& operator<< (ostream &os, ApproximateLatticePoint<ET,nfixed> const & approx_point) //output.
//{
//    assert(false); //does not work because we do not store length.
//    return os;
//}

template<class ET, int nfixed>
double approximate_scalar_product_d(ApproximateLatticePoint<ET,nfixed> const &arg1, ApproximateLatticePoint<ET,nfixed> const & arg2, int const dim);

template<class ET, int nfixed>
ApproximateLatticePoint<ET,nfixed> create_detached_approximation(ExactLatticePoint<ET,nfixed> const & exact_point)
{
    ApproximateLatticePoint<ET,nfixed> approx_lp(exact_point); return approx_lp;
}

namespace GaussSieve
{
//computes scalar products: mantissa variant only computes the scalar product of the mantissas (i.e. it lacks an 2^{exponent1+exponent2} - factor)
//                          vec version computes scalar products between two arrays.
template<class ET,int nfixed>   LatticeApproximationsNew::ApproximationNorm2Type compute_mantissa_sc_product(ApproximateLatticePoint<ET,nfixed> const & arg1, ApproximateLatticePoint<ET,nfixed> const & arg2, int const dim);
inline                                 LatticeApproximationsNew::ApproximationNorm2Type compute_vec_sc_product(LatticeApproximationsNew::ApproximationEntriesType const * const arg1, LatticeApproximationsNew::ApproximationEntriesType const * const arg2, int const dim);
//FIXME: templating by nfixed does not work here. Correct way is to make ApproxEntryType* a template.
//The error message is weird and makes me suspect compiler bugs (the issue is with template function signature not depending on template param, but it give namespace errors (GCC)
};

#endif // APPROX_LATTICE_POINT_H
