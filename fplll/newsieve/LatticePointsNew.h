/* New header file for classes storing lattice points. */

#ifndef LATTICE_VECTOR_CLASS_NEW_H
#define LATTICE_VECTOR_CLASS_NEW_H

#include "sieve_common.h" //needed (at least for convertions from MatrixRow (the header has to be revised);
#include <string>
#include <numeric>
//probably some includes are missing, it's just working due to order-of-inclusion, as everything is really just one big header...

//forward declarations:

/*
-------------------------------------------------------------------------------------------------------------------------
ExactLatticePoint
-------------------------------------------------------------------------------------------------------------------------
*/

/*
ExactLatticePoint<ET,nfixed> stores an (exact) lattice point, together with its (squared) length.
nfixed denotes a fixed dimension (at compile-time) of vectors, with -1 meaning runtime-determined.
nfixed is not currently implemented.
*/

template <class ET,int nfixed=-1> class ExactLatticePoint;

/*
Usage notes:
    While ExactLatticePoint is currently implemented as derived from NumVect, usage of NumVect - specific classes is not recommended, as this might be changed.
    (NumVect's design clashes with the style used throughout the sieving code)
    Warning: Usage of NumVect - functions does not update norm2 appropriately.

    ExactLatticePoints may be copied, moved and default-constructed.
    Use ExactLatticePoint(dim,fill) to create a vector of dimension dim, filled with fill.

public members:
    ET norm2:                   square of norm. Use read-only.

public member functions defined:

    normalize()                 Updates norm2, usually no need to call directly.
    access_vector()             Used to access the underlying NumVect, READ-ONLY
    access_norm2()              Used to access norm2, READ-ONLY.  point.norm2 is equivalent to point.access_norm2(),  except that it enforces const-ness is and provided for consistency.
    get_dim()                   Used to retrieve the dimension of the vector. May become deprecated.
    print_exact_point(os=cout)  Prints to the output stream os, default is cout. Output includes norm^2.

Important operators supported:
    [i]                         Access i'th coordinate (via NumVect). Accessing out-of-range may cause undefined behaviour.
    =                           Copies / moves
    <, <=, >, >=                Compares by length. Note that these comparison operators only satisfy very weak guarantees. In particular, A <= B, B<=A does not imply A==B.
    stream <<, stream >>        Stream input/output. Input is not implemented yet.
    +, - (binary/unary)         Arithmetic

Non-member functions:
    void scalar_mult (A, x)     A := A*x where A is an ExactLatticePoint and x is an ET - scalar. Overwrites A.
    ET exact_scalar_product(A,B)Returns the scalar product of A and B.

    TODO: Cleanup conversion from matrix_row, more arithmetic, don't inherit from NumVect
*/

/*
-------------------------------------------------------------------------------------------------------------------------
ApproximateLatticePoint
-------------------------------------------------------------------------------------------------------------------------
*/

/*
    Note: We only have an implementation for nfixed == -1 for now.
*/

/*
    ApproximateLatticePoint stores an approximation to a lattice point in the form
    2^exponent * (mantissa), where exponent is shared among coordinates and
    (mantissa) is a vector of type LatticeApproximations::ApproxType (should be approx 16-32 bits)
    We also store a copy of the norm^2 of the mantissa (so the real norm is 2^2exponent * norm_of_mantissa
    Note that is does not really depend much on the underlying type ET, which is only used to select the correct conversion routines and to know whether the exponent may be negative.
*/

template<class ET, int nfixed=-1> class ApproximateLatticePoint;      //This class represents a stand-alone approximation to a lattice point

/*
Usage notes:
    ApproximateLatticePoint is essentially a wrapper around some vector of int's (called mantissas) and an exponent
    Since the lattice points we use in our algorithm(s) have all the same dimension, we do not explicitly store the dimension for each such vector.
    A consequence is that copy - routines need to be explicitly supplied with the dimension. In particular, we have NO copy constructor.
    Note that even if we would store the dimension, we would opt to have no copy constructor, since we want to force the user to explicitly express intent to copy (by design).

    ApproximateLatticePoints may be moved and default-constructed. The is no real copy constructor.
    There is a "copy" constructor ApproximateLatticePoint(old, dim) taking dimension as argument.
    explicit ApproximateLatticePoint(exact_point) creates an approximation from an exact point.
    TODO: Construct from Compressed & Shadow

public member functions defined:
    get_norm2_mantissa()            gets the mantissa of the squared norm.
    get_norm2_exponent()            gets the exponent of the squared norm. Squared norm = 2^get_norm2_exponent * get_norm2_mantissa
    access_vectors_mantissa()       used to access the actual vector. Avoid using it, as return value type (currently pointer) might change. Might become private.
    get_vectors_exponent()          gets the exponent of the vector. The actual approx is 2^get_vectors_exponent * vectors_mantissa. Might become private.
    replace_by(other,dim)           A.replace_by(B,dim) causes A:=B.
    AppLatPo make_copy(dim)         creates a copy.  A = B.make_copy(dim) causes A:=B.
    get_norm2_double()              returns norm2 as a double.

non-member functions:
    create_detached_approximation(exact_point) creates an approximation from an exact point. Same as constructor, really.

    To compute scalar products, there is
    GaussSieve::compute_mantissa_sc_product(A,B,dim) to compute <A,B>   for ApproximateLatticePoints A,B
    GaussSieve::compute_vec_sc_product(A,B,dim)                         for the underlying vectors (don't use if you can avoid it!)
    approximate_scalar_product_d(A,B,dim)       this one returns a double.

TODO:   Arithmetic, convert to vectors of doubles?
NOTE:   See utility.h for algorithms that compute&compare scalar products (without converting to double)
        We do not have a class for mantissa + exponent, because we can just use IEEE floats. doubles have enough precision.
*/

/*
-------------------------------------------------------------------------------------------------------------------------
CompressedPoint
-------------------------------------------------------------------------------------------------------------------------
*/

//This class represents an approximate point together with some extra data (called Details) that allows to recover the exact point again.
//Note that access to the exact point may be slow (in multi-threading, it might involve locks or fences), so working with the approximation is prefered.

template<class ET, bool MT, int nfixed=-1> class CompressedPoint;
/*
Usage notes:
    CompressedPoints are what are stored in the main list. The main queue also returns these.
    (Note: THe main queue might always use MT==false)
    A compressed point is a lattice point in a "compressed" format that allows quick access to the approximation and possibly slow access to the exact point.
    CompressedPoints have NO copy constructor, because there are cases where we need deep and where we need shallow copies.

    CompressedPoints may be moved and default-constructed, but not directly copied/copy-assigned.
    CompressedPoints may be created from ExactLatticePoints or pointers to exact lattice points:
    explicit CompressedPoint(exact_lattice_point)
    explicit CompressedPoint(*exact_lattice_point)
    Important: The latter transfers ownership, whereas the first does not!

public member functions:
    access_approximation_r()            accesses the approximation.
    get_approximation_const_ptr()       returns const pointer to approximation
    get_approximation_ptr               returns point to appoximation
    access_exact_point_r()              accesses the exact point -- DO NOT USE THIS IN MULTI-THREADED, will be made non-public

    get_approx_norm2_mantissa()         returns a copy of the approximate norm's mantissa
    get_approx_norm2_exponent()         returns a copy of the approximate norm's exponent
    get_approx_norm2_d()                returns a copy of the approximate norm as a double
    get_approx_vectors_exponent()       returns a copy of the exponent of the approximate vector.
    access_approx_vectors_mantissa      allows access to the underlying mantissa of the vector. Avoid using

    get_exact_norm2                     gets a copy of the exact norm (of type ET)
    get_exact_point                     returns a copy of the exact point (of type ExactLatticePoint)
    deep_copy_compressed_point          returns a copy of the compressed point (expensive, since it copies the underlying details)

    is_zero()                           tells whether a point is zero. Looks at the approximation only.

NOTE:
    We have a shallow copy version to another type rather than just doing ref-counting.
    The latter would be much easier, but has the problem that updating the ref-counts requires thread-synchronisation.
    Also, we do not really want to restrict the nature of "details" too much, these could be some kind of globally managed object.


*/



namespace LatticeApproximationsNew //internal helper types for approximations etc. enclosed in namespace
{
inline signed int get_exponent (Z_NR<mpz_t> const & val)     {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //returns smallest t s.t. abs(val) < 2^t (i.e. the bit-length for integral types)
inline signed int get_exponent (Z_NR<long> const & val)      {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //for val==0, we return the most negative value that fits into an int.
inline signed int get_exponent (Z_NR<double> const & val)    {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //Note that the function is already implemented in Z_NR< > - types.
template<class ET> [[ deprecated ("Using badly supported type") ]] signed int get_exponent (ET const & val);                               //Non-Z_NR< > - cases. Note that non-templates take precedence.


//ApproxType do_approximate(X val, int delta) returns 2^(-delta)*val, converted to appropriate ApproxType.
//Implementation depends on X

template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<double> >::value, ET>::type const & val, signed int const delta);
template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<mpz_t > >::value, ET>::type const & val, signed int const delta);
template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<long  > >::value, ET>::type const & val, signed int const delta);

using ApproximationNorm2Type = int32_t; //determines bit-length of approximation, by making sure that n* (MAXAPPROX^2) must fit into it.
using ApproximationEntriesType = int32_t; //at least half the size of the above. Same size might work better for vectorization.

template<class ET> class MaybeRational; //Compile-time function of ET that tells whether ET might represent values between 0 and 1
template<class ET> class MaybeRational       {public: static bool constexpr val=true;};  //helper template that selects whether the underlying type *might* be able to represent values strictly between 0 and 1.
template<> class MaybeRational<Z_NR<long > > {public: static bool constexpr val=false;}; //this distinction just to serves to avoid correct, but needless approximations if the values are small enough to not need approximations in the first place.
template<> class MaybeRational<Z_NR<mpz_t> > {public: static bool constexpr val=false;}; //(otherwise, we would pad with zeros from the right(!), which is correct (but useless and hurts readability for debug outputs).
template<> class MaybeRational<Z_NR<double>> {public: static bool constexpr val=true;};
//template<class ET>

//inline ApproxTypeNorm2 compute_sc_prod(ApproxType const * const arg1, ApproxType const * const arg2, unsigned int len);

//template<class ET>
//inline bool Compare_Sc_Prod(ApproxLatticePoint<ET,false,-1> const & arg1, ApproxLatticePoint<ET,false,-1> const & arg2, ApproxTypeNorm2 abslimit, int limit_exp, int dim);
};

//non-member functions belonging to the various container classes for lattice points.

//namespace GaussSieve
//{
////computes scalar products: mantissa variant only computes the scalar product of the mantissas (i.e. it lacks an 2^{exponent1+exponent2} - factor)
////                          vec version computes scalar products between two arrays.
//template<class ET,int nfixed>   LatticeApproximationsNew::ApproximationNorm2Type compute_mantissa_sc_product(ApproximateLatticePoint<ET,nfixed> const & arg1, ApproximateLatticePoint<ET,nfixed> const & arg2, int const dim);
//                                LatticeApproximationsNew::ApproximationNorm2Type compute_vec_sc_products(LatticeApproximationsNew::ApproximationEntriesType const * const arg1, LatticeApproximationsNew::ApproximationEntriesType const * const arg2, int const dim);
////FIXME: templating by nfixed does not work here. Correct way is to make ApproxEntryType* a template.
////The error message is weird and makes me suspect compiler bugs (the issue is with template function signature not depending on template param, but it give namespace errors (GCC)
//};
/* Creates a stand-alone Approximation from an exact point */
//template<class ET, int nfixed> ApproximateLatticePoint<ET,nfixed> create_detached_approximation(ExactLatticePoint<ET,nfixed> const & exact_point);

//template <class ET,int nfixed=-1> void compute_exact_sc_product (ET &result, const ExactLatticePoint<ET,nfixed> &p1, const ExactLatticePoint<ET,nfixed> &p2);   //computes the scalar product of two exact lattice points.
//template <class ET,int nfixed=-1> ET exact_scalar_product(ExactLatticePoint<ET,nfixed> const &p1, ExactLatticePoint<ET,nfixed> const & p2);                     //same, but returns result instead of storing in first arg.
//template <class ET,int nfixed=-1> void scalar_mult (ExactLatticePoint<ET,nfixed> &A, ET const & multiple); //A = A*multiple
//template <class ET,int nfixed=-1> std::ostream & operator<<(std::ostream &os, ExactLatticePoint<ET,nfixed> const & exact_point); //printing
//template<class ET> istream & operator>>(istream &is, ExactLatticePoint & exact_point); //reading (may also be used by constructor from istream)
//FIXME: operator >> is NOT IMPLEMENTED YET
//template<class ET, int nfixed=-1> std::ostream& operator<< (std::ostream &os, ApproximateLatticePoint<ET,nfixed> const & approx_point); //output.
//TODO: Input



/*actual class declarations in separate files*/

#include "ExactLatticePoint.h"
#include "ApproximateLatticePoint.h"
#include "CompressedPoint.h"



//template<int nfixed> LatticeApproximationsNew::ApproximationNorm2Type GaussSieve::compute_vec_sc_products<-1>(LatticeApproximationsNew::ApproximationEntriesType const * const arg1, LatticeApproximationsNew::ApproximationEntriesType const * const arg2, int const dim)
//{
//    return std::inner_product(arg1, arg1+nfixed,arg2,0);
//}




#endif

//Old & deleted for now:

//FIXME: Enable or use a constructor
//template<class ET,int nfixed=-1> ExactLatticePoint<ET,nfixed> conv_matrixrow_to_lattice_point (MatrixRow<ET> const &row);
