#ifndef GAUSS_SIEVE_UTILITY_H
#define GAUSS_SIEVE_UTILITY_H

#include <type_traits>
#include "LatticePointsNew.cpp"


/**
This macro is used to test the presence of a (public) member typedef in a class
Args:   TypeToCheck - typename whose presence to check
        CheckerClassName - Name of the checker class
This macro emits a new template class definition with the name CheckerClassName.
TypeToCheck must not be void.

Usage:
CREATE_MEMBER_TYPEDEF_CHECK_CLASS(NiceTypedef, IsMyClassGood);
This creates the template class IsMyClassGood.

Then CheckerClassName<SomeSuspiciousClass>::value will be true if
SomeSuspicousClass::NiceTypedef exists, false otherwise

Note that the missing semicolon at the end is intentional.
The user needs to put it.
*/

// clang-format off

#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TypeToCheck, CheckerClassName)\
template<class ClassToCheck> class CheckerClassName                     \
{                                                                       \
private:                                                                \
    template<class Arg> static typename Arg:: TypeToCheck foo(int);     \
    template<class ...> void                              foo(...);     \
public:                                                                 \
    using value_t = std::integral_constant<bool,                        \
        !(std::is_void<decltype(foo< ClassToCheck >(0))>::value)>;      \
    static bool constexpr value = value_t::value;                       \
}

/**
Similar to the above, creates a checker template class that checks wether
TypeToCheck exists and is equal to TypeShouldBe
*/

#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(TypeToCheck,TypeShouldBe,CheckerClassName) \
template<class ClassToCheck> class CheckerClassName                                         \
{                                                                                           \
private:                                                                                    \
    template<class Arg> static typename Arg:: TypeToCheck foo(int);                         \
    template<class ...> void                              foo(...);                         \
public:                                                                                     \
    using value_t = std::integral_constant<bool, std::is_same< TypeShouldBe,                \
                decltype(foo< ClassToCheck > (0)) >::value>;                                \
                static bool constexpr value = value_t::value;                               \
}

// clang-format on


//class that ignores its argument. Can be used to optimize away unused parameters in function templates...
class IgnoreAnyArg{
    public:
    template<class T> constexpr IgnoreAnyArg(T val){};
    constexpr IgnoreAnyArg(){};
};

//same, but enforces the type of the ignored argument.
template<class T>
class IgnoreArg{
    public:
    inline constexpr IgnoreArg(T val){};
    IgnoreArg(){};
};

template<int nfixed=-1> class Dimension;

template<>
class Dimension<-1>{
    public:
    using IsFixed=false_type;
    Dimension(unsigned int const new_dim):dim(new_dim){};
    Dimension()=delete;
    inline operator unsigned int() const {return dim;};
    unsigned int dim;
};

template<int nfixed>
class Dimension{
    public:
    using IsFixed=true_type;
    Dimension(){};
    Dimension(IgnoreArg<unsigned int> new_dim){}; //assert(new_dim==nfixed);}
//    Dimension(unsigned int){};
    inline constexpr operator unsigned int() const {return nfixed;};
    static constexpr unsigned int dim = nfixed;
};

namespace GaussSieve //helper functions
{
inline    bool string_consume(istream &is, std::string const & str, bool elim_ws= true, bool verbose=true);   //helper function for dumping/reading
inline    Z_NR<mpz_t> compute_mink_bound(ZZ_mat<mpz_t> const & basis);                                        //computes a meaningful Minkowski bound for the length of the shortest vector

//    template<class ET>
//    inline bool Compare_Sc_Prod(ApproxLatticePoint<ET,false,-1> const & arg1, ApproxLatticePoint<ET,false,-1> const & arg2, ApproxTypeNorm2 abslimit, int limit_exp, int dim);

    //Compares | <arg1,arg2> | with 2^bound_exponent * bound_mantissa. Asserts that abslimit>=0
    //Returns true if  |<arg1,arg2>| < 2^bound_exponent*bound_mantissa
    //Returns false if |<arg1,arg2>| > 2^bound_exponent*bound_mantissa
    //In case of equality, the answer is not specified
    //Undefined behaviour may occur if abslimit < 0

    template<class ET, int nfixed> bool compare_abs_approx_scalar_product( ApproximateLatticePoint<ET,nfixed> const & arg1,
                                            ApproximateLatticePoint<ET,nfixed> const & arg2,
                                            LatticeApproximationsNew::ApproximationNorm2Type const bound_mantissa,
                                            int const bound_exponent,
                                            int const dim);

    //checks whether 2-reduction is possible, i.e. if there exists k s.t. p1 - k*p2 is shorter than p1.
    //If yes, returns true, otherwise false.
    //If the function returns yes, the value of k is stored in scalar.
    template<class ET,int nfixed> bool check2red_exact (ExactLatticePoint<ET,nfixed> const &p1, ExactLatticePoint<ET,nfixed> const &p2, ET &scalar); //only checks whether 2reduction is possible

    //performs 2-reduction, i.e. returns p1 - k*p2
    template<class ET,bool MT, int nfixed> CompressedPoint<ET,MT,nfixed> perform2red_exact_to_compressed(ExactLatticePoint<ET,nfixed> const &p1,ExactLatticePoint<ET,nfixed> const &p2, ET const & scalar);

}




 template<class ET,int nfixed> bool GaussSieve::compare_abs_approx_scalar_product( ApproximateLatticePoint<ET,nfixed> const & arg1,                                           ApproximateLatticePoint<ET,nfixed> const & arg2,
                                            LatticeApproximationsNew::ApproximationNorm2Type const bound_mantissa,
                                            int const bound_exponent,
                                            int const dim)
{
    int const rel = arg1.get_vectors_exponent() + arg2.get_vectors_exponent() - bound_exponent; //relative exponent between lhs and rhs of the expression we want to compare.
    LatticeApproximationsNew::ApproximationNorm2Type const abs_sc_prod_mantissa = abs(GaussSieve::compute_mantissa_sc_product(arg1,arg2,dim));
    if(rel > 0) //depending on whether the relative exponent is positive or not, we bit-shift either the lhs or the rhs.
                //This is to ensure we also shift right in order to avoid overflows.
    {
        return abs_sc_prod_mantissa > (bound_mantissa >> rel);
    }
    else
    {
        return (abs_sc_prod_mantissa >> -rel);
    }
}

template<class ET,int nfixed> bool GaussSieve::check2red_exact (ExactLatticePoint<ET,nfixed> const &p1, ExactLatticePoint<ET,nfixed> const &p2, ET &scalar) //only checks whether 2reduction is possible
{
    //assert(p1.norm2 >= p2.norm2); Not neccessarily true in multi-threaded case. -- Gotti
    ET abs_2scprod;
    scalar = 0;
    ET sc_prod=exact_scalar_product(p1, p2);    //stores scalar product of p1 and p2 into sc_prod.
    abs_2scprod.mul_ui(sc_prod,2);  //multiply result by two
    abs_2scprod.abs(abs_2scprod);   //and take the absolute value

    // check if |2 * <p1, p2>| <= |p2|^2. If yes, no reduction is possible
    if (abs_2scprod <= p2.norm2)
        return false;

    // compute the (integer) multiple for p1: mult = round(<p1, p2> / |p2|^2)
    FP_NR<double> mult, tmp; //maybe we can use another type
    mult.set_z(sc_prod); //conversions
    tmp.set_z(p2.norm2);


    mult.div(mult, tmp);
    mult.rnd(mult);
    scalar.set_f(mult); //converts mult to the type suitable for mult_const; This restricts the types of ET, as we rely on the Z_NR<foo> - routines here.
    return true;
}

template<class ET,bool MT, int nfixed> CompressedPoint<ET,MT,nfixed> GaussSieve::perform2red_exact_to_compressed(ExactLatticePoint<ET,nfixed> const &p1,ExactLatticePoint<ET,nfixed> const &p2, ET const & scalar)
{
    ExactLatticePoint<ET,nfixed> * res = new ExactLatticePoint<ET,nfixed>(p1);
    ET tmpnegscalar;
    tmpnegscalar.neg(scalar); //negscalar := -scalar.
    res->addmul(p2,std::move(tmpnegscalar)); //Z_NR< > does not allow just -scalar, we have to do stupid things.
    res->normalize();
    return static_cast<CompressedPoint<ET,MT,nfixed > >(res); //Note: This constructor transfers ownership, which is why we don't delete res.
}


//template<class ET>
//inline bool Compare_Sc_Prod(ApproxLatticePoint<ET,false,-1> const & arg1, ApproxLatticePoint<ET,false,-1> const & arg2, ApproxTypeNorm2 abslimit, int limit_exp, int dim)
//{
//    int rel = arg1.get_length_exponent() + arg2.get_length_exponent() - limit_exp;
//    ApproxTypeNorm2 sc = abs(compute_sc_prod(arg1.get_approx(),arg2.get_approx(),dim));
//    if(rel > 0)
//    {
//        return sc > (abslimit >> rel);
//    }
//    else
//    {
//        return (sc >> -rel) > abslimit;
//    }
//}



//    template<class ET>     bool check_perform_2red (LatticePoint<ET> &p1, const LatticePoint<ET> &p2); //2-reduces p1 with the help of p2.
                                                                                                       //p1 is overwritten, whereas p2 is const. Returns true if p1 actually changed.
//    template<class ET>     bool check2red_new (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET &scalar); //only checks whether 2reduction is possible
//    template<class ET>     LatticePoint<ET> perform2red (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET const & scalar); //replaces p1 by p1 - scalar * p2
//    template<class ET>     bool check3red(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, float px1, float px2, float x1x2, int & sgn1, int & sgn2);
//    template<class ET>     bool check3red_signs(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, int px1, int px2, int x1x2, int &sgn1, int &sgn2);
//    template<class ET>    LatticePoint<ET> perform3red (const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, const int & sgn1, const int &sgn2);



#endif // GAUSS_SIEVE_UTILITY_H
