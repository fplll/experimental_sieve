#ifndef GAUSS_SIEVE_UTILITY_CPP
#define GAUSS_SIEVE_UTILITY_CPP

#include "Utility.h"

//OLD FUNCTIONS:

namespace GaussSieve //helper functions
{
inline    bool string_consume(istream &is, std::string const & str, bool elim_ws= true, bool verbose=true);   //helper function for dumping/reading

    template<class ET>
    inline bool Compare_Sc_Prod(ApproxLatticePoint<ET,false,-1> const & arg1, ApproxLatticePoint<ET,false,-1> const & arg2, ApproxTypeNorm2 abslimit, int limit_exp, int dim);

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


template<class ET>
bool GaussSieve::check_perform_2red (LatticePoint<ET> &p1, const LatticePoint<ET> &p2)
{
    //assert(p1.norm2 >= p2.norm2);     Not neccessarily true in multi-threaded case. -- Gotti
    ET sc_prod, abs_2scprod, scalar;
    sc_product(sc_prod, p1, p2);
    abs_2scprod.mul_ui(sc_prod,2);
    abs_2scprod.abs(abs_2scprod);

    // check if |2 * <p1, p2>| <= |p2|^2. If yes, no reduction
    if (abs_2scprod <= p2.norm2)
        return false;

    // compute the (integer) multiple for p1: mult = round(<p1, p2> / |p2|^2)
    FP_NR<double> mult, tmp; //may be can use another type
    mult.set_z(sc_prod); //conversions
    tmp.set_z(p2.norm2);


    mult.div(mult, tmp);
    mult.rnd(mult);
    scalar.set_f(mult); //converts mult to the type suitable for mult_const;


    LatticePoint<ET> res(p2);
    scalar_mult(res, scalar);
    p1 = p1 - res;
    return true;
}

// separate check2Red and perform2Red

//if true, scalar is the multiple s.t. we reduce p1 = p1-sclar * p2;

template<class ET> bool GaussSieve::check2red_new (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET &scalar)
{

    //assert(p1.norm2 >= p2.norm2); Not neccessarily true in multi-threaded case. -- Gotti
    ET sc_prod, abs_2scprod;
    scalar = 0;
    sc_product(sc_prod, p1, p2);
    abs_2scprod.mul_ui(sc_prod,2);
    abs_2scprod.abs(abs_2scprod);

    // check if |2 * <p1, p2>| <= |p2|^2. If yes, no reduction
    if (abs_2scprod <= p2.norm2)
        return false;

    // compute the (integer) multiple for p1: mult = round(<p1, p2> / |p2|^2)
    FP_NR<double> mult, tmp; //may be can use another type
    mult.set_z(sc_prod); //conversions
    tmp.set_z(p2.norm2);


    mult.div(mult, tmp);
    mult.rnd(mult);
    scalar.set_f(mult); //converts mult to the type suitable for mult_const;
    return true;
}

// return res = p1 - scalar*p2;
template<class ET>
LatticePoint<ET> GaussSieve::perform2red (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET const & scalar)
{
    LatticePoint<ET> res(p2);
    scalar_mult(res, scalar);
    return (p1 - res);

}



//
//the first argument p is assumed to have the largest norm
// the function returns true if indeed || p \pm x1 \pm x2 || < || p ||
// The correct signs in front of x1, x2 are deduced from the sign os the corresp. inner-products px1, px2.
// The output is true if p can be reduced
//


template<class ET>
bool GaussSieve::check3red(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, float px1, float px2, float x1x2, int &sgn1, int &sgn2)
{
    //in case sgn(x1x2)*sgn(px1)*sgn(px2) = 0, we cannot produce a shorter p:
    //  either they are all positive (clearly, no combination of \pm can result in a shorter p
    //  or there are two inner-products with -1 and one point in common. Take the negative of the common vector, arrive to the first case

    //if (px1.sgn()*px2.sgn()*x1x2.sgn() == 1)
    //    return false;


    //TODO: look-up a sign-fnct for floats
    if (px1>0 && px2>0 && x1x2>0)
        return false;
    if (px1>0 && px2<0 && x1x2<0)
        return false;
    if (px2>0 && px1<0 && x1x2<0)
        return false;
    if (x1x2>0 && px1<0 && px2<0)
        return false;

    LatticePoint<ET> res(p);


    if (px1<0){
        res = res+x1;
        sgn1=1;
    }
    else {
        res = res-x1;
        sgn1=-1;
    }
    if(px2<0){
        res = res+x2;
        sgn2=1;
    }
    else{
        res = res-x2;
        sgn2=-1;
    }

    if(res.norm2>=p.norm2)
    {
        return false;
    }

    //cout <<"p: ";
    //p.printLatticePoint();
    //cout<< "res: ";
    //res.printLatticePoint();
    return true;
}

template<class ET>
bool GaussSieve::check3red_signs(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, int px1, int px2, int x1x2, int &sgn1, int &sgn2)
{
    //in case sgn(x1x2)*sgn(px1)*sgn(px2) = 0, we cannot produce a shorter p:
    //  either they are all positive (clearly, no combination of \pm can result in a shorter p
    //  or there are two inner-products with -1 and one point in common. Take the negative of the common vector, arrive to the first case

    //if (px1.sgn()*px2.sgn()*x1x2.sgn() == 1)
    //    return false;


    //TODO: look-up a sign-fnct for floats
    if (px1>0 && px2>0 && x1x2>0)
        return false;
    if (px1>0 && px2<0 && x1x2<0)
        return false;
    if (px2>0 && px1<0 && x1x2<0)
        return false;
    if (x1x2>0 && px1<0 && px2<0)
        return false;

    LatticePoint<ET> res(p);


    if (px1<0){
        res = res+x1;
        sgn1=1;
    }
    else {
        res = res-x1;
        sgn1=-1;
    }
    if(px2<0){
        res = res+x2;
        sgn2=1;
    }
    else{
        res = res-x2;
        sgn2=-1;
    }

    if(res.norm2>=p.norm2)
    {
        //cout << "res = " <<  res << " res.norm2 = " << res.norm2 <<  "p.norm2 = " << p.norm2 << endl;
        return false;
    }

    //cout <<"p: ";
    //p.printLatticePoint();
    //cout<< "res: ";
    //res.printLatticePoint();
    return true;
}

// return res = p +sgn1*x1+sgn2*x2; res is supposed to be shorter than p

template<class ET>
LatticePoint<ET> GaussSieve::perform3red (const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, const int & sgn1, const int &sgn2)
{
    LatticePoint<ET> res(p);
    if (sgn1<0)
        res = res-x1;
    else
        res = res+x1;

    if(sgn2<0)
        res = res-x2;
    else
        res = res+x2;

    //cout << "3 red is performed:" << endl;
    //cout <<"p: ";
    //p.printLatticePoint();
    //cout<< "res: ";
    //res.printLatticePoint();
    return res;
}

/*
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
*/

/*
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
*/


/*
template<class ET,bool MT, int nfixed> CompressedPoint<ET,MT,nfixed> GaussSieve::perform2red_exact_to_compressed(ExactLatticePoint<ET,nfixed> const &p1,ExactLatticePoint<ET,nfixed> const &p2, ET const & scalar)
{
    ExactLatticePoint<ET,nfixed> * res = new ExactLatticePoint<ET,nfixed>(p1);
    ET tmpnegscalar;
    tmpnegscalar.neg(scalar); //negscalar := -scalar.
    res->addmul(p2,std::move(tmpnegscalar)); //Z_NR< > does not allow just -scalar, we have to do stupid things.
    res->normalize();
    return static_cast<CompressedPoint<ET,MT,nfixed > >(res); //Note: This constructor transfers ownership, which is why we don't delete res.
}
*/



//    template<class ET>     bool check_perform_2red (LatticePoint<ET> &p1, const LatticePoint<ET> &p2); //2-reduces p1 with the help of p2.
                                                                                                       //p1 is overwritten, whereas p2 is const. Returns true if p1 actually changed.
//    template<class ET>     bool check2red_new (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET &scalar); //only checks whether 2reduction is possible
//    template<class ET>     LatticePoint<ET> perform2red (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET const & scalar); //replaces p1 by p1 - scalar * p2
//    template<class ET>     bool check3red(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, float px1, float px2, float x1x2, int & sgn1, int & sgn2);
//    template<class ET>     bool check3red_signs(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, int px1, int px2, int x1x2, int &sgn1, int &sgn2);
//    template<class ET>    LatticePoint<ET> perform3red (const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, const int & sgn1, const int &sgn2);



#endif
