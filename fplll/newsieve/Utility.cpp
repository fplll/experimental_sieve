#ifndef GAUSS_SIEVE_UTILITY_CPP
#define GAUSS_SIEVE_UTILITY_CPP

#include "Utility.h"

//template<class ET>
//bool GaussSieve::check_perform_2red (LatticePoint<ET> &p1, const LatticePoint<ET> &p2)
//{
//    //assert(p1.norm2 >= p2.norm2);     Not neccessarily true in multi-threaded case. -- Gotti
//    ET sc_prod, abs_2scprod, scalar;
//    sc_product(sc_prod, p1, p2);
//    abs_2scprod.mul_ui(sc_prod,2);
//    abs_2scprod.abs(abs_2scprod);
//
//    // check if |2 * <p1, p2>| <= |p2|^2. If yes, no reduction
//    if (abs_2scprod <= p2.norm2)
//        return false;
//
//    // compute the (integer) multiple for p1: mult = round(<p1, p2> / |p2|^2)
//    FP_NR<double> mult, tmp; //may be can use another type
//    mult.set_z(sc_prod); //conversions
//    tmp.set_z(p2.norm2);
//
//
//    mult.div(mult, tmp);
//    mult.rnd(mult);
//    scalar.set_f(mult); //converts mult to the type suitable for mult_const;
//
//
//    LatticePoint<ET> res(p2);
//    scalar_mult(res, scalar);
//    p1 = p1 - res;
//    return true;
//}

// separate check2Red and perform2Red

//if true, scalar is the multiple s.t. we reduce p1 = p1-sclar * p2;

//template<class ET> bool GaussSieve::check2red_new (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET &scalar)
//{
//
//    //assert(p1.norm2 >= p2.norm2); Not neccessarily true in multi-threaded case. -- Gotti
//    ET sc_prod, abs_2scprod;
//    scalar = 0;
//    sc_product(sc_prod, p1, p2);
//    abs_2scprod.mul_ui(sc_prod,2);
//    abs_2scprod.abs(abs_2scprod);
//
//    // check if |2 * <p1, p2>| <= |p2|^2. If yes, no reduction
//    if (abs_2scprod <= p2.norm2)
//        return false;
//
//    // compute the (integer) multiple for p1: mult = round(<p1, p2> / |p2|^2)
//    FP_NR<double> mult, tmp; //may be can use another type
//    mult.set_z(sc_prod); //conversions
//    tmp.set_z(p2.norm2);
//
//
//    mult.div(mult, tmp);
//    mult.rnd(mult);
//    scalar.set_f(mult); //converts mult to the type suitable for mult_const;
//    return true;
//}

// return res = p1 - scalar*p2;
//template<class ET>
//LatticePoint<ET> GaussSieve::perform2red (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET const & scalar)
//{
//    LatticePoint<ET> res(p2);
//    scalar_mult(res, scalar);
//    return (p1 - res);
//
//}
//


//
//the first argument p is assumed to have the largest norm
// the function returns true if indeed || p \pm x1 \pm x2 || < || p ||
// The correct signs in front of x1, x2 are deduced from the sign os the corresp. inner-products px1, px2.
// The output is true if p can be reduced
//


//template<class ET>
//bool GaussSieve::check3red(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, float px1, float px2, float x1x2, int &sgn1, int &sgn2)
//{
//    //in case sgn(x1x2)*sgn(px1)*sgn(px2) = 0, we cannot produce a shorter p:
//    //  either they are all positive (clearly, no combination of \pm can result in a shorter p
//    //  or there are two inner-products with -1 and one point in common. Take the negative of the common vector, arrive to the first case
//
//    //if (px1.sgn()*px2.sgn()*x1x2.sgn() == 1)
//    //    return false;
//
//
//    //TODO: look-up a sign-fnct for floats
//    if (px1>0 && px2>0 && x1x2>0)
//        return false;
//    if (px1>0 && px2<0 && x1x2<0)
//        return false;
//    if (px2>0 && px1<0 && x1x2<0)
//        return false;
//    if (x1x2>0 && px1<0 && px2<0)
//        return false;
//
//    LatticePoint<ET> res(p);
//
//
//    if (px1<0){
//        res = res+x1;
//        sgn1=1;
//    }
//    else {
//        res = res-x1;
//        sgn1=-1;
//    }
//    if(px2<0){
//        res = res+x2;
//        sgn2=1;
//    }
//    else{
//        res = res-x2;
//        sgn2=-1;
//    }
//
//    if(res.norm2>=p.norm2)
//    {
//        return false;
//    }
//
//    //cout <<"p: ";
//    //p.printLatticePoint();
//    //cout<< "res: ";
//    //res.printLatticePoint();
//    return true;
//}

//template<class ET>
//bool GaussSieve::check3red_signs(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, int px1, int px2, int x1x2, int &sgn1, int &sgn2)
//{
//    //in case sgn(x1x2)*sgn(px1)*sgn(px2) = 0, we cannot produce a shorter p:
//    //  either they are all positive (clearly, no combination of \pm can result in a shorter p
//    //  or there are two inner-products with -1 and one point in common. Take the negative of the common vector, arrive to the first case
//
//    //if (px1.sgn()*px2.sgn()*x1x2.sgn() == 1)
//    //    return false;
//
//
//    //TODO: look-up a sign-fnct for floats
//    if (px1>0 && px2>0 && x1x2>0)
//        return false;
//    if (px1>0 && px2<0 && x1x2<0)
//        return false;
//    if (px2>0 && px1<0 && x1x2<0)
//        return false;
//    if (x1x2>0 && px1<0 && px2<0)
//        return false;
//
//    LatticePoint<ET> res(p);
//
//
//    if (px1<0){
//        res = res+x1;
//        sgn1=1;
//    }
//    else {
//        res = res-x1;
//        sgn1=-1;
//    }
//    if(px2<0){
//        res = res+x2;
//        sgn2=1;
//    }
//    else{
//        res = res-x2;
//        sgn2=-1;
//    }
//
//    if(res.norm2>=p.norm2)
//    {
//        //cout << "res = " <<  res << " res.norm2 = " << res.norm2 <<  "p.norm2 = " << p.norm2 << endl;
//        return false;
//    }
//
//    //cout <<"p: ";
//    //p.printLatticePoint();
//    //cout<< "res: ";
//    //res.printLatticePoint();
//    return true;
//}

// return res = p +sgn1*x1+sgn2*x2; res is supposed to be shorter than p
//template<class ET>
//LatticePoint<ET> GaussSieve::perform3red (const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, const int & sgn1, const int &sgn2)
//{
//    LatticePoint<ET> res(p);
//    if (sgn1<0)
//        res = res-x1;
//    else
//        res = res+x1;
//
//    if(sgn2<0)
//        res = res-x2;
//    else
//        res = res+x2;
//
//    //cout << "3 red is performed:" << endl;
//    //cout <<"p: ";
//    //p.printLatticePoint();
//    //cout<< "res: ";
//    //res.printLatticePoint();
//    return res;
//
//}

//helper function for reading in from streams. Gobbles up str from the stream (and optionally whitespace before/after).
//If str is not on the stream, outputs an error.
//Note that whitespace inside str is OK, as long as it is not at the beginning or end.

inline bool GaussSieve::string_consume(istream &is, std::string const & str, bool elim_ws, bool verbose)
{
    unsigned int len = str.length();
    char *buf = new char[len+1];
    buf[len] = 0; //for error message.
    if (elim_ws)
    {
        is >> std::ws;
    }
    is.read(buf,len);
    if(is.gcount() != len)
    {
        if(verbose)
        {
            cerr << "Failure reading header: Expected to read" << str << endl;
            cerr << "Read only "<<is.gcount() << "bytes. String read was" << buf<<endl;
        }
        return false;
    }
    if(elim_ws)
    {
        is >> std::ws;
    }
    if(str.compare(0,len,buf,len)!=0)
    {
        if(verbose)
        {
            cerr << "Failure reading header: Expected to read" << str << endl;
            cerr << "Read instead:" << buf << endl;
        }
        return false;
    }
    return true;
}

//namespace GaussSieve{
inline Z_NR<mpz_t> GaussSieve::compute_mink_bound(ZZ_mat<mpz_t> const & basis)
{
    assert(basis.get_rows() == basis.get_cols()); //Note : Alg might work even otherwise. This assertion failure is just a reminder that this needs to be checked.
    //compute Gram-Schmidt-Orthogonalization.
    ZZ_mat<mpz_t> Empty_mat;
    ZZ_mat<mpz_t> basis2 = basis; //need to copy, as BGSO is not const-specified...
    MatGSO<Z_NR<mpz_t>, FP_NR<double>> BGSO(basis2, Empty_mat, Empty_mat, 0);
    BGSO.update_gso();

    FP_NR<double> entry;

    //for (int i=0; i<basis.get_rows(); i++)
    //{
// 	for (int j=0; j<basis.get_rows(); j++)
//        	//cout << BGSO.get_gram(entry, j, j) << endl;
//		cout << (BGSO.get_r(entry, j, j)) << "  " << log(BGSO.get_r(entry, j, j)) << endl;
//        cout << endl;
//   }

    // returns det(B)^{2/dim}

    FP_NR<double> root_det2 = BGSO.get_root_det (1, basis.get_rows());
    FP_NR<double> log_det2 = BGSO.get_log_det (1, basis.get_rows());
    //cout << "root_det2: " << root_det2 << endl;
    //cout << "log_det2: " << log_det2 << endl;

    //lambda_1^2 = n * det(B)^{2/n}

    FP_NR<double> MinkBound_double = 0.074 * root_det2 * static_cast<double> (basis.get_rows() ); //technically, we need to multiply by Hermite's constant in dim n here. We are at least missing a constant factor here.
    //DUE TO [KL79], the best know multiple (for the squared norm) whould be 1/(pi* exp(1)*2^{2*0.099} ~ 0.102) for n->infinity. Blichfeldt's bound: 1 / (pi*exp(1))=0.117.
    // Darmstadt's challenge suggests: 1.10 / (2*pi*exp(1)) = 0.0644;

    //cout << "after MinkBound_double is assigned... " << endl;
    Z_NR<mpz_t> Minkowski;
    Minkowski.set_f(MinkBound_double);
    cout << "Mink. bound = " << Minkowski << endl;
    return Minkowski;
}

//}
#endif
