#ifndef LATTICE_VECTOR_CLASS_NEW_CPP
#define LATTICE_VECTOR_CLASS_NEW_CPP

#include "LatticePointsNew.h"


#include "ExactLatticePoint.cpp"
#include "ApproximateLatticePoint.cpp"
#include "CompressedPoint.cpp"


//ApproxType do_approximate(X val, int delta) returns 2^(-delta)*val, converted to appropriate ApproxType.

template<class ApproxType, class ET>
ApproxType LatticeApproximationsNew::do_approximate( typename enable_if< is_same<ET, Z_NR<mpz_t> >::value , ET>::type const & val, signed int const delta)
//ApproxType do_approximate(Z_NR<mpz_t> const & val, signed int const delta )
{
    assert(delta>=0);
    Z_NR<mpz_t> temp;
    temp.div_2si (val,delta);
    return temp.get_si();
}

//Note that nr_Z_l 's div_2si uses >> and << operators on signed values.
//This would lead to undefined behaviour (as explicitly stated in nr_Z) in our case as of C++17 (although it would work with most compilers).
//(This is the single reason we don't default to Z_NR's capabilities here.


template<class ApproxType, class ET>
ApproxType LatticeApproximationsNew::do_approximate( typename enable_if< is_same<ET, Z_NR<long> >::value , ET>::type const & val, signed int const delta)
//ApproxType do_approximate(Z_NR<long> const & val, signed int const delta)
{
    assert(delta>=0);
    return ( val.get_data() / (2L << delta ) );
}

template<class ApproxType, class ET>
ApproxType LatticeApproximationsNew::do_approximate( typename enable_if< is_same<ET, Z_NR<double> >::value , ET>::type const & val, signed int const delta)
//Read : ApproxType do_approximate(Z_NR<double> const & val, signed int const delta)
{
    return (std::ldexp(val.get_data(),-delta) ); //Note : Conversion from double to integral type rounds towards zero. This is needed to prevent overflows.
}

/*get_exponent :
returns minimal r, s.t. |val| < 2^r (or INTMIN, if val = 0)
*/

template <class ET> //fallback version, should never be called anyway.
[[ deprecated ("Using badly supported type") ]] signed int LatticeApproximationsNew::get_exponent(ET const & val)
{

    if(val==0) return std::numeric_limits<signed int>::min();
    ET absval = abs(val);
    if(absval >= 1)
    {
        signed int res=0;
        while(absval>=1)
        {
            ++res;
            absval = absval /2;
        }
        return res;
    }
    else
    {
        signed int res=1;
        while(absval<1)
        {
            --res;
            absval = absval*2;
        }
        return res;
    }
}

//instantiate for appropriate types
template class ExactLatticePoint<Z_NR<long>>;
template class ExactLatticePoint<Z_NR<double>>;
template class ExactLatticePoint<Z_NR<mpz_t>>;
template class ApproximateLatticePoint<Z_NR<long>>;
template class ApproximateLatticePoint<Z_NR<double>>;
template class ApproximateLatticePoint<Z_NR<mpz_t>>;
template class CompressedPoint<Z_NR<long>,false>;
template class CompressedPoint<Z_NR<mpz_t>,false>;
template class CompressedPoint<Z_NR<double>,false>;

#endif
