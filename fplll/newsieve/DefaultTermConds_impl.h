#ifndef DEFAULT_TERM_CONDS_IMPL_H
#define DEFAULT_TERM_CONDS_IMPL_H

#include "DefaultTermConds.h"

inline fplll::Z_NR<mpz_t> GaussSieve::compute_mink_bound(fplll::ZZ_mat<mpz_t> const & basis)
{
    assert(basis.get_rows() == basis.get_cols()); //Note : Alg might work even otherwise. This assertion failure is just a reminder that this needs to be checked.
    //compute Gram-Schmidt-Orthogonalization.
    fplll::ZZ_mat<mpz_t> Empty_mat;
    fplll::ZZ_mat<mpz_t> basis2 = basis; //need to copy, as BGSO is not const-specified...
    fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>> BGSO(basis2, Empty_mat, Empty_mat, 0);
    BGSO.update_gso();

    fplll::FP_NR<double> entry;

    //for (int i=0; i<basis.get_rows(); i++)
    //{
// 	for (int j=0; j<basis.get_rows(); j++)
//        	//cout << BGSO.get_gram(entry, j, j) << endl;
//		cout << (BGSO.get_r(entry, j, j)) << "  " << log(BGSO.get_r(entry, j, j)) << endl;
//        cout << endl;
//   }

    // returns det(B)^{2/dim}

    fplll::FP_NR<double> root_det2 = BGSO.get_root_det (1, basis.get_rows());
    fplll::FP_NR<double> log_det2 = BGSO.get_log_det (1, basis.get_rows());
    //cout << "root_det2: " << root_det2 << endl;
    //cout << "log_det2: " << log_det2 << endl;

    //lambda_1^2 = n * det(B)^{2/n}

    fplll::FP_NR<double> MinkBound_double = 0.074 * root_det2 * static_cast<double> (basis.get_rows() ); //technically, we need to multiply by Hermite's constant in dim n here. We are at least missing a constant factor here.
    //DUE TO [KL79], the best know multiple (for the squared norm) whould be 1/(pi* exp(1)*2^{2*0.099} ~ 0.102) for n->infinity. Blichfeldt's bound: 1 / (pi*exp(1))=0.117.
    // Darmstadt's challenge suggests: 1.10 / (2*pi*exp(1)) = 0.0644;

    //cout << "after MinkBound_double is assigned... " << endl;
    fplll::Z_NR<mpz_t> Minkowski;
    Minkowski.set_f(MinkBound_double);
    std::cout << "Mink. bound = " << Minkowski << std::endl;
    return Minkowski;
}

template<class ET, bool MT, int nfixed>
inline int LengthTerminationCondition<ET,MT,nfixed>::check(Sieve<ET,MT,nfixed> * const sieve)
{
return (sieve -> get_best_length2()<=target_length)?1:0;
}

template<class ET, bool MT, int nfixed>
inline int MinkowskiTerminationCondition<ET,MT,nfixed>::check(Sieve<ET,MT,nfixed>  * const sieve)
{
return (sieve -> get_best_length2()<=target_length)?1:0;
}

template<class ET, bool MT, int nfixed>
inline void MinkowskiTerminationCondition<ET, MT, nfixed>::init(Sieve<ET,MT,nfixed> * const sieve)
{
target_length = GaussSieve::compute_mink_bound(sieve->get_original_basis());
}



#endif
