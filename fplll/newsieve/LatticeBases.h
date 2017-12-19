/**
  The purpose of this file is to provide glue code for the lattice basis processing capabilities
  of non-sieve fplll with the sieve part. This processing is mainly to compute GSOs.
  Some glue is needed because we do not neccessarily want to use Z_NR - types (exclusively)
  inside the sieving code, whereas the interface to compute GSOs is heavily intertwined with Z_NR.
*/

#ifndef GAUSS_SIEVE_LATTICE_BASES_H
#define GAUSS_SIEVE_LATTICE_BASES_H

#include "DefaultIncludes.h"
#include "PlainLatticePoint.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include "gmp.h"
#include "gmpxx.h"

namespace GaussSieve
{

/**
  LatticeBasisType is the type used to represent lattice bases.
  A basis of this type is stored inside the sieve.
  We support converting to an array of PlainLatticePoints and extracting GSO information.
*/

// NOTE: The last template parameter Enabled is a dummy parameter and is always void.
// Its purpose is to allow partial specializations that depend on properties of SieveTraits.
//
// Notably, you can do partial specializations
// template<class SieveTraits, bool MT>
// class SieveLatticeBasis<SieveTraits, MT, mystd::enable_if_t< CONDITION > >
// {
// ...
// };
// where CONDITION is a bool that depends on SieveTraits and MT.
// Since template argument deduction takes places before template argument substitution and SFINAE,
// this has the effect that the specialization is only considered iff CONDITION is true.

// This way, we can write various versions, depending on whether we receive Z_NR - input types or
// not. Currently, we only take ZZ_Mat - classes

template <class SieveTraits, bool MT, class Enabled = void> class SieveLatticeBasis;

// unusable default (by design)
// clang-format off
template<class SieveTraits, bool MT, class Enabled>
class SieveLatticeBasis
// clang-format on
{
  // This must never be instantiated.
  static_assert(std::is_void<Enabled>::value == false, "");  // this will always trigger!
  static_assert(SieveTraits::IsSieveTraitsClass::value, "Invalid Traits class");

public:
  SieveLatticeBasis(...) = delete;
};

/**
  Specialization for classes where SieveTraits::InputBasisType is recognized as a ZZMat-class.
  This is the only implementation we currently have.
  This essentially wraps around the fplll GSO class interfac(es).

  Note that the GSO class interfaces may not have been designed with multi-threading in mind.
  In particular, they perfom lazy evaluation and store values already computed.
  So even what should be read-only operations are probably not thread-safe.
  For that reason, we only use the GSO methods during creation of the SieveLatticeBasis object and
  store the results in our own classes. This is quite suboptimal and ugly, but it's not dominating
  the cost of the whole algorithm anyway.
*/

// TODO: Make GSO object local to constructor.

// clang-format off
template <class SieveTraits, bool MT>
class SieveLatticeBasis<SieveTraits, MT,  // next param selectively enables this specialization:
    mystd::enable_if_t<IsZZMatClass<typename SieveTraits::InputBasisType>::value> >
// clang-format on
{

private:
  // clang-format off
  // Transform the received InputBasisType into something appropriate...
  using InputBasisType     = typename SieveTraits::InputBasisType;

  // Note : InputET_NOZNR is not Z_NR - wrapped (because of the way ZZMat works...)
  using InputET_NOZNR      = typename IsZZMatClass<InputBasisType>::GetET;
  using InputET            = fplll::Z_NR<InputET_NOZNR>;

  // Same as ET_NOZNR, but ET_NOZNRFixed may be mpz_class instead of mpz_t
  using InputET_NOZNRFixed = PossiblyMpztToMpzClass<InputET_NOZNR>;
  using OutputET           = typename SieveTraits::LengthType;
  // clang-format on

  using DimensionType               = typename SieveTraits::DimensionType;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;

  // This one is hard-coded to a PlainLatticePoint
  using BasisVectorType = PlainLatticePoint<OutputET, SieveTraits::get_nfixed>;

  using GSOType = fplll::MatGSO<InputET, fplll::FP_NR<double>>;

public:
  // Note: We copy the basis into originial_basis.
  // This is because the GSO object actually uses a reference and modifies original_basis.
  // clang-format off
  explicit SieveLatticeBasis(InputBasisType const &input_basis,
                             GlobalStaticDataInitializer const &static_data)
      : original_basis(input_basis),
        ambient_dimension(input_basis.get_cols()),
        init_basis_vector_type(static_data),
        lattice_rank(input_basis.get_rows()),
        // u,u_inv intentionally uninitialized
        // will be initialized in call to
        // GSO(original_basis, u,u_inv, fplll::MatGSOInterfaceFlags::GSO_INT_GRAM),
        mu_matrix(lattice_rank, std::vector<double>(lattice_rank)),
        g_matrix(lattice_rank, std::vector<InputET_NOZNRFixed>(lattice_rank)),
        basis_vectors(nullptr)

  // clang-format on
  {
    fplll::Matrix<InputET> u, u_inv;  //, g;
    GSOType GSO(original_basis, u, u_inv, fplll::MatGSOInterfaceFlags::GSO_INT_GRAM);
    GSO.update_gso();  // todo: raise exception in case of error
    // We compute these at creation to simplify thread-safety issues.
    compute_mu_matrix(GSO);
    compute_g_matrix(GSO);

    // extract and convert the actual lattice vectors.

    basis_vectors = new BasisVectorType[lattice_rank];
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      basis_vectors[i] = make_from_znr_vector<BasisVectorType>(input_basis[i], ambient_dimension);
    }
    maxbistar2 = GSO.get_max_bstar().get_d();

    compute_minkowski_bound(GSO);

#ifdef PROGRESSIVE
    progressive_bounds.resize(lattice_rank);
    compute_progressive_bounds();
#endif
  }

  // clang-format off
  SieveLatticeBasis(SieveLatticeBasis const &)            = delete;
  SieveLatticeBasis(SieveLatticeBasis &&)                 = default;
  SieveLatticeBasis &operator=(SieveLatticeBasis const &) = delete;
  SieveLatticeBasis &operator=(SieveLatticeBasis &&)      = default;
  // clang-format on

  ~SieveLatticeBasis() { delete[] basis_vectors; }

  // Note: Const-correctness is strange wrt. the fplll::GSO classes.
  // We "patch" this up by marking GSO mutable.
  // TODO: Patch the MatGSO class upstream instead (use mutable for lazy evaluation...)

  // IMPORTANT: get_mu_matrix()[i][j] is only meaningful for j>i!

  // helper function that precomputes, converts and stores the whole mu_matrix:
  void compute_mu_matrix(GSOType &GSO)  // Note: GSO is not passed as const-ref.
  {
    // use GSO's capabilities and convert to non-FP_NR - type
    fplll::Matrix<fplll::FP_NR<double>> ZNR_mu = GSO.get_mu_matrix();
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        mu_matrix[i][j] = ZNR_mu[i][j].get_d();
      }
    }
    return;
  }

  void compute_g_matrix(GSOType &GSO)  // Note: GSO is not passed as const-ref.
  {
    fplll::Matrix<InputET> const gmatrix_GSO = GSO.get_g_matrix();  // class returned by GSO
    // convert fplll::Matrix<ET> to vector<vector<ET_NOZNRFixed>>
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        g_matrix[i][j] = static_cast<InputET_NOZNRFixed>(gmatrix_GSO(i, j).get_data());
      }
    }
  }

  std::vector<std::vector<double>> get_mu_matrix() const { return mu_matrix; }

  double get_maxbistar2() const { return maxbistar2; }

  std::vector<std::vector<InputET_NOZNRFixed>> get_g_matrix() const { return g_matrix; }

  // returns (i,j)th entry of mu. Assumes j>i
  double get_mu(uint_fast16_t i, uint_fast16_t j) const
  {
#ifdef DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    assert(j > i);
#endif
    return (mu_matrix)[i][j];
  }

  // returns (i,j)th entry of g. Assumes j>=i
  InputET_NOZNRFixed get_g(uint_fast16_t i, uint_fast16_t j) const
  {
#ifdef DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    assert(j >= i);
#endif  // DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    return (g_matrix)[i][j];
  }

  // Returns ith basis vector.
  BasisVectorType const &get_basis_vector(uint_fast16_t i) const
  {
    assert(i < lattice_rank);
    return basis_vectors[i];
  }

  // Computes a bound on (lambda_1)^2 for n-rank lattice:
  //    const * sqrt(n) * det(B)^{1/n}
  // DUE TO [KL79], the best know const (for the squared norm) is 1/(pi* exp(1)*2^{2*0.099} ~ 0.102)
  // Blichfeldt's bound: 1 / (pi*exp(1))=0.117.
  // Darmstadt's challenge suggests: 1.10 / (2*pi*exp(1)) = 0.0644;
  void compute_minkowski_bound(GSOType &GSO)
  {
    // returns det(B)^{2/dim}
    // fplll::FP_NR<double> root_det2 =  GSO.get_root_det(1, lattice_rank);
    double root_det     = GSO.get_root_det(1, lattice_rank).get_d();
    double mink_bound_d = 0.076 * root_det * static_cast<double>(lattice_rank);
    mink_bound          = static_cast<InputET_NOZNRFixed>(mink_bound_d);
    std::cout << "mink_bound is set to: " << mink_bound << std::endl;
  }

#ifdef PROGRESSIVE
  void compute_progressive_bounds()
  {
    //std::cout << "lattice_rank = " << lattice_rank << std::endl;
    progressive_bounds[0] = log( convert_to_double(get_g(0,0)) );
    double accumulate_sum = progressive_bounds[0];
    for (unsigned int i = 1; i<lattice_rank; ++i)
    {
      accumulate_sum+=log(convert_to_double (get_g(i,i)));
      progressive_bounds[i] = exp( (1. / (i+1.) ) * accumulate_sum );
      //sstd::cout << i <<" accumulate_prod "<< accumulate_sum << " progressive_bounds[i] = " << progressive_bounds[i] << std::endl;
    }
  }
#endif

  InputET_NOZNRFixed get_minkowski_bound() const { return mink_bound; }

private:
  InputBasisType original_basis;

public:
  DimensionType const ambient_dimension;
  StaticInitializer<BasisVectorType> const init_basis_vector_type;
  uint_fast16_t const lattice_rank;  // Technically, just number of vectors.
                                     // We don't verify linear independence ourselves.
                                     // (even though GSO computation does, probably)
#ifdef PROGRESSIVE
  std::vector<double> progressive_bounds; // progressive_bounds[i] stores GH^2 for
                                          // L[n_start+i] (to determine when we increase the rank)
#endif
private:
  //  fplll::Matrix<InputET> u, u_inv; //, g;
  //  fplll::MatGSO<InputET, fplll::FP_NR<double>> GSO;
  // precomputed on demand:
  // clang-format off
  std::vector<std::vector<double            >> mu_matrix;
  std::vector<std::vector<InputET_NOZNRFixed>> g_matrix;
  // clang-format on
  InputET_NOZNRFixed mink_bound;

  // Note: We use a dynamically allocated array here rather than std::vector.
  // The reason is that PlainLatticePoint<...> might contain static data that needs to be set
  // before we (even default-) construct any PlainLatticePoint.
  // This includes constructing PlainLatticePoints from the initializer list of SieveLatticeBasis.

  // TODO: This argument should no longer apply with our RAII style initializers.
  BasisVectorType *basis_vectors;
  double maxbistar2;
};  // end of class

}  // namespace GaussSieve

#endif
