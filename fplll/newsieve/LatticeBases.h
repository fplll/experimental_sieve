/**
  The purpose of this file is to provide glue code for the lattice basis processing capabilities
  of non-sieve fplll with the sieve part. This processing is mainly to compute GSOs.
  Some glue is needed because we do not neccessarily want to use Z_NR - types (exclusively)
  inside the sieving code, whereas the interface to compute GSOs is heavily intertwined with Z_NR.
*/

#ifndef GAUSS_SIEVE_LATTICE_BASES_H
#define GAUSS_SIEVE_LATTICE_BASES_H

#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include "gmp.h"
#include "gmpxx.h"
#include <mutex>
#include "DebugAll.h"
#include "PlainLatticePoint.h"
#include <iostream>

namespace GaussSieve
{

/**
  LatticeBasisType is the type used to represent lattice bases.
  We support converting to PlainLatticePoints and extracting GSO information.
*/

template<class SieveTraits, bool MT, bool Enabled=true> class SieveLatticeBasis;


// unusable default, until someone implements more general GSOs...
template<class SieveTraits, bool MT, bool Enabled>
class SieveLatticeBasis
{
  // This must never be instantiated anyway, but if SieveTraits is ill-formed, the SFINAE selections
  // below may fail.
  static_assert(Enabled,"");
  static_assert(SieveTraits::IsSieveTraitsClass::value,"Invalid Traits class");
  public: SieveLatticeBasis(...) = delete;
};

/**
  Specialization for classes where SieveTraits::InputBasisType is recognized as a ZZMat-class.
  This essentially wraps around the fplll GSO class interfac(es).

  Note that the GSO class interfaces have not been designed with multi-threading in mind.
  In particular, they perfom lazy evaluation and store values already computed.
  So even what should be read-only operations are probably not thread-safe.
  For that reason, we only use the GSO methods during creation of the SieveLatticeBasis object and
  store the results in our own classes. This is quite suboptimal and ugly, but it's not dominating
  the cost of the whole algorithm anyway.
*/

// TODO: Make GSO object local to constructor.

template<class SieveTraits, bool MT>
class SieveLatticeBasis< SieveTraits, MT,true> //TODO: LAST ARGUMENT: NOT CORRECT
{
  static_assert(IsZZMatClass<typename SieveTraits::InputBasisType>::value, ""); // should *never* happen.
  public:
  using InputBasisType = typename SieveTraits::InputBasisType;
  // Note : InputET_NOZNR is not Z_NR - wrapped (because of the way ZZMat works...)
  using InputET_NOZNR  = typename IsZZMatClass<InputBasisType>::GetET;
  using InputET        = fplll::Z_NR<InputET_NOZNR>;
  // Same as ET_NOZNR, but ET_NOZNRFixed may be mpz_class instead of mpz_t
  using InputET_NOZNRFixed = typename UnZNR<InputET>::type;

  using OutputET       = typename SieveTraits::EntryType;

  using DimensionType  = typename SieveTraits::DimensionType;

  // This one is hard-coded to PlainLatticePoint
  using BasisVectorType= PlainLatticePoint<OutputET,SieveTraits::get_nfixed>;

  using GSOType  = fplll::MatGSO<InputET, fplll::FP_NR<double>>;

  // Note: We copy the basis into originial_basis.
  // This is because the GSO object actually uses a reference and modifies original_basis.
  explicit SieveLatticeBasis(InputBasisType const &input_basis):
    original_basis(input_basis),
    ambient_dimension(input_basis.get_cols()),
    lattice_rank(input_basis.get_rows()),
    // u,u_inv intentionally uninitialized
//    GSO(original_basis, u,u_inv, fplll::MatGSOInterfaceFlags::GSO_INT_GRAM),
    mu_matrix(lattice_rank,std::vector<double>(lattice_rank) ),
    g_matrix (lattice_rank,std::vector<InputET_NOZNRFixed>(lattice_rank) ),
    basis_vectors(nullptr)
  {
    fplll::Matrix<InputET> u, u_inv; //, g;
    GSOType GSO(original_basis, u, u_inv, fplll::MatGSOInterfaceFlags::GSO_INT_GRAM);
    GSO.update_gso();  // todo: raise exception in case of error
    // We compute these at creation to simplify thread-safety issues.
    compute_mu_matrix(GSO);
    compute_g_matrix(GSO);

    // extract and convert the actual lattice vectors.

    bool s = BasisVectorType::class_init(MaybeFixed<SieveTraits::get_nfixed>{ambient_dimension});
    assert(s); // TODO: Clean up and throw exception instead.
    basis_vectors = new BasisVectorType[lattice_rank];
    for(uint_fast16_t i=0;i<lattice_rank;++i)
    {
      basis_vectors[i] = make_from_znr_vector<BasisVectorType>(input_basis[i], ambient_dimension);
    }
    maxbistar2=GSO.get_max_bstar().get_d();
  }


  SieveLatticeBasis(SieveLatticeBasis const &old) = delete;
  SieveLatticeBasis(SieveLatticeBasis &&old) = default;
  SieveLatticeBasis & operator= (SieveLatticeBasis const &old) = delete;
  SieveLatticeBasis & operator= (SieveLatticeBasis &&old) = default;


  ~SieveLatticeBasis()
  {
    delete[] basis_vectors;
    BasisVectorType::class_uninit();
  }


  // Note: Const-correctness is strange wrt. the fplll::GSO classes.
  // We "patch" this up by marking GSO mutable.
  // TODO: Patch the MatGSO class upstream instead (use mutable for lazy evaluation...)

  // IMPORTANT: get_mu_matrix()[i][j] is only meaningful for j>i!

  // helper function that precomputes, converts and stores the whole mu_matrix:
  void compute_mu_matrix(GSOType &GSO) // Note: GSO is not passed as const-ref.
  {
    // use GSO's capabilities and convert to non-FP_NR - type
    fplll::Matrix<fplll::FP_NR<double>> ZNR_mu = GSO.get_mu_matrix();
    for (uint_fast16_t i = 0;i < lattice_rank;++i)
    {
      for(uint_fast16_t j=0; j < lattice_rank;++j)
      {
        mu_matrix[i][j] = ZNR_mu[i][j].get_d();
      }
    }
    return;
  }

  void compute_g_matrix(GSOType &GSO) // Note: GSO is not passed as const-ref.
  {
    fplll::Matrix<InputET> const gmatrix_GSO = GSO.get_g_matrix(); // class returned by GSO
    // convert fplll::Matrix<ET> to vector<vector<ET_NOZNRFixed>>
    for(uint_fast16_t i=0; i< lattice_rank; ++i)
    {
      for(uint_fast16_t j=0; j<lattice_rank; ++j)
      {
          g_matrix[i][j] = static_cast<InputET_NOZNRFixed> ( gmatrix_GSO(i,j).get_data() );
      }
    }
  }

  std::vector<std::vector<double>> get_mu_matrix() const
  {
    return mu_matrix;
  }

  double get_maxbistar2() const
  {
    return maxbistar2;
  }

  std::vector<std::vector<  InputET_NOZNRFixed   > > get_g_matrix() const
  {
    return g_matrix;
  }

  // returns (i,j)th entry of mu. Assumes j>i
  double get_mu(uint_fast16_t i, uint_fast16_t j) const
  {
  #ifdef DEBUG_SIEVE_LOWERTRIANGULAR_MUG
  assert(j>i);
  #endif
    return (mu_matrix)[i][j];
  }

  // returns (i,j)th entry of g. Assumes j>=i
  InputET_NOZNRFixed get_g(uint_fast16_t i,uint_fast16_t j) const
  {
  #ifdef DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    assert(j>=i);
  #endif // DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    return (g_matrix)[i][j];
  }

  // Returns ith basis vector.
  BasisVectorType const & get_basis_vector(uint_fast16_t i) const
  {
    assert(i<lattice_rank);
    return basis_vectors[i];
  }

  InputET_NOZNRFixed get_minkowski_bound() const
  {
    return 0;
  }

  private:

  InputBasisType original_basis;
  public:
  DimensionType const ambient_dimension;
  uint_fast16_t const lattice_rank;      // Technically, just number of vectors.
                                  // We don't verify linear independence ourselves.
                                  // (even though GSO computation does, probably)
  private:
//  fplll::Matrix<InputET> u, u_inv; //, g;
//  fplll::MatGSO<InputET, fplll::FP_NR<double>> GSO;
  // precomputed on demand:
  std::vector<std::vector<double>> mu_matrix;
  std::vector<std::vector<InputET_NOZNRFixed>> g_matrix;
  // Note: We use a dynamically allocated array here rather than std::vector.
  // The reason is that PlainLatticePoint<...> might contain static data that needs to be set
  // before we (even default-) construct any PlainLatticePoint.
  // This includes constructing PlainLatticePoints from the initializer list of SieveLatticeBasis.
  BasisVectorType *basis_vectors;
  double maxbistar2;
};

} // namespace

#endif
