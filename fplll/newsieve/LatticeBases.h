/**
  The purpose of this file is to provide glue code for the lattice basis processing capabilities
  of non-sieve fplll with the sieve. This processing is mainly to compute GSOs.
  Some glue is needed because we do not neccessarily want to use Z_NR - types (exclusively)
  inside the sieving code, whereas the interface to compute GSO is heavily intertwined with Z_NR.
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


namespace GaussSieve
{

/**
  LatticeBasisType is the type used to represent lattice bases.
*/
template<class SieveTraits,bool MT, bool Enabled=true> class SieveLatticeBasis;


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

template<class SieveTraits, bool MT>
class SieveLatticeBasis< SieveTraits, MT,
                         IsZZMatClass<typename SieveTraits::InputBasisType>::value>
{
  static_assert(IsZZMatClass<typename SieveTraits::InputBasisType>::value, ""); // should *never* happen.
  public:
  // Note : ET_NOZNR is not Z_NR - wrapped.
  using ET_NOZNR       = typename IsZZMatClass<typename SieveTraits::InputBasisType>::GetET;
  using ET             = fplll::Z_NR<ET_NOZNR>;
  using InputBasisType = typename SieveTraits::InputBasisType;
  using DimensionType  = typename SieveTraits::DimensionType;

  // Note: We pass by value (because the GSO constructor may modify its arguments???)
  SieveLatticeBasis(InputBasisType input_basis):
  original_basis(input_basis),
  ambient_dimension(input_basis.get_cols()),
  lattice_rank(input_basis.get_rows()),
  // u,u_inv intentionally uninitialized
  GSO(input_basis, u,u_inv, fplll::MatGSOInterfaceFlags::GSO_INT_GRAM)
  {
    GSO.update_gso();  // todo: raise exception in case of error
  }


  // Note: Const-correctness is broken, really.
  // We "patch" this up by marking GSO mutable.
  // TODO: Patch the MatGSO class upstream instead (use mutable for lazy evaluation...)
  //

  std::vector<std::vector<double>> get_mu_matrix() const
  {
    std::vector< std::vector<double> > mu(lattice_rank, std::vector<double>(lattice_rank));
    fplll::Matrix<fplll::FP_NR<double>> ZNR_mu = GSO.get_mu_matrix();
    for (uint_fast16_t i = 0;i < lattice_rank;++i)
    {
      for(uint_fast16_t j=0; j < lattice_rank;++j)
      {
        mu[i][j] = ZNR_mu[i][j].get_d();
      }
    }
    return mu;
  }

  double get_maxbistar2() const
  {
    return GSO.get_max_bstar().get_d();
  }

  fplll::Matrix<ET> get_g_matrix() const
  {
    return GSO.get_g_matrix();
  }

  private:

  // Note that we are using unsigned int const types here, rather than
  // MaybeFixed<...>
  InputBasisType original_basis;
  DimensionType const ambient_dimension;
  uint_fast16_t const lattice_rank;      // Technically, just number of vectors.
                                  // We don't verify linear independence ourselves.
                                  // (even though GSO computation does, probably)

  fplll::Matrix<ET> u, u_inv; //, g;
  mutable fplll::MatGSO<ET, fplll::FP_NR<double>> GSO;

};

} // namespace

#endif
