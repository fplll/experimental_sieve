/**
  Very simple lattice point class that stores precomputed values for its length.
  Otherwise, just a wrapper around a vector, really.
*/

#ifndef EXACT_LATTICE_POINT_H
#define EXACT_LATTICE_POINT_H

#include "DebugAll.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include "assert.h"
#include <array>
#include <vector>

#define FOR_FIXED_DIM template <int X = nfixed, typename std::enable_if<X >= 0, int>::type = 0>
#define FOR_VARIABLE_DIM template <int X = nfixed, typename std::enable_if<X == -1, int>::type = 0>

namespace GaussSieve
{

// most simple Lattice Point that just wraps around a vector / array of ET's.
// dimension is static
template <class ET, int nfixed> class ExactLatticePoint;

template <class ET, int nfixed> class LatticePointTraits<ExactLatticePoint<ET, nfixed>>
{
public:
  using AuxDataType             = MaybeFixed<nfixed>;
  using ScalarProductReturnType = ET;
  using CoordinateVector        = std::true_type;
  using CoordinateAccess        = std::true_type;
  using AbsoluteCoos            = std::true_type;
  using CoordinateType          = ET;
  using CheapNorm2              = std::true_type;
};

template <class ET, int nfixed>
class ExactLatticePoint : public GeneralLatticePoint<ExactLatticePoint<ET, nfixed>>
{
public:
  using LatticePointTag         = std::true_type;
  using AuxDataType             = typename GetAuxDataType<ExactLatticePoint>::type;
  using ScalarProductReturnType = ET;
  using Container =
      typename std::conditional<nfixed >= 0, std::array<ET, nfixed>,  // if nfixed >= 0
                                std::vector<ET>>::type;               // if nfixed <0
  static void class_init(AuxDataType const aux_data)
  {
#ifdef DEBUG_SIEVE_LP_INIT
    class_initialized = true;
#endif
    dim = aux_data;
  };

  FOR_FIXED_DIM
  static constexpr MaybeFixed<nfixed> get_dim()
  {
    static_assert(X == nfixed);
    return MaybeFixed<nfixed>(nfixed);
  }

  FOR_VARIABLE_DIM
  static MaybeFixed<-1> get_dim()
  {
    static_assert(nfixed == -1);
    return dim;
  }

#ifdef DEBUG_SIEVE_LP_INIT
  FOR_FIXED_DIM
  explicit ExactLatticePoint()
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((ExactLatticePoint<ET, nfixed>::class_initialized));
  }

  FOR_VARIABLE_DIM
  explicit ExactLatticePoint() : data(static_cast<unsigned int>(get_dim()))
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((ExactLatticePoint<ET, nfixed>::class_initialized));
  }
#else
  FOR_FIXED_DIM
  explicit ExactLatticePoint() = default;

  FOR_VARIABLE_DIM
  explicit ExactLatticePoint() : data(static_cast<unsigned int>(get_dim())) {}
#endif

  FOR_FIXED_DIM
  explicit ExactLatticePoint(MaybeFixed<nfixed>)
  {
#ifdef DEBUG_SIEVE_LP_INIT
    assert((ExactLatticePoint<ET, nfixed>::class_initialized));
#endif
  };

  FOR_VARIABLE_DIM
  explicit ExactLatticePoint(MaybeFixed<nfixed> dim) : data(static_cast<unsigned int>(dim))
  {
    static_assert(nfixed == -1);
#ifdef DEBUG_SIEVE_LP_INIT
    assert((ExactLatticePoint<ET, -1>::class_initialized));
#endif
#ifdef DEBUG_SIEVE_LP_MATCHDIM
    assert(dim == get_dim());
#endif
  }

  ExactLatticePoint(ExactLatticePoint const &old) = delete;
  ExactLatticePoint(ExactLatticePoint &&old)      = default;
  ExactLatticePoint &operator=(ExactLatticePoint const &other) = delete;
  ExactLatticePoint &operator=(ExactLatticePoint &&other) = default;
  ET &operator[](unsigned int idx) { return data[idx]; };
  ET const &operator[](unsigned int idx) const { return data[idx]; };
  static std::string class_name() { return "Exact Lattice Point"; };

  void sanitize() { norm2 = compute_sc_product(*this, *this); }

  ET get_norm2() const { return norm2; }

private:
  static MaybeFixed<nfixed> dim;  // note that for nfixed != -1, this variable is actually unused.
#ifdef DEBUG_SIEVE_LP_INIT
  static bool class_initialized;
#endif  // DEBUG_SIEVE_LP_INIT
  Container data;
  ET norm2;
};

// initialize static data:

template <class ET, int nfixed>
MaybeFixed<nfixed> ExactLatticePoint<ET, nfixed>::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);

#ifdef DEBUG_SIEVE_LP_INIT
template <class ET, int nfixed> bool ExactLatticePoint<ET, nfixed>::class_initialized = false;
#endif  // DEBUG_SIEVE_LP_INIT

}  // end of namespace

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM

#endif
