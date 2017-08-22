#ifndef PLAIN_LATTICE_POINT_H
#define PLAIN_LATTICE_POINT_H

#include "DebugAll.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include "assert.h"
#include <array>
#include <vector>
#include <iostream>

namespace GaussSieve
{

// most simple Lattice Point that just wraps around a vector / array of ET's.
// dimension is static
template <class ET, int nfixed> class PlainLatticePoint;

template <class ET, int nfixed> class LatticePointTraits<PlainLatticePoint<ET, nfixed>>
{
public:
  using AuxDataType             = MaybeFixed<nfixed>;
  using ScalarProductReturnType = ET;
  using CoordinateVector        = std::true_type;
  using CoordinateAccess        = std::true_type;
  using AbsoluteCoos            = std::true_type;
  using CoordinateType          = ET;
};

template <class ET, int nfixed>
class PlainLatticePoint : public GeneralLatticePoint<PlainLatticePoint<ET, nfixed>>
{
  static_assert(nfixed >= 0);  // we have to specialize for nfixed==-1
public:
  using LatticePointTag         = std::true_type;
  using AuxDataType             = typename GetAuxDataType<PlainLatticePoint>::type;
  using ScalarProductReturnType = ET;
  using Container               = std::array<ET, nfixed>;
  static bool class_init(AuxDataType const aux_data)
  {
    static_assert(aux_data == dim);
// dim = aux_data;
#ifdef DEBUG_SIEVE_LP_INIT
    class_initialized = true;
#endif
    return true;
  };
  static constexpr MaybeFixed<nfixed> get_dim() { return dim; };

  static bool class_uninit() { return true; }

#ifdef DEBUG_SIEVE_LP_INIT
  explicit PlainLatticePoint()
  {
    // The extra () are needed, because assert is a macro and the argument contains a ","
    assert((PlainLatticePoint<ET, nfixed>::class_initialized));
  }
#else
  explicit PlainLatticePoint() = default;
#endif
  explicit PlainLatticePoint(MaybeFixed<nfixed>)
  {
#ifdef DEBUG_SIEVE_LP_INIT
    assert((PlainLatticePoint<ET, nfixed>::class_initialized));
#endif
  };
  PlainLatticePoint(PlainLatticePoint const &old) = delete;
  PlainLatticePoint(PlainLatticePoint &&old)      = default;
  PlainLatticePoint &operator=(PlainLatticePoint const &other) = delete;
  PlainLatticePoint &operator=(PlainLatticePoint &&other) = default;
  ET &operator[](unsigned int idx) { return data[idx]; };
  ET const &operator[](unsigned int idx) const { return data[idx]; };
  static std::string class_name() { return "Plain Lattice Point, fixed dim"; };

private:
  static constexpr MaybeFixed<nfixed> dim = MaybeFixed<nfixed>(nfixed);
#ifdef DEBUG_SIEVE_LP_INIT
  static bool class_initialized;
#endif  // DEBUG_SIEVE_LP_INIT
  Container data;
};

template <class ET>
class PlainLatticePoint<ET, -1> : public GeneralLatticePoint<PlainLatticePoint<ET, -1>>
{
public:
  using LatticePointTag         = std::true_type;
  using AuxDataType             = typename GetAuxDataType<PlainLatticePoint>::type;
  using ScalarProductReturnType = ET;
  using Container               = std::vector<ET>;
  using GeneralLatticePoint<PlainLatticePoint<ET, -1>>::get_dim;
  static bool class_init(AuxDataType const aux_data)
  {
    if( (user_counter>0) && (dim!=aux_data) )
    {
#ifdef DEBUG_SIEVE_LP_INIT
      assert(class_initialized);
      std::cerr << "Warning: Class initialization failed for " << class_name()
                << std::endl << std::flush;
      return false;
#endif
    }
    dim = aux_data;
#ifdef DEBUG_SIEVE_LP_INIT
    class_initialized = true;
#endif
    ++user_counter;
    return true;
  };

  static bool class_uninit()
  {
    assert(user_counter>0);
    --user_counter;
#ifdef DEBUG_SIEVE_LP_INIT
    assert(class_initialized);
    class_initialized=(user_counter>0);
#endif
    return (user_counter==0);
  }

  static MaybeFixed<-1> get_dim() { return dim; };

  explicit PlainLatticePoint() : data(static_cast<unsigned int>(get_dim()))
  {
#ifdef DEBUG_SIEVE_LP_INIT
    assert((PlainLatticePoint<ET, -1>::class_initialized));
#endif
  };

  explicit PlainLatticePoint(MaybeFixed<-1> const dim) : data(static_cast<unsigned int>(get_dim()))
  {
#ifdef DEBUG_SIEVE_LP_INIT
    assert((PlainLatticePoint<ET, -1>::class_initialized));
#endif
#ifdef DEBUG_SIEVE_LP_MATCHDIM
    assert(dim == get_dim());
#endif
  };

  PlainLatticePoint(PlainLatticePoint const &old) = delete;
  PlainLatticePoint(PlainLatticePoint &&old)      = default;
  PlainLatticePoint &operator=(PlainLatticePoint const &other) = delete;
  PlainLatticePoint &operator=(PlainLatticePoint &&other) = default;
  ET &operator[](unsigned int idx) { return data[idx]; };
  ET const &operator[](unsigned int idx) const { return data[idx]; };

  static std::string class_name() { return "Plain Lattice Point, variable dim"; };

private:
  static MaybeFixed<-1> dim;
  static unsigned int user_counter;
#ifdef DEBUG_SIEVE_LP_INIT
  static bool class_initialized;
#endif  // DEBUG_SIEVE_LP_INIT
  Container data;

public:
  //    std::ostream& write_to_stream(std::ostream &os) const { return os;};
};

template <class ET, int nfixed> constexpr MaybeFixed<nfixed> PlainLatticePoint<ET, nfixed>::dim;
// = MaybeFixed<nfixed>(nfixed<0?0:nfixed);
template <class ET> MaybeFixed<-1> PlainLatticePoint<ET, -1>::dim = MaybeFixed<-1>(0);

template <class ET> unsigned int PlainLatticePoint<ET,-1>::user_counter = 0;

#ifdef DEBUG_SIEVE_LP_INIT
template <class ET, int nfixed> bool PlainLatticePoint<ET, nfixed>::class_initialized = false;
template <class ET> bool PlainLatticePoint<ET, -1>::class_initialized = false;
#endif  // DEBUG_SIEVE_LP_INIT

}  // end of namespace

#endif
