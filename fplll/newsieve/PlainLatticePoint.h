#ifndef PLAIN_LATTICE_POINT_H
#define PLAIN_LATTICE_POINT_H

#include "DefaultIncludes.h"
#include "GlobalStaticData.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"

namespace GaussSieve
{

// most simple Lattice Point that just wraps around a vector / array of ET's.
// dimension is static
template <class ET, int nfixed> class PlainLatticePoint;
template <class ET, int nfixed> class ExactLatticePoint;   // for friend declaration
template <class ET, int nfixed> class HashedLatticePoint;  // for friend declaration

template <class ET, int nfixed> class LatticePointTraits<PlainLatticePoint<ET, nfixed>>
{
public:
  using Trait_ScalarProductStorageType = ET;
  using Trait_AbsoluteCoos             = std::true_type;
  using Trait_CoordinateType           = ET;
  using Trait_CheapNegate              = std::true_type;
  using Trait_InternalRepIsAbsolute    = std::true_type;
  using Trait_InternalRepByCoos        = std::true_type;
  using Trait_InternalRepLinear        = std::true_type;
  using Trait_InternalRep_RW           = std::true_type;
};

// for nfixed >=0 :
template <class ET, int nfixed>
class PlainLatticePoint : public GeneralLatticePoint<PlainLatticePoint<ET, nfixed>>
{
  static_assert(nfixed >= 0, "");  // we have to specialize for nfixed==-1
public:
  friend ExactLatticePoint<ET, nfixed>;
  friend StaticInitializer<PlainLatticePoint<ET, nfixed>>;

  using LatticePointTag = std::true_type;
  using Container       = std::array<ET, nfixed>;
  static constexpr MaybeFixed<nfixed> get_dim() { return dim; };

//  static bool class_uninit() { return true; }

#ifdef DEBUG_SIEVE_LP_INIT
  explicit PlainLatticePoint()
  {
    assert((StaticInitializer<PlainLatticePoint<ET, nfixed>>::is_initialized()));
  }
#else
  explicit PlainLatticePoint() = default;
#endif

  // Creates a lattice point with the stated dimension
  explicit PlainLatticePoint(MaybeFixed<nfixed>) : PlainLatticePoint(){};

  PlainLatticePoint(PlainLatticePoint const &old) = delete;
  PlainLatticePoint(PlainLatticePoint &&old)      = default;
  PlainLatticePoint &operator=(PlainLatticePoint const &other) = delete;
  PlainLatticePoint &operator=(PlainLatticePoint &&other) = default;
  ET &operator[](unsigned int idx) { return data[idx]; };
  ET const &operator[](unsigned int idx) const { return data[idx]; };
  static std::string class_name() { return "Plain Lattice Point, fixed dim"; };

private:
  static constexpr MaybeFixed<nfixed> dim = MaybeFixed<nfixed>(nfixed);

  Container data;
};
// template <class ET, int nfixed> constexpr MaybeFixed<nfixed> PlainLatticePoint<ET, nfixed>::dim;

/*******************************
  Specialization for nfixed==-1
********************************/

template <class ET>
class PlainLatticePoint<ET, -1> : public GeneralLatticePoint<PlainLatticePoint<ET, -1>>
{
public:
  friend ExactLatticePoint<ET, -1>;
  friend HashedLatticePoint<ET, -1>;
  friend StaticInitializer<PlainLatticePoint<ET, -1>>;

  using LatticePointTag = std::true_type;
  using Container       = std::vector<ET>;
  using GeneralLatticePoint<PlainLatticePoint<ET, -1>>::get_dim;

  static MaybeFixed<-1> get_dim() { return dim; };

  explicit PlainLatticePoint() : data(static_cast<unsigned int>(get_dim()))
  {
    assert((StaticInitializer<PlainLatticePoint<ET, -1>>::is_initialized()));
  };

  explicit PlainLatticePoint(MaybeFixed<-1> const dim) : data(static_cast<unsigned int>(get_dim()))
  {
    assert((StaticInitializer<PlainLatticePoint<ET, -1>>::is_initialized()));
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
  Container data;

  // public:
  //    std::ostream& write_to_stream(std::ostream &os) const { return os;};
};
template <class ET> MaybeFixed<-1> PlainLatticePoint<ET, -1>::dim = MaybeFixed<-1>(0);

/**************************
 Static Initializers
**************************/

// Static Initializer:
template <class ET, int nfixed>
class StaticInitializer<PlainLatticePoint<ET, nfixed>>
    : public DefaultStaticInitializer<PlainLatticePoint<ET, nfixed>>
{
  using Parent = DefaultStaticInitializer<PlainLatticePoint<ET, nfixed>>;
  static_assert(nfixed >= 0, "nfixed == -1 is specialized.");

public:
  template <class T, TEMPL_RESTRICT_DECL(IsArgForStaticInitializer<T>::value)>
  StaticInitializer(T const &initializer) : StaticInitializer(initializer.dim)
  {
  }

  StaticInitializer(MaybeFixed<nfixed> const new_dim)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing PlainLatticePoint class with nfixed = "
                                      << nfixed << "Counter is" << Parent::user_count)
  }

  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("DeInitializing PlainLatticePoint class with nfixed = "
                                      << nfixed << "Counter is" << Parent::user_count)
  }
};

template <class ET>
class StaticInitializer<PlainLatticePoint<ET, -1>>
    : public DefaultStaticInitializer<PlainLatticePoint<ET, -1>>
{
  using Parent = DefaultStaticInitializer<PlainLatticePoint<ET, -1>>;

public:
  template <class T, TEMPL_RESTRICT_DECL(IsArgForStaticInitializer<T>::value)>
  StaticInitializer(T const &initializer) : StaticInitializer(initializer.dim)
  {
  }

  StaticInitializer(MaybeFixed<-1> const new_dim)
  {
    assert(Parent::user_count > 0);
    if (Parent::user_count > 1)
    {
      assert((new_dim == PlainLatticePoint<ET, -1>::dim));
      // TODO: Throw exception!
    }
    else
    {
      PlainLatticePoint<ET, -1>::dim = new_dim;
    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing PlainLatticePoint with nfixed = -1"
                                      << "REALDIM = " << new_dim << "Counter is"
                                      << Parent::user_count)
  }

  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("DeInitializing PlainLatticePoint with nfixed = -1"
                                      << "Counter is" << Parent::user_count)
  }
};

}  // end of namespace GaussSieve

#endif
