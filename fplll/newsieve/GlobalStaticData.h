#ifndef GLOBAL_STATIC_DATA_H
#define GLOBAL_STATIC_DATA_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"

// This struct holds all data that is used to initialize the static data associated to our various classes.
// Note that it need not be a singleton itself

namespace GaussSieve{

namespace TraitHelpers
{
template<class T> using Predicate_StaticInitializerArg = typename T::StaticInitializerArgTag;
template<class T> using Predicate_DefaultStaticInitializer = typename T::HasDefaultStaticInitializer;
}
template<class T> using IsArgForStaticInitializer =
    mystd::is_detected<TraitHelpers::Predicate_StaticInitializerArg,T>;
template<class T> using IsStaticInitializerDefaulted =
    mystd::is_detected<TraitHelpers::Predicate_DefaultStaticInitializer,T>;


// forward declarations:
template<class T> class StaticInitializer;
template<class T> class DefaultStaticInitializer;
template<class DimensionType> struct StaticInitializerArg;

/**
  This is the default initializer, which does nothing apart from counting number of instances.
  Note that typically, we inherit from this, which is the reason why it is templated by T.
  Usage: This class shall only be used by StaticInitializer<T>
*/

template<class T>
class DefaultStaticInitializer
{
  friend StaticInitializer<T>;
  public:
#ifndef DEBUG_SIEVE_LP_INIT
  static bool constexpr is_initialized() { return true; } // may be overloaded
#else
  static bool is_initialized(){ return user_count > 0; }; // Does an object exist?
#endif
  static unsigned int get_user_count() { return user_count; }
  static unsigned int user_count; // counts the number of objects of this type that exist, essentially.
  private:
  explicit DefaultStaticInitializer(){ ++user_count; };

  template<class X,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<X>)>
  explicit DefaultStaticInitializer(X const &) : DefaultStaticInitializer(){}


  template<class X,TEMPL_RESTRICT_DECL2(std::is_integral<X>)>
  [[deprecated]]explicit DefaultStaticInitializer(X const &) : DefaultStaticInitializer(){}
  template<int nfixed, class IntType>
  [[deprecated]]explicit DefaultStaticInitializer(MaybeFixed<nfixed,IntType> const &) : DefaultStaticInitializer(){}

  ~DefaultStaticInitializer()
  {
    --user_count;
  }
};
// initialize static data this class:
template<class T> unsigned int DefaultStaticInitializer<T>::user_count = 0;

/**
  StaticInitializer is a RAII class that manages static initialization of Data.
  Usage:
  Create an object StaticInitializer<T> init_T(args);
  During the livetime of the object init_T, the class T is considered initialized.

  The exact implementation of StaticInitializer<T> depends on T.
  The implementation below is a default implementation that essentially does nothing.
  It may be specialized for certain T's.
  Safeguard: In order to enable the default implementation, the class T has to declare a
  public member typedef HasDefaultStaticInitializer, e.g.
  using HasDefaultStaticInitializer = std::true_type; otherwise, you get a compilation error.

  The type of args may be any integral class (deprecated) or
  any class that has a public typedef StaticInitializerArgTag.
  It is the job of the the specialization to ensure that arg has the required data.
  SieveTraits has a GlobalStaticData class supposed to be usable for all such args.
*/

// StaticInitializer<T> for classes T that have the IsStaticInitializerDefaulted Trait
template<class T>
class StaticInitializer : public DefaultStaticInitializer<T>
{
  static_assert(IsStaticInitializerDefaulted<T>::value,"Missing Static Initializer");
  explicit StaticInitializer() = default;
  template<class X,TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<X>)>
  explicit StaticInitializer(X const &) : StaticInitializer() {}
  template<class X,TEMPL_RESTRICT_DECL2(std::is_integral<X>)>
  [[deprecated]] explicit StaticInitializer(X const &) : StaticInitializer() {}
  template<int nfixed, class IntType>
  [[deprecated]] explicit StaticInitializer(MaybeFixed<nfixed,IntType> const &) : StaticInitializer() {}
};



/**
  StaticInitializerArg<DimensionType> encapsulates an argument of type
  DimensionType.
  It is meant to be an argument for static initializers.
*/

template<class DimensionType>
struct StaticInitializerArg
{
  using StaticInitializerArgTag = std::true_type;
  DimensionType const dim;
//  unsigned int const dim_int;
  constexpr StaticInitializerArg(DimensionType const &new_dim) : dim(new_dim) {}
};

} // namespace

#endif
