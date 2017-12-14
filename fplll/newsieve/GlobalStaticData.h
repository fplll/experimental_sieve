#ifndef GLOBAL_STATIC_DATA_H
#define GLOBAL_STATIC_DATA_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"

// clang-format changes applied on a case-by-case basis manually.

/**
  This class deals with static initializations of static members of classes.
  For reasons of efficiency, keep certain data as *static* members of our classes.
  In particular, the dimension of lattice points is kept as static data.
  ( This is to avoid passing the dimension as a parameter in every operation, which would prevent us
    from using operator+, etc. Alternatively, we could store the dimension inside the points, but
    that seems like a waste (although this is actually done because the containers we use inside
    our current lattice point implementations do that).
    The main reason has to do with multi-threading: If we rely on the value of the size of the
    containers inside the point, we have to be very careful about about synchronization,
    because data races for the size member will cause us to read/write to the wrong memory
    (as opposed to just getting some data wrong -- in particular, it prevents certain strategies of
    optimistic concurrency control)
    Furthermore, we might want to store a slow function of the dimension (which we do not need to
    recompute every time), e.g. some optimal sieve parameters that depend on the dimension...)
  Essentially, we have some global variables this way. While this has certain advantages (simpler
  code, some easier-to-write optimizations), one downside is that we need to care about their
  initialization. Furthermore, this choice implies that
  *** WE MIGHT NOT BE ABLE TO RUN TWO INSTANCES OF THE SIEVE IN PARALLEL***,
  because they will both need to access the same global static data.
  ( It's not that bad: actually, we can, provided the static paramters are the same for all
    instances (which is what happens if you use e.g. extreme pruning). Furthermore, we can
    parallelize inside the loop. In fact, since available memory is a shared constraint, running
    several sieves in parallel is probably a BAD IDEA anyway... )
  Note that classes that have such global data as static member variables *must* actually be class
  templates (or store their static data in another way).
  (The reason is that this makes them "inline", i.e. even if several cpp-files see a
  header with such static data members, they are recognized as the same object by the linker.
  This is ONLY true for members of class templates, not for non-template classes)

  In order to manage such static data, we use RAII - style wrappers, called static initializers.
  e.g. StaticInitializer<LatticePoint>'s constructor will perform the initialization and its
  destructor the (implicit) deinitialization. Essentially the "Is Initialized" state of the
  LatticePoint class is tied to the lifetime of a StaticInitializer<LatticePoint> object.
  (having such RAII wrappers ensures deinitialization and exception safety, guards against
   reinitialization and allows us to avoid some initialization order issues)

  Notably, these static initializers are kept as member objects of e.g. the main_sieve.
  Note that the order of member fields (which equals the order of construction) matters here:
  we keep these as early members and construct them in the initializer list. This way we can have
  e.g. LatticePoint members of main_sieve (if we would perform initialization in the constructor's
  body, we would be in trouble, because this is done after the (implicit) construction of all member
  objects)

  Our static initializers StaticInitializer<ForObject> keeps track of how many of them exist (on a
  For-Object basis), essentially performing some kind of reference counting.
  We allow reinitialization with the same data (which is a no-op). Reinitializing with different
  data is only possibly if the reference count is at 0. This ensures that all previous users are no
  longer using the global data.

  All StaticInitializers<ForClass> shall be derived from DefaultStaticInitializer<ForClass>.
  This parent class template takes care of reference counting.

  Note that template<class AnyClass> StaticInitializer is already defined as a default (which does
  nothing but inconsequetial reference counting). To define a custom static initializer for a given
  class, you need to specialize that template.

  For classes Y for which our generic code assumes that a static initializer is needed, you
  *have to* set a public member typedef HasDefaultStaticInitializer = std::true_type in Y to enable
  the unspecialized template.
  (Otherwise, some friendly static_asserts will remind you ;-). This is to ensure correct usage.)

  TODO: List classes Y for which we use StaticInitializers

  ****StaticInitializerArgs:****

  Note that StaticInitializer<ForObject>'s constructor might need some argument/data that depends on
  ForObject. Since we need to forward such data around (the queue has its own static_initializer
  objects of type StaticInitializer<LatticePointsUsedByQueue>, for example), but do not know exactly
  which data every class needs (and do not want to change the interface in too many places if we add
  something), we encapsulate the data in a class X.

  Technically, this is a concept (i.e. a predicate on a class). Which one we are using is determined
  by SieveTraits.

  A class X is a ArgForStaticInitializer iff it has a public member typedef
  using StaticInitializerArgTag = std::true_type;

  All static initializers for *any* class Y take *any* such ArgForStaticInitializer class X as an
  argument to their constructors StaticInitializer<Y>(X const &x)
  and will look for their required member objects (which depend on Y) such as x.dim or x.memory_pool
  or whatever. It is the task of SieveTraits to select a concrete ArgForStaticInitializer class X
  that fits every static initializer used.
  (The purpose of these rules is that if a static initializer for Y needs some new input argument,
   you can add a single member to X without having to change any interface of any other class but
   X and Y)

   We have and currently use StaticInitializerArg<DimensionType> that contains a member .dim
   of type DimensionType. It is used to pass the (ambient) dimension around.
*/

namespace GaussSieve
{

// forward declarations:
template <class T> class StaticInitializer;
template <class T> class DefaultStaticInitializer;
template <class DimensionType> struct StaticInitializerArg;


/**
  IsArgForStaticInitializer<T>    == std::true_type iff T sets the tag StaticInitializerArgTag.
  IsStaticInitializerDefaulted<T> == std::true_type iff T sets the tag HasDefaultStaticInitializer.

  See above for their meaning.
*/

namespace TraitHelpers  // for IsArgForStaticInitializer and IsStaticInitializerDefaulted
{
// SFINAE - Utilities for IsArgForStaticInitializer and IsStaticInitializerDefaulted
// T::StaticInitializerArgTag / T::HasDefaultStaticInitializer should equal std::true_type
// (or similar), hence the ::value.
// Note that we can get SFINAE - substitution failure in two ways:
// Either T::StaticInitializerArgTag / T::HasDefaultStaticInitializer might not exists or
// ::value might not be true and enable_if_t causes SFINAE.
// clang-format off
template <class T> using Predicate_StaticInitializerArg     = mystd::enable_if_t<T::StaticInitializerArgTag::value>;
template <class T> using Predicate_DefaultStaticInitializer = mystd::enable_if_t<T::HasDefaultStaticInitializer::value>;
// clang-format on
}  // end namespace TraitHelpers

// We use (my)std::is_detected to obtain the actual Trait checkers:
// clang-format off
template<class T> using IsArgForStaticInitializer =
    mystd::is_detected<TraitHelpers::Predicate_StaticInitializerArg, T>;
template<class T> using IsStaticInitializerDefaulted =
    mystd::is_detected<TraitHelpers::Predicate_DefaultStaticInitializer, T>;
// clang-format on

/**
  This is the default initializer, which does nothing apart from counting number of instances.
  StaticInitializer<T> inherits from this, which is the reason why it is templated by T.
  (otherwise the instance-count is not specific to a given T)
  Usage: This class shall only be used by StaticInitializer<T>
*/

// clang-format off
template<class T>
class DefaultStaticInitializer
// clang-format on
{
  friend StaticInitializer<T>;

public:
// TODO: Fix debug symbol usage here. -> only 1 symbol for all such initializations

// returns whether there exists any object of this class (hence whether the class it is supposed
// to to initialize is indeed initialized).
// This is supposed to be overloaded by classes that actually do perform some initialization.
// The version given here always returns true (unless in debug mode), because in cases when the
// initializer does nothing, we do not really need one.
// TODO: reconsider this definition (in non-debug mode)
#ifndef DEBUG_SIEVE_LP_INIT
  static bool constexpr is_initialized() { return true; }  // may be overloaded
#else
  static bool is_initialized() { return user_count > 0; };  // Does an object exist?
#endif
  // counts the number of objects of this type that exist, essentially.
  static unsigned int get_user_count() { return user_count; }

private:
  static unsigned int user_count;
  explicit DefaultStaticInitializer() noexcept { ++user_count; }
  // pointless, but anyway...
  // clang-format off
  explicit DefaultStaticInitializer(DefaultStaticInitializer const &) noexcept { ++user_count; }
  explicit DefaultStaticInitializer(DefaultStaticInitializer &&)      noexcept { ++user_count; }
  // clang-format on

  ~DefaultStaticInitializer() { --user_count; }
};
// initialize static data this class:
template <class T> unsigned int DefaultStaticInitializer<T>::user_count = 0;

/**
  StaticInitializer is a RAII class that manages static initialization of Data.
  Usage:
  Create an object StaticInitializer<T> init_T(args);
  During the livetime of the object init_T, the class T is considered initialized.

  The exact implementation of StaticInitializer<T> depends on T.
  The implementation below is a default implementation that essentially does nothing.
  It may be specialized for certain T's.
  Safeguard: In order to enable the default implementation, the class T has to declare a
  public member typedef HasDefaultStaticInitializer as std::true_type, e.g.
  using HasDefaultStaticInitializer = std::true_type; otherwise, you get a compilation error.

  The type of args may be any integral class (deprecated) or
  any class that has a public typedef StaticInitializerArgTag.
  It is the job of the the specialization to ensure that arg has the required data.
  SieveTraits has a GlobalStaticData class supposed to be usable for all such args.

  (see the explanation above)
*/

// StaticInitializer<T> for classes T that have the IsStaticInitializerDefaulted Trait
// does essentially nothing. Its parent class' constructor increments a counter that is never used.
// clang-format off
template <class T>
class StaticInitializer : public DefaultStaticInitializer<T>
{
  static_assert(IsStaticInitializerDefaulted<T>::value,
                "Missing Static Initializer or forgot to set HasDefaultStaticInitializer trait.");
  explicit StaticInitializer() = default;

  template <class X, TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<X>)>
  explicit StaticInitializer(X const &) noexcept
      : StaticInitializer() {}

  // old versions did not have StaticInitializerArg and were hard-wired to "Dimension"
  template <class X, TEMPL_RESTRICT_DECL2(std::is_integral<X>)>
  [[deprecated]] explicit StaticInitializer(X const &)
      : StaticInitializer() {}

  template <int nfixed, class IntType>
  [[deprecated]] explicit StaticInitializer(MaybeFixed<nfixed, IntType> const &)
      : StaticInitializer() {}
};
// clang-format on

/**
  StaticInitializerArg<DimensionType> encapsulates an argument of type
  DimensionType.
  It is meant to be an argument for static initializers, where dim denotes the ambient dimension in
  which the sieve operates.
*/

// clang-format off
template <class DimensionType>
struct StaticInitializerArg
// clang-format on
{
  using StaticInitializerArgTag = std::true_type;
  DimensionType const dim;
  //  unsigned int const dim_int;
  constexpr StaticInitializerArg(DimensionType const &new_dim) noexcept : dim(new_dim) {}
};

}  // namespace GaussSieve

#endif
