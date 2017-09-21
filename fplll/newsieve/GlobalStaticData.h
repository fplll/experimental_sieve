#ifndef GLOBAL_STATIC_DATA_H
#define GLOBAL_STATIC_DATA_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"

// This struct holds all data that is used to initialize the static data associated to our various classes.
// Note that it need not be a singleton itself

namespace GaussSieve{

CREATE_MEMBER_TYPEDEF_CHECK_CLASS(StaticDataInitializerTag,IsStaticDataInitializer);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS(HasDefaultStaticInitializer,IsStaticInitializerDefaulted);

// Initializer for static data.
// This is the default initializer, which does nothing.

template<class T> class StaticInitializer;
template<class T> class DefaultStaticInitializer;
template<class DimensionType> struct StaticDataInitializer;


template<class T>
class DefaultStaticInitializer
{
  public:
#ifndef DEBUG_SIEVE_LP_INIT
  static bool constexpr is_initialized() { return true; } // may be overloaded
#else
  static bool is_initialized(){ return user_count > 0; }; // Does an object exist?
#endif
  static unsigned int get_user_count() { return user_count; }
  static unsigned int user_count; // counts the number of objects of this type that exist, essentially.
  explicit DefaultStaticInitializer(){++user_count;};
  template<class X,TEMPL_RESTRICT_DECL(IsStaticDataInitializer<X>::value)>
  explicit DefaultStaticInitializer(X const &) : DefaultStaticInitializer(){}


  template<class X,TEMPL_RESTRICT_DECL(std::is_integral<X>::value)>
  [[deprecated]]explicit DefaultStaticInitializer(X const &) : DefaultStaticInitializer(){}
  template<int nfixed, class IntType>
  [[deprecated]]explicit DefaultStaticInitializer(MaybeFixed<nfixed,IntType> const &) : DefaultStaticInitializer(){}

  ~DefaultStaticInitializer()
  {
//    assert(init_count>0);
    --user_count;
  }
};
template<class T> unsigned int DefaultStaticInitializer<T>::user_count = 0;

template<class T>
class StaticInitializer : public DefaultStaticInitializer<T>
{
  static_assert(IsStaticInitializerDefaulted<T>::value,"Missing Static Initializer");
  explicit StaticInitializer() = default;
  template<class X,TEMPL_RESTRICT_DECL(IsStaticDataInitializer<X>::value)>
  explicit StaticInitializer(X const &) : StaticInitializer() {}
  template<class X,TEMPL_RESTRICT_DECL(std::is_integral<X>::value)>
  [[deprecated]] explicit StaticInitializer(X const &) : StaticInitializer() {}
  template<int nfixed, class IntType>
  [[deprecated]] explicit StaticInitializer(MaybeFixed<nfixed,IntType> const &) : StaticInitializer() {}
};



template<class DimensionType>
struct StaticDataInitializer
{
  using StaticDataInitializerTag = std::true_type;
  DimensionType const dim;
  unsigned int const dim_int;

  constexpr StaticDataInitializer(DimensionType const &new_dim) : dim(dim), dim_int(dim) {}
};





} // namespace

#endif
