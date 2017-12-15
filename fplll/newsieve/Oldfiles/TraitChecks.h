#ifndef TRAIT_CHECKS_H
#define TRAIT_CHECKS_H

#warning This file is deprecated

// This class contains macros used to detect the presence of typedefs inside classes.
// This used to use some complicated manual SFINAE constructions.
// Changed to use (indirectly) std::void_t and std::[experimental::]is_detected and friends.
// (enabled via compat.h)

// This means that this file should not be needed at all anymore and will be removed in the future.
// It currently is still in use only in automatic testing code and deactivated Lazy.h

/**
This macro is used to test the presence of a (public) member typedef in a class
Args:   TypeToCheck - typename whose presence to check
        CheckerClassName - Name of the checker class
This macro emits a new template class definition with the name CheckerClassName.
TypeToCheck must not be void (or another incomplete type)

Usage:
CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TypeToCheck, CheckerClassName);
This creates(!) the template class CheckerClassName.

Then CheckerClassName<SomeSuspiciousClass>::value will be true if
SomeSuspicousClass::TypeToCheck exists, false otherwise

The missing semicolon at the end of the macro is intentional.
The user needs to put it to emphasize that this is a declaration.
*/

// clang-format off
/* old version
#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TypeToCheck, CheckerClassName)                           \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> static void                      foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t =                                                                                \
        mystd::bool_constant< !(std::is_void<decltype(foo<ClassToCheck>(0))>::value)>;      \
    static bool constexpr value = value_t::value;                                                  \
    constexpr operator bool() const { return value; };                                             \
  }
*/

// New Version:
// CheckerClassName = is_detected<Op, T>; where Op = typename T::TypeToCheck;

#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TypeToCheck, CheckerClassName)                           \
/*creating a unique name for typename T::TypeToCheck */                                            \
namespace TraitDetectionHelpers                                                                    \
{                                                                                                  \
template<class T>                                                                                  \
using MemberTypedefCheck_##TypeToCheck##_##CheckerClassName = typename T::TypeToCheck;             \
}                                                                                                  \
template<class T>                                                                                  \
using CheckerClassName = mystd::is_detected<TraitDetectionHelpers::MemberTypedefCheck_##TypeToCheck##_##CheckerClassName,T>;

// clang-format on


/**
Similar to the above, creates a checker template class that checks wether
TypeToCheck exists and is equal to TypeShouldBe
*/

// clang-format off

#define CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(TypeToCheck, TypeShouldBe, CheckerClassName)      \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> static void                      foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t =                                                                                \
        mystd::bool_constant<                                                               \
                               std::is_same<TypeShouldBe, decltype(foo<ClassToCheck>(0))>::value>; \
    static bool constexpr value = value_t::value;                                                  \
    constexpr operator bool() const { return value; };                                             \
  }

// clang-format on

/**
  Checks whether TypeToCheck exists in TraitClass<ClassToCheck>.
*/

// clang-format off

#define CREATE_TRAIT_CHECK_CLASS(TraitClass, TypeToCheck, CheckerClassName)                    \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> static void                      foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t = mystd::bool_constant<                                                   \
        !(std::is_void<decltype(foo<TraitClass<ClassToCheck>>(0))>::value)>;                       \
    static bool constexpr value = value_t::value;                                                  \
    constexpr operator bool() const { return value; };                                             \
  }

// clang-format on

/**
  Checks whether TraitClass<ClassToCheck>::TypeToCheck exists and equals TypeShouldBe.
*/

// clang-format off

#define CREATE_TRAIT_EQUALS_CHECK(TraitClass, TypeToCheck, TypeShouldBe, CheckerClassName)         \
  template <class ClassToCheck> class CheckerClassName                                             \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class ...> static void                      foo(...);                                \
                                                                                                   \
  public:                                                                                          \
    using value_t =                                                                                \
        mystd::bool_constant<                                                               \
                  std::is_same<TypeShouldBe, decltype(foo<TraitClass<ClassToCheck>>(0))>::value>;  \
    static bool constexpr value = value_t::value;                                                  \
    constexpr operator bool() const { return value; };                                             \
  }

// clang-format on

/**
  This is used to obtain traits from a trait class, with default settings.
  Notably CheckerClassName<T> is equal to
    TraitClass<T>::TypeToCheck if this exists,
    DefaultType otherwise.
*/

#define MAKE_TRAIT_GETTER(TraitClass, TypeToCheck, DefaultType, CheckerClassName)                  \
namespace TraitGetterHelper{                                                                       \
  template <class ClassToCheck> class CheckerClassName##_Helper                                    \
  {                                                                                                \
  private:                                                                                         \
    template <class Arg> static typename Arg::TypeToCheck foo(int);                                \
    template <class...> static DefaultType foo(...);                                               \
                                                                                                   \
  public:                                                                                          \
    using type = decltype(foo<TraitClass<ClassToCheck>>(0));                                       \
  };                                                                                               \
} \
template<class ClassToCheck> \
using CheckerClassName = typename TraitGetterHelper::CheckerClassName##_Helper<ClassToCheck>::type

#endif
