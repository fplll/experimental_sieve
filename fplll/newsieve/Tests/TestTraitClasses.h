#ifndef TEST_TRAIT_CLASSES_H
#define TEST_TRAIT_CLASSES_H

#include <type_traits>
#include "../SieveUtility.h"
#include "../Typedefs.h"
#include "gmpxx.h"

class C1{};
class C2{};
class C3{};
class C4{};

template<class T> class TraitClass
{
  public:
  using Invalid = std::true_type;
};

template<>
class TraitClass<C1>
{
  public:
  using TraitA = std::true_type;
  using TraitB = int;
};

template<>
class TraitClass<C2>
{
  public:
  using TraitB = double;
};

template<>
class TraitClass<C3>
{
  // just defaults.
};


CREATE_MEMBER_TYPEDEF_CHECK_CLASS(Invalid, HasInvalidTrait);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TraitA, HasTraitA);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS(TraitB, HasTraitB);
CREATE_MEMBER_TYPEDEF_CHECK_CLASS_EQUALS(TraitB, int, BIsInt);
CREATE_TRAIT_CHECK_CLASS(TraitClass, Invalid, IsInvalid);
CREATE_TRAIT_CHECK_CLASS(TraitClass, TraitB, HasB);
CREATE_TRAIT_EQUALS_CHECK(TraitClass, TraitB, int, TraitBIsInt);
MAKE_TRAIT_GETTER(TraitClass, TraitB, bool, GetB);


bool test_trait_classes()
{
  static_assert(HasInvalidTrait<TraitClass<C4>>::value,"");
  static_assert(HasInvalidTrait<TraitClass<C4>>(),"");
  static_assert(HasInvalidTrait<TraitClass<C4>>{},"");
  static_assert(!HasInvalidTrait<TraitClass<C1>>::value,"");
  static_assert(IsInvalid<C4>::value,"");
  static_assert(!IsInvalid<C3>::value,"");
  static_assert(HasTraitA<TraitClass<C1>>::value,"");
  static_assert(!HasTraitA<TraitClass<C2>>::value,"");
  static_assert(BIsInt<TraitClass<C1>>::value,"");
  static_assert(!BIsInt<TraitClass<C2>>::value,"");
  static_assert(!BIsInt<TraitClass<C3>>::value,"");
  static_assert(IsInvalid<C4>{},"");
  static_assert(!IsInvalid<C1>{},"");
  static_assert(!IsInvalid<C3>{},"");
  static_assert(!HasB<C3>{},"");
  static_assert(!HasB<C4>{},"");
  static_assert(HasB<C1>{},"");
  static_assert(HasB<C2>{},"");
  static_assert(TraitBIsInt<C1>{},"");
  static_assert(!TraitBIsInt<C2>{},"");
  static_assert(!TraitBIsInt<C3>{},"");
  static_assert(std::is_same<int, GetB<C1>::type>::value,"");
  static_assert(std::is_same<double, GetB<C2>::type>::value,"");
  static_assert(std::is_same<bool, GetB<C3>::type>::value,"");

  using GaussSieve::IsZNRClass;
  using GaussSieve::UnZNR;
  using GaussSieve::AddZNR;
  using GaussSieve::FixZNR;
  using fplll::Z_NR;

  static_assert(!IsZNRClass<long>::value,"");
  static_assert(IsZNRClass<Z_NR<long>>{},"");
  static_assert(IsZNRClass<Z_NR<mpz_t>>(),"");
  static_assert(!IsZNRClass<mpz_class>::value);

  static_assert(std::is_same<long,UnZNR<Z_NR<long>>::type>::value,"" );
  static_assert(std::is_same<double,UnZNR<Z_NR<double>>::type>::value,"" );
  static_assert(std::is_same<mpz_class,UnZNR<Z_NR<mpz_t>>::type>::value,"" );
  static_assert(std::is_same<bool,UnZNR<bool>::type>::value,"" );
  static_assert(std::is_same<AddZNR<long>::type, Z_NR<long> >::value,"");
  static_assert(std::is_same<AddZNR<mpz_t>::type, Z_NR<mpz_t> >::value,"");
  static_assert(std::is_same<AddZNR<mpz_class>::type, Z_NR<mpz_t> >::value,"");

  static_assert(std::is_same<FixZNR<bool>::type,bool>::value,"" );
  static_assert(std::is_same<FixZNR<long>::type,long>::value,"" );
  static_assert(std::is_same<FixZNR<Z_NR<mpz_t>>::type,Z_NR<mpz_t>>::value,"" );
  static_assert(std::is_same<FixZNR<mpz_class>::type,mpz_t>::value,"" );
return true;
}

#endif
