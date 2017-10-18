/**
  This file provides generic implementation of lattice point classes
*/

#ifndef LATTICE_POINT_GENERIC_H
#define LATTICE_POINT_GENERIC_H

namespace GaussSieve{

/**********************
  operators
**********************/

/*********************
  comparison operators
**********************/

template<class LatP>
template<class LatP2, class Impl, TEMPL_RESTRICT_IMPL(
  IsALatticePoint<LatP2>::value && HasInternalRep<Impl>::value && HasInternalRep<LatP2>::value)>
inline bool GeneralLatticePoint<LatP>::operator==(LatP2 const &x2) const
{
  IMPL_IS_LATP;
  // This *might* actually not be an error. However, it is extremely likely.
  static_assert(std::is_same<typename GetCooType<LatP>::type,typename GetCooType<LatP2>::type>::value,
  "Different coordinate types. Probably an error.");

  DEBUG_TRACEGENERIC("Generically comparing " << LatP::class_name() "and" << LatP2::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = CREALTHIS->get_internal_rep_size();
  auto const dim2 = x2.get_internal_rep_size();
  assert(dim1 == dim2);
  #endif // DEBUG_SIEVE_LP_MATCHDIM

  CPP17CONSTEXPRIF(std::is_same<LatP,LatP2>::value && IsNorm2Cheap<LatP>::value)
  {
    if(CREALTHIS->get_norm2() != x2.get_norm2() )
    {
      return false;
    }
  }

  auto const dim = CREALTHIS->get_internal_rep_size();
  for(uint_fast16_t i=0;i<dim;++i)
  {
    if (CREALTHIS->get_internal_rep(i) != x2.get_internal_rep(i))
    {
      return false;
    }
  }
  return true;
}


template<class LatP>
inline bool GeneralLatticePoint<LatP>::operator<(LatP const &rhs) const
{
  DEBUG_TRACEGENERIC("Generically comparing < for" << LatP::class_name() )
  return CREALTHIS->get_norm2_exact() < rhs.get_norm2_exact();
}

template<class LatP>
inline bool GeneralLatticePoint<LatP>::operator>( LatP const &rhs) const
{
  DEBUG_TRACEGENERIC("Generically comparing > for" << LatP::class_name() )
  return CREALTHIS->get_norm2_exact() > rhs.get_norm2_exact();
}

template<class LatP>
inline bool GeneralLatticePoint<LatP>::operator<= ( LatP const &rhs ) const
{
  DEBUG_TRACEGENERIC("Generically comparing <= for" << LatP::class_name() )
  return CREALTHIS->get_norm2_exact() <= rhs.get_norm2_exact();
}

template<class LatP>
inline bool GeneralLatticePoint<LatP>::operator>= ( LatP const &rhs ) const
{
  DEBUG_TRACEGENERIC("Generically comparing >= for" << LatP::class_name() )
  return CREALTHIS->get_norm2_exact() >= rhs.get_norm2_exact();
}


/******************************************
  arithmetic operators modifying the point
******************************************/

template<class LatP>
template<class LatP2, class Impl, TEMPL_RESTRICT_IMPL(
  IsALatticePoint<LatP2>::value && IsRepLinear_RW<Impl>::value && IsRepLinear<LatP2>::value)>
inline LatP& GeneralLatticePoint<LatP>::operator+=(LatP2 const &x2)
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC( "generically adding" << LatP::class_name() << " and " << LatP2::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = CREALTHIS->get_internal_rep_size();
  auto const dim2 = x2.get_internal_rep_size();
  assert( dim1 == dim2 );
  auto const real_dim1 = CREALTHIS->get_dim();
  auto const real_dim2 = x2.get_dim();
  assert(real_dim1 == real_dim2);
  #endif
  auto const dim = CREALTHIS->get_internal_rep_size();
//  auto const real_dim = x1.get_dim();
//  LP NewLP(real_dim);
  for(uint_fast16_t i = 0; i < dim; ++i )
  {
    ( REALTHIS->get_internal_rep(i) ) += x2.get_internal_rep(i);
  }
  REALTHIS->sanitize();
  return *REALTHIS;
}

template<class LatP>
template<class LatP2, class Impl, TEMPL_RESTRICT_IMPL(
  IsALatticePoint<LatP2>::value && IsRepLinear_RW<Impl>::value && IsRepLinear<LatP2>::value)>
inline LatP& GeneralLatticePoint<LatP>::operator-=(LatP2 const &x2)
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC( "generically adding" << LatP::class_name() << " and " << LatP2::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = CREALTHIS->get_internal_rep_size();
  auto const dim2 = x2.get_internal_rep_size();
  assert( dim1 == dim2 );
  auto const real_dim1 = CREALTHIS->get_dim();
  auto const real_dim2 = x2.get_dim();
  assert(real_dim1 == real_dim2);
  #endif
  auto const dim = CREALTHIS->get_internal_rep_size();
//  auto const real_dim = x1.get_dim();
//  LP NewLP(real_dim);
  for(uint_fast16_t i = 0; i < dim; ++i )
  {
    REALTHIS->get_internal_rep(i) -= x2.get_internal_rep(i);
  }
  REALTHIS->sanitize();
  return *REALTHIS;
}

template<class LatP>
template<class Integer, class Impl, TEMPL_RESTRICT_IMPL(
  IsRepLinear_RW<Impl>::value && std::is_integral<Integer>::value)>
inline LatP& GeneralLatticePoint<LatP>::operator*=(Integer const multiplier)
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Generically scalar-multiplying for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_internal_rep_size();
  for(uint_fast16_t i=0;i<dim;++i)
  {
    REALTHIS->get_internal_rep(i) *= multiplier;
  }
  CPP17CONSTEXPRIF(IsNorm2Cheap<LatP>::value)
  {
    REALTHIS->sanitize(CREALTHIS->get_norm2_exact() * multiplier * multiplier );
  }
  else
  {
    REALTHIS->sanitize();
  }
  return *REALTHIS;
}

//inline LatP& operator*=(mpz_class const &multiplier);

/****************************************
  out-of-class definitions of +,-,* etc.
******************************************/

// This is actually tricky:
// Be aware that LP1 && LP && below are universal references, NOT rvalue references!

// The whole code relies on the fact that IsALatticePoint<LP> will only work
// for LP a lattice point, but not a lvalue reference type to it.
// (Note that the FOR_LATTICE_POINT* macros use IsALatticePoint)
// Consequently, the unversal references behave (almost) like rvalue references.

// binary +

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator+(LP1 const &x1, LP2 const &x2)
{
  LP1 NewLP(x1.make_copy());
  NewLP+=x2;
  return NewLP;
}

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator+(LP1 && x1, LP2 const &x2)
{
  auto tmp = std::move(x1);
  tmp+=x2;
  return tmp;
}

// We don't want to return LP2 here...
// If this causes trouble, the caller should change the order of arguments.
FOR_LATTICE_POINT_LP
LP operator+(LP const &x1, LP && x2)
{
  auto tmp = std::move(x2);
  tmp+=x1;
  return tmp;
}

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator+(LP1 &&x1, LP2 && x2)
{
  auto tmp = std::move(x1);
  tmp+=std::move(x2);
  return tmp;
//LP1 tmp = std::move(x1);
//return addval(tmp,std::move(x2));
}

// binary minus

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 const &x1, LP2 const &x2)
{
  LP1 NewLP(x1.make_copy());
  NewLP-=x2;
  return NewLP;
}

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 && x1, LP2 const &x2)
{
  LP1 tmp = std::move(x1);
  tmp-=x2;
  return tmp;
}

FOR_LATTICE_POINT_LP
LP operator-(LP const &x1, LP &&x2)
{
  static_assert(IsNorm2Cheap<LP>::value,"Improve this code");
  LP tmp = std::move(x2);
  tmp.make_negative();
  tmp+=x1;
  return tmp;
}

FOR_LATTICE_POINTS_LP1_LP2
LP1 operator-(LP1 && x1, LP2 && x2)
{
  LP1 tmp = std::move(x1);
  tmp-=std::move(x2);
  return tmp;
}



// unary minus

template<class LatP>
inline LatP GeneralLatticePoint<LatP>::operator-() &&
{
  LatP tmp = std::move(*REALTHIS);
  tmp.make_negative();
  return tmp;
}

// Note Integer is passed by value, even for mpz_classes. Not optimal...
template<class LP, class Integer, TEMPL_RESTRICT_DECL(
  IsALatticePoint<LP>::value &&
  (std::is_integral<Integer>::value || std::is_same<Integer,mpz_class>::value) )>
inline LP operator*(LP const &x1, Integer const multiplier)
{
//  assert(false);
  LP tmp = x1.make_copy();
  tmp*=multiplier;
  return tmp;
}

template<class LP, class Integer, TEMPL_RESTRICT_DECL(
  IsALatticePoint<LP>::value &&
  (std::is_integral<Integer>::value || std::is_same<Integer,mpz_class>::value) )>
inline LP operator*(LP &&x1, Integer const multiplier)
{
  LP tmp = std::move(x1);
  tmp*=multiplier;
  return tmp;
}


/*************************
I/O
*************************/

template<class LatP>
inline std::ostream& GeneralLatticePoint<LatP>::write_lp_to_stream(std::ostream &os, bool const include_norm2) const
{
  DEBUG_TRACEGENERIC("Using generic writer (absolute) for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_dim();
  os << "[ "; // makes spaces symmetric
  for (uint_fast16_t i =0; i<dim; ++i)
  {
    os << CREALTHIS->get_absolute_coo(i) << " ";
  }
  os << "]";
  if(include_norm2)
  {
    os <<", norm2= " << CREALTHIS->get_norm2();
  }
      // No endl here (this is the caller's job).
  return os;
}

template<class LatP>
template<class Impl, TEMPL_RESTRICT_IMPL(HasInternalRep<Impl>::value)>
inline std::ostream& GeneralLatticePoint<LatP>::write_lp_rep_to_stream(std::ostream &os) const
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using generic writer (internal rep) for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_internal_rep_size();
  os << "[ ";
  for (uint_fast16_t i=0; i<dim; ++i)
  {
    os << CREALTHIS->get_internal_rep(i) << " ";
  }
  os << "]";
  return os;
}


FOR_LATTICE_POINT_LP
std::istream& operator>>(std::istream &is, LP &LatP)
{
  return LatP.read_from_stream(is);
}

FOR_LATTICE_POINT_LP
std::ostream& operator<<(std::ostream &os, LP const &LatP)
{
  return LatP.write_lp_to_stream(os);
}


/********************************
Getter functions
********************************/

template<class LatP>
template<class Impl, TEMPL_RESTRICT_IMPL(HasInternalRep<Impl>::value)>
inline auto GeneralLatticePoint<LatP>::get_internal_rep_size() const -> decltype( std::declval<Impl>().get_dim() )
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Generically getting vec_size for" << LatP::class_name() )
  return CREALTHIS->get_dim();
}




template<class LatP>
inline typename GetScalarProductStorageType<LatP>::type GeneralLatticePoint<LatP>::get_norm2() const
{
  DEBUG_TRACEGENERIC("Generically computing norm2 for " << LatP::class_name() )
      // This function should not be called if IsNorm2Cheap is set,
      // since IsNorm2Cheap basically promises an overload.
  static_assert(IsNorm2Cheap<LatP>::value == false, "");
  return compute_sc_product(*(CREALTHIS),*(CREALTHIS) );
}


template<class LatP>
template<class Impl, TEMPL_RESTRICT_IMPL(IsRepLinear<Impl>::value)>
inline bool GeneralLatticePoint<LatP>::is_zero() const
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using (possibly inefficient) test for zero for " << LatP::class_name() )

  if (IsNorm2Cheap<LatP>::value) // constexpr if
  {
    return (CREALTHIS->get_norm2() == 0);
  }
  else
  {
    auto const dim = CREALTHIS->get_internal_rep_size();
    for (uint_fast16_t i=0;i<dim;++i)
    {
      if(CREALTHIS->get_internal_rep(i) != 0)
      {
        return false;
      }
    }
    return true;
  }
}

/**************************
  Initialization functions
***************************/

template<class LatP>
template<class Impl, TEMPL_RESTRICT_IMPL(IsRepLinear_RW<Impl>::value)>
inline void GeneralLatticePoint<LatP>::fill_with_zero()
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using generic fill with zero for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_internal_rep_size();
  for (uint_fast16_t i=0;i<dim;++i)
  {
    REALTHIS->get_internal_rep(i) = 0;
  }
  REALTHIS->sanitize(0);
}

template<class LatP>
template<class Impl, TEMPL_RESTRICT_IMPL(IsRepLinear_RW<Impl>::value)>
inline void GeneralLatticePoint<LatP>::make_negative()
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using generic negation function for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_internal_rep_size();
  for (uint_fast16_t i=0; i<dim; ++i)
  {
    REALTHIS->get_internal_rep(i) = - CREALTHIS->get_internal_rep(i);
  }
  CPP17CONSTEXPRIF( IsNegateCheap<Impl>::value) // constexpr if
  {
    return;
  }
  else
  {
    if(IsNorm2Cheap<Impl>::value)
    {
      REALTHIS->sanitize(CREALTHIS->get_norm2_full() );
    }
    else
    {
    REALTHIS->sanitize();
    }
  }
  return;
}



template<class LatP>
template<class Impl, TEMPL_RESTRICT_IMPL(IsRepLinear_RW<Impl>::value)>
inline LatP GeneralLatticePoint<LatP>::make_copy() const
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using generic copy for " << LatP::class_name() )
  auto const real_dim=CREALTHIS->get_dim(); // means ambient dimension.
  auto const dim = CREALTHIS->get_internal_rep_size(); // number of coordinates stored. May be rank.
  LatP NewLP(real_dim);
  for (uint_fast16_t i=0; i<dim; ++i)
  {
    NewLP.get_internal_rep(i) = CREALTHIS->get_internal_rep(i);
  }
  if (IsNorm2Cheap<LatP>::value)
  {
    NewLP.sanitize(CREALTHIS->get_norm2_full() );
  }
  else
  {
    NewLP.sanitize();
  }
  return NewLP;
}

/*********************
  Scalar product
*********************/



template<class LatP>
inline typename GeneralLatticePoint<LatP>::ScalarProductStorageType GeneralLatticePoint<LatP>::do_compute_sc_product(LatP const &x2) const
{
  DEBUG_TRACEGENERIC("Generically computing scalar product for" << LP::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = CREALTHIS->get_dim();
  auto const dim2 = x2.get_dim();
  assert(dim1 == dim2 );
  #endif // DEBUG_SIEVE_LP_MATCHDIM
  using ET = typename GetAbsoluteCooType<LatP>::type;
  auto const dim = CREALTHIS->get_dim();
  ET result = 0; // assumes that ET can be initialized from 0...
  for(uint_fast16_t i=0; i<dim; ++i)
  {
    result += CREALTHIS->get_absolute_coo(i) * x2.get_absolute_coo(i);
  }
  return static_cast<typename LatP::ScalarProductStorageType>(result);
}






} // end namespace

#endif
