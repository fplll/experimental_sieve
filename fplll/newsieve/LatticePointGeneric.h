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
template<class LatP2, TEMPL_RESTRICT_IMPL(IsALatticePoint<LatP2>::value && HasCoos<LatP>::value && HasCoos<LatP2>::value)>
inline bool GeneralLatticePoint<LatP>::operator==(LatP2 const &x2) const
{
  // This *might* actually not be an error. However, it is extremely likely.
  static_assert(std::is_same<typename GetCooType<LatP>::type,typename GetCooType<LatP2>::type>::value,"Different coordinate types. Probably an error.");

  DEBUG_TRACEGENERIC("Generically comparing " << LatP::class_name() "and" << LatP2::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = CREALTHIS->get_vec_size();
  auto const dim2 = x2.get_vec_size();
  assert(dim1 == dim2);
  #endif // DEBUG_SIEVE_LP_MATCHDIM

  if(std::is_same<LatP,LatP2>::value && IsNorm2Cheap<LatP>::value) // constexpr if
  {
    if(CREALTHIS->get_norm2() != x2.get_norm2() )
    {
      return false;
    }
  }

  auto const dim = CREALTHIS->get_vec_size();
  for(uint_fast16_t i=0;i<dim;++i)
  {
    if (CREALTHIS->operator[](i) != x2[i])
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
template<class LatP2, TEMPL_RESTRICT_IMPL(IsALatticePoint<LatP2>::value && IsCooVector<LatP>::value && HasCoos<LatP2>::value)>
inline LatP& GeneralLatticePoint<LatP>::operator+=(LatP2 const &x2)
{
  DEBUG_TRACEGENERIC( "generically adding" << LatP::class_name() << " and " << LatP2::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = CREALTHIS->get_vec_size();
  auto const dim2 = x2.get_vec_size();
  assert( dim1 == dim2 );
  auto const real_dim1 = CREALTHIS->get_dim();
  auto const real_dim2 = x2.get_dim();
  assert(real_dim1 == real_dim2);
  #endif
  auto const dim = CREALTHIS->get_vec_size();
//  auto const real_dim = x1.get_dim();
//  LP NewLP(real_dim);
  for(uint_fast16_t i = 0; i < dim; ++i )
  {
    ( REALTHIS->operator[](i) ) += x2[i];
  }
  REALTHIS->sanitize();
  return *REALTHIS;
}

template<class LatP>
template<class LatP2, TEMPL_RESTRICT_IMPL(IsALatticePoint<LatP2>::value && IsCooVector<LatP>::value && HasCoos<LatP2>::value)>
inline LatP& GeneralLatticePoint<LatP>::operator-=(LatP2 const &x2)
{
  DEBUG_TRACEGENERIC( "generically adding" << LatP::class_name() << " and " << LatP2::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = CREALTHIS->get_vec_size();
  auto const dim2 = x2.get_vec_size();
  assert( dim1 == dim2 );
  auto const real_dim1 = CREALTHIS->get_dim();
  auto const real_dim2 = x2.get_dim();
  assert(real_dim1 == real_dim2);
  #endif
  auto const dim = CREALTHIS->get_vec_size();
//  auto const real_dim = x1.get_dim();
//  LP NewLP(real_dim);
  for(uint_fast16_t i = 0; i < dim; ++i )
  {
    REALTHIS->operator[](i) -= x2[i];
  }
  REALTHIS->sanitize();
  return *REALTHIS;
}

template<class LatP>
template<class Integer,
  TEMPL_RESTRICT_IMPL(IsCooVector<LatP>::value && std::is_integral<Integer>::value)>
inline LatP& GeneralLatticePoint<LatP>::operator*=(Integer const multiplier)
{
  DEBUG_TRACEGENERIC("Generically scalar-multiplying for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_vec_size();
  for(uint_fast16_t i=0;i<dim;++i)
  {
    REALTHIS->operator[](i) *= multiplier;
  }
  if(IsNorm2Cheap<LatP>::value) // constexpr if
  {
    REALTHIS->sanitize(CREALTHIS->get_norm2() * multiplier * multiplier );
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
//LP1 tmp = std::move(x1);
//return addval(tmp,x2);
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


// unary minus

template<class LatP>
inline LatP GeneralLatticePoint<LatP>::operator-() &&
{
  LatP tmp = std::move(*REALTHIS);
  tmp.make_negative();
  return tmp;
}

template<class LP, class Integer,
  typename std::enable_if<IsALatticePoint<LP>::value &&
  (std::is_integral<Integer>::value || std::is_same<Integer,mpz_class>::value),
  int>::type=0>
LP operator*(LP const &x1, Integer const multiplier)
{
//  assert(false);
  LP tmp = x1.make_copy();
  tmp*=multiplier;
  return tmp;
}

template<class LP, class Integer,
  typename std::enable_if<IsALatticePoint<LP>::value &&
  (std::is_integral<Integer>::value || std::is_same<Integer,mpz_class>::value),
  int>::type=0>
LP operator*(LP &&x1, Integer const multiplier)
{
  LP tmp = std::move(x1);
  tmp*=multiplier;
  return tmp;
}




/*************************
I/O
*************************/

template<class LatP>
MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE_IMPL // This may be too strict.
inline std::ostream& GeneralLatticePoint<LatP>::write_to_stream(std::ostream &os, bool const include_norm2) const
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using generic writer for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_vec_size();
  os << "[ "; // makes spaces symmetric
  for (uint_fast16_t i =0; i<dim; ++i)
  {
    os << CREALTHIS->operator[](i) << " ";
  }
    os <<"]";
  if(include_norm2)
  {
    os <<", norm2= " << CREALTHIS->get_norm2();
  }
      // No endl here (this is the caller's job).
  return os;
}

/********************************
Getter functions
********************************/

template<class LatP>
MEMBER_ONLY_EXISTS_IF_COO_READ_IMPL
inline auto GeneralLatticePoint<LatP>::get_vec_size() const -> decltype( std::declval<Impl>().get_dim() )
{
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


/********

******/

template<class LatP>
MEMBER_ONLY_EXISTS_IF_COO_WRITE_IMPL
inline void GeneralLatticePoint<LatP>::fill_with_zero()
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using generic fill with zero for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_vec_size();
  for (uint_fast16_t i=0;i<dim;++i)
  {
    REALTHIS->operator[](i) = 0;
  }
  REALTHIS->sanitize(0);
}

template<class LatP>
MEMBER_ONLY_EXISTS_IF_COO_WRITE_IMPL
inline void GeneralLatticePoint<LatP>::make_negative()
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using generic negation function for " << LatP::class_name() )
  auto const dim = CREALTHIS->get_vec_size();
  for (uint_fast16_t i=0; i<dim; ++i)
  {
    REALTHIS->operator[](i) = - REALTHIS->operator[](i);
  }
  if( IsNegateCheap<Impl>::value) // constexpr if
  {
    return;
  }
  else
  {
    if(IsNorm2Cheap<Impl>::value)
    {
      REALTHIS->sanitize(CREALTHIS->get_norm2() );
    }
    else
    {
    REALTHIS->sanitize();
    }
  }
  return;
}

template<class LatP>
MEMBER_ONLY_EXISTS_IF_COO_READ_IMPL
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
    auto const dim = CREALTHIS->get_vec_size();
    for (uint_fast16_t i=0;i<dim;++i)
    {
      if(CREALTHIS->operator[](i) != 0)
      {
        return false;
      }
    }
    return true;
  }
}

template<class LatP>
MEMBER_ONLY_EXISTS_IF_COO_WRITE_IMPL
inline LatP GeneralLatticePoint<LatP>::make_copy() const
{
  IMPL_IS_LATP;
  DEBUG_TRACEGENERIC("Using generic copy for " << LatP::class_name() )
  auto const real_dim=CREALTHIS->get_dim(); // means ambient dimension.
  auto const dim = CREALTHIS->get_vec_size(); // number of coordinates stored. May be rank.
  LatP NewLP(real_dim);
  for (uint_fast16_t i=0; i<dim; ++i)
  {
    NewLP[i] = CREALTHIS->operator[](i);
  }
  if (IsNorm2Cheap<LatP>::value)
  {
    NewLP.sanitize(CREALTHIS->get_norm2() );
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

template<class LP, TEMPL_RESTRICT_IMPL(IsALatticePoint<LP>::value)>
typename LP::ScalarProductStorageType compute_sc_product(LP const &lp1, LP const &lp2)
{
  return lp1.do_compute_sc_product(lp2);
}


template<class LP, TEMPL_RESTRICT_IMPL(IsALatticePoint<LP>::value)>
typename LP::ScalarProductStorageType compute_sc_product_exact(LP const &lp1, LP const &lp2)
{
  return lp1.do_compute_sc_product_exact(lp2);
}

template<class LP, TEMPL_RESTRICT_IMPL(IsALatticePoint<LP>::value)>
typename LP::ScalarProductStorageType compute_sc_product_full(LP const &lp1, LP const &lp2)
{
  return lp2.do_compute_sc_product_full(lp2);
}

template<class LatP>
MEMBER_ONLY_EXISTS_IF_COOS_ABSOLUTE_IMPL
inline typename GeneralLatticePoint<LatP>::ScalarProductStorageType GeneralLatticePoint<LatP>::do_compute_sc_product(LatP const &x2) const
{
  DEBUG_TRACEGENERIC("Generically computing scalar product for" << LP::class_name() )
  #ifdef DEBUG_SIEVE_LP_MATCHDIM
  auto const dim1 = CREALTHIS->get_dim();
  auto const dim2 = x2.get_dim();
  assert(dim1 == dim2 );
  #endif // DEBUG_SIEVE_LP_MATCHDIM
  using ET = typename GetCooType<LatP>::type;
  auto const dim = CREALTHIS->get_dim();
  ET result = 0; // assumes that ET can be initialized from 0...
  for(uint_fast16_t i=0; i<dim; ++i)
  {
    result += (*CREALTHIS)[i] * x2[i];
  }
  return static_cast<typename LatP::ScalarProductStorageType>(result);
}






} // end namespace

#endif
