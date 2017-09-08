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
inline bool GeneralLatticePoint<LatP>::operator<(LatP const &rhs) const
{
  DEBUG_TRACEGENERIC("Generically comparing < for" << LatP::class_name() )
  return CREALTHIS->get_norm2() < rhs.get_norm2();
}

template<class LatP>
inline bool GeneralLatticePoint<LatP>::operator>( LatP const &rhs) const
{
  DEBUG_TRACEGENERIC("Generically comparing > for" << LatP::class_name() )
  return CREALTHIS->get_norm2() > rhs.get_norm2();
}

template<class LatP>
inline bool GeneralLatticePoint<LatP>::operator<= ( LatP const &rhs ) const
{
  DEBUG_TRACEGENERIC("Generically comparing <= for" << LatP::class_name() )
  return CREALTHIS->get_norm2() <= rhs.get_norm2();
}

template<class LatP>
inline bool GeneralLatticePoint<LatP>::operator>= ( LatP const &rhs ) const
{
  DEBUG_TRACEGENERIC("Generically comparing >= for" << LatP::class_name() )
  return CREALTHIS->get_norm2() >= rhs.get_norm2();
}

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
inline typename GeneralLatticePoint<LatP>::ScalarProductReturnType GeneralLatticePoint<LatP>::get_norm2() const
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
  for (uint_fast16_t i=0;i<dim;++i)
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

template<class LatP>
template<class Integer,
  TEMPL_RESTRICT_IMPL(IsCooVector<LatP>::value && std::is_integral<Integer>::value)>
inline void GeneralLatticePoint<LatP>::scalar_multiply(Integer const multiplier)
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
}







} // end namespace

#endif
