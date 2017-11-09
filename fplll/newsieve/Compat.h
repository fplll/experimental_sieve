// Definitions for C++ version compatibility

#ifndef SIEVE_GAUSS_COMPAT_H
#define SIEVE_GAUSS_COMPAT_H

#include<type_traits>
#include<utility>

#if __cpp_constexpr >= 201304
  #define CPP14CONSTEXPR constexpr
#else
  #define CPP14CONSTEXPR
#endif

#if __if_constexpr
  #define CPP17CONSTEXPRIF if constexpr
#else
  #define CPP17CONSTEXPRIF if
#endif

namespace GaussSieve
{
#if __cpp_lib_logical_traits >= 201510
// untested:

  template<class... Bs> using MyConjunction = std::conjunction<Bs...>;    //AND
  template<class... Bs> using MyDisjunction = std::disjunction<Bs...>;    //OR
  template<class B>     using MyNegation    = std::negation<B>;           //NOT
  template<class... Bs> using MyNAND = MyNegation<MyConjunction<Bs...>>;  //NAND
  template<class... Bs> using MyNOR  = MyNegation<MyDisjunction<Bs...>>;  //NOR
#else
// just implement std::conjunction and friends myself:
  template<class...> struct MyConjunction     : std::true_type{};
  template<class B1> struct MyConjunction<B1> : B1 {};
  template<class B1,class... Bs> struct MyConjunction<B1,Bs...>
    : std::conditional<static_cast<bool>(B1::value), MyConjunction<Bs...>,B1>::type {};

  template<class...> struct MyDisjunction     : std::false_type{};
  template<class B1> struct MyDisjunction<B1> : B1 {};
  template<class B1, class... Bs> struct MyDisjunction<B1,Bs...>
    : std::conditional<static_cast<bool>(B1::value), B1, MyDisjunction<Bs...> >::type {};

  template<class B> struct MyNegation : std::integral_constant<bool,!static_cast<bool>(B::value)>{};
  template<class... Bs> using MyNAND = MyNegation<MyConjunction<Bs...>>;
  template<class... Bs> using MyNOR  = MyNegation<MyDisjunction<Bs...>>;

#endif

#if __cpp_lib_integer_sequence >= 201304
  template<std::size_t... Ints> using MyIndexSeq = std::index_sequence<Ints...>;
  template<std::size_t N>  using MyMakeIndexSeq = std::make_index_sequence<N>;
  template<class... T> using MyIndexSequenceFor = std::index_sequence_for<T...>;
#else
  template<std::size_t... Ints> class MyIndexSeq {}; // not equivalent to the above, lacks size()
  namespace Helpers
  {
  // encapsulates a integer sequence 0,...,N-1, RestArgs
    template<std::size_t N, std::size_t... RestArgs> struct GenIndexSeq
      :GenIndexSeq <N-1,N-1, RestArgs...> {};
    template<std::size_t... RestArgs> struct GenIndexSeq<0,RestArgs...>{using type = MyIndexSeq<RestArgs...>;};
  }
  template<std::size_t N> using MyMakeIndexSeq = typename Helpers::GenIndexSeq<N>::type;
  template<class... T> using MyIndexSequenceFor = MyMakeIndexSeq<sizeof...(T)>;
#endif

namespace mystd
{
  template<class T> using decay_t = typename std::decay<T>::type;
  template<bool b>  using bool_constant = std::integral_constant<bool, b>;
  template<class... Bs> using conjunction = MyConjunction<Bs...>;    //AND
  template<class... Bs> using disjunction = MyDisjunction<Bs...>;    //OR
  template<class B>     using negation    = MyNegation<B>;           //NOT
}
} // end namespace

#endif
