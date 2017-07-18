#ifndef COMPRESSED_LATTICE_POINT_H
#define COMPRESSED_LATTICE_POINT_H

#include "ExactLatticePoint.h"

template<class ET, bool MT, int nfixed>
class CompressedPoint
{
    public:
    using DetailType = ExactLatticePoint<ET,nfixed>;
    CompressedPoint() : approx_data(), details(nullptr) {};                                                                 //creates a dummy, invalid point
    CompressedPoint(CompressedPoint && old) : approx_data(std::move(old.approx_data)), details(old.details) {old.details = nullptr;}; //move semantics
    CompressedPoint(CompressedPoint const & old) = delete; //no implicit copying, we have several explicit functions for that.
    CompressedPoint<ET,MT,nfixed>& operator=(CompressedPoint<ET,MT,nfixed> const &that) =delete;
    CompressedPoint<ET,MT,nfixed>& operator=(CompressedPoint<ET,MT,nfixed> && other){approx_data = std::move(other.approx_data);swap(details,other.details);return *this;};
    //TODO: This function should be called sparingly, as it makes an unneccessary copy.
    //[[deprecated("Avoid copying")]]
    explicit CompressedPoint(ExactLatticePoint<ET,nfixed> const & exact_point);     //creates a CP from an exact point. Makes a copy!
    explicit CompressedPoint(ExactLatticePoint<ET,nfixed> * const exact_point_ptr); //transfers ownership
    ~CompressedPoint();

    ApproximateLatticePoint<ET,nfixed>  const &                         access_approximation_r() const          {return approx_data;};
    ApproximateLatticePoint<ET,nfixed>  const*                          get_approximation_const_ptr() const     {return &approx_data;};
    ApproximateLatticePoint<ET,nfixed> *                                get_approximation_ptr()                 {return &approx_data;};
    ExactLatticePoint<ET,nfixed> const &                                access_exact_point_r() const            {return *details;};
    typename ApproximateLatticePoint<ET,nfixed>::ApproxTypeNorm2        get_approx_norm2_mantissa() const       {return approx_data.get_norm2_mantissa();};
    signed int                                                          get_approx_norm2_exponent() const       {return approx_data.get_norm2_exponent();};
    double                                                              get_approx_norm2_d() const              {return approx_data.get_norm2_d();};
    signed int                                                          get_approx_vectors_exponent() const     {return approx_data.get_vectors_exponent();};
    typename ApproximateLatticePoint<ET,nfixed>::ApproxEntryType*       access_approx_vectors_mantissa() const  {return approx_data.access_vectors_mantissa();};
    ET                                                                  get_exact_norm2() const                 {return details->access_norm2();};
    ExactLatticePoint<ET,nfixed>                                        get_exact_point() const;                 //returns a copy of the underlying exact point

    bool                                                                is_zero() const {return (details->access_norm2()==0);};
    //TODO: approx. and exact scalar products?

    bool operator< (CompressedPoint const &other ) const {return (this->access_exact_point_r() < other.access_exact_point_r() );};
    CompressedPoint<ET,MT,nfixed> deep_copy_compressed_point() const; //returns a deep copy

    protected:
    ApproximateLatticePoint<ET,nfixed> approx_data;
    DetailType * details;
};

template<class ET, bool MT, int nfixed> CompressedPoint<ET,MT,nfixed> CompressedPoint<ET,MT,nfixed>::deep_copy_compressed_point() const //returns a deep copy
{
    return static_cast<CompressedPoint<ET,MT,nfixed> >(access_exact_point_r() );
}


template<class ET, bool MT, int nfixed> CompressedPoint<ET,MT,nfixed>::CompressedPoint(ExactLatticePoint<ET,nfixed> const  & exact_point)     //creates a CP from an exact point. Makes a copy!
:approx_data(exact_point), details(nullptr)
{
    details = new DetailType (exact_point);
}


template<class ET, bool MT, int nfixed> CompressedPoint<ET,MT,nfixed>::CompressedPoint(ExactLatticePoint<ET,nfixed> * const exact_point_ptr)
:approx_data(*exact_point_ptr), details(exact_point_ptr)
{
}

template<class ET, bool MT, int nfixed> typename CompressedPoint<ET,MT,nfixed>::DetailType CompressedPoint<ET,MT,nfixed>::get_exact_point() const
{
    return * details;
}

template<class ET, bool MT, int nfixed> CompressedPoint<ET,MT,nfixed>::~CompressedPoint()
{
    delete details;
}

#endif
