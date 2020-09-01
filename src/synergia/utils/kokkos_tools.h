#ifndef KOKKOS_TOOLS_H
#define KOKKOS_TOOLS_H

#include "synergia/utils/multi_array_typedefs.h"

namespace kt
{
    template< class ScalarType, int N >
    struct array_type 
    {
        ScalarType arr[N];

        KOKKOS_INLINE_FUNCTION
        array_type() { init(); }

        KOKKOS_INLINE_FUNCTION
        array_type(const array_type & rhs) 
        { for (int i=0; i<N; i++) arr[i] = rhs.arr[i]; }

        KOKKOS_INLINE_FUNCTION  // initialize arr[] to 0
        void init() 
        { for (int i=0; i<N; i++) arr[i] = 0; }

        KOKKOS_INLINE_FUNCTION
        array_type& operator += (const array_type& src) 
        {
            for (int i=0; i<N; i++) arr[i] += src.arr[i];
            return *this;
        }

        KOKKOS_INLINE_FUNCTION
        void operator += (const volatile array_type& src) volatile 
        { for (int i=0; i<N; i++) arr[i] += src.arr[i]; }
    };

    template<class T, int N>
    struct SumArray 
    {
    public:
        //Required
        typedef SumArray reducer;
        typedef array_type<T,N> value_type;
        typedef Kokkos::View<value_type*, Kokkos::MemoryUnmanaged> result_view_type;

    private:
        value_type & value;

    public:

        KOKKOS_INLINE_FUNCTION
        SumArray(value_type& value_): value(value_) {}

        //Required
        KOKKOS_INLINE_FUNCTION
        void join(value_type& dest, const value_type& src) const 
        { dest += src; }

        KOKKOS_INLINE_FUNCTION
        void join(volatile value_type& dest, const volatile value_type& src) const 
        { dest += src; }

        KOKKOS_INLINE_FUNCTION
        void init( value_type& val) const 
        { val.init(); }

        KOKKOS_INLINE_FUNCTION
        value_type& reference() const 
        { return value; }

        KOKKOS_INLINE_FUNCTION
        result_view_type view() const 
        { return result_view_type(&value); }

        KOKKOS_INLINE_FUNCTION
        bool references_scalar() const 
        { return true; }
    };

    struct alg_zeroer
    {
        karray1d_dev arr;

        KOKKOS_INLINE_FUNCTION
        void operator() (const int i) const
        { arr(i) = 0.0; }
    };

    inline void 
    zero_karray(karray1d_dev const& arr)
    {
        alg_zeroer alg{arr};
        Kokkos::parallel_for(arr.extent(0), alg);
    }
}

#endif
