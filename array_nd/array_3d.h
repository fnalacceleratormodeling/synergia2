#ifndef HAVE_ARRAY_3D_H
#define HAVE_ARRAY_3D_H true

#include "array_nd.h"
#include "vector_helper.h"

template<class T>
class Array_3d : public Array_nd<T>
{
public:
    Array_3d();
    Array_3d(const int nx, const int ny, const int nz, T *data_ptr=0);
    Array_3d(const int nx, const int ny, const int nz,
        const int stride_x, const int stride_y, const int stride_z,
        T *data_ptr=0);
    Array_3d(const Array_nd<T> &array_nd);

    void reshape(const int nx, const int ny, const int nz, T *data_ptr=0);
    void reshape(const int nx, const int ny, const int nz,
        const int stride_x, const int stride_y, const int stride_z,
        T *data_ptr=0);

    T& at(const int i, const int j, const int k);
    T at(const int i, const int j, const int k) const;
    T& operator()(const int i, const int j, const int k);
    T operator()(const int i, const int j, const int k) const;

    inline int offset(const int i, const int j, const int k) const;
    inline bool bounds_check(const int i, const int j, const int k) const;

    Array_nd<T> slice(Range i_range, Range j_range, Range k_range, bool reduce=true);

    virtual ~Array_3d();
};

#include "array_3d.tcc"

#endif
