#ifndef HAVE_ARRAY_3D_H
#define HAVE_ARRAY_3D_H true

#include "array_nd.h"
#include "vector_helper.h"

template<class T>
class Array_3d : public Array_nd<T>
{
public:
    Array_3d();
    Array_3d(const int nx, const int ny, const int nz);
    Array_3d(const int nx, const int ny, const int nz, T *data_ptr);

    void reshape(const int nx, const int ny, const int nz);
    void reshape(const int nx, const int ny, const int nz, T *data_ptr);

    T& at(const int i, const int j, const int k);
    T at(const int i, const int j, const int k) const;
    T& operator()(const int i, const int j, const int k);
    T operator()(const int i, const int j, const int k) const;

    inline int offset(const int i, const int j, const int k) const;
    inline bool bounds_check(const int i, const int j, const int k) const;

    virtual ~Array_3d();
};

// Begin implementation.
#include <iostream>

template<class T>
Array_3d<T>::Array_3d() : Array_nd<T>()
{
    this->shape = vector3(0,0,0);
}

template<class T>
Array_3d<T>::Array_3d(const int nx, const int ny, const int nz) :
    Array_nd<T>()
{
    this->construct(vector3(nx,ny,nz),true);
}

template<class T>
Array_3d<T>::Array_3d(const int nx, const int ny, const int nz, T *data_ptr) :
    Array_nd<T>()
{
    this->data_ptr = data_ptr;
    this->construct(vector3(nx,ny,nz),false);
}

template<class T>
void
Array_3d<T>::reshape(const int nx, const int ny, const int nz)
{
    bool shape_changed = this->different_shape(vector3(nx,ny,nz));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_3d");
    }
    if (shape_changed) {
        this->construct(vector3(nx,ny,nz),true);
    }
}

template<class T>
void
Array_3d<T>::reshape(const int nx, const int ny, const int nz, T *data_ptr)
{
    bool shape_changed = this->different_shape(vector3(nx,ny,nz));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_3d");
    }
    this->data_ptr = data_ptr;
    if (shape_changed) {
        this->construct(vector3(nx,ny,nz),false);
    }
}

template<class T>
inline int
Array_3d<T>::offset(const int i, const int j, const int k) const
{
    return k + this->shape[1]*(j+this->shape[2]*i);
}

    
template<class T>
inline T&
Array_3d<T>::at(const int i, const int j, const int k)
{
    if (bounds_check(i,j,k)) {
        return this->data_ptr[offset(i,j,k)];
    } else {
        throw
        std::out_of_range("Bounds check failed");
    }
}

template<class T>
inline T
Array_3d<T>::at(const int i, const int j, const int k) const
{
    if (bounds_check(i,j,k)) {
        return this->data_ptr[offset(i,j,k)];
    } else {
        throw
        std::out_of_range("Bounds check failed");
    }
}

template<class T>
inline T&
Array_3d<T>::operator()(const int i, const int j, const int k)
{
    return this->data_ptr[offset(i,j,k)];
}

template<class T>
inline T
Array_3d<T>::operator()(const int i, const int j, const int k) const
{
    return this->data_ptr[offset(i,j,k)];
}

template<class T>
inline bool
Array_3d<T>::bounds_check(const int i, const int j, const int k) const
{
    if ( (i<0) || (i>this->shape[0]) ||
        (j<0) || (j>this->shape[1]) ||
        (k<0) || (k>this->shape[2]) ) {
        return false;
    }
    return true;
}

template<class T>
Array_3d<T>::~Array_3d()
{
}

#endif
