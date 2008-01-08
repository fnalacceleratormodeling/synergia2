#ifndef HAVE_ARRAY_1d_H
#define HAVE_ARRAY_1d_H true

#include "array_nd.h"
#include "vector_helper.h"

template<class T>
class Array_1d : public Array_nd<T>
{
public:
    Array_1d();
    Array_1d(const int n);
    Array_1d(const int n, const int stride);
    Array_1d(const int n, T *data_ptr);
    Array_1d(const int n, const int stride, T *data_ptr);
    Array_1d(const Array_nd<T> &array_nd);

    void reshape(const int n);
    void reshape(const int n, const int stride);
    void reshape(const int n, T *data_ptr);
    void reshape(const int n, const int stride, T *data_ptr);

    T& at(const int i);
    T at(const int i) const;
    T& operator()(const int i);
    T operator()(const int i) const;

    inline int offset(const int i) const;
    inline bool bounds_check(const int i) const;

    virtual ~Array_1d();
};

// Begin implementation.
#include <iostream>

template<class T>
Array_1d<T>::Array_1d() : Array_nd<T>()
{
    this->shape = vector1(0);
    this->strides = vector1(0);
}

template<class T>
Array_1d<T>::Array_1d(const int n) : Array_nd<T>()
{
    this->construct(vector1(n),true);
}

template<class T>
Array_1d<T>::Array_1d(const int n, const int stride) : Array_nd<T>()
{
    this->construct(vector1(n),vector1(stride),true);
}

template<class T>
Array_1d<T>::Array_1d(const int n, T *data_ptr) : Array_nd<T>()
{
    this->data_ptr = data_ptr;
    this->construct(vector1(n),false);
}

template<class T>
Array_1d<T>::Array_1d(const int n, const int stride, T *data_ptr) : Array_nd<T>()
{
    this->data_ptr = data_ptr;
    this->construct(vector1(n),vector1(stride),false);
}

template<class T>
Array_1d<T>::Array_1d(const Array_nd<T>& original)
{
    //~ std::cout << "converting from nd to 1d\n";
    if (original.get_rank() != 1) {
        throw 
            std::runtime_error("Attempt to convert Array_nd of rank !=1 to Array_1d");
    }
    copy_construct(original);
}

template<class T>
void
Array_1d<T>::reshape(const int n)
{
    bool shape_changed = this->different_shape(vector1(n));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_1d");
    }
    if (shape_changed) {
        this->construct(vector1(n),true);
    }
}

template<class T>
void
Array_1d<T>::reshape(const int n, const int stride)
{
    bool shape_changed = this->different_shape(vector1(n));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_1d");
    }
    if (shape_changed) {
        this->construct(vector1(n),vector1(stride),true);
    }
}

template<class T>
void
Array_1d<T>::reshape(const int n, T *data_ptr)
{
    bool shape_changed = this->different_shape(vector1(n));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_1d");
    }
    this->data_ptr = data_ptr;
    if (shape_changed) {
        this->construct(vector1(n),false);
    }
}

template<class T>
void
Array_1d<T>::reshape(const int n, const int stride, T *data_ptr)
{
    bool shape_changed = this->different_shape(vector1(n));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_1d");
    }
    this->data_ptr = data_ptr;
    if (shape_changed) {
        this->construct(vector1(n),vector1(stride),false);
    }
}

template<class T>
inline int
Array_1d<T>::offset(const int i) const
{
    return i*this->strides[0];
}

    
template<class T>
inline T&
Array_1d<T>::at(const int i)
{
    if (bounds_check(i)) {
        return this->data_ptr[offset(i)];
    } else {
        throw
        std::out_of_range("Bounds check failed");
    }
}

template<class T>
inline T
Array_1d<T>::at(const int i) const
{
    if (bounds_check(i)) {
        return this->data_ptr[offset(i)];
    } else {
        throw
        std::out_of_range("Bounds check failed");
    }
}

template<class T>
inline T&
Array_1d<T>::operator()(const int i)
{
    return this->data_ptr[offset(i)];
}

template<class T>
inline T
Array_1d<T>::operator()(const int i) const
{
    return this->data_ptr[offset(i)];
}

template<class T>
inline bool
Array_1d<T>::bounds_check(const int i) const
{
    if ((i<0) || (i>this->shape[0])) {
        return false;
    }
    return true;
}

template<class T>
Array_1d<T>::~Array_1d()
{
}

#endif
