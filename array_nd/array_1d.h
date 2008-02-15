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

#include "array_1d.tcc"

#endif
