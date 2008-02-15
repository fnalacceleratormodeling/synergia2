#ifndef HAVE_ARRAY_2D_H
#define HAVE_ARRAY_2D_H true

#include "array_nd.h"
#include "vector_helper.h"

template<class T>
class Array_2d : public Array_nd<T>
{
public:
    //~ class Iterator : public Array_nd<T>::Iterator
    //~ {
        //~ protected:
            //~ int step0,max0,current_index0;
            //~ int step1,max1,current_index1;
        //~ public:
            //~ Iterator(T* begin_ptr, std::vector<int> &shape,
                //~ std::vector<int> &strides);
             //~ void operator++();
    //~ };        
    Array_2d();
    Array_2d(const int nx, const int ny);
    Array_2d(const int nx, const int ny, 
        const int stride_x, const int stride_y);
    Array_2d(const int nx, const int ny, T *data_ptr);
    Array_2d(const int nx, const int ny, 
        const int stride_x, const int stride_y, T *data_ptr);
    Array_2d(const Array_nd<T> &array_nd);

    void reshape(const int nx, const int ny);
    void reshape(const int nx, const int ny, 
        const int stride_x, const int stride_y);
    void reshape(const int nx, const int ny, T *data_ptr);
    void reshape(const int nx, const int ny, 
        const int stride_x, const int stride_y, T *data_ptr);

    //~ Iterator begin();
    //~ Iterator end();

    T& at(const int i, const int j);
    T at(const int i, const int j) const;
    T& operator()(const int i, const int j);
    T operator()(const int i, const int j) const;

    inline int offset(const int i, const int j) const;
    inline bool bounds_check(const int i, const int j) const;

    virtual ~Array_2d();
};

//~ #include "array_2d_Iterator.tcc"
#include "array_2d.tcc"

#endif
