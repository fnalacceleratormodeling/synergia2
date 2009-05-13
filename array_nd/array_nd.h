// -*- C++ -*-
#ifndef HAVE_ARRAY_ND_H
#define HAVE_ARRAY_ND_H true

#include <vector>
#include <string>
#include <stdexcept>
#include <memory>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "range.h"

template<class T>
class Array_nd;

template<class T>
class Array_nd
{
protected:
    T *data_ptr, *end_ptr;
    bool own_data;
    std::vector<int> shape;
    std::vector<int> strides;
    bool shape_frozen;
    bool contiguous;
    int size;

    std::allocator<T> myallocator;

    void construct(const std::vector<int> shape,
                   const std::vector<int> strides, const bool allocate);
    void construct(const std::vector<int> shape, const bool allocate);
    void copy_construct(const Array_nd& original);
    bool different_shape(const std::vector<int> shape) const;
    std::vector<int>
    default_strides_from_shape(const std::vector<int> &shape);
    void recursive_print(const std::string name, const int which_index,
                         std::vector<int> indices) const;

public:
    class Iterator
    {
    private:
        void noncontig_plusplus();
    protected:
        bool contiguous;
        T* begin_ptr;
        T *end_ptr;
        T* current_ptr;
        std::vector<int> current_indices;
        std::vector<int> *shape;
        int index_start, index_end, index_step;
        std::vector<int> *strides;
    public:
        Iterator(T* begin_ptr, T* end_ptr, bool contiguous,
                 std::vector<int> &shape, std::vector<int> &strides);
        inline T& operator*();
        void operator++();
        inline bool operator==(const Iterator& other);
        inline bool operator!=(const Iterator& other);
        inline bool operator==(const T *ptr);
        inline bool operator!=(const T *ptr);
    };

    Array_nd();
    Array_nd(const std::vector<int> shape, T *data_ptr=0);
    Array_nd(const std::vector<int> shape,
             const std::vector<int> strides, T *data_ptr=0);
    Array_nd(const Array_nd& original);

    void copy();

    void reshape(const std::vector<int> shape, T *data_ptr=0);
    void reshape(const std::vector<int> shape,
                 const std::vector<int> strides, T *data_ptr=0);
    void freeze_shape();

    void set_all(const T value);
    void scale(const T factor);
    void add(const T constant);

    T& at(const std::vector<int> &indices);
    T at(const std::vector<int> &indices) const;
    T& operator()(const std::vector<int> &indices);
    T operator()(const std::vector<int> &indices) const;

    T* get_data_ptr() const;
    bool owns_data() const;
    Iterator begin();
    T* end();
    std::vector<int> get_shape() const;
    std::vector<int> get_strides() const;
    int get_rank() const;
    bool shape_is_frozen() const;
    inline int offset(const std::vector<int> &indices) const;
    inline bool bounds_check(const std::vector<int> &indices) const;

    Array_nd<T> slice(std::vector<Range>);

    void describe() const;
    void print(const std::string name) const;

    void write_to_fstream(std::ofstream& stream);
    void write_to_file(std::string filename);
    void read_from_fstream(std::ifstream& stream);
    void read_from_file(std::string filename);

    virtual ~Array_nd();
};

#include "array_nd_iterator.tcc"
#include "array_nd.tcc"

#endif
