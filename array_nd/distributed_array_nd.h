// -*- C++ -*-
#ifndef HAVE_DISTRIBUTED_ARRAY_ND_H
#define HAVE_DISTRIBUTED_ARRAY_ND_H true

#include <vector>
#include <string>
#include <stdexcept>
#include <memory>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "range.h"

#include "array_nd.h"

template<class T>
class Distributed_array_nd;

template<class T>
class Distributed_array_nd
{
protected:
	Array_nd<T> local_array;
    std::vector<int> global_min,global_max;
    std::vector<int> local_min,local_max;
    std::vector<int> global_guard_min,global_guard_max;
    std::vector<int> local_guard_min,local_guard_max;

    void construct(const std::vector<int> global_min,
    		const std::vector<int> global_max,
    		const std::vector<int> local_min,
    		const std::vector<int> local_max,
    		const std::vector<bool> periodic,
    		const int guard_size,
    		T *data_ptr=0);
    void copy_construct(const Distributed_array_nd& original);

public:
    Distributed_array_nd();
    Distributed_array_nd(const std::vector<int> global_min,
    		const std::vector<int> global_max,
    		const std::vector<int> local_min,
    		const std::vector<int> local_max,
    		const std::vector<bool> periodic,
    		const int guard_size,
    		T *data_ptr=0);
    Distributed_array_nd(const Distributed_array_nd& original);

    void reshape(const std::vector<int> global_min,
    		const std::vector<int> global_max,
			const std::vector<int> local_min,
			const std::vector<int> local_max,
			const std::vector<bool> periodic,
			const int guard_size,
			T *data_ptr=0);

    void set_all(const T value);
    void scale(const T factor);
    void add(const T constant);

    T& at(const std::vector<int> &indices);
    T at(const std::vector<int> &indices) const;
    T& operator()(const std::vector<int> &indices);
    T operator()(const std::vector<int> &indices) const;

    T* get_data_ptr() const;
    bool owns_data() const;
    std::vector<int> get_global_min() const;
    std::vector<int> get_global_max() const;
    std::vector<int> get_local_min() const;
    std::vector<int> get_local_max() const;
    std::vector<int> get_global_guard_min() const;
    std::vector<int> get_global_guard_max() const;
    std::vector<int> get_local_guard_min() const;
    std::vector<int> get_local_guard_max() const;
    std::vector<bool> get_periodic() const;
    int get_rank() const;
    inline int offset(const std::vector<int> &indices) const;
    inline bool bounds_check(const std::vector<int> &indices) const;

    Array_nd<T>& get_local_array();
    Array_nd<T>& get_local_array_wguards();

    void describe() const;
    virtual ~Distributed_array_nd();
};

#include "distributed_array_nd.tcc"

#endif
