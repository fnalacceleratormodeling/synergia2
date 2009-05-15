#include <iostream>
#include <iomanip>

template<class T>
Distributed_array_nd<T>::Distributed_array_nd()
{
}

template<class T>
void
Distributed_array_nd<T>::construct(const std::vector<int> global_min,
		const std::vector<int> global_max,
		const std::vector<int> local_min,
		const std::vector<int> local_max,
		const std::vector<bool> periodic,
		const int guard_size,
		T *data_ptr)
{
	if ((global_min.size() != global_max.size()) ||
			(global_min.size() != local_min.size()) ||
			(global_min.size() != local_max.size()) ||
			(global_min.size() != periodic.size())) {
		throw std::runtime_error("Distributed_array_nd must be constructed or reshaped with same lengths for global_min,global_max,local_min,local_max and periodic.");
	}
	this->global_min = global_min;
	this->global_max = global_max;
	this->local_min = local_min;
	this->local_max = local_max;
	this->periodic = periodic;
	std::vector<int> local_shape(global_min.size());
	for(int i=0; i<global_min.size(); ++i) {
		if (periodic[i]) {
			global_guard_min[i] = global_min[i] - guard_size;
			global_guard_max[i] = global_max[i] + guard_size;
		} else {
			global_guard_min[i] = global_min[i];
			global_guard_max[i] = global_max[i];
		}
		local_guard_min[i] = local_guard_min[i] - guard_size;
		if (local_guard_min[i] < global_guard_min[i]) {
			local_guard_min[i] = global_guard_min[i];
		}
		local_guard_max[i] = local_guard_max[i] + guard_size;
		if (local_guard_max[i] > global_guard_max[i]) {
			local_guard_max[i] = global_guard_max[i];
		}
		local_shape[i] = local_guard_max[i] - local_guard_min[i];
	}
	local_array.reshape(local_shape,data_ptr);
}

template<class T>
Distributed_array_nd<T>::Distributed_array_nd(const std::vector<int> global_min,
	const std::vector<int> global_max,
	const std::vector<int> local_min,
	const std::vector<int> local_max,
	const std::vector<bool> periodic,
	const int guard_size,
	T *data_ptr)
{
    construct(global_min,global_max,local_min,local_max,periodic,guard_size,data_ptr);
}

template<class T>
void
Distributed_array_nd<T>::copy_construct(const Distributed_array_nd& original)
{
	throw std::runtime_error("Distributed_array_nd copy constructor not yet implemented.");
//#if defined(DEBUG_ALL) || defined(DEBUG_ARRAY_ND_ALL) || defined(DEBUG_ARRAY_ND_COPY_CTOR)
//    std::cout << "calling Distributed_array_nd copy constructor:";
//#endif
//    data_ptr = original.data_ptr;
//    T* original_data_ptr;
//    if (original.own_data){
//    	original_data_ptr = 0;
//    } else {
//    	original_data_ptr = original.data_ptr;
//    }
//    construct(original.shape, original.strides, original_data_ptr);
}

template<class T>
Distributed_array_nd<T>::Distributed_array_nd(const Distributed_array_nd& original)
{
    copy_construct(original);
}

template<class T>
void
Distributed_array_nd<T>::reshape(const std::vector<int> global_min,
		const std::vector<int> global_max,
		const std::vector<int> local_min,
		const std::vector<int> local_max,
		const std::vector<bool> periodic,
		const int guard_size,
		T *data_ptr)
{
	construct(global_min,global_max,local_min,local_max,periodic,guard_size,data_ptr);
}

template<class T>
std::vector<int>
Distributed_array_nd<T>::get_global_min() const
{
    return global_min;
}

template<class T>
std::vector<int>
Distributed_array_nd<T>::get_global_max() const
{
    return global_max;
}

template<class T>
std::vector<int>
Distributed_array_nd<T>::get_local_min() const
{
    return local_min;
}

template<class T>
std::vector<int>
Distributed_array_nd<T>::get_local_max() const
{
    return local_max;
}

template<class T>
std::vector<int>
Distributed_array_nd<T>::get_global_guard_min() const
{
    return global_guard_min;
}

template<class T>
std::vector<int>
Distributed_array_nd<T>::get_global_guard_max() const
{
    return global_guard_max;
}

template<class T>
std::vector<int>
Distributed_array_nd<T>::get_local_guard_min() const
{
    return local_guard_min;
}

template<class T>
std::vector<int>
Distributed_array_nd<T>::get_local_guard_max() const
{
    return local_guard_max;
}

template<class T>
int
Distributed_array_nd<T>::get_rank() const
{
    return global_min.size();
}

template<class T>
void
Distributed_array_nd<T>::set_all(const T value)
{
	local_array.set_all(value);
}

template<class T>
void
Distributed_array_nd<T>::scale(const T factor)
{
	local_array.scale(factor);
}

template<class T>
void
Distributed_array_nd<T>::add(const T constant)
{
	local_array.add(constant);
}

template<class T>
inline T&
Distributed_array_nd<T>::at(const std::vector<int> &indices)
{
	std::vector<int> la_indices(global_min.size());
	for(int i=0; i<global_min.size(); ++i) {
		la_indices[i] = indices[i] - local_guard_min[i];
	}
	return local_array.at(la_indices);
}

template<class T>
inline T
Distributed_array_nd<T>::at(const std::vector<int> &indices) const
{
	std::vector<int> la_indices(global_min.size());
	for(int i=0; i<global_min.size(); ++i) {
		la_indices[i] = indices[i] - local_guard_min[i];
	}
	return local_array.at(la_indices);
}

template<class T>
inline T&
Distributed_array_nd<T>::operator()(const std::vector<int> &indices)
{
	std::vector<int> la_indices(global_min.size());
	for(int i=0; i<global_min.size(); ++i) {
		la_indices[i] = indices[i] - local_guard_min[i];
	}
	return local_array[la_indices];
}

template<class T>
inline T
Distributed_array_nd<T>::operator()(const std::vector<int> &indices) const
{
	std::vector<int> la_indices(global_min.size());
	for(int i=0; i<global_min.size(); ++i) {
		la_indices[i] = indices[i] - local_guard_min[i];
	}
	return local_array[la_indices];
}

template<class T>
inline int
Distributed_array_nd<T>::offset(const std::vector<int> &indices) const
{
	std::vector<int> la_indices(global_min.size());
	for(int i=0; i<global_min.size(); ++i) {
		la_indices[i] = indices[i] + local_guard_min[i];
	}
	return local_array.offset(la_indices);
}

template<class T>
inline bool
Distributed_array_nd<T>::bounds_check(const std::vector<int> &indices) const
{
	std::vector<int> la_indices(global_min.size());
	for(int i=0; i<global_min.size(); ++i) {
		la_indices[i] = indices[i] + local_guard_min[i];
	}
	return local_array.bounds_check(la_indices);
}

template<class T>
Array_nd<T>&
Distributed_array_nd<T>::get_local_array()
{
	std::vector<Range> ranges(global_min.size());
	for(int i=0; i<global_min.size(); ++i) {
		ranges[i].set(-local_guard_min[i],local_max[i]-local_guard_min[i]);
	}
	return local_array.slice(ranges);
}

template<class T>
Array_nd<T>&
Distributed_array_nd<T>::get_local_array_wguards()
{
	return local_array;
}

template<class T>
void
Distributed_array_nd<T>::describe() const
{
	std::cout << " would describe this array as... distributed\n";
//    std::cout << "array has order " << shape.size() << " ";
//    std::cout << "and dimensions (";
//    for (unsigned int i = 0; i < shape.size(); i++) {
//        if (i > 0) {
//            std::cout << ",";
//        }
//        std::cout << shape[i];
//    }
//    std::cout << ")\n";
}

template<class T>
Distributed_array_nd<T>::~Distributed_array_nd()
{
}
