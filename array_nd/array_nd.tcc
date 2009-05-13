#include <iostream>
#include <iomanip>

template<class T>
Array_nd<T>::Array_nd()
{
    shape_frozen = false;
    own_data = false;
    data_ptr = 0;
    end_ptr = 0;
}

template<class T>
void
Array_nd<T>::construct(const std::vector<int> shape,
                       const std::vector<int> strides, T *data_ptr)
{
    if (own_data) {
#if defined(DEBUG_ALL) || defined(DEBUG_ARRAY_ND_ALL) || defined(DEBUG_ARRAY_ND_DEALLOCATE)
        std::cout << "deallocating\n";
#endif
        myallocator.deallocate(this->data_ptr, size);
    }
    this->shape = shape;
    this->strides = strides;
    size = 1;
    std::vector<int> last_point(shape);
    for (unsigned int i = 0; i < shape.size() ; i++) {
        size *= shape[i];
        last_point[i] = shape[i] - 1;
    }
    if (offset(last_point) == size - 1) {
        contiguous = true;
    } else {
        contiguous = false;
    }
    if (data_ptr == 0) {
        this->data_ptr = myallocator.allocate(size);
        own_data = true;
    } else {
    	this->data_ptr = data_ptr;
        own_data = false;
    }
    end_ptr = this->data_ptr + offset(last_point) + 1;
}

template<class T>
void
Array_nd<T>::construct(const std::vector<int> shape, T *data_ptr)
{
    construct(shape, default_strides_from_shape(shape), data_ptr);
}

template<class T>
std::vector<int>
Array_nd<T>::default_strides_from_shape(const std::vector<int> &shape)
{
    std::vector<int> strides(shape);
    int dim = shape.size();
    strides.at(dim - 1) = 1;
    for (int i = dim - 2;i >= 0; --i) {
        strides.at(i) = strides.at(i + 1) * shape.at(i + 1);
    }
    return strides;
}

template<class T>
Array_nd<T>::Array_nd(const std::vector<int> shape, T *data_ptr)
{
    shape_frozen = false;
    own_data = false;
    construct(shape,data_ptr);
}

template<class T>
Array_nd<T>::Array_nd(const std::vector<int> shape,
                      const std::vector<int> strides, T *data_ptr)
{
    shape_frozen = false;
    own_data = false;
    construct(shape,strides,data_ptr);
}

template<class T>
void
Array_nd<T>::copy_construct(const Array_nd& original)
{
#if defined(DEBUG_ALL) || defined(DEBUG_ARRAY_ND_ALL) || defined(DEBUG_ARRAY_ND_COPY_CTOR)
    std::cout << "calling Array_nd copy constructor:";
#endif
    shape_frozen = false;
    own_data = false;
    data_ptr = original.data_ptr;
    T* original_data_ptr;
    if (original.own_data){
    	original_data_ptr = 0;
    } else {
    	original_data_ptr = original.data_ptr;
    }
    construct(original.shape, original.strides, original_data_ptr);
}

template<class T>
Array_nd<T>::Array_nd(const Array_nd& original)
{
    copy_construct(original);
}

template<class T>
void
Array_nd<T>::copy()
{
    if (not own_data) {
        T * tmp_ptr = data_ptr;
        data_ptr = myallocator.allocate(size);
        own_data = true;
        for (unsigned int i = 0; i < size; ++i) {
            data_ptr[i] = tmp_ptr[i];
        }
    }
}

template<class T>
bool
Array_nd<T>::different_shape(const std::vector<int> shape) const
{
    bool shape_changed = false;
    if (shape.size() == this->shape.size()) {
        for (unsigned int i = 0; i < shape.size(); ++i) {
            if (shape[i] != this->shape[i]) {
                shape_changed = true;
            }
        }
    } else {
        shape_changed = true;
    }
    return shape_changed;
}

template<class T>
void
Array_nd<T>::reshape(const std::vector<int> shape, T *data_ptr)
{
    reshape(shape, default_strides_from_shape(shape), data_ptr);
}

template<class T>
void
Array_nd<T>::reshape(const std::vector<int> shape,
                     const std::vector<int> strides, T *data_ptr)
{
    bool shape_changed = different_shape(shape);
    if (shape_changed && shape_frozen) {
        throw
        std::runtime_error("Attempt to change the shape of a frozen Array_nd");
    }
    if (shape_changed) {
        construct(shape, strides, data_ptr);
    } else {
    	if (data_ptr != 0) {
    		this->data_ptr = data_ptr;
    	}
    }
}

template<class T>
void
Array_nd<T>::freeze_shape()
{
    shape_frozen = true;
}

template<class T>
bool
Array_nd<T>::shape_is_frozen() const
{
    return shape_frozen;
}

template<class T>
std::vector<int>
Array_nd<T>::get_shape() const
{
    return shape;
}

template<class T>
int
Array_nd<T>::get_rank() const
{
    return shape.size();
}

template<class T>
std::vector<int>
Array_nd<T>::get_strides() const
{
    return strides;
}

template<class T>
void
Array_nd<T>::set_all(const T value)
{
    for (Iterator it = begin(); it != end(); ++it) {
        *it = value;
    }
}

template<class T>
void
Array_nd<T>::scale(const T factor)
{
    for (Iterator it = begin(); it != end(); ++it) {
        *it *= factor;
    }
}

template<class T>
void
Array_nd<T>::add(const T constant)
{
    for (Iterator it = begin(); it != end(); ++it) {
        *it += constant;
    }
}

template<class T>
inline int
Array_nd<T>::offset(const std::vector<int> &indices) const
{
    int val = 0;
    for (int i = 0; i < shape.size() ; ++i) {
        val += indices[i] * strides[i];
    }
    return val;
}


template<class T>
inline T&
Array_nd<T>::at(const std::vector<int> &indices)
{
    if (bounds_check(indices)) {
        return data_ptr[offset(indices)];
    } else {
        throw
        std::out_of_range("Bounds check failed");
    }
}

template<class T>
inline T
Array_nd<T>::at(const std::vector<int> &indices) const
{
    if (bounds_check(indices)) {
        return data_ptr[offset(indices)];
    } else {
        throw
        std::out_of_range("Bounds check failed");
    }
}

template<class T>
inline T&
Array_nd<T>::operator()(const std::vector<int> &indices)
{
    return data_ptr[offset(indices)];
}

template<class T>
inline T
Array_nd<T>::operator()(const std::vector<int> &indices) const
{
    return data_ptr[offset(indices)];
}

template<class T>
inline bool
Array_nd<T>::bounds_check(const std::vector<int> &indices) const
{
    if (shape.size() == 0) {
        return false;
    }
    for (unsigned int i = 0; i < shape.size(); ++i) {
        if ((indices[i] < 0) || (indices[i] >= shape[i])) {
            return false;
        }
    }
    return true;
}

template<class T>
T*
Array_nd<T>::get_data_ptr() const
{
    return data_ptr;
}

template<class T>
bool
Array_nd<T>::owns_data() const
{
    return own_data;
}

template<class T>
typename Array_nd<T>::Iterator
Array_nd<T>::begin()
{
    return Iterator(data_ptr, end_ptr, contiguous, shape, strides);
}

template<class T>
T*
Array_nd<T>::end()
{
    return end_ptr;
}

template<class T>
Array_nd<T>
Array_nd<T>::slice(std::vector<Range> ranges)
{
    if (ranges.size() != get_rank()) {
        throw std::runtime_error("Array_nd.slice(ranges) called with incorrect ranges length.");
    }
    int new_rank = 0;
    for (int i = 0; i < ranges.size(); ++i) {
        if (! ranges.at(i).is_unit_length()) {
            ++new_rank;
        }
    }
    std::vector<int> new_shape(new_rank);
    std::vector<int> new_strides(new_rank);
    std::vector<int> offset_indices(get_rank());
    int which_new = 0;
    for (int i = 0; i < ranges.size(); ++i) {
        int min = ranges.at(i).get_first(0, shape.at(i) - 1);
        int max = ranges.at(i).get_last(0, shape.at(i) - 1);
        if ( (min < 0) || (max >= shape.at(i)) ) {
            throw std::runtime_error("Array_nd.slice: requested range outside range of parent array.");
        }
        offset_indices.at(i) = min;
        if (! ranges.at(i).is_unit_length()) {
            new_shape.at(which_new) = (max - min) / ranges.at(i).get_step() + 1;
            new_strides.at(which_new) = strides.at(i) *
                                        ranges.at(i).get_step();
            ++which_new;
        }
    }
    return Array_nd<T>(new_shape, new_strides,
                            get_data_ptr() + offset(offset_indices));
}

template<class T>
void
Array_nd<T>::describe() const
{
    std::cout << "array has order " << shape.size() << " ";
    std::cout << "and dimensions (";
    for (unsigned int i = 0; i < shape.size(); i++) {
        if (i > 0) {
            std::cout << ",";
        }
        std::cout << shape[i];
    }
    std::cout << ")\n";
}

template<class T>
void
Array_nd<T>::recursive_print(const std::string name,
                             const int which_index, std::vector<int> indices) const
{
    if (which_index == 2) {
        std::cout << name << "(:,:";
        for (unsigned int i = 2; i < shape.size(); i++) {
            std::cout << "," << indices[i];
        }
        std::cout << ")" << std::endl;
        for (int i = 0; i < shape[0]; i++) {
            indices[0] = i;
            for (int j = 0; j < shape[1]; j++) {
                indices[1] = j;
                std::cout << std::setw(12);
                std::cout << this->operator()(indices);
            }
            std::cout << std::endl;
        }
    } else {
        for (int i = 0; i < shape[which_index-1]; i++) {
            indices[which_index-1] = i;
            recursive_print(name, which_index - 1, indices);
            std::cout << std::endl;
        }
    }
}

template<class T>
void
Array_nd<T>::print(const std::string name) const
{
    if (shape.size() == 0) {
        std::cout << name << ": empty\n";
    } else if (shape.size() == 1) {
        std::cout << name << ":\n";
        std::vector<int> index(1);
        for (int i = 0; i < shape[0]; ++i) {
            index[0] = i;
            std::cout << std::setw(12);
            std::cout << this->operator()(index) << std::endl;
        }
    } else {
        std::vector<int> indices(shape.size());
        recursive_print(name, shape.size(), indices);
    }
}

template<class T>
void
Array_nd<T>::write_to_fstream(std::ofstream& stream)
{
    for (unsigned int i = 0; i < shape.size(); ++i) {
        stream << shape[i];
        if (i < shape.size() - 1) {
            stream << " ";
        }
    }
    stream << std::endl;

    for (unsigned int i = 0; i < shape.size(); ++i) {
        stream << strides[i];
        if (i < strides.size() - 1) {
            stream << " ";
        }
    }
    stream << std::endl;

    for (unsigned int i = 0; i < size; ++i ) {
        stream << data_ptr[i] << std::endl;
    }
}

template<class T>
void
Array_nd<T>::write_to_file(std::string filename)
{
    std::ofstream stream(filename.c_str());
    write_to_fstream(stream);
    stream.close();
}

template<class T>
std::vector<T>
read_vector_line(std::vector<T>& v, std::ifstream& stream)
{
    const int buffer_length(1000);
    char buffer[1000];
    stream.getline(buffer, buffer_length);
    std::stringstream sstream(buffer);
    while (! sstream.eof()) {
        T val;
        sstream >> val;
        v.push_back(val);
    }
    return v;
}

template<class T>
void
Array_nd<T>::read_from_fstream(std::ifstream& stream)
{
    std::vector<int> read_shape;
    read_vector_line(read_shape, stream);
    std::vector<int> read_strides;
    read_vector_line(read_strides, stream);
    int size = 1;
    for (int i = 0; i < read_shape.size(); ++i) {
        size *= read_shape[i];
    }
    reshape(shape, strides);

    for (unsigned int i = 0; i < size; ++i ) {
        stream >> data_ptr[i];
    }
}

template<class T>
void
Array_nd<T>::read_from_file(std::string filename)
{
    std::ifstream stream(filename.c_str());
    read_from_fstream(stream);
    stream.close();
}


template<class T>
Array_nd<T>::~Array_nd()
{
    if (own_data) {
        myallocator.deallocate(data_ptr, size);
    }
}
