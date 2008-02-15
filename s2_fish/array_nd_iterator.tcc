#include <iostream>

template<class T>
Array_nd<T>::Iterator::Iterator(T* begin_ptr, T* end_ptr, bool contiguous,
                                std::vector<int> &shape, std::vector<int> &strides)
{
    current_ptr = begin_ptr;
    this->contiguous = contiguous;
    this->begin_ptr = begin_ptr;
    this->end_ptr = end_ptr;
    this->shape = &shape;
    if (strides[shape.size()-1] < strides[0]) {
        index_start = shape.size() - 1;
        index_end = -1;
        index_step = -1;
    } else {
        index_start = 0;
        index_end = shape.size();
        index_step = + 1;
    }
    this->strides = &strides;
    for (int i = 0; i < shape.size(); ++i) {
        current_indices.push_back(0);
    }
}

template<class T>
inline T&
Array_nd<T>::Iterator::operator*()
{
    return *current_ptr;
}

template<class T>
inline void
Array_nd<T>::Iterator::noncontig_plusplus()
{
    int index = index_start;
    bool roll_over = true;
    while (roll_over) {
        ++current_indices[index];
        if (current_indices[index] >= (*shape)[index]) {
            current_indices[index] = 0;
            this->current_ptr -= ((*shape)[index] - 1) * (*strides)[index];
            index += index_step;
            if (index == index_end) {
                this->current_ptr = this->end_ptr;
                roll_over = false;
            }
        } else {
            this->current_ptr += (*strides)[index];
            roll_over = false;
        }
    }
}

template<class T>
inline void
Array_nd<T>::Iterator::operator++()
{
    if (contiguous) {
        ++current_ptr;
    } else {
        noncontig_plusplus();
    }
}

template<class T>
inline bool
Array_nd<T>::Iterator::operator==(const Array_nd<T>::Iterator& other)
{
    return (other.current_ptr == current_ptr);
}

template<class T>
inline bool
Array_nd<T>::Iterator::operator!=(const Array_nd<T>::Iterator& other)
{
    return (other.current_ptr != current_ptr);
}

template<class T>
inline bool
Array_nd<T>::Iterator::operator==(const T *ptr)
{
    return (ptr == current_ptr);
}

template<class T>
inline bool
Array_nd<T>::Iterator::operator!=(const T *ptr)
{
    return (ptr != current_ptr);
}
