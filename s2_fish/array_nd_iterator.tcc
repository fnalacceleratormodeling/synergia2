#include <iostream>

template<class T>
Array_nd<T>::iterator::iterator(T* begin_ptr, std::vector<int> &shape,
                std::vector<int> &strides)
{
    this->begin_ptr = begin_ptr;
    this->shape = &shape;
    if (strides[shape.size()-1] < strides[0]) {
        index_start = shape.size()-1;
        index_end = -1;
        index_step = -1;
    } else {
        index_start = 0;
        index_end = shape.size();
        index_step = +1;
    } 
    this->strides = &strides;
    for(int i=0; i<shape.size(); ++i) {
        current_indices.push_back(0);
        //~ std::cout << "strides["<<i<<"] = " << strides[i] << std::endl;
    }
    current_ptr = begin_ptr;
}

template<class T>
inline T&
Array_nd<T>::iterator::operator*()
{
    return *current_ptr;
}

//~ template<class T>
//~ inline void
//~ Array_nd<T>::iterator::operator++()
//~ {
    //~ current_ptr += (*strides)[0];
    //~ ++current_indices[0];
    //~ if (current_indices[0] == (*shape)[0]) {
        //~ current_indices[0] = 0;
        //~ current_ptr += (*strides)[1] - (*shape)[0] * (*strides)[0];
        //~ ++current_indices[1];
        //~ if (current_indices[1] == (*shape)[1]) {
            //~ current_ptr = 0;
        //~ }
    //~ }
//~ }

template<class T>
inline void
Array_nd<T>::iterator::operator++()
{
    int index = index_start;
    bool roll_over = true;
    while(roll_over) {
        ++current_indices[index];
        if (current_indices[index] >= (*shape)[index]) {
            current_indices[index] = 0;
            current_ptr -= ((*shape)[index] - 1) * (*strides)[index];
            index += index_step;
            if (index == index_end) {
                current_ptr = 0;
                roll_over=false;
            }
        } else {
            current_ptr += (*strides)[index];
            roll_over = false;
        }
    }
}

//~ template<class T>
//~ inline void
//~ Array_nd<T>::iterator::operator++()
//~ {
    //~ count++;
    //~ if (count > max_count) {
        //~ current_ptr = 0;
    //~ } else {
        //~ current_ptr ++;
    //~ }
//~ }

template<class T>
inline bool
Array_nd<T>::iterator::operator==(const Array_nd<T>::iterator& other)
{
    return (other.current_ptr == current_ptr);
}

template<class T>
inline bool
Array_nd<T>::iterator::operator!=(const Array_nd<T>::iterator& other)
{
    return (other.current_ptr != current_ptr);
}

