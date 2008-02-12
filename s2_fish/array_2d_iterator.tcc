template<class T>
Array_2d<T>::iterator::iterator(T* begin_ptr, std::vector<int> &shape,
                std::vector<int> &strides) : 
    Array_nd<T>::iterator::iterator(begin_ptr,shape,strides)
{
    step0 = (*this->strides)[0];
    max0 = (*this->shape)[0];
    current_index0 = this->current_indices[0];
    step1 = (*this->strides)[1] - step0*max0;
    max1 = (*this->shape)[1];
    current_index1 = this->current_indices[1];
}

template<class T>
void
Array_2d<T>::iterator::operator++()
{
    this->current_ptr += step0;
    ++current_index0;
    if (current_index0== max0) {
        current_index0 = 0;
        this->current_ptr += step1;
        ++current_index1;
        if(current_index1 == max1) {
            this->current_ptr = 0;
        }
    }
}
