template<class T>
Array_1d<T>::Array_1d() : Array_nd<T>()
{
    this->shape = vector1(0);
    this->strides = vector1(0);
}

template<class T>
Array_1d<T>::Array_1d(const int n, T *data_ptr) : Array_nd<T>()
{
	bool allocate;
	if (data_ptr == 0) {
		allocate = true;
	} else {
		allocate = false;
	    this->data_ptr = data_ptr;
	}
    this->construct(vector1(n),allocate);
}

template<class T>
Array_1d<T>::Array_1d(const int n, const int stride, T *data_ptr) : Array_nd<T>()
{
	bool allocate;
	if (data_ptr == 0) {
		allocate = true;
	} else {
		allocate = false;
	    this->data_ptr = data_ptr;
	}
    this->construct(vector1(n),vector1(stride),allocate);
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
Array_1d<T>::reshape(const int n, T *data_ptr)
{
    bool shape_changed = this->different_shape(vector1(n));
    if (shape_changed && this->shape_frozen) {
        throw
        std::out_of_range("Attempt to change the shape of a frozen Array_1d");
    }
    bool allocate;
    if (data_ptr == 0) {
    	allocate = true;
    } else {
    	allocate = false;
        this->data_ptr = data_ptr;
    }
    if (shape_changed) {
        this->construct(vector1(n),allocate);
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
    bool allocate;
    if (data_ptr == 0) {
    	allocate = true;
    } else {
    	allocate = false;
        this->data_ptr = data_ptr;
    }
    if (shape_changed) {
        this->construct(vector1(n),vector1(stride),allocate);
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
inline int
Array_1d<T>::get_length() const
{
    return this->shape[0];
}

template<class T>
inline bool
Array_1d<T>::bounds_check(const int i) const
{
    if ((i<0) || (i>=this->shape[0])) {
        return false;
    }
    return true;
}

template<class T>
Array_1d<T>::~Array_1d()
{
}
